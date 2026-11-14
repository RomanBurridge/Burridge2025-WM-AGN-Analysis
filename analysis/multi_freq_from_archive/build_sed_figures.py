#!/usr/bin/env python3
"""Build the six manuscript continuum SED figures without network access.

The historical figures were made by ``Beam_Date_Flux_Graphs.ipynb``.  This
builder keeps its shared log-space axes, beam-size colour encoding, upper-limit
triangles, angular-scale inset, and output dimensions while replacing notebook
state and ``eval`` with a deterministic command-line workflow.
"""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import math
import shutil
import tempfile
from dataclasses import dataclass
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.cm as cmx  # noqa: E402
import matplotlib.colors  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT = SCRIPT_DIR / "fitsumfiles/final/fitsummary_final_withextras.txt"
DEFAULT_OUTPUT_DIR = SCRIPT_DIR / "SpectraGraphs"
DEFAULT_BHS_PARAMS = SCRIPT_DIR.parent / "CountingBHS_and_Spectra/BHS_params.csv"

TARGET_SOURCES = (
    "NGC4258",
    "NGC1194",
    "NGC1068",
    "NGC3079",
    "NGC4945",
    "Circinus",
)

# The notebook used these ten sources to establish common plot bounds.  Keep
# that cohort so the six manuscript panels retain their established axes.
AXIS_SOURCES = TARGET_SOURCES + ("NGC2273", "NGC4388", "NGC2960", "UGC3789")

# One arcsecond at one megaparsec, expressed in kiloparsecs.  The small-angle
# conversion is sufficiently exact at these distances and matches the values
# reported in the manuscript.
KPC_PER_ARCSEC_PER_MPC = math.radians(1.0 / 3600.0) * 1000.0

EXPECTED_SOURCE_ROWS = {
    "NGC4258": 82,
    "NGC1194": 9,
    "NGC1068": 88,
    "NGC3079": 79,
    "NGC4945": 45,
    "Circinus": 53,
}


@dataclass(frozen=True)
class ContinuumPoint:
    source: str
    log_frequency_ghz: float
    log_flux_density_mjy: float
    log_beam_arcsec: float
    upper_limit: bool


@dataclass(frozen=True)
class SharedAxes:
    frequency_limits: tuple[float, float]
    frequency_ticks: tuple[float, ...]
    flux_limits: tuple[float, float]
    flux_ticks: tuple[float, ...]
    beam_color_limits: tuple[float, float]
    beam_color_ticks: tuple[float, ...]
    point_color_limits: tuple[float, float]


@dataclass(frozen=True)
class AngularScale:
    distance_mpc: float
    distance_unc_mpc: float | None
    scale_kpc_per_arcsec: float
    scale_unc_kpc_per_arcsec: float | None
    display_value: str
    display_uncertainty: str | None


def _parse_positive(value: object, label: str, line_number: int) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Line {line_number}: invalid {label}: {value!r}") from exc
    if not math.isfinite(number) or number <= 0:
        raise ValueError(f"Line {line_number}: {label} must be positive: {value!r}")
    return number


def _parse_beam(value: object, line_number: int) -> float:
    text = str(value).strip()
    if "-" not in text:
        return _parse_positive(text, "beam", line_number)
    pieces = text.split("-")
    if len(pieces) != 2:
        raise ValueError(f"Line {line_number}: invalid beam range: {value!r}")
    low = _parse_positive(pieces[0], "beam-range lower bound", line_number)
    high = _parse_positive(pieces[1], "beam-range upper bound", line_number)
    return (low + high) / 2.0


def load_continuum_points(path: Path) -> list[ContinuumPoint]:
    """Read the canonical literal-list table and return validated log points."""

    points: list[ContinuumPoint] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, text in enumerate(handle, start=1):
            if not text.strip():
                continue
            try:
                row = ast.literal_eval(text)
            except (SyntaxError, ValueError) as exc:
                raise ValueError(f"Line {line_number}: invalid literal row") from exc
            if not isinstance(row, list) or len(row) < 6:
                raise ValueError(f"Line {line_number}: expected at least six fields")
            source = str(row[0]).strip()
            if source not in AXIS_SOURCES:
                continue
            flux_text = str(row[3]).strip()
            frequency = _parse_positive(row[2], "frequency", line_number)
            flux = _parse_positive(flux_text.rstrip("*"), "flux", line_number)
            beam = _parse_beam(row[5], line_number)
            points.append(
                ContinuumPoint(
                    source=source,
                    log_frequency_ghz=math.log10(frequency),
                    log_flux_density_mjy=math.log10(flux),
                    log_beam_arcsec=math.log10(beam),
                    upper_limit=flux_text.endswith("*"),
                )
            )

    counts = {source: 0 for source in TARGET_SOURCES}
    for point in points:
        if point.source in counts:
            counts[point.source] += 1
    missing = [source for source, count in counts.items() if count == 0]
    if missing:
        raise ValueError(f"Canonical table has no rows for: {', '.join(missing)}")
    wrong_counts = {
        source: (counts[source], expected)
        for source, expected in EXPECTED_SOURCE_ROWS.items()
        if counts[source] != expected
    }
    if wrong_counts:
        details = ", ".join(
            f"{source}={actual} (expected {expected})"
            for source, (actual, expected) in wrong_counts.items()
        )
        raise ValueError(f"Unexpected manuscript-source row counts: {details}")
    return points


def _optional_positive(value: str, label: str, source: str) -> float | None:
    text = value.strip()
    if text == "NA":
        return None
    try:
        number = float(text)
    except ValueError as exc:
        raise ValueError(f"{source}: invalid {label}: {value!r}") from exc
    if not math.isfinite(number) or number <= 0:
        raise ValueError(f"{source}: {label} must be positive: {value!r}")
    return number


def _decimal_places_for_sigfigs(value: float, significant_digits: int = 2) -> int:
    decimal_value = Decimal(str(value))
    if not decimal_value.is_finite() or decimal_value <= 0:
        raise ValueError(f"Cannot format a non-positive uncertainty: {value!r}")
    return max(0, significant_digits - 1 - decimal_value.adjusted())


def _format_fixed_half_up(value: float, decimal_places: int) -> str:
    decimal_value = Decimal(str(value))
    quantum = Decimal("1").scaleb(-decimal_places)
    rounded = decimal_value.quantize(quantum, rounding=ROUND_HALF_UP)
    return f"{rounded:.{decimal_places}f}"


def _format_angular_scale(
    scale: float, uncertainty: float | None
) -> tuple[str, str | None]:
    if uncertainty is not None:
        decimals = _decimal_places_for_sigfigs(uncertainty)
        return (
            _format_fixed_half_up(scale, decimals),
            _format_fixed_half_up(uncertainty, decimals),
        )
    decimals = _decimal_places_for_sigfigs(scale, significant_digits=3)
    return _format_fixed_half_up(scale, decimals), None


def load_angular_scales(path: Path) -> dict[str, AngularScale]:
    """Derive each inset from the canonical adopted distance table.

    ``distance_unc_mpc`` is the canonical quadrature total of the numerical
    distance components.  The SED inset has room for one uncertainty, so that
    total is used rather than silently selecting only one component.  Sources
    without a numerical distance uncertainty retain a central scale only.
    """

    scales: dict[str, AngularScale] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"source_key", "distance_mpc", "distance_unc_mpc"}
        missing = required.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(
                f"BHS parameter table is missing fields: {', '.join(sorted(missing))}"
            )
        for row in reader:
            source = row["source_key"].strip()
            if source not in TARGET_SOURCES:
                continue
            if source in scales:
                raise ValueError(f"Duplicate BHS distance row for {source}")
            distance = _optional_positive(row["distance_mpc"], "distance", source)
            if distance is None:
                raise ValueError(f"{source}: adopted distance is unavailable")
            distance_unc = _optional_positive(
                row["distance_unc_mpc"], "distance uncertainty", source
            )
            scale = distance * KPC_PER_ARCSEC_PER_MPC
            scale_unc = (
                None
                if distance_unc is None
                else distance_unc * KPC_PER_ARCSEC_PER_MPC
            )
            display_value, display_uncertainty = _format_angular_scale(
                scale, scale_unc
            )
            scales[source] = AngularScale(
                distance_mpc=distance,
                distance_unc_mpc=distance_unc,
                scale_kpc_per_arcsec=scale,
                scale_unc_kpc_per_arcsec=scale_unc,
                display_value=display_value,
                display_uncertainty=display_uncertainty,
            )

    missing_sources = sorted(set(TARGET_SOURCES).difference(scales))
    if missing_sources:
        raise ValueError(
            "BHS parameter table has no distance row for: "
            + ", ".join(missing_sources)
        )
    return scales


def _round_down_half(number: float) -> float:
    return math.floor(number * 2.0) / 2.0


def _round_up_half(number: float) -> float:
    return math.ceil(number * 2.0) / 2.0


def _float_range(start: float, stop: float, step: float) -> tuple[float, ...]:
    return tuple(float(value) for value in np.arange(start, stop, step))


def shared_axes(points: Sequence[ContinuumPoint]) -> SharedAxes:
    """Reproduce the notebook's shared half-decade bounds and five-percent pad."""

    if not points:
        raise ValueError("Cannot establish SED axes without continuum points")

    freq_min = _round_down_half(min(point.log_frequency_ghz for point in points))
    freq_max = _round_up_half(max(point.log_frequency_ghz for point in points))
    flux_min = _round_down_half(min(point.log_flux_density_mjy for point in points))
    flux_max = _round_up_half(max(point.log_flux_density_mjy for point in points))
    beam_min = _round_down_half(min(point.log_beam_arcsec for point in points))
    beam_max = _round_up_half(max(point.log_beam_arcsec for point in points))

    freq_pad = (freq_max - freq_min) * 0.05
    flux_pad = (flux_max - flux_min) * 0.05
    beam_pad = (beam_max - beam_min) * 0.05

    return SharedAxes(
        frequency_limits=(freq_min - freq_pad, freq_max + freq_pad),
        frequency_ticks=_float_range(freq_min, freq_max + 0.5, 0.5),
        flux_limits=(flux_min - flux_pad, flux_max + flux_pad),
        flux_ticks=_float_range(flux_min, flux_max + 0.5, 0.5),
        beam_color_limits=(beam_min, beam_max),
        beam_color_ticks=_float_range(beam_min, beam_max + 1.0, 1.0),
        # The historical notebook used padded beam bounds for point colours.
        point_color_limits=(beam_min - beam_pad, beam_max + beam_pad),
    )


def _notebook_colormap() -> LinearSegmentedColormap:
    jet = plt.get_cmap("jet")
    colors = jet(np.linspace(0, 1, 256))
    colors = colors[np.logical_not((colors[:, 0] > 0.9) & (colors[:, 1] > 0.9))]
    return LinearSegmentedColormap.from_list("jet_no_yellow_deterministic", colors)


def _scale_text(scale: AngularScale) -> str:
    if scale.display_uncertainty is None:
        return rf"$1^{{\prime\prime}} = {scale.display_value}\ \mathrm{{kpc}}$"
    return (
        rf"$1^{{\prime\prime}} = {scale.display_value}\,\pm\,"
        rf"{scale.display_uncertainty}\ "
        rf"\mathrm{{kpc}}$"
    )


def render_source(
    source: str,
    points: Sequence[ContinuumPoint],
    axes: SharedAxes,
    angular_scale: AngularScale,
    output_path: Path,
) -> None:
    """Render one source panel with the established notebook visual language."""

    source_points = [point for point in points if point.source == source]
    if not source_points:
        raise ValueError(f"No rows available for {source}")

    cmap = _notebook_colormap()
    point_norm = matplotlib.colors.Normalize(
        vmin=axes.point_color_limits[0], vmax=axes.point_color_limits[1]
    )
    point_map = cmx.ScalarMappable(norm=point_norm, cmap=cmap)

    detections = [point for point in source_points if not point.upper_limit]
    limits = [point for point in source_points if point.upper_limit]

    fig, ax = plt.subplots(dpi=300)
    if limits:
        ax.scatter(
            [point.log_frequency_ghz for point in limits],
            [point.log_flux_density_mjy for point in limits],
            c=[point_map.to_rgba(point.log_beam_arcsec) for point in limits],
            marker="v",
            s=10,
        )
    if detections:
        ax.scatter(
            [point.log_frequency_ghz for point in detections],
            [point.log_flux_density_mjy for point in detections],
            c=[point_map.to_rgba(point.log_beam_arcsec) for point in detections],
            s=10,
        )

    colorbar_norm = matplotlib.colors.Normalize(
        vmin=axes.beam_color_limits[0], vmax=axes.beam_color_limits[1]
    )
    colorbar = fig.colorbar(
        cmx.ScalarMappable(norm=colorbar_norm, cmap=cmap),
        ax=ax,
        ticks=axes.beam_color_ticks,
    )
    colorbar.set_label(r"log$_{10}$ Beam size [arcsec]")
    minor_locations = np.arange(
        axes.beam_color_limits[0], axes.beam_color_limits[1], 0.1
    )
    major_locations = set(float(value) for value in colorbar.get_ticks())
    colorbar.ax.yaxis.set_ticks(
        [float(value) for value in minor_locations if float(value) not in major_locations],
        minor=True,
    )

    ax.set_xlim(axes.frequency_limits)
    ax.set_xticks(axes.frequency_ticks)
    ax.set_xticklabels([f"{value:.2f}" for value in axes.frequency_ticks])
    ax.set_ylim(axes.flux_limits)
    ax.set_yticks(axes.flux_ticks)
    ax.set_yticklabels([f"{value:.2f}" for value in axes.flux_ticks])
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.grid(True)

    ax.text(
        0.98,
        0.96,
        _scale_text(angular_scale),
        transform=ax.transAxes,
        horizontalalignment="right",
        verticalalignment="top",
        fontsize=8.5,
        bbox={
            "boxstyle": "round,pad=0.25",
            "facecolor": "white",
            "edgecolor": "0.55",
            "linewidth": 0.5,
            "alpha": 0.9,
        },
    )
    ax.set_xlabel(r"log$_{10}$ Frequency [GHz]")
    ax.set_ylabel(r"log$_{10}$ Flux Density [mJy]")
    ax.set_title(f"{source} Spectra")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
        metadata={"Software": "AE8A42 deterministic SED builder"},
    )
    plt.close(fig)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_figures(
    input_path: Path,
    output_dir: Path,
    bhs_params_path: Path = DEFAULT_BHS_PARAMS,
) -> dict[str, object]:
    """Build all six PNGs and their deterministic provenance manifest."""

    points = load_continuum_points(input_path)
    angular_scales = load_angular_scales(bhs_params_path)
    axes = shared_axes(points)
    output_dir.mkdir(parents=True, exist_ok=True)

    artifacts: dict[str, dict[str, object]] = {}
    with matplotlib.rc_context(rc=matplotlib.rcParamsDefault):
        for source in TARGET_SOURCES:
            output = output_dir / f"{source}_Spectra.png"
            render_source(source, points, axes, angular_scales[source], output)
            source_points = [point for point in points if point.source == source]
            scale = angular_scales[source]
            artifacts[source] = {
                "filename": output.name,
                "rows": len(source_points),
                "detections": sum(not point.upper_limit for point in source_points),
                "upper_limits": sum(point.upper_limit for point in source_points),
                "distance_mpc": scale.distance_mpc,
                "distance_unc_mpc": scale.distance_unc_mpc,
                "scale_kpc_per_arcsec": scale.scale_kpc_per_arcsec,
                "scale_unc_kpc_per_arcsec": scale.scale_unc_kpc_per_arcsec,
                "scale_display": _scale_text(scale),
                "sha256": sha256(output),
            }

    manifest: dict[str, object] = {
        "schema_version": 2,
        "source_table": "fitsumfiles/final/fitsummary_final_withextras.txt",
        "distance_table": "analysis/CountingBHS_and_Spectra/BHS_params.csv",
        "axis_sources": list(AXIS_SOURCES),
        "frequency_limits_log10_ghz": list(axes.frequency_limits),
        "flux_density_limits_log10_mjy": list(axes.flux_limits),
        "beam_color_limits_log10_arcsec": list(axes.beam_color_limits),
        "artifacts": artifacts,
    }
    manifest_path = output_dir / "sed_figure_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def sync_revision_images(output_dir: Path, revision_dir: Path) -> None:
    """Copy the exact checked-in PNG bytes to the manuscript work directory."""

    if not revision_dir.is_dir():
        raise FileNotFoundError(f"Revision directory does not exist: {revision_dir}")
    for source in TARGET_SOURCES:
        shutil.copyfile(
            output_dir / f"{source}_Spectra.png",
            revision_dir / f"{source}_Spectra.png",
        )


def _assert_exact_artifacts(
    generated_dir: Path, expected_dir: Path, label: str
) -> None:
    for source in TARGET_SOURCES:
        filename = f"{source}_Spectra.png"
        generated = generated_dir / filename
        expected = expected_dir / filename
        if not expected.is_file():
            raise RuntimeError(f"Missing {label} artifact: {expected}")
        if generated.read_bytes() != expected.read_bytes():
            raise RuntimeError(f"Stale {label} artifact: {expected}")


def check_figures(
    input_path: Path,
    output_dir: Path,
    revision_dir: Path | None,
    bhs_params_path: Path = DEFAULT_BHS_PARAMS,
) -> None:
    """Rebuild in a temporary directory and require byte-identical artifacts."""

    with tempfile.TemporaryDirectory(prefix="ae8a42-sed-check-") as temp_dir:
        generated_dir = Path(temp_dir)
        build_figures(input_path, generated_dir, bhs_params_path)
        _assert_exact_artifacts(generated_dir, output_dir, "repository")
        expected_manifest = output_dir / "sed_figure_manifest.json"
        generated_manifest = generated_dir / "sed_figure_manifest.json"
        if not expected_manifest.is_file():
            raise RuntimeError(f"Missing repository manifest: {expected_manifest}")
        if generated_manifest.read_bytes() != expected_manifest.read_bytes():
            raise RuntimeError(f"Stale repository manifest: {expected_manifest}")
        if revision_dir is not None:
            _assert_exact_artifacts(generated_dir, revision_dir, "revision")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument(
        "--bhs-params",
        type=Path,
        default=DEFAULT_BHS_PARAMS,
        help="canonical table supplying adopted distances for the angular-scale insets",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--revision-dir",
        type=Path,
        help="also synchronize the six PNGs into this manuscript directory",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="rebuild in a temporary directory and require byte-identical outputs",
    )
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    input_path = args.input.resolve()
    bhs_params_path = args.bhs_params.resolve()
    output_dir = args.output_dir.resolve()
    revision_dir = args.revision_dir.resolve() if args.revision_dir else None

    if args.check:
        check_figures(input_path, output_dir, revision_dir, bhs_params_path)
        print(f"verified six deterministic SED PNGs in {output_dir}")
        if revision_dir is not None:
            print(f"verified exact revision copies in {revision_dir}")
        return 0

    manifest = build_figures(input_path, output_dir, bhs_params_path)
    if revision_dir is not None:
        sync_revision_images(output_dir, revision_dir)
    for source in TARGET_SOURCES:
        artifact = manifest["artifacts"][source]  # type: ignore[index]
        print(
            f"{source}: {artifact['rows']} rows, "
            f"sha256={artifact['sha256']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
