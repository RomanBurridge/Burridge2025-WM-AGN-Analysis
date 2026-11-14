#!/usr/bin/env python3
"""Build the publication black-hole-shadow-size versus 230 GHz flux figure.

This standalone builder replaces the stateful plotting cells in
``flux_size_graph.ipynb``.  It reads only the audited shadow-parameter table
and the selected 230 GHz flux table.  The accepted publication figure omits
error bars for clarity; the audited input tables retain the flux-density and
BHS uncertainty values.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402


SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY_ROOT = SCRIPT_DIR.parents[1]
DEFAULT_BHS_PATH = (
    REPOSITORY_ROOT / "analysis" / "CountingBHS_and_Spectra" / "BHS_params.csv"
)
DEFAULT_FLUX_PATH = SCRIPT_DIR / "machinetables" / "S230_extrapolation_comparison.csv"
DEFAULT_OUTPUT_DIR = SCRIPT_DIR

SPEED_OF_LIGHT = 299_792_458.0  # m s^-1, exact
RAD_TO_UAS = 180.0 / math.pi * 3_600.0 * 1_000_000.0

BHS_REQUIRED_FIELDS = {
    "source_key",
    "source_label",
    "size_uas",
    "size_stat_unc_minus_uas",
    "size_stat_unc_plus_uas",
    "size_total_unc_minus_uas",
    "size_total_unc_plus_uas",
    "size_stat_component_status",
    "size_systematic_status",
    "size_total_status",
}
FLUX_REQUIRED_FIELDS = {
    "Source",
    "flux_status",
    "S_230_mJy",
    "sigma_S_230_mJy",
}

EXPECTED_MISSING_FORMAL_SOURCES = {"NGC4945", "NGC3079"}
LABELLED_SOURCES = {
    "NGC4258",
    "NGC1194",
    "NGC3079",
    "NGC4945",
    "Circinus",
    "NGC1068",
}
LABEL_OFFSETS = {
    "NGC1068": (0.0, 0.08),
    "NGC4945": (0.0, 0.08),
    "Circinus": (0.0, -0.16),
}
FLUX_THRESHOLDS_MJY = (0.1, 1.0, 10.0, 100.0)
SIZE_STAT_COMPONENT_STATUSES = {
    "paired_chain_quantile_interval",
    "reported_undecomposed_uncertainty",
    "fit_plus_absolute_position_mixed_component",
    "fixed_distance_formal_fit_interval",
    "angular_radius_only_component",
    "no_formal_uncertainty",
    "not_available",
}
NGC4388_PRIMARY_SIZE_UAS = 0.04535053125888249
NGC4388_PRIMARY_COMPONENT_UAS = 0.0010797745537829164


@dataclass(frozen=True)
class PlotPoint:
    """One joined source record in logarithmic plotting coordinates."""

    source_key: str
    display_name: str
    flux_status: str
    size_uas: float
    flux_mjy: float
    log_size: float
    log_flux: float
    size_log_unc_minus: float | None
    size_log_unc_plus: float | None
    flux_log_unc_minus: float | None
    flux_log_unc_plus: float | None
    uncertainty_basis: str
    size_stat_component_status: str
    size_systematic_status: str
    size_total_status: str


@dataclass(frozen=True)
class PlotData:
    points: tuple[PlotPoint, ...]
    skipped_without_size: tuple[str, ...]


@dataclass(frozen=True)
class FigureBuildResult:
    pdf_path: Path
    png_path: Path
    data: PlotData


@dataclass(frozen=True)
class ReferencePoint:
    label: str
    size_uas: float
    flux_mjy: float
    size_log_unc_minus: float
    size_log_unc_plus: float
    flux_log_unc_minus: float
    flux_log_unc_plus: float


def _normalize_source(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]", "", value).upper()


def _optional_float(value: str | None) -> float | None:
    if value is None or value.strip() in {"", "NA"}:
        return None
    parsed = float(value)
    if not math.isfinite(parsed):
        raise ValueError(f"non-finite numeric value: {value!r}")
    return parsed


def _require_complete_pair(
    lower: float | None,
    upper: float | None,
    *,
    source: str,
    component: str,
) -> None:
    if (lower is None) != (upper is None):
        raise ValueError(
            f"{source}: {component} uncertainty must provide both minus and plus values"
        )


def exact_log10_uncertainties(
    value: float,
    uncertainty_minus: float,
    uncertainty_plus: float,
) -> tuple[float, float]:
    """Map asymmetric linear bounds into exact asymmetric log10 bounds."""

    if value <= 0.0:
        raise ValueError("central value must be positive")
    if uncertainty_minus < 0.0 or uncertainty_plus < 0.0:
        raise ValueError("uncertainties must be non-negative")
    lower_value = value - uncertainty_minus
    upper_value = value + uncertainty_plus
    if lower_value <= 0.0:
        raise ValueError("lower uncertainty bound must remain positive in log space")
    return (
        math.log10(value) - math.log10(lower_value),
        math.log10(upper_value) - math.log10(value),
    )


def _read_csv(path: Path, required_fields: set[str]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])
        missing = sorted(required_fields - fields)
        if missing:
            raise ValueError(f"{path} is missing required columns: {missing}")
        rows = list(reader)
    if any(None in row for row in rows):
        raise ValueError(f"{path} contains one or more rows with extra columns")
    return rows


def load_plot_data(
    bhs_path: Path = DEFAULT_BHS_PATH,
    flux_path: Path = DEFAULT_FLUX_PATH,
    *,
    enforce_publication_sample: bool = True,
) -> PlotData:
    """Join the audited size and flux tables and select plotted uncertainties."""

    bhs_rows = _read_csv(Path(bhs_path), BHS_REQUIRED_FIELDS)
    flux_rows = _read_csv(Path(flux_path), FLUX_REQUIRED_FIELDS)

    by_source: dict[str, dict[str, str]] = {}
    for row in bhs_rows:
        normalized = _normalize_source(row["source_key"])
        if normalized in by_source:
            raise ValueError(f"duplicate BHS source after normalization: {row['source_key']}")
        by_source[normalized] = row

    points: list[PlotPoint] = []
    skipped_without_size: list[str] = []
    seen_flux_sources: set[str] = set()
    for flux_row in flux_rows:
        display_name = flux_row["Source"]
        normalized = _normalize_source(display_name)
        if normalized in seen_flux_sources:
            raise ValueError(f"duplicate flux source after normalization: {display_name}")
        seen_flux_sources.add(normalized)

        if normalized not in by_source:
            raise ValueError(f"flux source has no BHS record: {display_name}")
        bhs_row = by_source[normalized]
        source_key = bhs_row["source_key"]
        if bhs_row["size_stat_component_status"] not in SIZE_STAT_COMPONENT_STATUSES:
            raise ValueError(
                f"{source_key}: unsupported size statistical component status"
            )
        size = _optional_float(bhs_row["size_uas"])
        if size is None:
            skipped_without_size.append(source_key)
            continue
        if size <= 0.0:
            raise ValueError(f"{source_key}: size_uas must be positive")

        stat_minus = _optional_float(bhs_row["size_stat_unc_minus_uas"])
        stat_plus = _optional_float(bhs_row["size_stat_unc_plus_uas"])
        total_minus = _optional_float(bhs_row["size_total_unc_minus_uas"])
        total_plus = _optional_float(bhs_row["size_total_unc_plus_uas"])
        _require_complete_pair(
            stat_minus,
            stat_plus,
            source=source_key,
            component="statistical/formal",
        )
        _require_complete_pair(
            total_minus,
            total_plus,
            source=source_key,
            component="total",
        )

        if source_key == "NGC4388":
            if (
                bhs_row["size_stat_component_status"]
                != "fit_plus_absolute_position_mixed_component"
                or bhs_row["size_systematic_status"] != "unquantified_model"
                or bhs_row["size_total_status"]
                != "incomplete_unquantified_systematic"
                or not math.isclose(
                    size,
                    NGC4388_PRIMARY_SIZE_UAS,
                    rel_tol=0.0,
                    abs_tol=1e-17,
                )
                or stat_minus is None
                or stat_plus is None
                or not math.isclose(
                    stat_minus,
                    NGC4388_PRIMARY_COMPONENT_UAS,
                    rel_tol=0.0,
                    abs_tol=1e-18,
                )
                or not math.isclose(
                    stat_plus,
                    NGC4388_PRIMARY_COMPONENT_UAS,
                    rel_tol=0.0,
                    abs_tol=1e-18,
                )
            ):
                raise ValueError(
                    "NGC4388 reviewed primary size/component semantics do not match"
                )

        total_status = bhs_row["size_total_status"]
        if total_status in {
            "quadrature_statistical_plus_systematic",
            "quadrature_mixed_fit_position_plus_systematic",
        }:
            if total_minus is None or total_plus is None:
                raise ValueError(f"{source_key}: calibrated total lacks numeric bounds")
            uncertainty_basis = "total"
            error_minus, error_plus = exact_log10_uncertainties(
                size, total_minus, total_plus
            )
        elif total_status == "quadrature_for_display_model_envelope_not_probability_interval":
            if total_minus is None or total_plus is None:
                raise ValueError(f"{source_key}: model span lacks numeric bounds")
            uncertainty_basis = "model_span"
            error_minus, error_plus = exact_log10_uncertainties(
                size, total_minus, total_plus
            )
        elif total_status in {
            "formal_component_only",
            "reported_undecomposed_uncertainty_component_only",
        }:
            if stat_minus is None or stat_plus is None:
                raise ValueError(f"{source_key}: component-only status lacks numeric bounds")
            uncertainty_basis = "component_only"
            error_minus, error_plus = exact_log10_uncertainties(
                size, stat_minus, stat_plus
            )
        elif total_status == "incomplete_unquantified_systematic":
            if stat_minus is None and stat_plus is None:
                if (
                    source_key != "NGC3079"
                    or bhs_row["size_stat_component_status"]
                    != "no_formal_uncertainty"
                    or bhs_row["size_systematic_status"] != "unquantified_model"
                ):
                    raise ValueError(
                        f"{source_key}: incomplete status has no retained component"
                    )
                uncertainty_basis = "missing"
                error_minus = None
                error_plus = None
            else:
                if stat_minus is None or stat_plus is None:
                    raise ValueError(
                        f"{source_key}: component-only status lacks numeric bounds"
                    )
                uncertainty_basis = "component_only"
                error_minus, error_plus = exact_log10_uncertainties(
                    size, stat_minus, stat_plus
                )
        elif total_status == "no_formal_uncertainty":
            uncertainty_basis = "missing"
            error_minus = None
            error_plus = None
        else:
            raise ValueError(
                f"{source_key}: unsupported size_total_status {total_status!r}"
            )

        flux = _optional_float(flux_row["S_230_mJy"])
        if flux is None or flux <= 0.0:
            raise ValueError(f"{display_name}: S_230_mJy must be positive")
        flux_status = flux_row["flux_status"]
        if flux_status not in {"detection", "upper_limit"}:
            raise ValueError(f"{display_name}: invalid flux_status {flux_status!r}")
        flux_unc = _optional_float(flux_row["sigma_S_230_mJy"])
        if flux_unc is not None and flux_unc < 0.0:
            raise ValueError(f"{display_name}: sigma_S_230_mJy must be non-negative")
        if flux_unc is None:
            flux_error_minus = None
            flux_error_plus = None
        else:
            flux_error_minus, flux_error_plus = exact_log10_uncertainties(
                flux,
                flux_unc,
                flux_unc,
            )

        points.append(
            PlotPoint(
                source_key=source_key,
                display_name=display_name,
                flux_status=flux_status,
                size_uas=size,
                flux_mjy=flux,
                log_size=math.log10(size),
                log_flux=math.log10(flux),
                size_log_unc_minus=error_minus,
                size_log_unc_plus=error_plus,
                flux_log_unc_minus=flux_error_minus,
                flux_log_unc_plus=flux_error_plus,
                uncertainty_basis=uncertainty_basis,
                size_stat_component_status=bhs_row["size_stat_component_status"],
                size_systematic_status=bhs_row["size_systematic_status"],
                size_total_status=bhs_row["size_total_status"],
            )
        )

    data = PlotData(tuple(points), tuple(skipped_without_size))
    if enforce_publication_sample:
        missing_formal = {
            point.source_key
            for point in data.points
            if point.uncertainty_basis == "missing"
        }
        if missing_formal != EXPECTED_MISSING_FORMAL_SOURCES:
            raise ValueError(
                "publication sample must retain explicit missing formal uncertainties for "
                f"{sorted(EXPECTED_MISSING_FORMAL_SOURCES)}; found {sorted(missing_formal)}"
            )
    return data


def _reference_point(
    *,
    label: str,
    size_uas: float,
    size_unc_minus_uas: float,
    size_unc_plus_uas: float,
    flux_mjy: float,
    flux_unc_mjy: float,
) -> ReferencePoint:
    log_minus, log_plus = exact_log10_uncertainties(
        size_uas,
        size_unc_minus_uas,
        size_unc_plus_uas,
    )
    flux_log_minus, flux_log_plus = exact_log10_uncertainties(
        flux_mjy,
        flux_unc_mjy,
        flux_unc_mjy,
    )
    return ReferencePoint(
        label=label,
        size_uas=size_uas,
        flux_mjy=flux_mjy,
        size_log_unc_minus=log_minus,
        size_log_unc_plus=log_plus,
        flux_log_unc_minus=flux_log_minus,
        flux_log_unc_plus=flux_log_plus,
    )


def reference_points() -> tuple[ReferencePoint, ReferencePoint]:
    """Return the direct EHT comparison sizes printed in manuscript Table 1."""

    m87 = _reference_point(
        label="M87*",
        size_uas=42.0,
        size_unc_minus_uas=3.0,
        size_unc_plus_uas=3.0,
        flux_mjy=1.571 * 1_000.0,
        flux_unc_mjy=0.261 * 1_000.0,
    )
    sgr_a = _reference_point(
        label="Sgr A*",
        size_uas=51.8,
        size_unc_minus_uas=2.3,
        size_unc_plus_uas=2.3,
        flux_mjy=3.22 * 1_000.0,
        flux_unc_mjy=0.62 * 1_000.0,
    )
    return m87, sgr_a


def baseline_resolutions_uas() -> tuple[tuple[str, float], ...]:
    """Return the 230 GHz nominal angular resolution of each prior baseline."""

    earth_diameter_km = 12_756.0
    earth_radius_km = earth_diameter_km / 2.0
    baselines_km = (
        ("Moon", 382_500.0 + earth_radius_km),
        ("Earth", earth_diameter_km),
        ("L2", 1.5e6 + earth_radius_km),
        # BHEX's 20,200 km value is orbital altitude.  The manuscript's
        # illustrative marker uses one Earth radius plus that altitude as a
        # characteristic projected space-ground baseline (~20 Glambda).
        ("BHEX", 20_200.0 + earth_radius_km),
        ("L4/L5", 150.0e6 + earth_radius_km),
    )
    wavelength_m = SPEED_OF_LIGHT / (230.0e9)
    return tuple(
        (
            label,
            wavelength_m / (baseline_km * 1_000.0) * RAD_TO_UAS,
        )
        for label, baseline_km in baselines_km
    )


def _round_down_half(value: float) -> float:
    return math.floor(value * 2.0) / 2.0


def _round_up_half(value: float) -> float:
    return math.ceil(value * 2.0) / 2.0


def _asymmetric_errorbar(lower: float, upper: float) -> np.ndarray:
    return np.asarray([[lower], [upper]], dtype=float)


def _plot_source_point(ax: plt.Axes, point: PlotPoint) -> None:
    color = "blue" if point.flux_status == "detection" else "red"
    marker = "o" if point.flux_status == "detection" else "v"
    ax.errorbar(
        point.log_size,
        point.log_flux,
        color=color,
        markersize=4.5,
        marker=marker,
        yerr=None,
        linestyle="none",
        elinewidth=0.9,
        capsize=2,
        capthick=0.9,
        zorder=3,
    )

    if point.display_name in LABELLED_SOURCES:
        dx, dy = LABEL_OFFSETS.get(point.display_name, (0.0, 0.03))
        ax.text(
            point.log_size + dx,
            point.log_flux + dy,
            point.display_name,
            fontsize=9,
            ha="center",
        )


def _marker_legend(ax: plt.Axes) -> None:
    handles = [
        Line2D(
            [0],
            [0],
            color="blue",
            marker="o",
            markersize=4.5,
            linestyle="none",
            label="Detection",
        ),
        Line2D(
            [0],
            [0],
            color="red",
            marker="v",
            markersize=4.5,
            linestyle="none",
            label="Upper limit",
        ),
        Line2D(
            [0],
            [0],
            color="green",
            marker="o",
            markersize=4.5,
            linestyle="none",
            label="M87* and Sgr A*",
        ),
    ]
    ax.legend(
        handles=handles,
        loc="lower left",
        frameon=False,
        fontsize=8,
        ncol=3,
        columnspacing=1.4,
        handletextpad=0.4,
    )


def _plot_reference(ax: plt.Axes, reference: ReferencePoint) -> None:
    log_size = math.log10(reference.size_uas)
    log_flux = math.log10(reference.flux_mjy)
    ax.errorbar(
        log_size,
        log_flux,
        color="green",
        markersize=4.5,
        marker="o",
        yerr=None,
        linestyle="none",
        elinewidth=0.9,
        capsize=2,
        capthick=0.9,
        zorder=3,
    )
    dx = 0.21 if reference.label == "M87*" else 0.31
    ax.text(log_size + dx, log_flux - 0.035, reference.label, ha="center")


def make_figure(data: PlotData) -> plt.Figure:
    if not data.points:
        raise ValueError("no sources have both a shadow size and 230 GHz flux")

    references = reference_points()
    baselines = baseline_resolutions_uas()

    with matplotlib.rc_context(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans"],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    ):
        figure, ax = plt.subplots(figsize=(10, 6), dpi=300)

        for point in data.points:
            _plot_source_point(ax, point)
        for reference in references:
            _plot_reference(ax, reference)

        min_flux_value = min(point.log_flux for point in data.points)
        max_flux_value = max(
            [point.log_flux for point in data.points]
            + [math.log10(reference.flux_mjy) for reference in references]
        )
        min_size_value = min(point.log_size for point in data.points)
        max_size_value = max(
            [point.log_size for point in data.points]
            + [math.log10(reference.size_uas) for reference in references]
        )

        min_flux_tick = _round_down_half(min_flux_value)
        max_flux_tick = _round_up_half(max_flux_value)
        min_size_tick = _round_down_half(min_size_value)
        max_size_tick = _round_up_half(max_size_value)
        baseline_logs = [math.log10(resolution) for _, resolution in baselines]
        min_size_tick = min(min_size_tick, _round_down_half(min(baseline_logs)))
        max_size_tick = max(max_size_tick, _round_up_half(max(baseline_logs)))

        flux_padding = (max_flux_tick - min_flux_tick) * 0.05
        size_padding = (max_size_tick - min_size_tick) * 0.05
        ax.set_ylim(min_flux_tick - flux_padding, max_flux_tick + flux_padding)
        ax.set_xlim(min_size_tick - size_padding, max_size_tick + size_padding)

        x_ticks = np.arange(min_size_tick, max_size_tick + 0.5, 0.5)
        y_ticks = np.arange(min_flux_tick, max_flux_tick + 0.5, 0.5)
        ax.set_xticks(x_ticks, [f"{value:.1f}" for value in x_ticks])
        ax.set_yticks(y_ticks, [f"{value:.1f}" for value in y_ticks])
        ax.xaxis.set_minor_locator(MultipleLocator(0.1))
        ax.yaxis.set_minor_locator(MultipleLocator(0.1))

        for label, resolution in baselines:
            log_resolution = math.log10(resolution)
            ax.axvline(
                log_resolution,
                color="0.35",
                linestyle="-",
                linewidth=0.8,
                zorder=1,
            )
            label_x = log_resolution - 0.1
            label_y = max_flux_tick
            horizontal_alignment = "center"
            if label == "L4/L5":
                label_x = log_resolution + 0.04
                label_y = max_flux_tick - 0.08
                horizontal_alignment = "left"
            ax.text(
                label_x,
                label_y,
                label,
                va="center",
                ha=horizontal_alignment,
                color="0.25",
                rotation=90,
            )

        for threshold in FLUX_THRESHOLDS_MJY:
            ax.axhline(
                math.log10(threshold),
                color="0.82",
                linestyle="-",
                linewidth=0.6,
                zorder=0,
            )

        _marker_legend(ax)
        ax.set_xlabel(r"$\log_{10}(\theta_{\rm BHS}/\mu{\rm as})$")
        ax.set_ylabel(r"$\log_{10}(S_{230}/{\rm mJy})$")
        ax.grid(False)
        return figure


def build_figure(
    bhs_path: Path = DEFAULT_BHS_PATH,
    flux_path: Path = DEFAULT_FLUX_PATH,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> FigureBuildResult:
    """Build deterministic PDF and PNG figure products."""

    data = load_plot_data(Path(bhs_path), Path(flux_path))
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = output_dir / "BHS_Flux_Graph.pdf"
    png_path = output_dir / "BHS_Flux_Graph.png"
    figure = make_figure(data)

    fixed_time = datetime(2000, 1, 1, tzinfo=timezone.utc)
    figure.savefig(
        pdf_path,
        bbox_inches="tight",
        dpi=300,
        metadata={
            "Title": "Predicted BHS angular size versus nuclear-region S230",
            "Creator": "build_bhs_flux_figure.py",
            "Producer": "Matplotlib",
            "CreationDate": fixed_time,
            "ModDate": fixed_time,
        },
    )
    figure.savefig(
        png_path,
        bbox_inches="tight",
        dpi=300,
        metadata={"Software": "build_bhs_flux_figure.py"},
    )
    plt.close(figure)
    return FigureBuildResult(pdf_path=pdf_path, png_path=png_path, data=data)


def _status_summary(points: Iterable[PlotPoint]) -> str:
    counts = {
        basis: 0
        for basis in ("total", "component_only", "model_span", "missing")
    }
    for point in points:
        counts[point.uncertainty_basis] += 1
    return ", ".join(f"{basis}={count}" for basis, count in counts.items())


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bhs-params", type=Path, default=DEFAULT_BHS_PATH)
    parser.add_argument("--flux-table", type=Path, default=DEFAULT_FLUX_PATH)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    result = build_figure(args.bhs_params, args.flux_table, args.output_dir)
    print(f"Wrote {result.pdf_path}")
    print(f"Wrote {result.png_path}")
    print(f"Plotted {len(result.data.points)} sources ({_status_summary(result.data.points)})")
    if result.data.skipped_without_size:
        print(
            "Skipped sources without a robust shadow size: "
            + ", ".join(result.data.skipped_without_size)
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
