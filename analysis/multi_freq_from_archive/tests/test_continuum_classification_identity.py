import ast
import csv
import hashlib
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


PACKAGE_DIR = Path(__file__).resolve().parents[1]
if str(PACKAGE_DIR) not in sys.path:
    sys.path.insert(0, str(PACKAGE_DIR))

import continuum_classification_identity as identity  # noqa: E402
from continuum_noise_inventory import GENERIC_FINAL_NOISE_AUDIT_FIELDS  # noqa: E402
from continuum_product_additions import (  # noqa: E402
    ADDITION_POLICY_FIELDS,
    ADDITION_POLICY_VERSION,
    ADMISSION_DECISION,
)
from continuum_boundary_rejections import (  # noqa: E402
    BOUNDARY_REJECTION_EVIDENCE_FIELDS,
    BOUNDARY_REJECTION_EVIDENCE_NAME,
    BOUNDARY_REJECTION_EVIDENCE_VERSION,
    BOUNDARY_REJECTION_POLICY_FIELDS,
    BOUNDARY_REJECTION_POLICY_VERSION,
    EXPECTED_CLASSIFICATION,
    EXPECTED_ROW_DISPOSITION,
    REJECTION_DECISION,
    REJECTION_GATE,
)


class ContinuumClassificationIdentityTests(unittest.TestCase):
    BOUNDARY_SPECS = (
        ("BoundaryA", "VLA/BoundaryA/a.fits", "1.50", "1985/01/01", "91.0"),
        ("BoundaryB", "VLA/BoundaryB/b.fits", "5.00", "1986/02/02", "95.0"),
        ("BoundaryC", "VLA/BoundaryC/c.fits", "15.00", "1987/03/03", "93.0"),
    )

    @staticmethod
    def write_rows(path: Path, rows: list[list[object]]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("".join(f"{row!r}\n" for row in rows), encoding="utf-8")

    @staticmethod
    def write_csv(
        path: Path,
        fields: tuple[str, ...],
        rows: list[dict[str, object]],
    ) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    @staticmethod
    def write_test_fits(path: Path, ra_deg: float, dec_deg: float) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(f"synthetic FITS fixture {ra_deg} {dec_deg}\n".encode())

    @staticmethod
    def write_header_fits(path: Path, date_obs: str, frequency_hz: float) -> None:
        def card(keyword: str, value: str) -> bytes:
            return f"{keyword:<8}= {value}".ljust(80).encode("ascii")

        cards = [
            card("SIMPLE", "                    T"),
            card("BITPIX", "                    8"),
            card("NAXIS", "                    0"),
            card("DATE-OBS", f"'{date_obs}'"),
            card("CTYPE3", "'FREQ'"),
            card("CUNIT3", "'Hz'"),
            card("CRVAL3", f"{frequency_hz:20.12E}"),
            b"END".ljust(80),
        ]
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"".join(cards).ljust(2880, b" "))

    @staticmethod
    def sha256(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    def make_fixture(self, root: Path) -> dict[str, Path]:
        workspace = root / "release"
        analysis = workspace / "analysis/multi_freq_from_archive"
        products = (
            "ALMA/Test/product_000.fits",
            "VLA/Test/product_001.fits",
        )
        centers = ((150.0, 2.0), (151.0, 3.0))
        for product, (ra_deg, dec_deg) in zip(products, centers):
            self.write_test_fits(analysis / product, ra_deg, dec_deg)

        canonical = root / "canonical.txt"
        generated = root / "generated.txt"
        exclusions = root / "exclusions.csv"
        audit = root / "audit.csv"
        policy = analysis / "continuum_classification_identity_policy.csv"
        canonical_rows = [
            ["TestA", "ALMA", "230.00", "0.30*", "2020/01/01_1", "1.0", "NA", products[0]],
            ["TestB", "VLA", "5.00", "0.20", "2020/01/02_1", "1.0", "4.000", products[1]],
        ]
        generated_rows = [
            ["TestA", "ALMA", "230.00", "0.20", "2020/01/01_1", "1.0", "4.000", products[0]],
            ["TestB", "VLA", "5.00", "0.20", "2020/01/02_1", "1.0", "4.000", products[1]],
        ]
        self.write_rows(canonical, canonical_rows)
        self.write_rows(generated, generated_rows)
        self.write_csv(exclusions, ("product",), [])

        audit_rows = []
        for index, product in enumerate(products):
            row = {field: "NA" for field in GENERIC_FINAL_NOISE_AUDIT_FIELDS}
            row.update(
                {
                    "source": f"Test{'AB'[index]}",
                    "telescope": ("ALMA", "VLA")[index],
                    "frequency_ghz": ("230.00", "5.00")[index],
                    "product": product,
                    "forced_policy": "no",
                    "aperture_center_x_pixel": "10.0",
                    "aperture_center_y_pixel": "10.0",
                    "fit_center_reference_frame": "ICRS",
                    "fit_center_ra_deg": f"{centers[index][0]:.8f}",
                    "fit_center_dec_deg": f"{centers[index][1]:.8f}",
                    "beam_major_arcsec": "1.0",
                    "flux_mjy": "0.20",
                    "rms_mjy_per_beam": ("0.04", "0.10")[index],
                    "scaled_mad_mjy_per_beam": "0.05",
                    "classification": "detection",
                }
            )
            audit_rows.append(row)
        self.write_csv(audit, GENERIC_FINAL_NOISE_AUDIT_FIELDS, audit_rows)

        policy_rows = []
        for index, (product, (ra_deg, dec_deg)) in enumerate(zip(products, centers)):
            policy_rows.append(
                {
                    "schema_version": "1",
                    "policy_version": identity.POLICY_VERSION,
                    "review_id": f"identity-review-{index + 1:02d}",
                    "trigger": (
                        "new_detection"
                        if index == 0
                        else "rms_mad_threshold_crossing"
                    ),
                    "source": f"Test{'AB'[index]}",
                    "telescope": ("ALMA", "VLA")[index],
                    "frequency_ghz": ("230.00", "5.00")[index],
                    "product": product,
                    "product_sha256": self.sha256(analysis / product),
                    "reference_role": "published_nuclear_position",
                    "reference_frame": "ICRS",
                    "reference_ra_deg": f"{ra_deg:.8f}",
                    "reference_dec_deg": f"{dec_deg:.8f}",
                    "reference_evidence": "doi:10.0000/example; table 1",
                    "identity_evidence_class": "a_priori_position_match",
                    "component_identity_status": "confirmed_target_component",
                    "component_label": "nuclear_component",
                    "measurement_semantics": "compact_component_flux",
                    "background_morphology": "noise_like",
                    "expected_scaled_mad_classification": "detection",
                    "review_decision": "accept_scaled_mad_classification",
                    "reviewer": "continuum_identity_audit",
                    "review_date": "2026-08-19",
                    "rationale": "Fresh fitted centroid matches the reviewed nuclear position.",
                }
            )
        self.write_csv(policy, identity.POLICY_FIELDS, policy_rows)
        return {
            "workspace": workspace,
            "canonical": canonical,
            "generated": generated,
            "exclusions": exclusions,
            "audit": audit,
            "policy": policy,
        }

    def add_boundary_contract(self, paths: dict[str, Path]) -> Path:
        analysis = paths["workspace"] / "analysis/multi_freq_from_archive"
        evidence_path = analysis / BOUNDARY_REJECTION_EVIDENCE_NAME
        evidence_rows: list[dict[str, str]] = []
        for index, (_, product, _, _, percent) in enumerate(
            self.BOUNDARY_SPECS, start=1
        ):
            model_hash = f"{index}" * 64
            evidence_rows.append(
                {
                    "schema_version": "1",
                    "evidence_version": BOUNDARY_REJECTION_EVIDENCE_VERSION,
                    "rejection_id": f"CBR-{index:03d}",
                    "product": product,
                    "independent_source_file": "synthetic/excluded_results.csv",
                    "independent_source_sha256": "a" * 64,
                    "independent_source_record_index": str(index),
                    "fit_converged": "true",
                    "boundary_result": "true",
                    "rerun_outcome": "on_boundary",
                    "historical_boundary_diagnostic": "true",
                    "boundary_highest_percent": percent,
                    "archived_final_model_exists": "true",
                    "archived_final_model_flux_mjy": str(index),
                    "rerun_model_flux_mjy": str(index),
                    "model_flux_delta_mjy": "0.0",
                    "retry3_model_payload_sha256": model_hash,
                    "authoritative7_model_payload_sha256": model_hash,
                    "model_payload_size_bytes": str(100 + index),
                    "input_metadata_unchanged": "true",
                    "status": "ok",
                }
            )
        self.write_csv(
            evidence_path,
            BOUNDARY_REJECTION_EVIDENCE_FIELDS,
            evidence_rows,
        )
        evidence_sha = self.sha256(evidence_path)

        with paths["audit"].open(encoding="utf-8", newline="") as handle:
            audit = list(csv.DictReader(handle))
        policy_rows: list[dict[str, str]] = []
        for index, (source, product, frequency, date, percent) in enumerate(
            self.BOUNDARY_SPECS, start=1
        ):
            product_path = analysis / product
            self.write_header_fits(
                product_path,
                date.replace("/", "-") + "T00:00:00",
                float(frequency) * 1e9,
            )
            model_hash = f"{index}" * 64
            policy_rows.append(
                {
                    "schema_version": "1",
                    "policy_version": BOUNDARY_REJECTION_POLICY_VERSION,
                    "rejection_id": f"CBR-{index:03d}",
                    "source": source,
                    "telescope": "VLA",
                    "expected_frequency_ghz": frequency,
                    "expected_observation_date": date,
                    "product": product,
                    "product_sha256": self.sha256(product_path),
                    "size_bytes": str(product_path.stat().st_size),
                    "independent_evidence_record": (
                        f"{BOUNDARY_REJECTION_EVIDENCE_NAME}#CBR-{index:03d}"
                    ),
                    "independent_evidence_sha256": evidence_sha,
                    "historical_boundary_highest_percent": percent,
                    "historical_model_payload_sha256": model_hash,
                    "expected_row_disposition": EXPECTED_ROW_DISPOSITION,
                    "expected_classification": EXPECTED_CLASSIFICATION,
                    "rejection_gate": REJECTION_GATE,
                    "rejection_decision": REJECTION_DECISION,
                    "reviewer": "identity-test",
                    "review_date": "2026-08-20",
                    "historical_evidence": "synthetic exact boundary reproduction",
                    "rationale": "Synthetic reviewed boundary rejection.",
                }
            )
            record = {field: "NA" for field in GENERIC_FINAL_NOISE_AUDIT_FIELDS}
            record.update(
                {
                    "source": source,
                    "telescope": "VLA",
                    "frequency_ghz": frequency,
                    "product": product,
                    "forced_policy": "no",
                    "classification": EXPECTED_CLASSIFICATION,
                    "row_disposition": EXPECTED_ROW_DISPOSITION,
                }
            )
            audit.append(record)
        audit.sort(key=lambda row: row["product"])
        self.write_csv(paths["audit"], GENERIC_FINAL_NOISE_AUDIT_FIELDS, audit)
        policy_path = analysis / "continuum_boundary_rejection_policy.csv"
        self.write_csv(
            policy_path,
            BOUNDARY_REJECTION_POLICY_FIELDS,
            policy_rows,
        )
        return policy_path

    def validate(
        self,
        paths: dict[str, Path],
        boundary_policy: Path | None = None,
    ) -> dict[str, object]:
        with patch.multiple(
            identity,
            EXPECTED_CANONICAL_ROWS=2,
            EXPECTED_GENERATED_ROWS=2,
            EXPECTED_EXCLUSIONS=0,
        ):
            return identity.validate_classification_identity_policy(
                workspace=paths["workspace"],
                canonical_table=paths["canonical"],
                generated_table=paths["generated"],
                exclusion_policy=paths["exclusions"],
                generic_noise_audit=paths["audit"],
                identity_policy=paths["policy"],
                boundary_rejection_policy=boundary_policy,
            )

    def test_exact_fresh_required_inventory_and_centroids_pass(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            paths = self.make_fixture(Path(temporary))
            result = self.validate(paths)
        self.assertEqual(result["status"], "exact_required_identity_reviews_resolved")
        self.assertEqual(result["required_review_rows"], 2)
        self.assertEqual(result["new_detection_rows"], 1)
        self.assertEqual(result["rms_mad_threshold_crossing_rows"], 1)
        self.assertEqual(len(result["validated_reviews"]), 2)
        for record in result["validated_reviews"]:
            self.assertLess(float(record["separation_bmaj"]), 1e-6)

    def test_exact_boundary_diagnostics_do_not_create_identity_reviews(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            paths = self.make_fixture(Path(temporary))
            boundary_policy = self.add_boundary_contract(paths)
            result = self.validate(paths, boundary_policy)
        self.assertEqual(result["required_review_rows"], 2)
        self.assertEqual(len(result["validated_reviews"]), 2)

    def test_unreviewed_or_stale_boundary_audit_rows_fail_closed(self) -> None:
        cases = (
            "unreviewed extra",
            "missing reviewed row",
            "wrong disposition",
            "wrong classification",
        )
        for mutation in cases:
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as temporary:
                paths = self.make_fixture(Path(temporary))
                boundary_policy = self.add_boundary_contract(paths)
                with paths["audit"].open(encoding="utf-8", newline="") as handle:
                    audit = list(csv.DictReader(handle))
                if mutation == "unreviewed extra":
                    extra = dict(audit[0])
                    extra["product"] = "VLA/Unreviewed/extra.fits"
                    audit.append(extra)
                elif mutation == "missing reviewed row":
                    audit = [
                        row
                        for row in audit
                        if row["product"] != self.BOUNDARY_SPECS[0][1]
                    ]
                else:
                    boundary_row = next(
                        row
                        for row in audit
                        if row["product"] == self.BOUNDARY_SPECS[0][1]
                    )
                    if mutation == "wrong disposition":
                        boundary_row["row_disposition"] = "generic_table_retained"
                    else:
                        boundary_row["classification"] = "ordinary_upper_limit"
                audit.sort(key=lambda row: row["product"])
                self.write_csv(
                    paths["audit"], GENERIC_FINAL_NOISE_AUDIT_FIELDS, audit
                )
                with self.assertRaises(RuntimeError):
                    self.validate(paths, boundary_policy)

    def test_policy_mutations_fail_closed(self) -> None:
        mutations = {
            "missing row": lambda rows: rows.pop(),
            "bad hash": lambda rows: rows[0].__setitem__("product_sha256", "0" * 64),
            "stale class": lambda rows: rows[0].__setitem__(
                "expected_scaled_mad_classification", "ordinary_upper_limit"
            ),
            "far reference": lambda rows: rows[0].__setitem__(
                "reference_ra_deg", "150.01000000"
            ),
            "duplicate review id": lambda rows: rows[1].__setitem__(
                "review_id", rows[0]["review_id"]
            ),
        }
        for label, mutate in mutations.items():
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                paths = self.make_fixture(Path(temporary))
                with paths["policy"].open(encoding="utf-8", newline="") as handle:
                    rows = list(csv.DictReader(handle))
                mutate(rows)
                self.write_csv(paths["policy"], identity.POLICY_FIELDS, rows)
                with self.assertRaises(RuntimeError):
                    self.validate(paths)

    def test_audit_identity_and_native_frame_mutations_fail_closed(self) -> None:
        mutations = {
            "stale audit source": lambda rows: rows[0].__setitem__("source", "Wrong"),
            "invalid native frame": lambda rows: rows[0].__setitem__(
                "fit_center_reference_frame", "FK5"
            ),
            "invalid native RA": lambda rows: rows[0].__setitem__(
                "fit_center_ra_deg", "360"
            ),
            "stale classification": lambda rows: rows[0].__setitem__(
                "classification", "ordinary_upper_limit"
            ),
        }
        for label, mutate in mutations.items():
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                paths = self.make_fixture(Path(temporary))
                with paths["audit"].open(encoding="utf-8", newline="") as handle:
                    rows = list(csv.DictReader(handle))
                mutate(rows)
                self.write_csv(paths["audit"], GENERIC_FINAL_NOISE_AUDIT_FIELDS, rows)
                with self.assertRaises(RuntimeError):
                    self.validate(paths)

        with tempfile.TemporaryDirectory() as temporary:
            paths = self.make_fixture(Path(temporary))
            with paths["policy"].open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            rows[0]["reference_frame"] = "J2000"
            self.write_csv(paths["policy"], identity.POLICY_FIELDS, rows)
            with self.assertRaisesRegex(RuntimeError, "frame differs"):
                self.validate(paths)

    def test_lost_detection_without_explicit_review_schema_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            paths = self.make_fixture(Path(temporary))
            with paths["generated"].open(encoding="utf-8") as handle:
                rows = [ast.literal_eval(line) for line in handle if line.strip()]
            rows[1][3] = "0.15*"
            rows[1][6] = "NA"
            self.write_rows(paths["generated"], rows)
            with paths["audit"].open(encoding="utf-8", newline="") as handle:
                audit = list(csv.DictReader(handle))
            audit[1]["flux_mjy"] = "0.10"
            audit[1]["classification"] = "ordinary_upper_limit"
            self.write_csv(paths["audit"], GENERIC_FINAL_NOISE_AUDIT_FIELDS, audit)
            with self.assertRaisesRegex(RuntimeError, "lost-detection"):
                self.validate(paths)

    def test_detected_addition_requires_identity_but_added_upper_limit_does_not(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            paths = self.make_fixture(Path(temporary))
            analysis = paths["workspace"] / "analysis/multi_freq_from_archive"
            addition_specs = [
                ("ZDetect", "ALMA/ZDetect/detect.fits", "500.00", "2021/02/03_1", True),
                ("ZLimit", "ALMA/ZLimit/limit.fits", "510.00", "2022/03/04_1", False),
            ]
            addition_policy_rows = []
            with paths["generated"].open(encoding="utf-8") as handle:
                generated = [ast.literal_eval(line) for line in handle if line.strip()]
            with paths["audit"].open(encoding="utf-8", newline="") as handle:
                audit = list(csv.DictReader(handle))
            with paths["policy"].open(encoding="utf-8", newline="") as handle:
                policy = list(csv.DictReader(handle))
            for index, (source, product, frequency, date_label, detected) in enumerate(addition_specs):
                self.write_header_fits(
                    analysis / product,
                    date_label.split("_")[0].replace("/", "-") + "T00:00:00",
                    float(frequency) * 1e9,
                )
                addition_policy_rows.append(
                    {
                        "schema_version": "1",
                        "policy_version": ADDITION_POLICY_VERSION,
                        "admission_id": f"ADD-{index + 1:03d}",
                        "source": source,
                        "telescope": "ALMA",
                        "expected_frequency_ghz": frequency,
                        "expected_date_label": date_label,
                        "product": product,
                        "product_sha256": self.sha256(analysis / product),
                        "expected_run_disposition": "publication_row_required",
                        "admission_decision": ADMISSION_DECISION,
                        "reviewer": "identity-test",
                        "review_date": "2026-08-20",
                        "rationale": "Synthetic admitted product.",
                    }
                )
                generated.append(
                    [
                        source,
                        "ALMA",
                        frequency,
                        "0.20" if detected else "0.15*",
                        date_label,
                        "1.0",
                        "4.000" if detected else "NA",
                        product,
                    ]
                )
                record = {field: "NA" for field in GENERIC_FINAL_NOISE_AUDIT_FIELDS}
                record.update(
                    {
                        "source": source,
                        "telescope": "ALMA",
                        "frequency_ghz": frequency,
                        "product": product,
                        "forced_policy": "no",
                        "fit_center_reference_frame": "ICRS",
                        "fit_center_ra_deg": "152.0",
                        "fit_center_dec_deg": "4.0",
                        "beam_major_arcsec": "1.0",
                        "flux_mjy": "0.20" if detected else "0.10",
                        "rms_mjy_per_beam": "0.10",
                        "scaled_mad_mjy_per_beam": "0.05",
                        "classification": "detection" if detected else "ordinary_upper_limit",
                    }
                )
                audit.append(record)
                if detected:
                    reviewed = dict(policy[0])
                    reviewed.update(
                        {
                            "review_id": "identity-review-addition",
                            "trigger": "new_product_admission",
                            "source": source,
                            "telescope": "ALMA",
                            "frequency_ghz": frequency,
                            "product": product,
                            "product_sha256": self.sha256(analysis / product),
                            "reference_ra_deg": "152.0",
                            "reference_dec_deg": "4.0",
                        }
                    )
                    policy.append(reviewed)
            write_additions = analysis / "continuum_product_addition_policy.csv"
            self.write_csv(write_additions, ADDITION_POLICY_FIELDS, addition_policy_rows)
            self.write_rows(paths["generated"], generated)
            self.write_csv(paths["audit"], GENERIC_FINAL_NOISE_AUDIT_FIELDS, audit)
            policy.sort(key=lambda row: row["product"])
            self.write_csv(paths["policy"], identity.POLICY_FIELDS, policy)

            with patch.multiple(
                identity,
                EXPECTED_CANONICAL_ROWS=2,
                EXPECTED_GENERATED_ROWS=2,
                EXPECTED_EXCLUSIONS=0,
            ):
                result = identity.validate_classification_identity_policy(
                    workspace=paths["workspace"],
                    canonical_table=paths["canonical"],
                    generated_table=paths["generated"],
                    exclusion_policy=paths["exclusions"],
                    generic_noise_audit=paths["audit"],
                    identity_policy=paths["policy"],
                    addition_policy=write_additions,
                )
            self.assertEqual(result["new_product_admission_rows"], 1)
            self.assertEqual(result["admitted_product_upper_limit_rows"], 1)

    def test_ngc4945_identity_records_use_the_fixed_source18_photometry(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "ngc4945_nuclear_photometry.csv"
            row = {field: "NA" for field in identity.NGC4945_PHOTOMETRY_FIELDS}
            row.update(
                {
                    "source": "NGC4945",
                    "telescope": "ALMA",
                    "frequency_ghz": "225.07",
                    "product": "ALMA/NGC4945/source18.fits",
                    "coordinate_frame": "ICRS",
                    "component_identity": "Emig2020_source18_AGN_maser_nucleus",
                    "centroid_policy": "fixed_xy",
                    "fixed_ra": "13:05:27.454",
                    "fixed_dec": "-49:28:04.944",
                    "beam_major_arcsec": "0.50",
                    "beam_minor_arcsec": "0.40",
                    "beam_pa_deg": "12.0",
                    "noise_aperture_half_width_x_pixels": "5",
                    "noise_aperture_half_width_y_pixels": "6",
                    "pixel_width_arcsec": "0.02",
                    "pixel_height_arcsec": "0.02",
                    "unrounded_flux_mjy": "62.54",
                    "new_snr": "45.42",
                    "row_disposition": "retained_source18_replacement",
                    "noise_guard_status": "passed",
                    "status": "detection",
                }
            )
            self.write_csv(path, identity.NGC4945_PHOTOMETRY_FIELDS, [row])

            with patch.object(identity, "EXPECTED_NGC4945_ROWS", 1):
                records = identity.load_ngc4945_source18_identity_records(path)

            record = records["ALMA/NGC4945/source18.fits"]
            self.assertEqual(record["classification"], "detection")
            self.assertEqual(
                record["measurement_scope"],
                "publication_fixed_source18_source_excluded",
            )
            self.assertAlmostEqual(float(record["fit_center_ra_deg"]), 196.3643917)
            self.assertAlmostEqual(float(record["fit_center_dec_deg"]), -49.46804)

    def test_archive_directory_source_normalizes_to_fitter_output_source(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            paths = self.make_fixture(Path(temporary))
            analysis = paths["workspace"] / "analysis/multi_freq_from_archive"
            product = "ALMA/CircinusGalaxy/addition.fits"
            self.write_header_fits(
                analysis / product,
                "2021-04-02T00:00:00",
                349.63e9,
            )
            addition_rows = [
                {
                    "schema_version": "1",
                    "policy_version": ADDITION_POLICY_VERSION,
                    "admission_id": "ADD-CIRCINUS",
                    "source": "CircinusGalaxy",
                    "telescope": "ALMA",
                    "expected_frequency_ghz": "349.63",
                    "expected_date_label": "2021/04/02_1",
                    "product": product,
                    "product_sha256": self.sha256(analysis / product),
                    "expected_run_disposition": "publication_row_required",
                    "admission_decision": ADMISSION_DECISION,
                    "reviewer": "identity-test",
                    "review_date": "2026-08-20",
                    "rationale": "Synthetic Circinus directory-label regression.",
                }
            ]
            additions = analysis / "continuum_product_addition_policy.csv"
            self.write_csv(additions, ADDITION_POLICY_FIELDS, addition_rows)

            with paths["generated"].open(encoding="utf-8") as handle:
                generated = [ast.literal_eval(line) for line in handle if line.strip()]
            generated.append(
                [
                    "Circinus",
                    "ALMA",
                    "349.63",
                    "1.00",
                    "2021/04/02_1",
                    "0.20",
                    "10.000",
                    product,
                ]
            )
            self.write_rows(paths["generated"], generated)

            with paths["audit"].open(encoding="utf-8", newline="") as handle:
                audit = list(csv.DictReader(handle))
            record = {field: "NA" for field in GENERIC_FINAL_NOISE_AUDIT_FIELDS}
            record.update(
                {
                    "source": "Circinus",
                    "telescope": "ALMA",
                    "frequency_ghz": "349.63",
                    "product": product,
                    "forced_policy": "no",
                    "fit_center_reference_frame": "ICRS",
                    "fit_center_ra_deg": "213.29145",
                    "fit_center_dec_deg": "-65.33922",
                    "beam_major_arcsec": "0.22",
                    "flux_mjy": "1.00",
                    "rms_mjy_per_beam": "0.20",
                    "scaled_mad_mjy_per_beam": "0.10",
                    "classification": "detection",
                }
            )
            audit.append(record)
            self.write_csv(paths["audit"], GENERIC_FINAL_NOISE_AUDIT_FIELDS, audit)

            with paths["policy"].open(encoding="utf-8", newline="") as handle:
                policy = list(csv.DictReader(handle))
            reviewed = dict(policy[0])
            reviewed.update(
                {
                    "review_id": "identity-review-circinus-addition",
                    "trigger": "new_product_admission",
                    "source": "Circinus",
                    "telescope": "ALMA",
                    "frequency_ghz": "349.63",
                    "product": product,
                    "product_sha256": self.sha256(analysis / product),
                    "reference_ra_deg": "213.29145",
                    "reference_dec_deg": "-65.33922",
                }
            )
            policy.append(reviewed)
            policy.sort(key=lambda row: row["product"])
            self.write_csv(paths["policy"], identity.POLICY_FIELDS, policy)

            with patch.multiple(
                identity,
                EXPECTED_CANONICAL_ROWS=2,
                EXPECTED_GENERATED_ROWS=2,
                EXPECTED_EXCLUSIONS=0,
            ):
                result = identity.validate_classification_identity_policy(
                    workspace=paths["workspace"],
                    canonical_table=paths["canonical"],
                    generated_table=paths["generated"],
                    exclusion_policy=paths["exclusions"],
                    generic_noise_audit=paths["audit"],
                    identity_policy=paths["policy"],
                    addition_policy=additions,
                )
            self.assertEqual(result["new_product_admission_rows"], 1)


if __name__ == "__main__":
    unittest.main()
