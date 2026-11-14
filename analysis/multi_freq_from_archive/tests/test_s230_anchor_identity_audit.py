from __future__ import annotations

import csv
import unittest
from pathlib import Path


ANALYSIS_DIR = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = ANALYSIS_DIR.parents[1]
AUDIT_PATH = ANALYSIS_DIR / "s230_anchor_identity_audit.csv"
S230_PATH = ANALYSIS_DIR / "machinetables" / "S230_extrapolation_comparison.csv"

EXPECTED_COLUMNS = [
    "source",
    "selected_flux_status",
    "selected_observation_id",
    "selected_reference",
    "identity_disposition",
    "component_identity",
    "evidence_basis",
    "canonical_component_reference",
    "rejected_component_reference",
    "required_handling",
]

EXPECTED_DISPOSITIONS = {
    "Circinus": "nuclear",
    "ESO558-G009": "nuclear",
    "IC2560": "ambiguous",
    "MRK1029": "upper_limit",
    "NGC1068": "nuclear",
    "NGC1194": "nuclear",
    "NGC1320": "ambiguous",
    "NGC2273": "conservative_upper_limit",
    "NGC2960": "upper_limit",
    "NGC3079": "blended",
    "NGC3393": "nuclear",
    "NGC4258": "blended",
    "NGC4388": "nuclear",
    "NGC4945": "nuclear",
    "NGC5765B": "nuclear",
    "NGC6264": "upper_limit",
    "NGC6323": "upper_limit",
    "UGC3789": "upper_limit",
    "UGC6093": "upper_limit",
}


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def resolve_record_reference(reference: str) -> tuple[Path, str]:
    relative_path, separator, record_id = reference.partition("#")
    if not separator or not relative_path or not record_id:
        raise ValueError(f"Invalid record reference: {reference!r}")
    return REPOSITORY_ROOT / relative_path, record_id


class S230AnchorIdentityAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.audit_fields, cls.audit_rows = read_csv(AUDIT_PATH)
        _, cls.s230_rows = read_csv(S230_PATH)
        cls.audit_by_source = {row["source"]: row for row in cls.audit_rows}
        cls.s230_by_source = {row["Source"]: row for row in cls.s230_rows}

    def test_audit_is_complete_unique_and_uses_the_fixed_schema(self) -> None:
        self.assertEqual(self.audit_fields, EXPECTED_COLUMNS)
        self.assertEqual(len(self.audit_rows), 19)
        self.assertEqual(len(self.audit_by_source), 19)
        self.assertEqual(set(self.audit_by_source), set(EXPECTED_DISPOSITIONS))
        self.assertEqual(
            {source: row["identity_disposition"] for source, row in self.audit_by_source.items()},
            EXPECTED_DISPOSITIONS,
        )

    def test_selected_observations_match_the_s230_product_exactly(self) -> None:
        self.assertEqual(set(self.audit_by_source), set(self.s230_by_source))
        for source, audit_row in self.audit_by_source.items():
            with self.subTest(source=source):
                s230_row = self.s230_by_source[source]
                self.assertEqual(audit_row["selected_flux_status"], s230_row["flux_status"])
                self.assertEqual(audit_row["selected_observation_id"], s230_row["nu0_observation_id"])
                self.assertEqual(audit_row["selected_reference"], s230_row["nu0_reference"])
                self.assertNotIn("/home/", audit_row["selected_reference"])

    def test_upper_limits_and_detection_dispositions_cannot_be_confused(self) -> None:
        for source, row in self.audit_by_source.items():
            with self.subTest(source=source):
                if row["identity_disposition"] in {"upper_limit", "conservative_upper_limit"}:
                    self.assertEqual(row["selected_flux_status"], "upper_limit")
                    self.assertEqual(row["component_identity"], "no_component_adopted")
                    self.assertIn("censored_upper_limit", row["required_handling"])
                else:
                    self.assertEqual(row["selected_flux_status"], "detection")

        ngc2273 = self.audit_by_source["NGC2273"]
        self.assertEqual(ngc2273["identity_disposition"], "conservative_upper_limit")
        self.assertIn("forced_nitercal_0", ngc2273["evidence_basis"])
        self.assertIn("conservative_beam_region_bound", ngc2273["evidence_basis"])

    def test_material_detection_caveats_are_explicit(self) -> None:
        self.assertEqual(
            {
                source
                for source, row in self.audit_by_source.items()
                if row["identity_disposition"] == "ambiguous"
            },
            {"IC2560", "NGC1320"},
        )
        self.assertEqual(
            {
                source
                for source, row in self.audit_by_source.items()
                if row["identity_disposition"] == "blended"
            },
            {"NGC3079", "NGC4258"},
        )
        for source in ("NGC3079", "NGC4258"):
            self.assertEqual(
                self.audit_by_source[source]["required_handling"],
                "do_not_describe_as_pure_compact_core",
            )

    def test_eso558_g009_uses_the_reviewed_x600_detection(self) -> None:
        row = self.audit_by_source["ESO558-G009"]
        self.assertEqual(row["selected_flux_status"], "detection")
        self.assertEqual(
            row["selected_observation_id"],
            "ESO558-G009|ALMA|2014-12-13|313.11GHz|"
            "member.uid-A001-X122-X600.calibrated-final-cont-image-ESO558-G009.pbcor.fits",
        )
        self.assertEqual(
            row["selected_reference"],
            "ALMA/ESO558-G009/"
            "member.uid___A001_X122_X600.calibrated_final_cont_image-ESO558-G009.pbcor.fits",
        )
        self.assertEqual(row["identity_disposition"], "nuclear")
        self.assertEqual(row["component_identity"], "nuclear_radio_continuum")
        self.assertEqual(row["required_handling"], "use_selected_observation")

    def test_ngc4945_uses_source18_and_preserves_the_source17_rejection(self) -> None:
        audit_row = self.audit_by_source["NGC4945"]
        self.assertEqual(audit_row["component_identity"], "Emig2020_source18_AGN_maser_nucleus")

        canonical_path, canonical_id = resolve_record_reference(
            audit_row["canonical_component_reference"]
        )
        rejected_path, rejected_id = resolve_record_reference(
            audit_row["rejected_component_reference"]
        )
        self.assertEqual(canonical_path, rejected_path)

        _, component_rows = read_csv(canonical_path)
        components = {row["record_id"]: row for row in component_rows}
        canonical = components[canonical_id]
        rejected = components[rejected_id]

        self.assertEqual(canonical_id, "accepted-source18-fixed-canonical")
        self.assertEqual(canonical["decision"], "accepted")
        self.assertEqual(canonical["included_in_canonical_table"], "yes")
        self.assertEqual(canonical["component_identity"], "Emig2020_source18_AGN_maser_nucleus")
        self.assertEqual(rejected_id, "rejected-source17-archived-anchor")
        self.assertEqual(rejected["decision"], "rejected")
        self.assertEqual(rejected["included_in_canonical_table"], "no")
        self.assertEqual(rejected["component_identity"], "Emig2020_source17_SSC_SNR_associated")

        s230_row = self.s230_by_source["NGC4945"]
        self.assertAlmostEqual(
            float(s230_row["S_nu0_mJy"]),
            round(float(canonical["flux_mjy"]), 2),
            places=12,
        )
        self.assertNotEqual(
            float(s230_row["S_nu0_mJy"]),
            round(float(rejected["flux_mjy"]), 2),
        )


if __name__ == "__main__":
    unittest.main()
