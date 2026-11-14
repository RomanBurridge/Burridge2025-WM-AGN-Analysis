import ast
import csv
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path


PACKAGE_DIR = Path(__file__).resolve().parents[1]
if str(PACKAGE_DIR) not in sys.path:
    sys.path.insert(0, str(PACKAGE_DIR))

from chosen_frequency_anchors import (  # noqa: E402
    SELECTION_RULE,
    angular_separation_arcsec,
    check_materialized_manifest,
    compare_reference_identities,
    consolidate,
    path_sha256,
    read_manifest,
    read_component_policy,
    read_reference,
    review_selected_anchors,
    select_anchors,
    sexagesimal_dec_to_deg,
    sexagesimal_ra_to_deg,
    tree_sha256,
    write_manifest,
)


REFERENCE = PACKAGE_DIR / "chosen_frequency_anchor_reference.csv"
COMPONENT_POLICY = PACKAGE_DIR / "chosen_frequency_component_policy.csv"
HISTORICAL_FINAL_STAGE = PACKAGE_DIR / "historical_final_stage_validation.json"


def write_policy(path: Path, rows: list[dict[str, str]]) -> None:
    fields = (
        "source",
        "selected_reference",
        "identity_disposition",
        "component_identity",
        "evidence_basis",
        "required_handling",
        "reviewed_frame",
        "reviewed_ra_deg",
        "reviewed_dec_deg",
        "reviewed_coordinate_basis",
        "max_separation_arcsec",
        "generic_seed_action",
    )
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def make_record(
    root: Path,
    *,
    source: str = "Example",
    product: str = "ALMA/Example/example.fits",
    model_relative: str = "candidate/Example/example.fits.fit",
    major: float = 0.4,
    minor: float = 0.25,
) -> dict[str, object]:
    input_path = root / product
    input_path.parent.mkdir(parents=True, exist_ok=True)
    input_path.write_bytes(product.encode())
    model = root / model_relative
    model.mkdir(parents=True, exist_ok=True)
    (model / "table.dat").write_bytes((source + product).encode())
    return {
        "source": source,
        "telescope": "ALMA",
        "frequency_ghz": 230.0,
        "beam_major_arcsec": major,
        "beam_minor_arcsec": minor,
        "beam_pa_deg": 10.0,
        "beam_geometric_arcsec": math.sqrt(major * minor),
        "anchor_frame": "ICRS",
        "anchor_ra_deg": 12.5,
        "anchor_dec_deg": -33.0,
        "icrs_ra": "00h50m00s",
        "icrs_dec": "-33d00m00s",
        "j2000_ra": "00h50m00s",
        "j2000_dec": "-33d00m00s",
        "input_product": product,
        "input_sha256": path_sha256(input_path),
        "model_name": Path(model_relative).name,
        "model_path": model_relative,
        "model_tree_sha256": tree_sha256(model),
        "run_context": "230316_2",
        "selection_rule": SELECTION_RULE,
        "scientific_disposition": "generic_anchor_seed",
    }


class ChosenFrequencyAnchorTests(unittest.TestCase):
    def test_minimum_beam_and_lexical_tie_break_are_deterministic(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            broad = make_record(
                root,
                product="ALMA/Example/z.fits",
                model_relative="models/z.fit",
                major=0.5,
                minor=0.4,
            )
            lexical_later = make_record(
                root,
                product="ALMA/Example/b.fits",
                model_relative="models/b.fit",
                major=0.2,
                minor=0.2,
            )
            lexical_first = make_record(
                root,
                product="ALMA/Example/a.fits",
                model_relative="models/a.fit",
                major=0.2,
                minor=0.2,
            )
            selected = select_anchors([broad, lexical_later, lexical_first])
            self.assertEqual(len(selected), 1)
            self.assertEqual(selected[0]["input_product"], "ALMA/Example/a.fits")

    def test_consolidation_copies_model_and_rewrites_canonical_path(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            candidate = make_record(root)
            candidate_manifest = root / "chosen_frequency/candidates/run.csv"
            write_manifest(candidate_manifest, [candidate])
            policy_path = root / "component_policy.csv"
            write_policy(
                policy_path,
                [
                    {
                        "source": "Example",
                        "selected_reference": candidate["input_product"],
                        "identity_disposition": "nuclear",
                        "component_identity": "reviewed_nucleus",
                        "evidence_basis": "synthetic_test_astrometry",
                        "required_handling": "use_selected_observation",
                        "reviewed_frame": "ICRS",
                        "reviewed_ra_deg": "12.5",
                        "reviewed_dec_deg": "-33.0",
                        "reviewed_coordinate_basis": "synthetic_truth",
                        "max_separation_arcsec": "0.5",
                        "generic_seed_action": "require_coordinate_match",
                    }
                ],
            )
            output_manifest = root / "chosen_frequency/chosen_frequency_anchors.csv"
            output_models = root / "images/chosenfreq"
            selected = consolidate(
                [candidate_manifest],
                analysis_root=root,
                output_manifest=output_manifest,
                models_output=output_models,
                component_policy=policy_path,
            )
            self.assertEqual(len(selected), 1)
            copied = output_models / "Example/example.fits.fit"
            self.assertTrue(copied.is_dir())
            self.assertEqual(tree_sha256(copied), candidate["model_tree_sha256"])
            loaded = read_manifest(output_manifest)
            self.assertEqual(
                loaded[0]["model_path"], "images/chosenfreq/Example/example.fits.fit"
            )
            self.assertEqual(loaded[0]["scientific_disposition"], "reviewed_nuclear_anchor")
            self.assertEqual(loaded[0]["component_separation_arcsec"], "0")
            self.assertEqual(check_materialized_manifest(output_manifest, root), loaded)

    def test_invalid_frequency_and_parent_path_fail_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            record = make_record(Path(temporary))
            record["frequency_ghz"] = 199.99
            with self.assertRaisesRegex(ValueError, "outside 200--400"):
                select_anchors([record])
            record["frequency_ghz"] = 230
            record["model_path"] = "../escape.fit"
            with self.assertRaisesRegex(ValueError, "portable relative path"):
                select_anchors([record])

    def test_historical_reference_has_nine_identities_and_rejects_ngc4945_seed(self):
        rows = read_reference(REFERENCE)
        self.assertEqual(len(rows), 9)
        self.assertEqual(compare_reference_identities(rows, rows), [])
        ngc4945 = next(row for row in rows if row["source"] == "NGC4945")
        self.assertEqual(ngc4945["scientific_status"], "rejected_component_seed")
        self.assertIn("source18", ngc4945["scientific_note"])
        final_rows = [
            ast.literal_eval(line)
            for line in (
                PACKAGE_DIR
                / "fitsumfiles/historical/"
                "fitsummary_fixtochosen_sorted_409_reference.txt"
            ).read_text().splitlines()
            if line.strip()
        ]
        self.assertEqual(len(final_rows), 409)
        for row in rows:
            reference_row = final_rows[int(row["reference_409_row"]) - 1]
            self.assertEqual(reference_row[0], row["source"])
            self.assertEqual(str(reference_row[7]).strip(), row["input_product"])

    def test_reference_check_preserves_nine_and_allows_only_policy_bound_addition(self):
        historical = read_reference(REFERENCE)
        added = {
            "source": "NewAnchor",
            "telescope": "ALMA",
            "input_product": "ALMA/NewAnchor/new.fits",
            "model_name": "new.fits.fit",
        }
        actual = [*historical, added]
        component = {
            "NewAnchor": {"selected_reference": "ALMA/NewAnchor/new.fits"}
        }
        self.assertEqual(
            compare_reference_identities(
                actual,
                historical,
                reviewed_addition_products={"ALMA/NewAnchor/new.fits"},
                component_policy=component,
            ),
            [],
        )
        wrong_historical = [dict(row) for row in actual]
        wrong_historical[0]["input_product"] = "ALMA/CircinusGalaxy/replacement.fits"
        self.assertTrue(
            compare_reference_identities(
                wrong_historical,
                historical,
                reviewed_addition_products={"ALMA/NewAnchor/new.fits"},
                component_policy=component,
            )
        )
        self.assertTrue(
            compare_reference_identities(
                actual,
                historical,
                reviewed_addition_products=set(),
                component_policy=component,
            )
        )

    def test_fainter_nucleus_brighter_off_nuclear_minimum_beam_fails_closed(self):
        """A bright off-nuclear fit in the smallest-beam image fails closed."""

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            reviewed_nucleus = make_record(
                root,
                source="SyntheticAGN",
                product="ALMA/SyntheticAGN/fainter_nucleus.fits",
                model_relative="models/nucleus.fit",
                major=0.4,
                minor=0.3,
            )
            brighter_off_nuclear = make_record(
                root,
                source="SyntheticAGN",
                product="ALMA/SyntheticAGN/high_resolution_two_components.fits",
                model_relative="models/off_nuclear.fit",
                major=0.1,
                minor=0.08,
            )
            # Simulate the free Gaussian fit locking onto an 8 mJy off-nuclear
            # component instead of the reviewed 1 mJy nucleus in this image.
            synthetic_peak_mjy = {"reviewed_nucleus": 1.0, "fitted_component": 8.0}
            self.assertGreater(
                synthetic_peak_mjy["fitted_component"],
                synthetic_peak_mjy["reviewed_nucleus"],
            )
            brighter_off_nuclear["anchor_ra_deg"] = 12.502
            generic_winner = select_anchors([reviewed_nucleus, brighter_off_nuclear])
            self.assertEqual(
                generic_winner[0]["input_product"],
                "ALMA/SyntheticAGN/high_resolution_two_components.fits",
            )
            policy_path = root / "policy.csv"
            write_policy(
                policy_path,
                [
                    {
                        "source": "SyntheticAGN",
                        "selected_reference": brighter_off_nuclear["input_product"],
                        "identity_disposition": "nuclear",
                        "component_identity": "fainter_reviewed_nucleus",
                        "evidence_basis": "synthetic_maser_astrometry",
                        "required_handling": "use_selected_observation",
                        "reviewed_frame": "ICRS",
                        "reviewed_ra_deg": "12.5",
                        "reviewed_dec_deg": "-33.0",
                        "reviewed_coordinate_basis": "synthetic_truth",
                        "max_separation_arcsec": "0.5",
                        "generic_seed_action": "require_coordinate_match",
                    }
                ],
            )
            with self.assertRaisesRegex(ValueError, "component coordinate mismatch"):
                review_selected_anchors(
                    generic_winner,
                    read_component_policy(policy_path),
                )

    def test_checked_in_policy_preserves_historical_components_and_rejects_source17(self):
        references = {row["source"]: row for row in read_reference(REFERENCE)}
        policy = read_component_policy(COMPONENT_POLICY)
        task28_replacements = {
            "NGC1068": (
                "ALMA/NGC1068/member.uid___A001_X2d20_X3f3f."
                "NGC_1068_sci.spw25_27_29_31.cont.I.pbcor.fits"
            ),
            "NGC1194": (
                "ALMA/NGC1194/member.uid___A001_X3827_X338."
                "NGC1194_sci.spw17_19_25_27.cont.I.pbcor.fits"
            ),
        }
        # New-source component policies may be staged before their exact
        # product-admission row.  Runtime anchor validation allows such a
        # source only when both policies bind the same selected product.
        self.assertTrue(set(references).issubset(policy))
        for source, reference in references.items():
            decision = policy[source]
            if source in task28_replacements:
                self.assertNotEqual(
                    reference["input_product"], decision["selected_reference"]
                )
                self.assertEqual(
                    decision["selected_reference"], task28_replacements[source]
                )
                self.assertIn(
                    "task28_smaller_beam_selection",
                    decision["reviewed_coordinate_basis"],
                )
            else:
                self.assertEqual(
                    reference["input_product"], decision["selected_reference"]
                )
            self.assertEqual(reference["anchor_frame"], decision["reviewed_frame"])
            separation = angular_separation_arcsec(
                float(reference["anchor_ra_deg"]),
                float(reference["anchor_dec_deg"]),
                float(decision["reviewed_ra_deg"]),
                float(decision["reviewed_dec_deg"]),
            )
            if source == "NGC4945":
                self.assertGreater(separation, float(decision["max_separation_arcsec"]))
                self.assertEqual(
                    decision["generic_seed_action"],
                    "reject_generic_seed_and_rerun_fixed_reviewed_nucleus",
                )
            else:
                self.assertLessEqual(
                    separation, float(decision["max_separation_arcsec"])
                )

    def test_angular_separation_handles_ra_wrap(self):
        self.assertAlmostEqual(
            angular_separation_arcsec(359.999, 0.0, 0.001, 0.0),
            7.2,
            places=8,
        )

    def test_historical_final_stage_evidence_is_scoped_and_hash_locked(self):
        evidence = json.loads(HISTORICAL_FINAL_STAGE.read_text())
        self.assertEqual(evidence["outputs"]["accepted_nonempty_rows"], 409)
        self.assertEqual(evidence["outputs"]["excluded_archive_paths"], 46)
        self.assertEqual(evidence["outputs"]["scientific_comparison"], "PASS")
        self.assertEqual(
            evidence["anchor_identity_reference"],
            "chosen_frequency_anchor_reference.csv",
        )
        self.assertIn("historical fourth stage", evidence["scientific_limit"])
        self.assertIn("NGC 4945 source-18", evidence["scientific_limit"])
        for digest in (
            *evidence["source_evidence_sha256"].values(),
            evidence["outputs"]["generated_sorted_summary_sha256"],
            evidence["outputs"]["generated_bad_paths_sha256"],
        ):
            self.assertEqual(len(digest), 64)
            self.assertEqual(set(digest).difference("0123456789abcdef"), set())

    def test_sexagesimal_coordinate_conversion(self):
        self.assertAlmostEqual(sexagesimal_ra_to_deg("10h16m18.704s"), 154.0779333333)
        self.assertAlmostEqual(sexagesimal_dec_to_deg("-33d33m49.729s"), -33.5638136111)


if __name__ == "__main__":
    unittest.main()
