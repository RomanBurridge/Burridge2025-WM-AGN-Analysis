from __future__ import annotations

import ast
import sys
import types
import unittest
from pathlib import Path


PACKAGE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE_DIR))

from celestial_frames import (  # noqa: E402
    casa_sexagesimal_declination,
    colon_sexagesimal_declination,
    casa_direction_for_reference_code,
    coordinate_pair_for_reference_code,
    coordinate_candidate_attempt_limit,
    discard_failed_coordinate_candidate_or_raise,
    normalize_celestial_reference_code,
)


ICRS_PAIR = ("10h16m18.700s", "-33d33m50.10s")
J2000_PAIR = ("10h16m18.710s", "-33d33m49.92s")


class CelestialFrameRoutingTests(unittest.TestCase):
    @staticmethod
    def _load_getbestrez(**overrides):
        """Load only getbestrez without executing the fitter's CLI body."""

        fitter_path = PACKAGE_DIR / "get_continuum_params.py"
        parsed = ast.parse(fitter_path.read_text(), filename=str(fitter_path))
        function = next(
            node
            for node in parsed.body
            if isinstance(node, ast.FunctionDef) and node.name == "getbestrez"
        )
        namespace = {
            "ICRS": "ICRS",
            "coordinate_candidate_attempt_limit": coordinate_candidate_attempt_limit,
            "discard_failed_coordinate_candidate_or_raise": (
                discard_failed_coordinate_candidate_or_raise
            ),
        }
        namespace.update(overrides)
        exec(
            compile(
                ast.fix_missing_locations(ast.Module(body=[function], type_ignores=[])),
                str(fitter_path),
                "exec",
            ),
            namespace,
        )
        return namespace["getbestrez"]

    def test_declination_normalization_preserves_fractional_seconds(self) -> None:
        expected_colon = "-065:20:20.952367"
        expected_casa = "-065d20m20.952367s"
        for raw_declination in (
            "-065:20:20.952367",
            "-065d20m20.952367s",
            "-065.20.20.952367",
        ):
            with self.subTest(raw_declination=raw_declination):
                self.assertEqual(
                    colon_sexagesimal_declination(raw_declination), expected_colon
                )
                self.assertEqual(
                    casa_sexagesimal_declination(raw_declination), expected_casa
                )

    def test_malformed_declinations_fail_closed(self) -> None:
        for raw_declination in ("", "-065:20", "-065d20m", "not-a-coordinate"):
            with self.subTest(raw_declination=raw_declination):
                with self.assertRaises(ValueError):
                    colon_sexagesimal_declination(raw_declination)

    def test_coordinate_candidate_retries_are_finite_and_make_progress(self) -> None:
        candidates = ["source", [["first.fit", "second.fit"], [], []]]
        self.assertEqual(
            coordinate_candidate_attempt_limit(candidates, choosefile=False), 2
        )

        first_error = ValueError("bad first coordinate")
        discard_failed_coordinate_candidate_or_raise(
            candidates,
            "first.fit",
            choosefile=False,
            error=first_error,
        )
        self.assertIsNone(candidates[1][0][0])
        self.assertEqual(candidates[1][0][1], "second.fit")

        with self.assertRaisesRegex(RuntimeError, "not present"):
            discard_failed_coordinate_candidate_or_raise(
                candidates,
                "missing.fit",
                choosefile=False,
                error=ValueError("bad missing coordinate"),
            )

    def test_chosen_coordinate_failure_raises_instead_of_retrying(self) -> None:
        candidates = ["source", [["chosen.fit"], [], []]]
        self.assertEqual(
            coordinate_candidate_attempt_limit(candidates, choosefile=True), 1
        )
        with self.assertRaisesRegex(RuntimeError, "cannot be retried"):
            discard_failed_coordinate_candidate_or_raise(
                candidates,
                "chosen.fit",
                choosefile=True,
                error=ValueError("bad chosen coordinate"),
            )
        self.assertEqual(candidates[1][0], ["chosen.fit"])

    def test_getbestrez_chosen_file_failure_attempts_only_once(self) -> None:
        imstat_calls = []

        def imstat(path, region=None):
            imstat_calls.append((path, region))
            return {"maxposf": "14:13:09.986, -65.20.20.970, 93GHz, I"}

        def fail_conversion(*_args):
            raise ValueError("bad converted coordinate")

        getbestrez = self._load_getbestrez(
            au=types.SimpleNamespace(
                J2000ToICRS=lambda *_args, **_kwargs: (
                    "14:13:09.985, -065:20:20.988"
                )
            ),
            sortfiles=lambda *_args: (None, None),
            image_celestial_frame=lambda _path: "ICRS",
            imstat=imstat,
            coherent_coordinate_pairs=fail_conversion,
        )
        chosen_file = "chosen.fit"
        with self.assertRaisesRegex(RuntimeError, chosen_file) as caught:
            getbestrez(
                ["source", [[chosen_file], [], []]],
                ["14:13:09.964", "-65:20:20.742"],
                "ALMA",
                ["ALMA", "VLA", "SMA"],
                1,
                chosen_file,
            )
        self.assertIsInstance(caught.exception.__cause__, ValueError)
        self.assertEqual(imstat_calls, [(chosen_file, None)])

    def test_getbestrez_automatic_selection_discards_failure_then_succeeds(self) -> None:
        candidates = ["source", [["bad.fit", "good.fit"], [], []]]
        imstat_calls = []

        def sortfiles(getfiles, _alltele, _tele):
            return next(
                (candidate, "ALMA")
                for candidate in getfiles[1][0]
                if candidate is not None
            )

        def imstat(path, region=None):
            imstat_calls.append((path, region))
            return {"maxposf": "14:13:09.986, -65.20.20.970, 93GHz, I"}

        def convert(_ra, _dec, _frame):
            if imstat_calls[-1][0] == "bad.fit":
                raise ValueError("bad coordinate")
            return (
                ["14h13m09.986s", "-65d20m20.970s"],
                ["14h13m09.986s", "-65d20m20.952s"],
            )

        getbestrez = self._load_getbestrez(
            au=types.SimpleNamespace(
                J2000ToICRS=lambda *_args, **_kwargs: (
                    "14:13:09.985, -065:20:20.988"
                )
            ),
            sortfiles=sortfiles,
            image_celestial_frame=lambda _path: "ICRS",
            imstat=imstat,
            coherent_coordinate_pairs=convert,
        )
        result = getbestrez(
            candidates,
            ["14:13:09.964", "-65:20:20.742"],
            "NA",
            ["ALMA", "VLA", "SMA"],
            0,
            None,
        )
        self.assertEqual(result[2], "good.fit")
        self.assertIsNone(candidates[1][0][0])
        self.assertEqual([call[0] for call in imstat_calls], ["bad.fit", "good.fit"])

    def test_casa_referencecode_sequence_routes_icrs(self) -> None:
        self.assertEqual(normalize_celestial_reference_code(["ICRS"]), "ICRS")
        self.assertEqual(
            coordinate_pair_for_reference_code(
                ["ICRS"], icrs=ICRS_PAIR, j2000=J2000_PAIR
            ),
            ICRS_PAIR,
        )

    def test_fits_fk5_and_casa_j2000_route_the_coherent_j2000_pair(self) -> None:
        for metadata in ("FK5", "J2000", ["J2000"], "FK5_J2000"):
            with self.subTest(metadata=metadata):
                self.assertEqual(
                    coordinate_pair_for_reference_code(
                        metadata, icrs=ICRS_PAIR, j2000=J2000_PAIR
                    ),
                    J2000_PAIR,
                )

    def test_direction_label_and_values_use_the_same_frame(self) -> None:
        self.assertEqual(
            casa_direction_for_reference_code(
                "FK5", icrs=ICRS_PAIR, j2000=J2000_PAIR
            ),
            f"J2000 {J2000_PAIR[0]} {J2000_PAIR[1]}",
        )
        self.assertEqual(
            casa_direction_for_reference_code(
                "ICRS", icrs=ICRS_PAIR, j2000=J2000_PAIR
            ),
            f"ICRS {ICRS_PAIR[0]} {ICRS_PAIR[1]}",
        )

    def test_unknown_missing_or_conflicting_metadata_fails_closed(self) -> None:
        for metadata in (None, "", [], "B1950", ["ICRS", "J2000"]):
            with self.subTest(metadata=metadata):
                with self.assertRaises(ValueError):
                    normalize_celestial_reference_code(metadata)

    def test_coordinate_pairs_must_be_complete(self) -> None:
        with self.assertRaises(ValueError):
            coordinate_pair_for_reference_code(
                "ICRS", icrs=("only-ra",), j2000=J2000_PAIR
            )
        with self.assertRaises(ValueError):
            coordinate_pair_for_reference_code(
                "ICRS", icrs=("", "dec"), j2000=J2000_PAIR
            )


if __name__ == "__main__":
    unittest.main()
