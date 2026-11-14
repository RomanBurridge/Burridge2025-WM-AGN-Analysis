from __future__ import annotations

import copy
import hashlib
import io
import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from analysis.multi_freq_from_archive import frozen_fitter_inputs as frozen


PUBLIC_NVAS_EXTRA = (
    "VLA/NGC1068/14.9I0.13_AC0467_1999SEP08_1_123.U42.6S.imfits"
)


class FrozenManifestTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.rows = frozen.load_frozen_fitter_inputs()

    def test_manifest_hash_and_exact_counts(self) -> None:
        self.assertEqual(
            frozen.manifest_sha256(), frozen.EXPECTED_MANIFEST_SHA256
        )
        self.assertEqual(len(self.rows), 550)
        self.assertEqual(
            sum(row["discovery_member"] == "yes" for row in self.rows), 544
        )
        self.assertEqual(
            sum(row["discovery_member"] == "no" for row in self.rows), 6
        )
        self.assertEqual(
            sum(
                row["input_role"] == "archive_fitter_input"
                and row["telescope"] == "ALMA"
                for row in self.rows
            ),
            222,
        )
        self.assertEqual(
            sum(
                row["input_role"] == "archive_fitter_input"
                and row["telescope"] == "VLA"
                for row in self.rows
            ),
            307,
        )
        self.assertEqual(
            sum(row["input_role"] == "sma_discovery_input" for row in self.rows),
            15,
        )
        self.assertEqual(
            sum(
                row["input_role"] == "sma_selected_clean_support"
                for row in self.rows
            ),
            6,
        )

    def test_manifest_bytes_are_frozen_not_only_row_counts(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            altered = Path(temp) / "frozen_fitter_inputs.csv"
            altered.write_bytes(frozen.DEFAULT_MANIFEST.read_bytes() + b"\n")
            with self.assertRaisesRegex(ValueError, "manifest hash mismatch"):
                frozen.load_frozen_fitter_inputs(altered)

    def test_historical_307_vla_whitelist_excludes_public_341_superset(self) -> None:
        vla_paths = {
            row["relative_path"] for row in self.rows if row["telescope"] == "VLA"
        }
        self.assertEqual(len(vla_paths), 307)
        self.assertNotIn(PUBLIC_NVAS_EXTRA, vla_paths)

    def test_recovery_routes_are_public_safe_and_complete(self) -> None:
        serialized = frozen.DEFAULT_MANIFEST.read_text(encoding="utf-8")
        self.assertNotIn("/home/", serialized)
        self.assertNotIn("existing_destination", serialized)
        direct = [
            row for row in self.rows if row["recovery_kind"] == "direct_public_file"
        ]
        derived = [
            row
            for row in self.rows
            if row["recovery_kind"] == "derive_from_public_parent"
        ]
        sma = [
            row
            for row in self.rows
            if row["recovery_kind"] == "sma_archive_reduction_product"
        ]
        self.assertEqual((len(direct), len(derived), len(sma)), (524, 5, 21))
        for row in direct:
            self.assertTrue(row["recovery_url"].startswith("https://"))
            self.assertEqual(row["recovery_sha256"], row["expected_sha256"])
        for row in derived:
            self.assertIn(".IQUV.", row["recovery_url"])
            self.assertNotEqual(
                Path(row["relative_path"]).name,
                Path(row["recovery_url"]).name,
            )
            self.assertEqual(
                row["transformation"], "extract_singleton_stokes_i_from_iquv"
            )
            self.assertRegex(row["recovery_sha256"], r"^[0-9a-f]{64}$")

    def test_selected_clean_support_is_three_image_residual_pairs(self) -> None:
        support = [
            row
            for row in self.rows
            if row["input_role"] == "sma_selected_clean_support"
        ]
        grouped: dict[str, set[str]] = {}
        for row in support:
            grouped.setdefault(row["source"], set()).add(
                Path(row["relative_path"]).suffix
            )
        self.assertEqual(
            grouped,
            {
                "NGC1068": {".image", ".residual"},
                "NGC3079": {".image", ".residual"},
                "NGC4258": {".image", ".residual"},
            },
        )

    def test_manifest_validation_fails_on_extra_vla_identity(self) -> None:
        rows = copy.deepcopy(self.rows)
        target = next(row for row in rows if row["telescope"] == "VLA")
        target["relative_path"] = PUBLIC_NVAS_EXTRA
        with self.assertRaisesRegex(ValueError, "source token|archive path|duplicate|recovery"):
            frozen.validate_manifest_rows(rows)

    def test_manifest_validation_fails_on_duplicate_or_unsafe_path(self) -> None:
        duplicate = copy.deepcopy(self.rows)
        duplicate[1]["relative_path"] = duplicate[0]["relative_path"]
        with self.assertRaisesRegex(ValueError, "duplicate path"):
            frozen.validate_manifest_rows(duplicate)
        unsafe = copy.deepcopy(self.rows)
        unsafe[0]["relative_path"] = "ALMA/../outside.fits"
        with self.assertRaisesRegex(ValueError, "unsafe or non-normalized"):
            frozen.validate_manifest_rows(unsafe)

    def test_derived_row_cannot_use_nonexistent_output_url(self) -> None:
        rows = copy.deepcopy(self.rows)
        row = next(
            row
            for row in rows
            if row["recovery_kind"] == "derive_from_public_parent"
        )
        row["recovery_url"] = (
            "https://almascience.eso.org/dataPortal/" + Path(row["relative_path"]).name
        )
        with self.assertRaisesRegex(ValueError, "IQUV parent"):
            frozen.validate_manifest_rows(rows)


class FrozenHashAndSelectionTests(unittest.TestCase):
    def test_workflow_tree_hash_matches_frozen_algorithm_fixture(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "empty").mkdir()
            (root / "nested").mkdir()
            (root / "a.txt").write_bytes(b"alpha")
            (root / "nested/b.bin").write_bytes(b"beta")
            (root / "nested/link").symlink_to("../a.txt")
            self.assertEqual(
                frozen.sha256_workflow_tree(root),
                "de4e9a8be088925c85369bfe33ba709c49fe5a511177423517cb4f6a9ef7f3ff",
            )

    def test_validate_and_select_returns_exact_whitelists(self) -> None:
        rows = frozen.load_frozen_fitter_inputs()
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            for row in rows:
                product = frozen._product_path(root, row)
                if row["sha256_algorithm"] == "sha256-file-v1":
                    product.parent.mkdir(parents=True, exist_ok=True)
                    product.touch()
                else:
                    product.mkdir(parents=True, exist_ok=True)
            with mock.patch.object(
                frozen, "_verify_product", side_effect=lambda path, row: row["expected_sha256"]
            ) as verifier:
                selection = frozen.validate_and_select_frozen_inputs(
                    root, workers=2
                )
            self.assertEqual(verifier.call_count, 550)
            self.assertEqual(selection.manifest_sha256, frozen.EXPECTED_MANIFEST_SHA256)
            self.assertEqual(len(selection.discovery_paths), 544)
            self.assertEqual(len(selection.selected_clean_paths), 6)
            selected_text = {str(path.relative_to(root)) for path in selection.discovery_paths}
            self.assertNotIn(
                "analysis/multi_freq_from_archive/" + PUBLIC_NVAS_EXTRA,
                selected_text,
            )


class FakeResponse(io.BytesIO):
    def __init__(self, payload: bytes, url: str) -> None:
        super().__init__(payload)
        self._url = url

    def __enter__(self) -> "FakeResponse":
        return self

    def __exit__(self, *_args: object) -> None:
        self.close()

    def geturl(self) -> str:
        return self._url


class FrozenRecoveryTests(unittest.TestCase):
    def test_download_does_not_replace_a_concurrently_created_destination(self) -> None:
        payload = b"approved archive bytes"
        url = "https://www.vla.nrao.edu/example.imfits"
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp).resolve()
            destination = root / "example.imfits"

            class ConcurrentResponse(FakeResponse):
                def __enter__(self):
                    destination.write_bytes(b"concurrent user data")
                    return self

            with mock.patch.object(frozen.urllib.request, "urlopen",
                                   return_value=ConcurrentResponse(payload, url)):
                with self.assertRaises(FileExistsError):
                    frozen._download_verified(url, hashlib.sha256(payload).hexdigest(),
                                              destination, 10, root)
            self.assertEqual(destination.read_bytes(), b"concurrent user data")

    def test_failed_download_cleanup_never_removes_a_foreign_partial(self) -> None:
        url = "https://www.vla.nrao.edu/example.imfits"
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp).resolve()
            destination = root / "example.imfits"
            foreign = root / "example.imfits.part"

            class ForeignPartialResponse(FakeResponse):
                def __enter__(self):
                    foreign.write_bytes(b"another download owns this")
                    return self

            with mock.patch.object(frozen.urllib.request, "urlopen",
                                   return_value=ForeignPartialResponse(b"bad", url)):
                with self.assertRaisesRegex(ValueError, "download hash mismatch"):
                    frozen._download_verified(url, hashlib.sha256(b"good").hexdigest(),
                                              destination, 10, root)
            self.assertEqual(foreign.read_bytes(), b"another download owns this")

    def test_derived_install_does_not_replace_a_concurrent_destination(self) -> None:
        payload = b"approved Stokes-I"
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp).resolve()
            parent, destination = root / "parent.fits", root / "output.fits"
            parent.write_bytes(b"parent")

            def concurrent_extractor(_parent, partial):
                destination.write_bytes(b"concurrent user data")
                partial.write_bytes(payload)

            with self.assertRaises(FileExistsError):
                frozen._derive_stokes_i_verified(
                    parent, hashlib.sha256(b"parent").hexdigest(), destination,
                    hashlib.sha256(payload).hexdigest(), root, concurrent_extractor)
            self.assertEqual(destination.read_bytes(), b"concurrent user data")

    def test_failed_extraction_cleanup_preserves_a_foreign_partial(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp).resolve()
            parent, destination = root / "parent.fits", root / "output.fits"
            parent.write_bytes(b"parent")
            foreign = root / "output.fits.stokes-i.part"

            def failed_extractor(_parent, partial):
                foreign.write_bytes(b"another extraction owns this")
                partial.write_bytes(b"partial output")
                raise RuntimeError("synthetic extraction failure")

            with self.assertRaisesRegex(RuntimeError, "synthetic extraction failure"):
                frozen._derive_stokes_i_verified(
                    parent, hashlib.sha256(b"parent").hexdigest(), destination,
                    hashlib.sha256(b"unused").hexdigest(), root, failed_extractor)
            self.assertEqual(foreign.read_bytes(), b"another extraction owns this")

    def test_verified_download_is_no_clobber(self) -> None:
        payload = b"archive bytes"
        expected = hashlib.sha256(payload).hexdigest()
        url = "https://www.vla.nrao.edu/example.imfits"
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp).resolve()
            destination = root / "analysis/multi_freq_from_archive/VLA/S/example.imfits"
            with mock.patch.object(
                frozen.urllib.request,
                "urlopen",
                return_value=FakeResponse(payload, url),
            ):
                self.assertEqual(
                    frozen._download_verified(url, expected, destination, 10, root),
                    "downloaded_verified",
                )
            self.assertEqual(destination.read_bytes(), payload)
            with mock.patch.object(frozen.urllib.request, "urlopen") as opener:
                self.assertEqual(
                    frozen._download_verified(url, expected, destination, 10, root),
                    "already_verified",
                )
            opener.assert_not_called()

    def test_failed_download_never_installs_bad_bytes(self) -> None:
        url = "https://www.vla.nrao.edu/example.imfits"
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp).resolve()
            destination = root / "analysis/multi_freq_from_archive/VLA/S/example.imfits"
            with mock.patch.object(
                frozen.urllib.request,
                "urlopen",
                return_value=FakeResponse(b"wrong", url),
            ):
                with self.assertRaisesRegex(ValueError, "download hash mismatch"):
                    frozen._download_verified(
                        url, hashlib.sha256(b"right").hexdigest(), destination, 10, root
                    )
            self.assertFalse(destination.exists())
            self.assertFalse(destination.with_name(destination.name + ".part").exists())

    def test_nonmatching_direct_output_is_never_overwritten(self) -> None:
        url = "https://www.vla.nrao.edu/example.imfits"
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp).resolve()
            destination = root / "analysis/multi_freq_from_archive/VLA/S/example.imfits"
            destination.parent.mkdir(parents=True)
            destination.write_bytes(b"preexisting user data")
            with mock.patch.object(frozen.urllib.request, "urlopen") as opener:
                with self.assertRaisesRegex(FileExistsError, "nonmatching file"):
                    frozen._download_verified(
                        url,
                        hashlib.sha256(b"expected archive bytes").hexdigest(),
                        destination,
                        10,
                        root,
                    )
            opener.assert_not_called()
            self.assertEqual(destination.read_bytes(), b"preexisting user data")

    def test_symlinked_destination_is_rejected_before_network(self) -> None:
        url = "https://www.vla.nrao.edu/example.imfits"
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp).resolve()
            outside = root / "outside"
            outside.mkdir()
            parent = root / "analysis"
            parent.symlink_to(outside, target_is_directory=True)
            destination = parent / "example.imfits"
            with mock.patch.object(frozen.urllib.request, "urlopen") as opener:
                with self.assertRaisesRegex(ValueError, "symlinked recovery parent"):
                    frozen._download_verified(
                        url, hashlib.sha256(b"x").hexdigest(), destination, 10, root
                    )
            opener.assert_not_called()

    def test_derived_output_is_exact_hash_verified_and_no_clobber(self) -> None:
        parent_payload = b"reviewed IQUV parent"
        output_payload = b"deterministic Stokes-I product"
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp).resolve()
            parent = root / ".frozen_input_parents/S/parent.IQUV.fits"
            destination = (
                root
                / "analysis/multi_freq_from_archive/ALMA/S/product.I.pbcor.fits"
            )
            parent.parent.mkdir(parents=True)
            parent.write_bytes(parent_payload)

            def extractor(_parent: Path, output: Path) -> None:
                output.write_bytes(output_payload)

            self.assertEqual(
                frozen._derive_stokes_i_verified(
                    parent,
                    hashlib.sha256(parent_payload).hexdigest(),
                    destination,
                    hashlib.sha256(output_payload).hexdigest(),
                    root,
                    extractor,
                ),
                "derived_verified",
            )
            self.assertEqual(destination.read_bytes(), output_payload)
            unused_extractor = mock.Mock()
            self.assertEqual(
                frozen._derive_stokes_i_verified(
                    parent,
                    hashlib.sha256(parent_payload).hexdigest(),
                    destination,
                    hashlib.sha256(output_payload).hexdigest(),
                    root,
                    unused_extractor,
                ),
                "already_verified_derived",
            )
            unused_extractor.assert_not_called()

    def test_nonmatching_derived_output_is_never_overwritten(self) -> None:
        parent_payload = b"reviewed IQUV parent"
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp).resolve()
            parent = root / ".frozen_input_parents/S/parent.IQUV.fits"
            destination = (
                root
                / "analysis/multi_freq_from_archive/ALMA/S/product.I.pbcor.fits"
            )
            parent.parent.mkdir(parents=True)
            destination.parent.mkdir(parents=True)
            parent.write_bytes(parent_payload)
            destination.write_bytes(b"preexisting user data")
            extractor = mock.Mock()
            with self.assertRaisesRegex(FileExistsError, "nonmatching derived"):
                frozen._derive_stokes_i_verified(
                    parent,
                    hashlib.sha256(parent_payload).hexdigest(),
                    destination,
                    hashlib.sha256(b"expected output").hexdigest(),
                    root,
                    extractor,
                )
            extractor.assert_not_called()
            self.assertEqual(destination.read_bytes(), b"preexisting user data")

    def test_bad_derived_parent_hash_is_rejected_before_extraction(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp).resolve()
            parent = root / ".frozen_input_parents/S/parent.IQUV.fits"
            destination = (
                root
                / "analysis/multi_freq_from_archive/ALMA/S/product.I.pbcor.fits"
            )
            parent.parent.mkdir(parents=True)
            parent.write_bytes(b"wrong parent")
            extractor = mock.Mock()
            with self.assertRaisesRegex(ValueError, "IQUV parent hash mismatch"):
                frozen._derive_stokes_i_verified(
                    parent,
                    hashlib.sha256(b"expected parent").hexdigest(),
                    destination,
                    hashlib.sha256(b"expected output").hexdigest(),
                    root,
                    extractor,
                )
            extractor.assert_not_called()
            self.assertFalse(destination.exists())


if __name__ == "__main__":
    unittest.main()
