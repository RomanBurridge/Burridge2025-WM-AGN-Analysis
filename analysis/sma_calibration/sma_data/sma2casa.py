#!/usr/bin/env python3
"""Convert an already-binned SMA MIR directory to a NEW Measurement Set.

The archive request, not this wrapper, supplies spectral averaging. MIR read
and MS-write numerical options retain pyuvdata's defaults. The former unused
UVFITS export has been removed: it did not feed the MS and could overwrite
and then delete a user's existing UVFITS file.

A failed conversion leaves any partial MS in place for diagnosis. Reruns
must use a new output path, not silently overwrite or delete that evidence.
"""
from __future__ import annotations

import argparse
from pathlib import Path


def process_folder(folder, *, output_ms=None, uvdata_factory=None):
    """Convert a real MIR directory without changing or deleting user files.

    The documented basename invocation keeps its legacy output location beside
    this script. Absolute/relative input paths are now resolved normally;
    ``output_ms`` may explicitly select another new output path.
    ``uvdata_factory`` is dependency injection for file-safety tests only.
    """
    requested_input = Path(folder).expanduser()
    if requested_input.is_symlink() or not requested_input.is_dir():
        raise ValueError(f"MIR input must be a real directory: {requested_input}")
    source = requested_input.resolve(strict=True)
    requested_output = (
        Path(output_ms).expanduser() if output_ms is not None
        else Path(__file__).resolve().parent / f"{source.name}_master.ms"
    )
    if requested_output.is_symlink() or requested_output.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {requested_output}")
    destination = requested_output.resolve()
    if destination == source or source in destination.parents:
        raise ValueError("The MS output must not be the MIR input or inside it")
    if destination.suffix != ".ms":
        raise ValueError("The new Measurement Set output must end in .ms")
    if not destination.parent.is_dir():
        raise ValueError(f"MS output parent directory does not exist: {destination.parent}")

    if uvdata_factory is None:
        try:
            from casacore import tables  # noqa: F401: fail before loading MIR
            from pyuvdata import UVData
        except ImportError as exc:
            raise RuntimeError(
                "MIR conversion requires pyuvdata and its optional casa/MS writer "
                "dependency (python-casacore); install the documented MIR extra"
            ) from exc
        uvdata_factory = UVData
    uv = uvdata_factory()
    uv.read_mir(str(source))
    # pyuvdata's own no-clobber guard is also retained in addition to our
    # early check. No temporary UVFITS or cleanup/deletion path is needed.
    uv.write_ms(str(destination), clobber=False)
    return destination


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path, help="Already-binned MIR input directory")
    parser.add_argument("--output-ms", type=Path, help="New .ms output path; must not already exist")
    args = parser.parse_args(argv)
    output = process_folder(args.folder, output_ms=args.output_ms)
    print(f"Wrote new Measurement Set: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
