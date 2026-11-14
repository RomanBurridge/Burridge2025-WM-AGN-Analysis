#!/usr/bin/env python3
"""Fail-fast checks for the two released SMA calibration recipes.

This reads serialized spectral-window metadata only. It does not certify an
earlier raw-MIR conversion, apply averaging, rescale weights, change flags or
establish the absolute flux scale.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


SUPPORTED_TRACKS = {
    220131: {"native_channels": 16384, "spectral_windows": 12, "bandwidth_hz": 2288000000.0},
    230316: {"native_channels": 16384, "spectral_windows": 12, "bandwidth_hz": 2288000000.0},
}
CHANNEL_COLUMNS = ("CHAN_FREQ", "CHAN_WIDTH", "EFFECTIVE_BW", "RESOLUTION")


def validate_spectral_windows(records, *, date, averaging=128):
    """Validate count, spacing, bandwidth and the declared averaging layout."""
    try:
        expected = SUPPORTED_TRACKS[int(date)]
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"No released SMA channel-layout contract for date {date}") from exc
    if (not isinstance(averaging, int) or isinstance(averaging, bool)
            or averaging <= 0 or expected["native_channels"] % averaging):
        raise ValueError("Averaging must be a positive integer divisor of native channel count")
    expected_channels = expected["native_channels"] // averaging
    if len(records) != expected["spectral_windows"]:
        raise ValueError(f"Expected {expected['spectral_windows']} spectral windows, found {len(records)}")
    summaries = []
    for index, row in enumerate(records):
        channels = int(row["NUM_CHAN"])
        if channels != expected_channels:
            raise ValueError(f"SPW {index}: expected {expected_channels} already-binned channels, found {channels}")
        arrays = {}
        for column in CHANNEL_COLUMNS:
            values = [float(value) for value in row[column]]
            if len(values) != channels or not all(math.isfinite(value) for value in values):
                raise ValueError(f"SPW {index}: invalid {column} shape/nonfinite value")
            if column in ("CHAN_FREQ", "EFFECTIVE_BW", "RESOLUTION") and not all(value > 0 for value in values):
                raise ValueError(f"SPW {index}: {column} must be positive")
            if column == "CHAN_WIDTH" and not all(value != 0 for value in values):
                raise ValueError(f"SPW {index}: CHAN_WIDTH must be nonzero")
            arrays[column] = values
        bandwidth = float(row["TOTAL_BANDWIDTH"])
        if not math.isfinite(bandwidth) or not math.isclose(
                abs(bandwidth), expected["bandwidth_hz"], rel_tol=1e-8, abs_tol=1.0):
            raise ValueError(f"SPW {index}: total bandwidth differs from the released track contract")
        width = expected["bandwidth_hz"] / expected_channels
        for column in ("CHAN_WIDTH", "EFFECTIVE_BW", "RESOLUTION"):
            if not all(math.isclose(abs(value), width, rel_tol=1e-8, abs_tol=1.0)
                       for value in arrays[column]):
                raise ValueError(f"SPW {index}: {column} is inconsistent with averaging/bandwidth")
        frequencies = arrays["CHAN_FREQ"]
        differences = [b - a for a, b in zip(frequencies, frequencies[1:])]
        if (differences and (not all(value > 0 for value in differences)
                             and not all(value < 0 for value in differences))):
            raise ValueError(f"SPW {index}: channel frequencies are not monotonic")
        if not all(math.isclose(abs(value), width, rel_tol=1e-8, abs_tol=1.0)
                   for value in differences):
            raise ValueError(f"SPW {index}: channel spacing differs from CHAN_WIDTH")
        summaries.append({"spw": index, "channels": channels,
                          "channel_width_hz": width, "bandwidth_hz": bandwidth,
                          "first_frequency_hz": frequencies[0],
                          "last_frequency_hz": frequencies[-1]})
    return {"date": int(date), "archive_averaging": averaging,
            "native_channels_per_spw": expected["native_channels"],
            "metadata_valid": True, "spectral_windows": summaries,
            "scope": "Serialized channel-layout check only; no averaging, rescaling, flags or calibration changed"}


def check_measurement_set(path, *, date, averaging=128, table_factory=None):
    """Read SPECTRAL_WINDOW in read-only mode and validate before calibration."""
    ms = Path(path)
    spectral_window = ms / "SPECTRAL_WINDOW"
    if ms.is_symlink() or spectral_window.is_symlink() or not spectral_window.is_dir():
        raise ValueError(f"A real Measurement Set SPECTRAL_WINDOW directory is required: {spectral_window}")
    if table_factory is None:
        from casatools import table
        table_factory = table
    tool = table_factory()
    opened = False
    try:
        tool.open(str(spectral_window), nomodify=True)
        opened = True
        records = []
        for row in range(tool.nrows()):
            record = {column: list(tool.getcell(column, row)) for column in CHANNEL_COLUMNS}
            record["NUM_CHAN"] = int(tool.getcell("NUM_CHAN", row))
            record["TOTAL_BANDWIDTH"] = float(tool.getcell("TOTAL_BANDWIDTH", row))
            records.append(record)
    finally:
        if opened:
            tool.close()
    result = validate_spectral_windows(records, date=date, averaging=averaging)
    result["measurement_set"] = str(ms.resolve())
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("measurement_set", type=Path)
    parser.add_argument("--date", type=int, choices=sorted(SUPPORTED_TRACKS), required=True)
    parser.add_argument("--averaging", type=int, default=128)
    args = parser.parse_args(argv)
    print(json.dumps(check_measurement_set(args.measurement_set, date=args.date,
                                          averaging=args.averaging), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
