"""Compatibility entry point for the deterministic BHS table builder.

Historically this script executed every code cell in ``BH_Parameters.ipynb``.
That made stale exploratory cells part of the publication pipeline.  The
audited workflow now calls one ordinary Python builder with explicit input
tables; the notebook is only a readable diagnostic companion.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BUILDER_PATH = ROOT / "analysis" / "CountingBHS_and_Spectra" / "build_bhs_params.py"
OUTPUT_PATH = ROOT / "analysis" / "CountingBHS_and_Spectra" / "BHS_params.csv"


def main() -> None:
    subprocess.run(
        [sys.executable, str(BUILDER_PATH), "--require-final"],
        cwd=ROOT,
        check=True,
    )

    if not OUTPUT_PATH.exists():
        raise FileNotFoundError(f"Expected BHS export was not created: {OUTPUT_PATH}")

    print(f"Refreshed {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
