from __future__ import annotations

import json
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
NOTEBOOK_PATH = ROOT / "BH_Parameters.ipynb"
OUTPUT_PATH = ROOT / "analysis" / "CountingBHS_and_Spectra" / "BHS_params.csv"


def main() -> None:
    os.chdir(ROOT)
    nb = json.loads(NOTEBOOK_PATH.read_text())
    namespace = {"__name__": "__main__"}

    for cell in nb.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        source = "".join(cell.get("source", []))
        if not source.strip():
            continue
        exec(compile(source, str(NOTEBOOK_PATH), "exec"), namespace)

    if not OUTPUT_PATH.exists():
        raise FileNotFoundError(f"Expected BHS export was not created: {OUTPUT_PATH}")

    print(f"Refreshed {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
