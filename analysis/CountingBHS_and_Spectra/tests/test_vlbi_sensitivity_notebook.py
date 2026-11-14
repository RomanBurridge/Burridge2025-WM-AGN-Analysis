from __future__ import annotations

import ast
import json
import math
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
NOTEBOOK_PATH = REPOSITORY_ROOT / "VLBI_sensitivity.ipynb"


def astrometric_source() -> str:
    notebook = json.loads(NOTEBOOK_PATH.read_text(encoding="utf-8"))
    cells = [
        cell
        for cell in notebook["cells"]
        if cell.get("cell_type") == "code"
        and "#astrometric precision" in "".join(cell.get("source", []))
    ]
    if len(cells) != 1:
        raise AssertionError("expected exactly one astrometric-precision cell")
    return "".join(cells[0]["source"])


def scalar_assignments(source: str) -> dict[str, float]:
    assignments: dict[str, float] = {}
    for node in ast.walk(ast.parse(source, filename=NOTEBOOK_PATH.name)):
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        target = node.targets[0]
        if not isinstance(target, ast.Name):
            continue
        try:
            value = ast.literal_eval(node.value)
        except (ValueError, TypeError):
            continue
        if isinstance(value, (int, float)):
            assignments[target.id] = float(value)
    return assignments


class VlbiSensitivityNotebookTests(unittest.TestCase):
    def test_accepted_astrometric_precision_benchmark_is_current(self) -> None:
        source = astrometric_source()
        assignments = scalar_assignments(source)

        self.assertEqual(assignments["prec"], 0.075)
        self.assertEqual(assignments["precunc"], 0.007)
        self.assertNotIn("0.086", source)

        factor = (4.0 / math.pi) ** 0.25 / math.sqrt(8.0 * math.log(2.0))
        self.assertAlmostEqual(factor * 0.18 / assignments["prec"], 1.0826324679225185)
        self.assertAlmostEqual(factor * 0.69 / assignments["prec"], 4.150091127036321)

    def test_notebook_is_unexecuted_source_only(self) -> None:
        notebook = json.loads(NOTEBOOK_PATH.read_text(encoding="utf-8"))
        self.assertEqual(notebook["nbformat"], 4)
        for cell in notebook["cells"]:
            if cell.get("cell_type") != "code":
                continue
            self.assertIsNone(cell.get("execution_count"))
            self.assertEqual(cell.get("outputs"), [])


if __name__ == "__main__":
    unittest.main()
