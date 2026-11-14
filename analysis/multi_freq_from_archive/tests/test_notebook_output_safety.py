from __future__ import annotations

import ast
import json
import unittest
import warnings
from collections import defaultdict
from pathlib import Path


ANALYSIS_DIR = Path(__file__).resolve().parents[1]
NOTEBOOK_OUTPUT_ROOTS = {
    "Beam_Date_Flux_Graphs.ipynb": "SpectraGraphs_notebook_diagnostic",
    "Variability_Search.ipynb": "VariabilityGraphs_notebook_diagnostic",
    "get_contams.ipynb": "ContaminationGraphs_notebook_diagnostic",
    "flux_size_graph.ipynb": "FluxSizeGraphs_notebook_diagnostic",
}
MUTATING_CALLS = {
    "os.makedirs",
    "os.mkdir",
    "os.remove",
    "os.rmdir",
    "os.unlink",
    "plt.savefig",
    "pylab.savefig",
    "shutil.rmtree",
}


def _call_name(node: ast.AST) -> str | None:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        parent = _call_name(node.value)
        return f"{parent}.{node.attr}" if parent else node.attr
    return None


def _write_mode(call: ast.Call) -> str | None:
    mode: ast.AST | None = call.args[1] if len(call.args) > 1 else None
    for keyword in call.keywords:
        if keyword.arg == "mode":
            mode = keyword.value
    if isinstance(mode, ast.Constant) and isinstance(mode.value, str):
        return mode.value
    return None


def _assignments(tree: ast.AST) -> dict[str, list[ast.AST]]:
    values: dict[str, list[ast.AST]] = defaultdict(list)
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name):
                    values[target.id].append(node.value)
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            if node.value is not None:
                values[node.target.id].append(node.value)
    return values


def _literal_prefix(node: ast.AST) -> str:
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    if isinstance(node, ast.JoinedStr):
        prefix = ""
        for value in node.values:
            if not isinstance(value, ast.Constant) or not isinstance(value.value, str):
                break
            prefix += value.value
        return prefix
    return ""


def _is_diagnostic_target(
    node: ast.AST,
    diagnostic_root: str,
    assignments: dict[str, list[ast.AST]],
    resolving: frozenset[str] = frozenset(),
) -> bool:
    prefix = _literal_prefix(node)
    if prefix:
        return prefix == diagnostic_root or prefix.startswith(diagnostic_root + "/")

    if isinstance(node, ast.Name):
        if node.id in resolving:
            return False
        values = assignments.get(node.id, [])
        return bool(values) and all(
            _is_diagnostic_target(
                value,
                diagnostic_root,
                assignments,
                resolving | {node.id},
            )
            for value in values
        )

    # Derived names such as ``contam_table.replace('.tex', '_new.tex')``
    # retain the safety of the rooted input path.
    if (
        isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr == "replace"
    ):
        return _is_diagnostic_target(
            node.func.value,
            diagnostic_root,
            assignments,
            resolving,
        )

    return False


class NotebookOutputSafetyTests(unittest.TestCase):
    def test_notebooks_are_valid_json_and_mutate_only_diagnostic_roots(self) -> None:
        for notebook_name, diagnostic_root in NOTEBOOK_OUTPUT_ROOTS.items():
            with self.subTest(notebook=notebook_name):
                notebook = json.loads(
                    (ANALYSIS_DIR / notebook_name).read_text(encoding="utf-8")
                )
                source = "\n".join(
                    "".join(cell.get("source", []))
                    for cell in notebook["cells"]
                    if cell.get("cell_type") == "code"
                )
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", SyntaxWarning)
                    tree = ast.parse(source, filename=notebook_name)
                assignments = _assignments(tree)
                checked_targets = 0

                for node in ast.walk(tree):
                    if not isinstance(node, ast.Call):
                        continue
                    call_name = _call_name(node.func)
                    targets: list[ast.AST] = []
                    if call_name in MUTATING_CALLS and node.args:
                        targets.append(node.args[0])
                    elif call_name == "open" and node.args:
                        mode = _write_mode(node)
                        if mode and any(flag in mode for flag in "wax+"):
                            targets.append(node.args[0])

                    targets.extend(
                        keyword.value
                        for keyword in node.keywords
                        if keyword.arg == "output_csv_path"
                    )

                    for target in targets:
                        checked_targets += 1
                        rendered = ast.unparse(target)
                        self.assertTrue(
                            _is_diagnostic_target(
                                target,
                                diagnostic_root,
                                assignments,
                            ),
                            f"{notebook_name}:{node.lineno} writes outside "
                            f"{diagnostic_root}: {rendered}",
                        )

                self.assertGreater(checked_targets, 0, notebook_name)
                self.assertIn("HISTORICAL EXPLORATORY NOTEBOOK", source)
                self.assertIn(diagnostic_root, source)


if __name__ == "__main__":
    unittest.main()
