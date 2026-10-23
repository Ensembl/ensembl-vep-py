# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""The layering rule that keeps the render package reusable."""

from __future__ import annotations

import ast
from pathlib import Path

RENDER_PACKAGE = Path(__file__).resolve().parents[2] / "ensembl" / "vep" / "render"


def _imported_modules(path: Path) -> set[str]:
    """Return the module names one file imports, at module scope or inside a function."""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    imported: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                imported.add(alias.name)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module)
    return imported


class TestRenderIsDomainFree:
    """render/ must import nothing else in ensembl.vep."""

    def test_no_module_reaches_outside_the_render_package(self) -> None:
        """It is what lets any subsystem draw without pulling the domain in behind it.

        The builders that read facts, rules and traces live beside those instead, and hand a plain
        Graph across.
        """
        offenders: dict[str, set[str]] = {}
        for path in RENDER_PACKAGE.rglob("*.py"):
            outside = set()
            for module in _imported_modules(path):
                if module.startswith("ensembl.vep.") and not module.startswith("ensembl.vep.render"):
                    outside.add(module)
            if outside:
                offenders[str(path.relative_to(RENDER_PACKAGE))] = outside
        assert not offenders, f"render/ imports the rest of the package: {offenders}"

    def test_the_package_imports_without_the_rest_of_ensembl_vep(self) -> None:
        """A stronger check than reading the imports: it must actually load on its own."""
        # pylint: disable=import-outside-toplevel
        import subprocess
        import sys

        result = subprocess.run(
            [sys.executable, "-c", "import ensembl.vep.render, ensembl.vep.render.graph"],
            capture_output=True,
            check=False,
        )
        assert result.returncode == 0, result.stderr.decode()
