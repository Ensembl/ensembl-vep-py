# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Generate the derivation graphs shown in the "Explaining a consequence call" tutorial.

Runs ``vep explain --graph`` on the tutorial's overlap over the example data, writing the full
HTML page, then copies the SVG of a few rules out of that page so the tutorial can show them as
images. Run it from the repository root.

Usage: ``python docs/_scripts/explaindoc.py [OUTPUT_DIR]``
"""

from __future__ import annotations

import argparse
import re
import sys
import tempfile
from collections.abc import Sequence
from pathlib import Path

from ensembl.vep.cli import main as vep_main

EXAMPLES = Path("examples/homo_sapiens")
VARIANT = "22:17188416:rs2231495"
TRANSCRIPT = "ENST00000262607"
PAGE_NAME = "derivation.html"
# One rule per outcome: emitted, fired but pruned, and not fired.
FIGURE_TERMS = ("missense_variant", "protein_altering_variant", "synonymous_variant")

CONFIG = f"""\
[source]
input_file = "{EXAMPLES / "homo_sapiens_GRCh38.vcf"}"
gtf = "{EXAMPLES / "genes.gtf.bgz"}"
fasta = "{EXAMPLES / "unmasked.fa.bgz"}"
"""


def _rule_svg(page: str, term: str) -> str:
    """Return the SVG drawn for one rule in a derivation page."""
    match = re.search(rf"<h2>{re.escape(term)}</h2>.*?(<svg\b.*?</svg>)", page, re.DOTALL)
    if match is None:
        raise LookupError(f"the derivation page has no graph for {term}")
    return match.group(1)


def main(argv: Sequence[str] | None = None) -> int:
    """Write the derivation page and the rule figures."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("output_dir", nargs="?", default="docs/auto/static/explain", type=Path)
    args = parser.parse_args(argv)

    output_dir: Path = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    page_path = output_dir / PAGE_NAME
    page_path.unlink(missing_ok=True)
    with tempfile.TemporaryDirectory() as work_dir:
        config_path = Path(work_dir) / "explain.toml"
        config_path.write_text(CONFIG, encoding="utf-8")
        exit_code = vep_main(
            [
                "explain",
                "--config",
                str(config_path),
                "--variant",
                VARIANT,
                "--transcript",
                TRANSCRIPT,
                "--graph",
                "--output_file",
                str(page_path),
            ]
        )
    if exit_code:
        return int(exit_code)

    page = page_path.read_text(encoding="utf-8")
    for term in FIGURE_TERMS:
        (output_dir / f"{term}.svg").write_text(_rule_svg(page, term) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
