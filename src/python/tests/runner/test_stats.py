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
"""Tests for the stats runner."""

from __future__ import annotations

from pathlib import Path

import pytest

import ensembl.vep.annotators  # noqa: F401  pylint: disable=unused-import

# ^ registers the annotators the profile runs and the facts the catalogue compiles against.
from ensembl.vep.cli import main
from ensembl.vep.runner.base import ExitCode

ANNOTATED_VCF = """##fileformat=VCFv4.2
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	rs1	A	T	.	.	CSQ=T|missense_variant|MODERATE
1	200	rs2	AC	A	.	.	CSQ=-|frameshift_variant|HIGH
2	300	rs3	G	A	.	.	CSQ=A|intron_variant|MODIFIER
"""


@pytest.fixture(name="annotated")
def fixture_annotated(tmp_path: Path) -> Path:
    """A small annotated VCF on disk."""
    path = tmp_path / "annotated.vcf"
    path.write_text(ANNOTATED_VCF, encoding="utf-8")
    return path


class TestSummary:
    """``vep stats --input_file``."""

    def test_it_summarises_an_annotated_file(
        self, annotated: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """The headline counts and the charts come out of the file alone."""
        assert main(["stats", "--input_file", str(annotated), "--output_format", "text"]) == ExitCode.OK
        out = capsys.readouterr().out
        assert "3 variants" in out
        assert "missense_variant" in out
        assert "Coding consequences" in out

    def test_html_and_json_render(self, annotated: Path, capsys: pytest.CaptureFixture[str]) -> None:
        """Every declared format works over the same summary."""
        assert main(["stats", "--input_file", str(annotated), "--output_format", "html"]) == ExitCode.OK
        assert "<svg" in capsys.readouterr().out
        assert main(["stats", "--input_file", str(annotated), "--output_format", "json"]) == ExitCode.OK
        assert '"Run summary"' in capsys.readouterr().out

    def test_it_can_write_to_a_file(self, annotated: Path, tmp_path: Path) -> None:
        """--output_file takes the report instead of stdout."""
        destination = tmp_path / "summary.html"
        assert (
            main(["stats", "--input_file", str(annotated), "--output_file", str(destination)]) == ExitCode.OK
        )
        assert "Run summary" in destination.read_text(encoding="utf-8")

    def test_a_missing_file_is_a_failure(self, tmp_path: Path) -> None:
        """The command ran and could not do its job."""
        assert main(["stats", "--input_file", str(tmp_path / "absent.vcf")]) == ExitCode.FAILURE


class TestProfile:
    """``vep stats --profile``."""

    def test_it_profiles_a_corpus(self, capsys: pytest.CaptureFixture[str]) -> None:
        """A fact table and the cost-table provenance come out."""
        assert (
            main(
                [
                    "stats",
                    "--profile",
                    "--kind",
                    "synthetic.default",
                    "--fuzz",
                    "0",
                    "--output_format",
                    "text",
                ]
            )
            == ExitCode.OK
        )
        out = capsys.readouterr().out
        assert "Fact profile: consequence.columnar-staged" in out
        assert "declared tiers" in out
        assert "overlap_cds" in out


class TestUsage:
    """What has to be given."""

    def test_nothing_is_a_usage_error(self) -> None:
        """Neither an input nor --profile means there is nothing to report on."""
        assert main(["stats"]) == ExitCode.USAGE

    def test_both_is_a_usage_error(self, annotated: Path) -> None:
        """The two modes read different things."""
        assert main(["stats", "--profile", "--input_file", str(annotated)]) == ExitCode.USAGE
