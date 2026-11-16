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
"""Tests for the Perl oracle comparison."""

from __future__ import annotations

import os
import re
import shutil
from pathlib import Path

import pytest

import ensembl.vep.annotators  # noqa: F401  pylint: disable=unused-import

# ^ registers the annotators and the facts they read.
from ensembl.vep.annotators.factory import build_annotator
from ensembl.vep.bench.oracle import (
    CASE_NAMES,
    covered_cases,
    covered_through_subsumption,
    oracle_cases,
    oracle_transcript,
    perl_terms,
    perl_unavailable_reason,
    python_terms,
)

PERL_SCRIPT = Path(__file__).resolve().parents[4] / "bench" / "variation_effect.pl"
PERL_MODULES_ENV = "ENSEMBL_VEP_PERL_MODULES"
"""Set to ensembl-variation's modules directory to run the Perl side in the test suite."""


class TestCases:
    """The case list, and its alignment with the Perl script."""

    def test_the_transcript_builds_with_a_derived_translation(self) -> None:
        """The builder this replaced no longer ran; this one derives the translation like the corpus."""
        transcript = oracle_transcript()
        assert transcript.translation is not None
        assert transcript.cds_length == 69

    def test_the_cases_are_named_and_ordered_as_the_perl_script_builds_them(self) -> None:
        """Drift between the two lists would pair each Perl result with the wrong case."""
        perl_names = re.findall(r"name => '(\w+)'", PERL_SCRIPT.read_text(encoding="utf-8"))
        assert tuple(perl_names) == CASE_NAMES
        assert tuple(case.name for case in oracle_cases()) == CASE_NAMES

    def test_each_case_emits_the_term_its_position_calls_for(self) -> None:
        """What linear emits after subsumption.

        Both intronic cases lie inside ``intron.start + 2 .. intron.end - 2``, so ``intron_variant``
        comes out beside the splice term, as VEP reports it: neither term suppresses the other.
        """
        emitted = python_terms(build_annotator(annotator="consequence.linear"), oracle_cases())
        terms = dict(zip(CASE_NAMES, emitted, strict=True))
        assert terms["coding_exact"] == ("missense_variant",)
        assert terms["intron_primary"] == ("splice_polypyrimidine_tract_variant", "intron_variant")
        assert terms["intron_secondary"] == ("splice_donor_region_variant", "intron_variant")
        assert terms["upstream_near"] == ("upstream_gene_variant",)
        assert terms["downstream_far"] == ("downstream_gene_variant",)
        assert terms["no_hit"] == ()


class TestPerlSide:
    """Whether and how the Perl side runs."""

    def test_missing_pieces_are_reported_as_reasons(self, tmp_path: Path) -> None:
        """Each missing prerequisite says what is missing."""
        assert "--perl_modules" in (perl_unavailable_reason(None, PERL_SCRIPT) or "")
        assert "not found" in (perl_unavailable_reason(tmp_path / "absent", PERL_SCRIPT) or "")
        assert "not found" in (perl_unavailable_reason(tmp_path, tmp_path / "absent.pl") or "")

    def test_coverage_counts_cases_whose_perl_terms_we_cover(self) -> None:
        """A superset covers; a missing term does not."""
        assert covered_cases([("a",), ("b",), ()], [("a", "x"), ("c",), ()]) == 2

    def test_suppressed_terms_count_as_covered_through_a_chain(self) -> None:
        """A Perl term our emitted term suppresses, even two steps down, is not a disagreement."""
        subsumes = {
            "missense_variant": ("protein_altering_variant",),
            "protein_altering_variant": ("coding_sequence_variant",),
        }
        perl = [("coding_sequence_variant",), ("intron_variant",)]
        ours = [("missense_variant",), ("upstream_gene_variant",)]
        assert covered_cases(perl, ours) == 0
        assert covered_through_subsumption(perl, ours, subsumes) == 1

    @pytest.mark.skipif(
        shutil.which("perl") is None or not os.environ.get(PERL_MODULES_ENV),
        reason=f"needs perl and {PERL_MODULES_ENV}",
    )
    def test_the_perl_oracle_reports_one_term_list_per_case(self) -> None:
        """Run only where the Perl side is set up."""
        modules = Path(os.environ[PERL_MODULES_ENV])
        assert len(perl_terms(modules, PERL_SCRIPT, len(CASE_NAMES))) == len(CASE_NAMES)
