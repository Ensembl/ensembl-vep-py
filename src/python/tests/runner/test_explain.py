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
"""Tests for the explain runner."""

from __future__ import annotations

from pathlib import Path

import pytest

import ensembl.vep.annotators  # noqa: F401  pylint: disable=unused-import

# ^ registers the annotators explain runs and the facts they read.
from ensembl.vep.annotators.consequence.rules.catalogue import load_consequence_definitions
from ensembl.vep.cli import main
from ensembl.vep.corpus import build_corpus
from ensembl.vep.corpus.describe import describe_variant
from ensembl.vep.corpus.spec import CorpusSpec
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.runner.base import ExitCode

SMALL_CORPUS = ["--kind", "synthetic.default", "--fuzz", "5"]
"""The smallest corpus that still contains every boundary overlap, so these stay quick."""

CODING_OVERLAP = "ENST0001_cds3_del1"
"""A coding deletion: several rules fire and subsumption decides between them."""


def _corpus_overlap(name: str) -> TranscriptVariationAllele:
    """Return the corpus overlap yielded under a name, for a test that needs its locus."""
    corpus = build_corpus(CorpusSpec(kind="synthetic.default", fuzz=5))
    return next(overlap for overlap_name, overlap in corpus.named_overlaps if overlap_name == name)


def _explain(*arguments: str) -> ExitCode:
    """Run ``vep explain`` with the small corpus already selected."""
    return main(["explain", *SMALL_CORPUS, *arguments])


_TRANSCRIPT_GFF3 = """\
chr1\tensembl\tgene\t1000\t5000\t.\t+\t.\tID=gene:G1;biotype=protein_coding
chr1\tensembl\tmRNA\t1000\t5000\t.\t+\t.\tID=transcript:ENST001;Parent=gene:G1;biotype=protein_coding
chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tID=exon:E1;Parent=transcript:ENST001
chr1\tensembl\tCDS\t1100\t2000\t.\t+\t0\tID=cds:P1;Parent=transcript:ENST001
"""


class TestSelection:
    """Which overlap gets explained, and what happens when that is not answerable."""

    def test_listing_shows_the_selectable_overlaps(self, capsys: pytest.CaptureFixture[str]) -> None:
        """--list is how a caller finds a name to pass to --name."""
        assert _explain("--list") == ExitCode.OK
        assert CODING_OVERLAP in capsys.readouterr().out

    def test_selecting_by_name_explains_that_overlap(self, capsys: pytest.CaptureFixture[str]) -> None:
        """The rendered derivation names the overlap it explains."""
        assert _explain("--name", CODING_OVERLAP) == ExitCode.OK
        assert "ENST0001" in capsys.readouterr().out

    def test_selecting_by_index_explains_the_overlap_at_that_position(self) -> None:
        """--index selects from the source's own order, for a corpus whose names repeat."""
        assert _explain("--index", "0") == ExitCode.OK

    def test_an_unknown_name_is_a_usage_error(self) -> None:
        """Nothing was explained, and the caller mistyped rather than the command failing."""
        assert _explain("--name", "no_such_overlap") == ExitCode.USAGE

    def test_selecting_nothing_is_a_usage_error(self) -> None:
        """An overlap has to be named: there is no sensible default one to explain."""
        assert _explain() == ExitCode.USAGE

    def test_selecting_by_variant_explains_that_overlap(self, capsys: pytest.CaptureFixture[str]) -> None:
        """--variant is the selector a real input will use, and it works on the corpus too."""
        overlap = _corpus_overlap(CODING_OVERLAP)
        locus, alleles = describe_variant(overlap.overlapping_feature)
        assert _explain("--variant", f"{locus}:{alleles}", "--transcript", "ENST0001") == ExitCode.OK
        assert "transcript    ENST0001" in capsys.readouterr().out

    def test_a_malformed_variant_is_a_usage_error(self) -> None:
        """The selector's shape is CHR:POS:REF>ALT; anything else is the caller's mistake."""
        assert _explain("--variant", "1-12345-A-T") == ExitCode.USAGE

    def test_an_unknown_variant_is_a_usage_error(self) -> None:
        """A well-formed selector that matches nothing is reported, not explained."""
        assert _explain("--variant", "chrSynthetic:1:A>T", "--transcript", "ENST0009") == ExitCode.USAGE

    def test_an_index_past_the_end_is_a_usage_error(self) -> None:
        """Wrapping round or clamping would explain the wrong overlap without saying so."""
        assert _explain("--index", str(10**6)) == ExitCode.USAGE

    def test_an_ambiguous_name_is_a_usage_error_naming_the_features(
        self, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """A name identifies the variant, so the heavy locus pairs one name with many features."""
        corpus = build_corpus(CorpusSpec(kind="synthetic.heavy", step=600))
        ambiguous_name, _overlap = next(iter(corpus.named_overlaps))

        code = main(
            [
                "explain",
                "--kind",
                "synthetic.heavy",
                "--step",
                "600",
                "--name",
                ambiguous_name,
            ]
        )

        assert code == ExitCode.USAGE
        assert "disambiguate" in capsys.readouterr().err

    def test_naming_a_transcript_disambiguates(self) -> None:
        """Which is why selecting one overlap has to be able to name a feature too."""
        corpus = build_corpus(CorpusSpec(kind="synthetic.heavy", step=600))
        ambiguous_name, first = next(iter(corpus.named_overlaps))

        code = main(
            [
                "explain",
                "--kind",
                "synthetic.heavy",
                "--step",
                "600",
                "--name",
                ambiguous_name,
                "--transcript",
                first.feature.stable_id,
            ]
        )

        assert code == ExitCode.OK

    def test_a_variant_and_a_name_must_both_match(self) -> None:
        """Giving both is a stricter question, not a choice between them."""
        assert _explain("--name", CODING_OVERLAP, "--variant", "chrSynthetic:1:A>T") == ExitCode.USAGE

    def test_a_variant_file_without_transcripts_is_a_usage_error(self) -> None:
        """A variant file is only half the input: pairing it needs --gff or --gtf as well. Omitting
        them is the caller's mistake, the same kind as a malformed selector."""
        assert main(["explain", "--input_file", "variants.vcf", "--name", CODING_OVERLAP]) == (ExitCode.USAGE)

    def test_a_real_input_file_can_be_explained(self, tmp_path: Path) -> None:
        """Explaining a variant read from a file, rather than generated.

        This works only because the runner forwards --gff/--fasta alongside --input_file, and
        because an empty allele is spelled the same way on both sides of a selector comparison.
        Both are easy to regress without anything noticing, so the mode is pinned here.
        """
        variants = tmp_path / "variants.vcf"
        variants.write_text(
            "##fileformat=VCFv4.3\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t1500\trs_coding\tA\tC\t.\tPASS\t.\n",
            encoding="utf-8",
        )
        transcripts = tmp_path / "transcripts.gff3"
        transcripts.write_text(
            "chr1\tensembl\tgene\t1000\t5000\t.\t+\t.\tID=gene:G1;biotype=protein_coding\n"
            "chr1\tensembl\tmRNA\t1000\t5000\t.\t+\t.\tID=transcript:ENST001;Parent=gene:G1;"
            "biotype=protein_coding\n"
            "chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tID=exon:E1;Parent=transcript:ENST001\n"
            "chr1\tensembl\tCDS\t1100\t2000\t.\t+\t0\tID=cds:P1;Parent=transcript:ENST001\n",
            encoding="utf-8",
        )
        fasta = tmp_path / "ref.fa"
        fasta.write_text(">chr1\n" + "ACGT" * 2000 + "\n", encoding="utf-8")

        code = main(
            [
                "explain",
                "--input_file",
                str(variants),
                "--gff",
                str(transcripts),
                "--fasta",
                str(fasta),
                "--name",
                "rs_coding",
            ]
        )

        assert code == ExitCode.OK

    @pytest.mark.parametrize(
        "selector",
        [
            pytest.param(["--variant", "chr1:1500:A>C", "--allele_num", "2"], id="alleles-already-named"),
            pytest.param(["--name", CODING_OVERLAP, "--allele_num", "2"], id="no-variant"),
        ],
    )
    def test_a_misplaced_allele_number_is_a_usage_error(self, selector: list[str]) -> None:
        """--allele_num only chooses among the alleles of a CHR:POS:ID variant."""
        assert _explain(*selector) == ExitCode.USAGE

    def test_a_real_input_variant_can_be_selected_by_id(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """CHR:POS:ID reaches the same variant by its VCF ID and position."""
        variants = tmp_path / "variants.vcf"
        variants.write_text(
            "##fileformat=VCFv4.3\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t1500\trs_coding\tA\tC,G\t.\tPASS\t.\n",
            encoding="utf-8",
        )
        transcripts = tmp_path / "transcripts.gff3"
        transcripts.write_text(_TRANSCRIPT_GFF3, encoding="utf-8")
        fasta = tmp_path / "ref.fa"
        fasta.write_text(">chr1\n" + "ACGT" * 2000 + "\n", encoding="utf-8")

        code = main(
            [
                "explain",
                "--input_file",
                str(variants),
                "--gff",
                str(transcripts),
                "--fasta",
                str(fasta),
                "--variant",
                "chr1:1500:rs_coding",
                "--allele_num",
                "2",
            ]
        )

        assert code == ExitCode.OK
        assert "emitted:" in capsys.readouterr().out

    def test_an_unreadable_variant_file_is_reported_not_raised(self, tmp_path: Path) -> None:
        """The readers raise rather than calling sys.exit, which leaves the exit code to the runner.

        A file source reads lazily, so a missing input surfaces well after the source is built. A
        traceback escaping the CLI here would be the regression this guards.
        """
        transcripts = tmp_path / "transcripts.gff3"
        transcripts.write_text(
            "chr1\tensembl\tgene\t1000\t5000\t.\t+\t.\tID=gene:G1\n"
            "chr1\tensembl\tmRNA\t1000\t5000\t.\t+\t.\tID=transcript:ENST001;Parent=gene:G1\n"
            "chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tID=exon:E1;Parent=transcript:ENST001\n",
            encoding="utf-8",
        )
        arguments = [
            "explain",
            "--input_file",
            str(tmp_path / "nope.vcf"),
            "--gff",
            str(transcripts),
            "--name",
            CODING_OVERLAP,
        ]

        assert main(arguments) == ExitCode.USAGE


class TestExplanation:
    """What a rendered derivation contains."""

    def test_the_derivation_names_facts_rules_and_the_result(
        self, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """The whole chain is the point: which facts, which verdicts, what was emitted."""
        assert _explain("--name", CODING_OVERLAP) == ExitCode.OK
        rendered = capsys.readouterr().out
        assert "facts:" in rendered
        assert "rules:" in rendered
        assert "pruned:" in rendered
        assert "emitted:" in rendered

    def test_the_derivation_opens_with_a_summary_of_the_overlap(
        self, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """What was annotated comes before how: variant, transcript, codons, the whole CDS."""
        assert _explain("--name", CODING_OVERLAP) == ExitCode.OK
        rendered = capsys.readouterr().out
        assert rendered.startswith("summary:\n")
        assert rendered.index("summary:") < rendered.index("  codons") < rendered.index("facts:")
        assert "  full CDS" in rendered

    def test_a_comparison_opens_with_the_summary_too(self, capsys: pytest.CaptureFixture[str]) -> None:
        """The overlap is the same for both annotators, so it is summarised once, first."""
        exit_code = _explain(
            "--name", CODING_OVERLAP, "--compare", "consequence.linear,consequence.columnar-staged"
        )
        assert exit_code == ExitCode.OK
        rendered = capsys.readouterr().out
        assert rendered.index("summary:") < rendered.index("emitted terms AGREE")

    def test_a_graph_source_carries_no_summary(self, capsys: pytest.CaptureFixture[str]) -> None:
        """Mermaid and dot output stay parseable graph source."""
        assert _explain("--name", CODING_OVERLAP, "--graph", "--output_format", "mermaid") == ExitCode.OK
        assert "summary:" not in capsys.readouterr().out

    def test_a_chosen_annotator_is_named_in_the_derivation(self, capsys: pytest.CaptureFixture[str]) -> None:
        """Any pure-Python strategy can be explained, not only the default one."""
        assert _explain("--name", CODING_OVERLAP, "--annotator", "consequence.staged") == ExitCode.OK
        assert "annotator: consequence.staged" in capsys.readouterr().out

    def test_the_derivation_can_be_written_to_a_file(self, tmp_path: Path) -> None:
        """--output_file is what makes a derivation reviewable outside the terminal."""
        destination = tmp_path / "derivation.txt"
        assert _explain("--name", CODING_OVERLAP, "--output_file", str(destination)) == ExitCode.OK
        assert "emitted:" in destination.read_text(encoding="utf-8")

    def test_a_graph_is_rendered_instead_of_the_text(self, capsys: pytest.CaptureFixture[str]) -> None:
        """--graph draws the same derivation, for the cases a table of rules does not show."""
        assert _explain("--name", CODING_OVERLAP, "--graph", "--output_format", "mermaid") == ExitCode.OK
        assert "flowchart" in capsys.readouterr().out

    def test_svg_is_a_linked_html_explanation(self, capsys: pytest.CaptureFixture[str]) -> None:
        """One linked term tree at a time remains readable while keeping every outcome available."""
        assert _explain("--name", CODING_OVERLAP, "--graph", "--output_format", "svg") == ExitCode.OK
        output = capsys.readouterr().out
        document = output[output.index("<!DOCTYPE html>") :]
        assert 'href="#stop_gained"' in document
        assert 'id="stop_gained"' in document
        assert 'href="#transcript_ablation"' in document
        assert 'href="#intergenic_variant"' in document
        assert 'href="#summary"' in document
        assert document.count('<section class="vep-section"') == len(load_consequence_definitions()) + 1

    def test_term_links_are_grouped_by_outcome(self, capsys: pytest.CaptureFixture[str]) -> None:
        """The result comes first, followed by suppressed matches and rules that were blocked."""
        assert _explain("--name", CODING_OVERLAP, "--graph", "--output_format", "svg") == ExitCode.OK
        output = capsys.readouterr().out
        emitted = output.index('<li class="vep-menu-group">Emitted</li>')
        pruned = output.index('<li class="vep-menu-group">Pruned</li>')
        not_fired = output.index('<li class="vep-menu-group">Not fired</li>')
        assert emitted < pruned < not_fired
        assert emitted < output.index('href="#frameshift_variant"') < pruned
        assert pruned < output.index('href="#protein_altering_variant"') < not_fired
        assert not_fired < output.index('href="#stop_gained"')

    def test_condition_labels_leave_truth_to_colour(self, capsys: pytest.CaptureFixture[str]) -> None:
        """The rule remains readable as rule syntax rather than a list of value annotations."""
        assert _explain("--name", CODING_OVERLAP, "--graph", "--output_format", "svg") == ExitCode.OK
        output = capsys.readouterr().out
        assert "NOT is_structural" in output
        assert "is_structural = false" not in output
        assert "✓" not in output

    def test_svg_is_written_as_an_html_document(self, tmp_path: Path) -> None:
        """The links and all term sections remain usable outside the terminal."""
        destination = tmp_path / "explanation.html"
        exit_code = _explain(
            "--name",
            CODING_OVERLAP,
            "--graph",
            "--output_format",
            "svg",
            "--output_file",
            str(destination),
        )
        assert exit_code == ExitCode.OK
        assert destination.read_text(encoding="utf-8").startswith("<!DOCTYPE html>")


class TestCompare:
    """Diffing two annotators over one overlap."""

    def test_two_agreeing_annotators_exit_ok(self, capsys: pytest.CaptureFixture[str]) -> None:
        """Every annotator is byte-identical by conformance, so agreement is the normal outcome."""
        exit_code = _explain(
            "--name", CODING_OVERLAP, "--compare", "consequence.linear,consequence.columnar-staged"
        )
        assert exit_code == ExitCode.OK
        assert "emitted terms AGREE" in capsys.readouterr().out

    def test_a_comparison_draws_both_derivations_when_they_differ(self, tmp_path: Path) -> None:
        """Two strategies observe different facts, so each side's proofs are drawn separately.

        The menu groups carry the annotator's name because the same term appears on both sides.
        """
        page = tmp_path / "comparison.html"
        exit_code = _explain(
            "--name",
            CODING_OVERLAP,
            "--compare",
            "consequence.linear,consequence.staged",
            "--graph",
            "--output_file",
            str(page),
        )

        assert exit_code == ExitCode.OK
        rendered = page.read_text(encoding="utf-8")
        assert "their derivations differ" in rendered
        assert "consequence.linear \u00b7 Emitted" in rendered
        assert "consequence.staged \u00b7 Emitted" in rendered

    def test_a_comparison_draws_one_set_when_the_derivations_are_identical(self, tmp_path: Path) -> None:
        """A second identical copy of the same drawing says nothing, so one set stands for both."""
        page = tmp_path / "comparison.html"
        exit_code = _explain(
            "--name",
            CODING_OVERLAP,
            "--compare",
            "consequence.staged,consequence.columnar-staged",
            "--graph",
            "--output_file",
            str(page),
        )

        assert exit_code == ExitCode.OK
        rendered = page.read_text(encoding="utf-8")
        assert "their derivations are identical" in rendered
        assert "consequence.staged \u00b7 Emitted" not in rendered

    def test_a_comparison_page_carries_the_report_it_draws(self, tmp_path: Path) -> None:
        """The verdict is the point of a comparison, so the drawing does not replace it."""
        page = tmp_path / "comparison.html"
        _explain(
            "--name",
            CODING_OVERLAP,
            "--compare",
            "consequence.linear,consequence.staged",
            "--graph",
            "--output_file",
            str(page),
        )

        rendered = page.read_text(encoding="utf-8")
        assert 'id="comparison"' in rendered
        assert "emitted terms AGREE" in rendered

    def test_a_comparison_graph_source_emits_every_set_and_logs_the_verdict(
        self, capsys: pytest.CaptureFixture[str], caplog: pytest.LogCaptureFixture
    ) -> None:
        """A source format must stay parseable, so the report goes to the log, not to stdout."""
        with caplog.at_level("INFO"):
            exit_code = _explain(
                "--name",
                CODING_OVERLAP,
                "--compare",
                "consequence.linear,consequence.staged",
                "--graph",
                "--output_format",
                "dot",
            )

        assert exit_code == ExitCode.OK
        rendered = capsys.readouterr().out
        assert rendered.lstrip().startswith("digraph")
        assert "emitted terms AGREE" not in rendered
        assert "emitted terms AGREE" in caplog.text

    def test_one_annotator_is_a_usage_error(self) -> None:
        """A comparison needs two sides."""
        assert _explain("--name", CODING_OVERLAP, "--compare", "consequence.linear") == ExitCode.USAGE

    def test_a_diverging_pair_exits_failure(self) -> None:
        """A comparison that finds the two annotators disagree is a failure the caller can act on.

        No registered pair disagrees -- that is the conformance invariant -- so this pins the
        agreeing case's exit code instead, which is the half that can be asserted honestly here.
        A native rust.* annotator, the one side that cannot trace, is not registered in this build,
        so the diff_terms fallback is covered by tests/trace/test_explain.py rather than from here.
        """
        exit_code = _explain("--name", CODING_OVERLAP, "--compare", "consequence.staged,consequence.columnar")
        assert exit_code == ExitCode.OK
