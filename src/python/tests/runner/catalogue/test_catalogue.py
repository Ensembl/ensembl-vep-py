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
"""Tests for the catalogue runners."""

from __future__ import annotations

import json
import re
import shutil
from pathlib import Path

import pytest

import ensembl.vep.annotators  # noqa: F401  pylint: disable=unused-import

# ^ registers the facts and annotators the catalogue verbs read.
from ensembl.vep.annotators.consequence.rules.catalogue import catalogue_path
from ensembl.vep.cli import main
from ensembl.vep.runner.base import ExitCode


@pytest.fixture(name="catalogue_copy")
def fixture_catalogue_copy(tmp_path: Path) -> Path:
    """A writable copy of the packaged catalogue."""
    destination = tmp_path / "catalogue.json"
    shutil.copy(catalogue_path(), destination)
    return destination


def _terms(path: Path) -> set[str]:
    """Return the term names a catalogue file declares."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    return {entry["term"] for entry in payload["consequences"]}


class TestValidate:
    """The build gate."""

    def test_the_packaged_catalogue_passes(self) -> None:
        """A non-zero exit here is a real failure, not a known-bad baseline."""
        assert main(["catalogue", "validate"]) == ExitCode.OK

    def test_findings_are_reported(self, caplog: pytest.LogCaptureFixture) -> None:
        """The rank inversions are intentional and must stay visible as warnings."""
        with caplog.at_level("WARNING"):
            main(["catalogue", "validate"])
        assert "rank inversion" in caplog.text


class TestShowTables:
    """Listing terms or facts as a table."""

    def test_terms_are_the_default_subject(self, capsys: pytest.CaptureFixture) -> None:
        """Including the two with no ``when``, which the annotators skip."""
        assert main(["catalogue", "show", "--output_format", "text"]) == ExitCode.OK
        assert "intron_variant" in capsys.readouterr().out

    def test_the_table_on_stdout_carries_no_log_line(self, capsys: pytest.CaptureFixture) -> None:
        """Logs go to stderr, so ``vep catalogue show > terms.md`` redirects the table alone."""
        assert main(["catalogue", "show"]) == ExitCode.OK
        captured = capsys.readouterr()
        first_line = captured.out.splitlines()[0]
        assert first_line.startswith("|")
        assert "[INFO]" not in captured.out
        assert "[INFO]" in captured.err

    def test_facts_can_be_listed_instead(self, capsys: pytest.CaptureFixture) -> None:
        """The registry behind the terms, with each fact's declared tier."""
        assert main(["catalogue", "show", "--subject", "facts", "--output_format", "markdown"]) == ExitCode.OK
        output = capsys.readouterr().out
        assert "overlap_cds" in output
        assert "MAPPING" in output

    @pytest.mark.parametrize("table_format", ["markdown", "confluence", "text", "html"])
    def test_every_table_format_renders(self, table_format: str, capsys: pytest.CaptureFixture) -> None:
        """One renderer serves them all, so a format cannot be silently unsupported."""
        exit_code = main(["catalogue", "show", "--subject", "facts", "--output_format", table_format])
        assert exit_code == ExitCode.OK
        assert capsys.readouterr().out.strip()

    def test_a_graph_format_needs_the_graph_flag(self, caplog: pytest.LogCaptureFixture) -> None:
        """One vocabulary covers both kinds, so the combination is what is checked."""
        with caplog.at_level("ERROR"):
            assert main(["catalogue", "show", "--output_format", "svg"]) == ExitCode.USAGE
        assert "pass --graph" in caplog.text


class TestShowSorting:
    """Row order, and what happens when a key does not apply."""

    def test_terms_sort_by_name(self, capsys: pytest.CaptureFixture) -> None:
        """The order that used to be accepted for terms but silently ignored."""
        assert main(["catalogue", "show", "--output_format", "text", "--sort", "name"]) == ExitCode.OK
        terms = _column(capsys.readouterr().out, 1)
        assert terms == sorted(terms)

    def test_terms_sort_by_impact(self, capsys: pytest.CaptureFixture) -> None:
        """Impact is a term's property, so it is offered for terms."""
        assert main(["catalogue", "show", "--output_format", "text", "--sort", "impact"]) == ExitCode.OK
        impacts = _column(capsys.readouterr().out, 3)
        assert impacts == sorted(impacts)

    def test_facts_sort_by_name(self, capsys: pytest.CaptureFixture) -> None:
        """``name`` is the one key that applies to both subjects."""
        exit_code = main(
            ["catalogue", "show", "--subject", "facts", "--output_format", "text", "--sort", "name"]
        )
        assert exit_code == ExitCode.OK
        names = _column(capsys.readouterr().out, 0)
        assert names == sorted(names)

    def test_facts_sort_by_cost_tier(self, capsys: pytest.CaptureFixture) -> None:
        """Cost tier is a fact's property, so it is offered for facts."""
        exit_code = main(
            ["catalogue", "show", "--subject", "facts", "--output_format", "text", "--sort", "cost_tier"]
        )
        assert exit_code == ExitCode.OK
        assert capsys.readouterr().out.strip()

    def test_terms_default_to_catalogue_rank(self, capsys: pytest.CaptureFixture) -> None:
        """The natural order when no key is given."""
        assert main(["catalogue", "show", "--output_format", "text"]) == ExitCode.OK
        ranks = [int(value) for value in _column(capsys.readouterr().out, 0)]
        assert ranks == sorted(ranks)

    @pytest.mark.parametrize("subject,sort", [("facts", "impact"), ("terms", "cost_tier")])
    def test_a_key_that_does_not_apply_falls_back(
        self, subject: str, sort: str, caplog: pytest.LogCaptureFixture
    ) -> None:
        """Reported as a warning, then the subject's default order is used."""
        with caplog.at_level("WARNING"):
            exit_code = main(
                ["catalogue", "show", "--subject", subject, "--output_format", "text", "--sort", sort]
            )
        assert exit_code == ExitCode.OK
        assert "does not apply" in caplog.text

    def test_tier_is_no_longer_a_sort_key(self) -> None:
        """It meant the fact cost tier only, and was offered for terms where it did nothing."""
        with pytest.raises(SystemExit):
            main(["catalogue", "show", "--sort", "tier"])


class TestShowSelection:
    """Showing one term or fact rather than all of them."""

    def test_a_single_term_is_shown(self, capsys: pytest.CaptureFixture) -> None:
        """One row, for the term asked for."""
        assert (
            main(["catalogue", "show", "--output_format", "text", "--select", "intron_variant"])
            == ExitCode.OK
        )
        output = capsys.readouterr().out
        assert "intron_variant" in output
        assert "transcript_ablation" not in output

    def test_a_single_fact_is_shown(self, capsys: pytest.CaptureFixture) -> None:
        """The same selector works for the other subject."""
        exit_code = main(
            ["catalogue", "show", "--subject", "facts", "--output_format", "text", "--select", "overlap_cds"]
        )
        assert exit_code == ExitCode.OK
        output = capsys.readouterr().out
        assert "overlap_cds" in output
        assert "biotype_miRNA" not in output

    def test_an_unknown_term_is_a_usage_error(self, caplog: pytest.LogCaptureFixture) -> None:
        """Silently printing an empty table would look like a rendering fault."""
        with caplog.at_level("ERROR"):
            assert (
                main(["catalogue", "show", "--output_format", "text", "--select", "nope"]) == ExitCode.USAGE
            )
        assert "no catalogue term named" in caplog.text

    def test_an_unknown_fact_is_a_usage_error(self, caplog: pytest.LogCaptureFixture) -> None:
        """The counterpart for the fact subject."""
        with caplog.at_level("ERROR"):
            exit_code = main(
                ["catalogue", "show", "--subject", "facts", "--output_format", "text", "--select", "nope"]
            )
        assert exit_code == ExitCode.USAGE
        assert "no registered fact named" in caplog.text

    def test_selecting_in_graph_mode_keeps_the_neighbourhood(self, capsys: pytest.CaptureFixture) -> None:
        """A whole-registry diagram is unreadable; one fact and its edges is the useful view."""
        exit_code = main(
            [
                "catalogue",
                "show",
                "--subject",
                "facts",
                "--graph",
                "--output_format",
                "mermaid",
                "--select",
                "overlap_cds",
            ]
        )
        assert exit_code == ExitCode.OK
        output = capsys.readouterr().out
        assert "overlap_cds" in output
        assert "is_upstream" not in output

    def test_selecting_an_unknown_node_in_graph_mode_is_a_usage_error(self) -> None:
        """Same contract as the table."""
        assert (
            main(["catalogue", "show", "--graph", "--output_format", "dot", "--select", "nope"])
            == ExitCode.USAGE
        )

    def test_selecting_an_unknown_fact_in_graph_mode_is_a_usage_error(self) -> None:
        """The counterpart for the fact subject, whose selector keeps a neighbourhood."""
        exit_code = main(
            [
                "catalogue",
                "show",
                "--subject",
                "facts",
                "--graph",
                "--output_format",
                "dot",
                "--select",
                "nope",
            ]
        )
        assert exit_code == ExitCode.USAGE


class TestShowGraphs:
    """Drawing the same subjects as a graph."""

    @pytest.mark.parametrize("graph_format", ["svg", "dot", "mermaid"])
    def test_each_graph_format_renders(self, graph_format: str, capsys: pytest.CaptureFixture) -> None:
        """The three emitters the render layer offers."""
        assert main(["catalogue", "show", "--graph", "--output_format", graph_format]) == ExitCode.OK
        assert capsys.readouterr().out.strip()

    def test_the_subject_selects_which_graph(self, capsys: pytest.CaptureFixture) -> None:
        """Terms draw one decision tree per rule; facts draw the one dependency graph."""
        assert main(["catalogue", "show", "--graph", "--output_format", "dot"]) == ExitCode.OK
        assert 'digraph "intron_variant"' in capsys.readouterr().out
        assert (
            main(["catalogue", "show", "--subject", "facts", "--graph", "--output_format", "dot"])
            == ExitCode.OK
        )
        assert "Fact dependencies" in capsys.readouterr().out

    def test_unselected_svg_terms_form_an_html_catalogue(self, capsys: pytest.CaptureFixture) -> None:
        """Many standalone SVG roots are wrapped in one navigable document."""
        assert main(["catalogue", "show", "--graph", "--output_format", "svg"]) == ExitCode.OK
        output = capsys.readouterr().out
        document = output[output.index("<!DOCTYPE html>") :]
        assert 'href="#stop_gained"' in document
        assert 'id="stop_gained"' in document
        targets = set(re.findall(r'id="([^"]+)"', document))
        links = set(re.findall(r'href="#([^"]+)"', document))
        assert links <= targets

    def test_the_catalogue_overview_is_linked_last(self, capsys: pytest.CaptureFixture) -> None:
        """The final section summarises subsumption and links back to term details."""
        assert main(["catalogue", "show", "--graph", "--output_format", "svg"]) == ExitCode.OK
        output = capsys.readouterr().out
        overview = output.index('id="catalogue-overview"')
        assert 'href="#catalogue-overview"' in output
        assert '<a href="#stop_gained">' in output[overview:]
        assert overview > output.index('id="sequence_variant"')

    def test_fallback_terms_have_explanatory_sections(self, capsys: pytest.CaptureFixture) -> None:
        """The catalogue includes terms intentionally lacking a rule expression."""
        assert main(["catalogue", "show", "--graph", "--output_format", "svg"]) == ExitCode.OK
        output = capsys.readouterr().out
        intergenic = output.index('id="intergenic_variant"')
        assert "fallback term has no" in output[intergenic:]

    def test_terms_are_drawn_one_diagram_per_rule(self, capsys: pytest.CaptureFixture) -> None:
        """Merging them would show neither formula clearly."""
        # pylint: disable=import-outside-toplevel
        from ensembl.vep.annotators.consequence.rules.rulesets import build_fact_ruleset

        assert main(["catalogue", "show", "--graph", "--output_format", "mermaid"]) == ExitCode.OK
        drawn = capsys.readouterr().out.count("flowchart LR")
        assert drawn == len(build_fact_ruleset().rules)

    def test_selecting_a_term_draws_only_its_tree(self, capsys: pytest.CaptureFixture) -> None:
        """The whole catalogue at once is unreadable; one formula is the useful view."""
        exit_code = main(
            ["catalogue", "show", "--graph", "--output_format", "mermaid", "--select", "intron_variant"]
        )
        assert exit_code == ExitCode.OK
        output = capsys.readouterr().out
        assert output.count("flowchart LR") == 1
        assert "overlap_intron" in output

    def test_selecting_a_term_as_svg_is_standalone(self, capsys: pytest.CaptureFixture) -> None:
        """A selected expression remains usable as an SVG artefact."""
        exit_code = main(
            ["catalogue", "show", "--graph", "--output_format", "svg", "--select", "stop_gained"]
        )
        assert exit_code == ExitCode.OK
        output = capsys.readouterr().out
        svg = output[output.index("<svg") :]
        assert "<!DOCTYPE html>" not in svg
        assert "NOT is_structural" in svg

    def test_giving_a_graph_format_selects_graph_mode(self, capsys: pytest.CaptureFixture) -> None:
        """One option carries both the decision and the format."""
        assert main(["catalogue", "show", "--graph", "--output_format", "mermaid"]) == ExitCode.OK
        assert "flowchart LR" in capsys.readouterr().out

    def test_a_table_format_cannot_be_drawn(self, caplog: pytest.LogCaptureFixture) -> None:
        """A diagram in a table format would mean guessing which was meant."""
        with caplog.at_level("ERROR"):
            assert main(["catalogue", "show", "--graph", "--output_format", "text"]) == ExitCode.USAGE
        assert "--graph draws" in caplog.text

    def test_the_graph_flag_alone_takes_the_default_diagram_format(self) -> None:
        """--graph on its own draws, rather than demanding a format as well."""
        assert main(["catalogue", "show", "--graph"]) == ExitCode.OK

    def test_a_selected_svg_is_written_to_a_file(self, tmp_path: Path) -> None:
        """A selected graph is a standalone SVG artefact."""
        destination = tmp_path / "graph.svg"
        exit_code = main(
            [
                "catalogue",
                "show",
                "--graph",
                "--output_format",
                "svg",
                "--select",
                "stop_gained",
                "--output_file",
                str(destination),
            ]
        )
        assert exit_code == ExitCode.OK
        assert destination.read_text(encoding="utf-8").startswith("<svg")

    def test_the_catalogue_document_is_written_to_a_file(self, tmp_path: Path) -> None:
        """Unselected term graphs are one self-contained HTML document."""
        destination = tmp_path / "catalogue.html"
        exit_code = main(
            ["catalogue", "show", "--graph", "--output_format", "svg", "--output_file", str(destination)]
        )
        assert exit_code == ExitCode.OK
        assert destination.read_text(encoding="utf-8").startswith("<!DOCTYPE html>")


class TestCosts:
    """Measuring what each fact costs."""

    def test_it_reports_declared_against_measured(self, capsys: pytest.CaptureFixture) -> None:
        """The smallest corpus that still reaches every tier."""
        exit_code = main(["catalogue", "costs", "--fuzz", "0", "--step", "600", "--output_format", "text"])
        assert exit_code == ExitCode.OK
        output = capsys.readouterr().out
        assert "Measured ns" in output
        assert "overlap_cds" in output

    def test_it_writes_a_table_when_asked(self, tmp_path: Path) -> None:
        """The measured table is what apply_fact_cost_overrides reads back."""
        destination = tmp_path / "costs.json"
        exit_code = main(
            [
                "catalogue",
                "costs",
                "--fuzz",
                "0",
                "--step",
                "600",
                "--output_format",
                "text",
                "--output_file",
                str(destination),
            ]
        )
        assert exit_code == ExitCode.OK
        assert json.loads(destination.read_text(encoding="utf-8"))

    def test_it_leaves_the_registry_intact(self) -> None:
        """Profiling swaps every fact's compute for a shim and must restore them all."""
        # pylint: disable=import-outside-toplevel
        from ensembl.vep.annotators.consequence.facts.base import REGISTRY

        main(["catalogue", "costs", "--fuzz", "0", "--step", "600", "--output_format", "text"])
        for registered_fact in REGISTRY.values():
            assert registered_fact.compute.__name__ != "timed_compute"


class TestValidateOptions:
    """The optional validation tiers."""

    def test_verify_requires_runs_over_a_corpus(self, caplog: pytest.LogCaptureFixture) -> None:
        """The check that a declared requires actually holds on real overlaps."""
        with caplog.at_level("INFO"):
            main(["catalogue", "validate", "--verify_requires", "--fuzz", "0", "--step", "600"])
        assert "verify_requires over" in caplog.text

    def test_an_unavailable_solver_is_skipped_rather_than_failing(
        self, caplog: pytest.LogCaptureFixture, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """z3 is an optional extra, so its absence must not break the gate."""
        # pylint: disable=import-outside-toplevel
        from ensembl.vep.annotators.consequence.rules import validation

        def unavailable(_ruleset: object) -> list[str]:
            raise RuntimeError("z3 is not installed")

        monkeypatch.setattr(validation.SolverChecks, "findings", unavailable)
        with caplog.at_level("WARNING"):
            exit_code = main(["catalogue", "validate", "--solver", "z3"])
        assert exit_code == ExitCode.OK
        assert "solver skipped" in caplog.text


class TestAdd:
    """Adding a term through the command line."""

    def test_a_term_is_written(self, catalogue_copy: Path) -> None:
        """The end-to-end path, including validation before the write."""
        exit_code = main(
            [
                "catalogue",
                "add",
                "--catalogue",
                str(catalogue_copy),
                "--term",
                "sample_term",
                "--accession",
                "SO:0009999",
                "--impact",
                "MODIFIER",
                "--when",
                '"is_snv"',
                "--rank",
                "20",
            ]
        )
        assert exit_code == ExitCode.OK
        assert "sample_term" in _terms(catalogue_copy)

    def test_a_dry_run_writes_nothing(self, catalogue_copy: Path) -> None:
        """The change is reported and validated, but the file is untouched."""
        before = catalogue_copy.read_text(encoding="utf-8")
        exit_code = main(
            [
                "catalogue",
                "add",
                "--catalogue",
                str(catalogue_copy),
                "--term",
                "sample_term",
                "--accession",
                "SO:0009999",
                "--impact",
                "MODIFIER",
                "--when",
                '"is_snv"',
                "--dry_run",
            ]
        )
        assert exit_code == ExitCode.OK
        assert catalogue_copy.read_text(encoding="utf-8") == before

    def test_no_term_is_a_usage_error(self, catalogue_copy: Path) -> None:
        """There is nothing to add without one."""
        assert main(["catalogue", "add", "--catalogue", str(catalogue_copy)]) == ExitCode.USAGE

    def test_a_term_with_no_when_is_refused(self, catalogue_copy: Path) -> None:
        """A formula is never guessed from a term name."""
        exit_code = main(
            [
                "catalogue",
                "add",
                "--catalogue",
                str(catalogue_copy),
                "--term",
                "sample_term",
                "--accession",
                "SO:0009999",
                "--impact",
                "MODIFIER",
            ]
        )
        assert exit_code == ExitCode.FAILURE
        assert "sample_term" not in _terms(catalogue_copy)

    def test_a_bad_accession_is_a_usage_error(self, catalogue_copy: Path) -> None:
        """Caught before anything is written."""
        exit_code = main(
            [
                "catalogue",
                "add",
                "--catalogue",
                str(catalogue_copy),
                "--term",
                "sample_term",
                "--accession",
                "nope",
                "--impact",
                "MODIFIER",
                "--when",
                '"is_snv"',
            ]
        )
        assert exit_code == ExitCode.USAGE


class TestRemove:
    """Removing a term through the command line."""

    def test_a_term_is_removed(self, catalogue_copy: Path) -> None:
        """And cleaned out of every other entry's subsumes."""
        exit_code = main(
            ["catalogue", "remove", "--catalogue", str(catalogue_copy), "--term", "feature_elongation"]
        )
        assert exit_code == ExitCode.OK
        assert "feature_elongation" not in _terms(catalogue_copy)

    def test_an_unknown_term_is_a_usage_error(self, catalogue_copy: Path) -> None:
        """A typo must not silently succeed."""
        exit_code = main(["catalogue", "remove", "--catalogue", str(catalogue_copy), "--term", "nope"])
        assert exit_code == ExitCode.USAGE


class TestWriteSafety:
    """An installed vep must not rewrite its own wheel."""

    def test_a_write_without_a_destination_is_refused(self) -> None:
        """Neither --catalogue nor --catalogue_overwrite was given."""
        assert main(["catalogue", "remove", "--term", "feature_elongation"]) == ExitCode.USAGE

    def test_the_packaged_catalogue_is_untouched_by_a_refused_write(self) -> None:
        """The guard has to hold even when the command fails."""
        before = catalogue_path().read_text(encoding="utf-8")
        main(["catalogue", "remove", "--term", "feature_elongation"])
        assert catalogue_path().read_text(encoding="utf-8") == before


def _column(text: str, index: int) -> list[str]:
    """Return one column of a fixed-width table, skipping the header and its rule."""
    values = []
    for line in text.splitlines():
        if not line.strip() or line.lstrip().startswith("-") or "  " not in line:
            continue
        fields = [field for field in line.split("  ") if field.strip()]
        if len(fields) > index and not fields[0].startswith(("Rank", "Fact")):
            values.append(fields[index].strip())
    return values
