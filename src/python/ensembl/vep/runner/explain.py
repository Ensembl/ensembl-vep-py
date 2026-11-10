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
"""End-to-end overlap consequence explanations."""

from __future__ import annotations

import logging
from collections.abc import Sequence
from html import escape
from typing import ClassVar

from ensembl.vep.annotators.consequence.graphs import derivation_graphs
from ensembl.vep.annotators.consequence.rules.catalogue import load_consequence_definitions
from ensembl.vep.config import Config
from ensembl.vep.config.report import resolve_report_format
from ensembl.vep.explain.annotate import NO_TRACE, annotate_with_trace
from ensembl.vep.explain.compare import diff_terms, diff_traces
from ensembl.vep.explain.selection import OverlapSelection
from ensembl.vep.explain.summary import format_summary, summary_html
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.render.graph import Graph, render_dot, render_mermaid, render_svg, tree
from ensembl.vep.render.page import Section, render_page
from ensembl.vep.runner.base import ExitCode
from ensembl.vep.runner.registry import register_runner
from ensembl.vep.trace import AnnotationTrace

__all__ = ["Runner"]

logger = logging.getLogger(__name__)

GRAPH_EMITTERS = {"dot": render_dot, "mermaid": render_mermaid}


TERM_GROUPS = {"emitted": "Emitted", "pruned": "Pruned", "": "Not fired"}
TERM_GROUP_ORDER = {"emitted": 0, "pruned": 1, "": 2}


@register_runner("explain")
class Runner:
    """Runner for explaining one traced overlap annotation.

    Attributes:
        name (ClassVar[str]): The registry name.
        sections (ClassVar[tuple[str, ...]]): The configuration sections read.
        description (ClassVar[str]): One-line description for the subcommand list.
    """

    name: ClassVar[str] = "explain"
    sections: ClassVar[tuple[str, ...]] = ("logging", "source", "corpus", "annotation", "explain")
    description: ClassVar[str] = "explain one overlap's derivation, fact by fact and rule by rule"

    def __init__(self, config: Config) -> None:
        """Hold the configuration and the selection it describes."""
        self.config = config
        self._selection = OverlapSelection(config=config)

    def run(self) -> ExitCode:
        """Select an overlap, annotate it with tracing on, and render the derivation.

        Returns:
            ExitCode: OK once a derivation was rendered, or the code saying why not.
        """
        try:
            return self._explain()
        except ValueError as unusable:
            # Building the sources needs transcripts to pair against, and a selector has to have a
            # shape. Both are the caller's mistake, so a usage error rather than a failure.
            logger.error("%s", unusable)
            return ExitCode.USAGE
        except LookupError as unknown:
            # Nothing matched the selectors, several did, or an index ran past the end.
            logger.error("%s", unknown)
            return ExitCode.USAGE
        except OSError as unreadable:
            # The sources read lazily, so an unreadable input surfaces here rather than at build
            # time. The readers raise instead of exiting, which leaves the exit code to this runner.
            logger.error("%s", unreadable)
            return ExitCode.USAGE
        finally:
            self._selection.close()

    def _explain(self) -> ExitCode:
        """Render the derivation the options ask for.

        Returns:
            ExitCode: OK once a derivation was rendered, or the code saying why not.

        Raises:
            LookupError: If the selectors identify no overlap, or more than one.
            ValueError: If the configured sources or the selectors cannot be used.
        """
        explain = self.config.explain
        # Resolved before anything is selected: an unrenderable request is wrong whatever matched.
        resolve_report_format(explain.output_format, graph=explain.graph, table_default="text")
        if explain.list_overlaps:
            self._emit(self._selection.listing())
            return ExitCode.OK

        if explain.name is None and explain.index is None and explain.variant is None:
            logger.error(
                "select an overlap with --name, --variant or --index; "
                "'vep corpus describe' or --list shows what is available"
            )
            return ExitCode.USAGE

        selected_name, overlap = self._selection.select()
        logger.info("explaining %s over %s", selected_name, overlap.feature.stable_id)
        if explain.compare:
            return self._compare(overlap)
        return self._explain_overlap(overlap)

    def _explain_overlap(self, overlap: TranscriptVariationAllele) -> ExitCode:
        """Render one annotator's derivation of one overlap, after a summary of the overlap itself.

        The graph source formats (dot, mermaid) are emitted bare, so that they still parse.
        """
        annotator_name = self.config.annotation.annotator
        terms, trace = annotate_with_trace(self.config, annotator_name, overlap)
        if trace is None:
            logger.error("%s did not record a trace: %s", annotator_name, NO_TRACE)
            logger.info("emitted terms: %s", terms)
            return ExitCode.FAILURE

        report_format = resolve_report_format(
            self.config.explain.output_format, graph=self.config.explain.graph, table_default="text"
        )
        if report_format == "text":
            self._emit(f"{format_summary(overlap)}\n\n{trace.format_explanation()}")
        else:
            self._emit(self._graphs([(annotator_name, trace)], report_format, overlap))
        return ExitCode.OK

    def _compare(self, overlap: TranscriptVariationAllele) -> ExitCode:
        """Run two annotators over one overlap and report where their derivations part company.

        A graph format draws the comparison rather than only describing it. The report still
        carries the verdict, so for the source formats -- which must stay parseable -- it is
        logged rather than printed.
        """
        compare = self.config.explain.compare
        if len(compare) != 2:
            logger.error("--compare takes exactly two annotator names, given as NAME,NAME")
            return ExitCode.USAGE

        name_a, name_b = compare
        terms_a, trace_a = annotate_with_trace(self.config, name_a, overlap)
        terms_b, trace_b = annotate_with_trace(self.config, name_b, overlap)

        if trace_a is not None and trace_b is not None:
            report = diff_traces(trace_a, trace_b)
        elif trace_a is not None:
            logger.warning("%s recorded no trace: %s", name_b, NO_TRACE)
            report = diff_terms(trace_a, terms_b, name_b)
        elif trace_b is not None:
            logger.warning("%s recorded no trace: %s", name_a, NO_TRACE)
            report = diff_terms(trace_b, terms_a, name_a)
        else:
            logger.error("neither %s nor %s records a trace: %s", name_a, name_b, NO_TRACE)
            return ExitCode.FAILURE

        report_format = resolve_report_format(
            self.config.explain.output_format, graph=self.config.explain.graph, table_default="text"
        )
        if report_format == "text":
            self._emit(f"{format_summary(overlap)}\n\n{report}")
        else:
            traced = [
                (name, trace) for name, trace in ((name_a, trace_a), (name_b, trace_b)) if trace is not None
            ]
            if report_format != "svg":
                logger.info("%s", report)
            self._emit(self._graphs(traced, report_format, overlap, report=report))
        return ExitCode.OK if terms_a == terms_b else ExitCode.FAILURE

    @staticmethod
    def _graphs(
        traces: Sequence[tuple[str, AnnotationTrace]],
        graph_format: str,
        overlap: TranscriptVariationAllele,
        report: str = "",
    ) -> str:
        """Render every term's result as source graphs or one navigable HTML document.

        One annotator draws one set. Two draw one set when they produce the same drawing and the
        same terms -- a second identical copy says nothing -- and a set each when they do not,
        which is the divergence a comparison is asking about. Agreement is judged on the drawings
        rather than on the emitted terms: every annotator agrees on the terms by construction (the
        conformance suite proves it), while the strategies observe different facts to get there,
        and that difference is the whole subject of the picture.

        Args:
            traces (Sequence[tuple[str, AnnotationTrace]]): Each annotator's name and derivation.
            graph_format (str): ``svg`` for the HTML document, else a graph source format.
            overlap (TranscriptVariationAllele): The overlap being explained.
            report (str): A comparison report to carry into the document, for ``svg``.

        Returns:
            str: The rendering.
        """
        drawn = [(name, derivation_graphs(trace), trace) for name, trace in traces]
        compared = tuple(name for name, _trace in traces) if len(traces) > 1 else ()
        if (
            len(drawn) == 2
            and drawn[0][1] == drawn[1][1]
            and drawn[0][2].emitted_terms == drawn[1][2].emitted_terms
        ):
            drawn = drawn[:1]

        if graph_format == "svg":
            return Runner._graph_page(drawn, overlap, report=report, compared=compared)

        emitter = GRAPH_EMITTERS[graph_format]
        return "\n\n".join(emitter(graph) for _name, graphs, _trace in drawn for graph in graphs)

    @staticmethod
    def _graph_page(
        drawn: Sequence[tuple[str, tuple[Graph, ...], AnnotationTrace]],
        overlap: TranscriptVariationAllele,
        report: str = "",
        compared: Sequence[str] = (),
    ) -> str:
        """Return the overlap's summary, then every term proof in outcome order, from a side menu.

        Two sets are kept apart in the menu, each group prefixed by its annotator, so the same
        term from each side is two entries rather than one anchor fighting over two drawings.

        Args:
            drawn (Sequence[tuple[str, tuple[Graph, ...], AnnotationTrace]]): Each set's annotator
                name, graphs and derivation. One set, or one per annotator.
            overlap (TranscriptVariationAllele): The overlap being explained.
            report (str): A comparison report to show above the proofs, if any.
            compared (Sequence[str]): The annotators being compared, empty for a plain explanation.

        Returns:
            str: The HTML document.
        """
        definitions = load_consequence_definitions()
        ranks = {definition.term: definition.rank for definition in definitions}
        several = len(drawn) > 1
        entries: list[tuple[int, int, int, Section]] = []

        for position, (name, graphs, trace) in enumerate(drawn):

            def section(term: str, blocks: tuple[str, ...], state: str, name: str = name) -> Section:
                """Return one term's menu entry, namespaced by annotator when there are two sets."""
                return Section(
                    identifier=f"{name}:{term}" if several else term,
                    title=term,
                    blocks=blocks,
                    group=f"{name} · {TERM_GROUPS[state]}" if several else TERM_GROUPS[state],
                )

            for graph in graphs:
                state = graph.node(graph.title).state
                entries.append(
                    (
                        position,
                        TERM_GROUP_ORDER[state],
                        ranks[graph.title],
                        section(graph.title, (render_svg(graph, tree(graph)),), state),
                    )
                )

            drawn_terms = {graph.title for graph in graphs}
            for definition in definitions:
                if definition.term in drawn_terms:
                    continue
                state = "emitted" if definition.term in trace.emitted_terms else ""
                explanation = (
                    "This fallback term was emitted because no expression term survived."
                    if state == "emitted"
                    else "This fallback term has no expression to fire for this overlap."
                )
                entries.append(
                    (
                        position,
                        TERM_GROUP_ORDER[state],
                        definition.rank,
                        section(definition.term, (f"<p>{explanation}</p>",), state),
                    )
                )

        entries.sort(key=lambda entry: entry[:3])
        sections = tuple(entry[3] for entry in entries)
        trace = drawn[0][2]
        nodes = "expression nodes: blue true, red false, white not observed"
        if compared and several:
            subtitle = (
                f"{trace.overlap_label} · comparing {' vs '.join(compared)} · their derivations "
                f"differ, so each side's proofs are listed separately · {nodes}"
            )
        elif compared:
            subtitle = (
                f"{trace.overlap_label} · comparing {' vs '.join(compared)} · their derivations "
                f"are identical, so one set of proofs stands for both · {nodes}"
            )
        else:
            counts = {group: 0 for group in TERM_GROUPS.values()}
            for section_entry in sections:
                counts[section_entry.group] += 1
            subtitle = (
                f"{trace.overlap_label} · annotator: {trace.annotator_name} · "
                f"{counts['Emitted']} emitted · {counts['Pruned']} pruned · "
                f"{counts['Not fired']} not fired · {nodes}"
            )

        leading: list[Section] = [
            Section(identifier="summary", title="Summary", blocks=(summary_html(overlap),))
        ]
        if report:
            leading.append(
                Section(
                    identifier="comparison",
                    title="Comparison",
                    blocks=(f"<pre>{escape(report)}</pre>",),
                )
            )
        return render_page("VEP consequence explanation", (*leading, *sections), subtitle=subtitle)

    def _emit(self, rendered: str) -> None:
        """Write a rendering to the configured destination."""
        destination = self.config.explain.output_file
        if destination is None:
            print(rendered)
            return
        destination.write_text(rendered, encoding="utf-8")
        logger.info("wrote %s", destination)
