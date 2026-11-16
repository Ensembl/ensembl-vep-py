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
"""Explain, end to end, why an annotator did or did not emit a term for one overlap."""

from __future__ import annotations

import logging
from collections.abc import Iterator
from pathlib import Path
from typing import ClassVar

from ensembl.vep.annotators.consequence.graphs import derivation_graphs
from ensembl.vep.annotators.consequence.rules.catalogue import load_consequence_definitions
from ensembl.vep.annotators.factory import build_annotator
from ensembl.vep.config import Config
from ensembl.vep.config.report import resolve_report_format
from ensembl.vep.corpus.base import Corpus
from ensembl.vep.corpus.registry import build_corpus
from ensembl.vep.corpus.spec import CorpusSpec
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.output.parser import parser_for
from ensembl.vep.render.graph import Graph, render_dot, render_mermaid, render_svg, tree
from ensembl.vep.render.page import Section, render_page
from ensembl.vep.runner.base import ExitCode
from ensembl.vep.runner.registry import register_runner
from ensembl.vep.sources.build import TranscriptSource
from ensembl.vep.sources.file.vcf.variation import VCFVariationSource
from ensembl.vep.runner.explain.overlaps import (
    NamedOverlap,
    VariantSelector,
    build_explain_sources,
    find_named_overlaps,
    named_overlaps,
    names_feature,
)
from ensembl.vep.trace import AnnotationTrace, format_explanation, tracing_enabled
from ensembl.vep.trace.explain import diff_terms, diff_traces
from ensembl.vep.trace.summary import format_summary, summary_html

__all__ = ["Runner"]

logger = logging.getLogger(__name__)

GRAPH_EMITTERS = {"dot": render_dot, "mermaid": render_mermaid}


def _plain_allele(value: str | None) -> str:
    """Return an allele with either spelling of "empty" collapsed to one.

    A corpus spells an empty allele ``""`` and a VCF reader spells it ``-``, so matching an overlap
    to a record has to treat the two as the same allele.
    """
    return "" if value is None or value == "-" else value.upper()


TERM_GROUPS = {"emitted": "Emitted", "pruned": "Pruned", "": "Not fired"}
TERM_GROUP_ORDER = {"emitted": 0, "pruned": 1, "": 2}

NO_TRACE = (
    "the native rust.* annotators never record a trace. That is not a gap: every annotator is "
    "proven byte-identical to consequence.linear by the conformance suite, so the same overlap "
    "run through consequence.linear has the same derivation to show."
)


@register_runner("explain")
class Runner:
    """Run one overlap through an annotator with tracing on, and report the whole derivation.

    Attributes:
        name (ClassVar[str]): The registry name.
        sections (ClassVar[tuple[str, ...]]): The configuration sections read.
        description (ClassVar[str]): One-line description for the subcommand list.
    """

    name: ClassVar[str] = "explain"
    sections: ClassVar[tuple[str, ...]] = ("logging", "source", "corpus", "annotation", "explain")
    description: ClassVar[str] = "explain one overlap's derivation, fact by fact and rule by rule"

    def __init__(self, config: Config) -> None:
        """Hold the configuration."""
        self.config = config
        self._sources: tuple[VCFVariationSource, TranscriptSource] | None = None
        self._built_corpus: Corpus | None = None

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
            if self._built_corpus is not None:
                self._built_corpus.close()

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
            self._emit(self._listing())
            return ExitCode.OK

        if explain.name is None and explain.index is None and explain.variant is None:
            logger.error(
                "select an overlap with --name, --variant or --index; "
                "'vep corpus describe' or --list shows what is available"
            )
            return ExitCode.USAGE

        selected_name, overlap = self._select()
        logger.info("explaining %s over %s", selected_name, overlap.feature.stable_id)
        if explain.compare:
            return self._compare(overlap)
        if explain.against_input:
            return self._against_input(overlap)
        return self._explain_overlap(overlap)

    def _select(self) -> NamedOverlap:
        """Return the single overlap the selectors identify.

        A name or a variant identifies the variant, not the overlap, so one that pairs with several
        features needs ``--transcript`` as well. Saying which stable ids matched is the whole value
        of the error, so the message names them.

        Returns:
            NamedOverlap: The selected overlap and the name it came under.

        Raises:
            LookupError: If nothing matches, if several do and no feature was named, or if the
                index is past the end.
            ValueError: If the variant selector is malformed, or an allele number was given
                without a ``CHR:POS:ID`` variant.
        """
        explain = self.config.explain
        if explain.index is not None:
            for position, named_overlap in enumerate(self._named_overlaps()):
                if position == explain.index:
                    return named_overlap
            raise LookupError(f"index {explain.index} is past the end of {self._description()}")

        if explain.allele_num is not None and explain.variant is None:
            raise ValueError("an allele number selects from a CHR:POS:ID variant; give the variant too")
        selector = (
            VariantSelector.parse(explain.variant, allele_num=explain.allele_num)
            if explain.variant is not None
            else None
        )
        selector_text = explain.name if explain.name is not None else explain.variant
        feature = explain.transcript

        matches: list[NamedOverlap] = []
        for candidate_name, overlap in self._candidate_overlaps(explain.name, selector):
            if explain.name is not None and candidate_name != explain.name:
                continue
            if selector is not None and not selector.matches(overlap.overlapping_feature):
                continue
            if feature is not None and not names_feature(overlap.feature.stable_id, feature):
                continue
            matches.append((candidate_name, overlap))
            if feature is not None:
                # A named feature makes the pair unique, so nothing later can make this ambiguous.
                break

        if not matches:
            suffix = f" for feature {feature!r}" if feature is not None else ""
            raise LookupError(f"no overlap matching {selector_text!r} in {self._description()}{suffix}")
        if len(matches) > 1:
            stable_ids = sorted({overlap.feature.stable_id for _name, overlap in matches})
            raise LookupError(
                f"{len(matches)} overlaps match {selector_text!r} (features: {', '.join(stable_ids)}); "
                "disambiguate by naming one"
            )
        return matches[0]

    def _candidate_overlaps(
        self, name: str | None, selector: VariantSelector | None
    ) -> Iterator[NamedOverlap]:
        """Yield the overlaps that could match, pairing no more of the input than it takes.

        A real input is searched for the named variant alone; a generated corpus is walked, there
        being nothing to query.

        Args:
            name (str | None): The variant name to look for.
            selector (VariantSelector | None): The variant's locus and alleles, or its identifier.

        Returns:
            Iterator[NamedOverlap]: Every overlap that still has to be filtered by the caller.
        """
        if self.config.source.input_file is None:
            return iter(self._corpus().named_overlaps)
        variation_source, transcript_source = self._file_sources()
        return find_named_overlaps(
            variation_source,
            transcript_source,
            upstream_distance=self.config.source.upstream_distance,
            downstream_distance=self.config.source.downstream_distance,
            name=name,
            selector=selector,
        )

    def _named_overlaps(self) -> Iterator[NamedOverlap]:
        """Yield every selectable overlap, in the order ``--index`` counts."""
        if self.config.source.input_file is None:
            return iter(self._corpus().named_overlaps)
        variation_source, transcript_source = self._file_sources()
        return named_overlaps(
            variation_source,
            transcript_source,
            upstream_distance=self.config.source.upstream_distance,
            downstream_distance=self.config.source.downstream_distance,
        )

    def _file_sources(self) -> tuple[VCFVariationSource, TranscriptSource]:
        """Return the input and transcript sources the options name, built once per run."""
        if self._sources is None:
            self._sources = build_explain_sources(self.config.source)
        return self._sources

    def _corpus(self) -> Corpus:
        """Return the generated corpus the options describe."""
        if self._built_corpus is None:
            self._built_corpus = build_corpus(CorpusSpec.from_options(self.config.corpus))
        return self._built_corpus

    def _description(self) -> str:
        """Describe where the overlaps come from, for an error or a log line."""
        source = self.config.source
        if source.input_file is None:
            return self._corpus().description
        features = source.gff if source.gff is not None else source.gtf
        feature_name = "no transcripts" if features is None else Path(features).name
        distance = ",".join(str(value) for value in source.distance)
        return f"{Path(source.input_file).name} against {feature_name} (distance {distance})"

    def _listing(self) -> str:
        """Return every selectable overlap, one per line, as index, name and feature."""
        lines = []
        for index, (name, overlap) in enumerate(self._named_overlaps()):
            lines.append(f"[{index:>5}] {name:<42} {overlap.feature.stable_id}")
        return "\n".join(lines)

    def _annotate(
        self, annotator_name: str, overlap: TranscriptVariationAllele
    ) -> tuple[list[str], AnnotationTrace | None]:
        """Run one annotator over one overlap with tracing on.

        Args:
            annotator_name (str): A registered annotator name.
            overlap (TranscriptVariationAllele): The overlap to annotate.

        Returns:
            tuple[list[str], AnnotationTrace | None]: The emitted terms, and the derivation -- None
                for a native rust.* annotator, which never traces.
        """
        annotation = self.config.annotation
        annotator = build_annotator(
            annotator=annotator_name,
            batch_size=annotation.batch_size,
            early_subsumption_skip=annotation.early_subsumption_skip,
            zone_map=annotation.zone_map,
            zone_map_bytes=annotation.zone_map_bytes,
            codon_cache=annotation.codon_cache,
            codon_cache_entries=annotation.codon_cache_entries,
            mapper_cache=annotation.mapper_cache,
            mapper_cache_bytes=annotation.mapper_cache_bytes,
            readthrough=annotation.readthrough,
            fact_cost_override=annotation.fact_cost_override,
            upstream_distance=self.config.source.upstream_distance,
            downstream_distance=self.config.source.downstream_distance,
        )
        with tracing_enabled(quiet=self.config.logging.effective_level().upper() != "DEBUG"):
            consequences = annotator.annotate_feature_overlap(overlap)
        return [consequence.term for consequence in consequences], getattr(annotator, "last_trace", None)

    def _explain_overlap(self, overlap: TranscriptVariationAllele) -> ExitCode:
        """Render one annotator's derivation of one overlap, after a summary of the overlap itself.

        The graph source formats (dot, mermaid) are emitted bare, so that they still parse.
        """
        annotator_name = self.config.annotation.annotator
        terms, trace = self._annotate(annotator_name, overlap)
        if trace is None:
            logger.error("%s did not record a trace: %s", annotator_name, NO_TRACE)
            logger.info("emitted terms: %s", terms)
            return ExitCode.FAILURE

        report_format = resolve_report_format(
            self.config.explain.output_format, graph=self.config.explain.graph, table_default="text"
        )
        if report_format == "text":
            self._emit(f"{format_summary(overlap)}\n\n{format_explanation(trace)}")
        else:
            self._emit(self._graphs(trace, report_format, overlap))
        return ExitCode.OK

    @staticmethod
    def _graphs(trace: AnnotationTrace, graph_format: str, overlap: TranscriptVariationAllele) -> str:
        """Render every term's result as source graphs or one navigable HTML document."""
        graphs = derivation_graphs(trace)
        if graph_format == "svg":
            return Runner._graph_page(trace, graphs, overlap)

        emitter = GRAPH_EMITTERS[graph_format]
        return "\n\n".join(emitter(graph) for graph in graphs)

    @staticmethod
    def _graph_page(
        trace: AnnotationTrace, graphs: tuple[Graph, ...], overlap: TranscriptVariationAllele
    ) -> str:
        """Return the overlap's summary, then every term proof in outcome order, from a side menu."""
        definitions = load_consequence_definitions()
        ranks = {definition.term: definition.rank for definition in definitions}
        entries: list[tuple[int, int, Section]] = []

        for graph in graphs:
            state = graph.node(graph.title).state
            group = TERM_GROUPS[state]
            entries.append(
                (
                    TERM_GROUP_ORDER[state],
                    ranks[graph.title],
                    Section(
                        identifier=graph.title,
                        title=graph.title,
                        blocks=(render_svg(graph, tree(graph)),),
                        group=group,
                    ),
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
                    TERM_GROUP_ORDER[state],
                    definition.rank,
                    Section(
                        identifier=definition.term,
                        title=definition.term,
                        blocks=(f"<p>{explanation}</p>",),
                        group=TERM_GROUPS[state],
                    ),
                )
            )

        entries.sort(key=lambda entry: entry[:2])
        sections = tuple(entry[2] for entry in entries)
        counts = {group: 0 for group in TERM_GROUPS.values()}
        for section in sections:
            counts[section.group] += 1
        subtitle = (
            f"{trace.overlap_label} · annotator: {trace.annotator_name} · "
            f"{counts['Emitted']} emitted · {counts['Pruned']} pruned · "
            f"{counts['Not fired']} not fired · expression nodes: blue true, red false, "
            "white not observed"
        )
        summary = Section(identifier="summary", title="Summary", blocks=(summary_html(overlap),))
        return render_page("VEP consequence explanation", (summary, *sections), subtitle=subtitle)

    def _against_input(self, overlap: TranscriptVariationAllele) -> ExitCode:
        """Diff our derivation of one overlap against the annotation the input file already carries.

        The point is not the set difference but the reason: for a term the file has and we do not,
        ``diff_terms`` reports from our own ledger whether the rule was false, or fired and was
        subsumed, or is a term the catalogue does not model at all.

        A caveat worth knowing rather than hiding: the reader seeds terms per *allele*, because the
        CSQ layout it parses keeps the allele and consequence columns and drops the feature. So the
        expected terms are everything the file recorded for this allele, across every feature it
        annotated -- directly comparable when one feature overlaps, a superset when several do.

        Args:
            overlap (TranscriptVariationAllele): The overlap to explain.

        Returns:
            ExitCode: OK once the diff was written, or the code saying why not.
        """
        input_file = self.config.source.input_file
        if input_file is None:
            logger.error("--against_input diffs against an annotated input; give --input_file")
            return ExitCode.USAGE

        annotator_name = self.config.annotation.annotator
        terms, trace = self._annotate(annotator_name, overlap)
        if trace is None:
            logger.error("%s did not record a trace: %s", annotator_name, NO_TRACE)
            logger.info("emitted terms: %s", terms)
            return ExitCode.FAILURE

        expected = self._recorded_terms(overlap)
        if expected is None:
            logger.error("no record in %s matches the selected variant", input_file)
            return ExitCode.USAGE

        self._emit(f"{format_summary(overlap)}\n\n{diff_terms(trace, expected, f'{input_file.name} CSQ')}")
        return ExitCode.OK

    def _recorded_terms(self, overlap: TranscriptVariationAllele) -> list[str] | None:
        """Return the terms the input file recorded for this overlap's allele, or None if absent.

        Args:
            overlap (TranscriptVariationAllele): The overlap whose allele to look up.

        Returns:
            list[str] | None: The recorded terms, or None when no record matches the allele.
        """
        source = self.config.source
        assert source.input_file is not None
        allele = overlap.overlapping_feature
        wanted = (
            allele.slice.seq_region_name if allele.slice is not None else None,
            allele.start,
            _plain_allele(allele.reference_allele),
            _plain_allele(allele.alternate_allele),
        )

        reader = parser_for(source.input_file, source.format)
        for variation in reader.parse(source.input_file):
            for recorded in variation.variation_alleles or ():
                found = (
                    recorded.slice.seq_region_name if recorded.slice is not None else None,
                    recorded.start,
                    _plain_allele(recorded.reference_allele),
                    _plain_allele(recorded.alternate_allele),
                )
                if found == wanted:
                    return [consequence.term for consequence in (recorded.overlap_consequences or ())]
        return None

    def _compare(self, overlap: TranscriptVariationAllele) -> ExitCode:
        """Run two annotators over one overlap and report where their derivations part company."""
        compare = self.config.explain.compare
        if len(compare) != 2:
            logger.error("--compare takes exactly two annotator names, given as NAME,NAME")
            return ExitCode.USAGE

        name_a, name_b = compare
        terms_a, trace_a = self._annotate(name_a, overlap)
        terms_b, trace_b = self._annotate(name_b, overlap)

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

        self._emit(f"{format_summary(overlap)}\n\n{report}")
        return ExitCode.OK if terms_a == terms_b else ExitCode.FAILURE

    def _emit(self, rendered: str) -> None:
        """Write a rendering to the configured destination."""
        destination = self.config.explain.output_file
        if destination is None:
            print(rendered)
            return
        destination.write_text(rendered, encoding="utf-8")
        logger.info("wrote %s", destination)
