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
"""Completed annotation traces and text rendering.

An ``AnnotationTrace`` is immutable and identifies both its overlap and annotator. Reconstruction
builds it after annotation from the fact memo, rule results, and the small amount of live state in
``TraceRecorder``. Optional event collections use empty tuples so consumers can iterate without
special cases.
"""

from __future__ import annotations

import logging
from collections.abc import Iterator
from contextlib import contextmanager

from attrs import field, frozen

from ensembl.vep.common import CostTier
from ensembl.vep.common.constants import UNKNOWN_VALUE
from ensembl.vep.trace.events import FactEvent, FactOrigin, PruneEvent, RuleEvent

__all__ = [
    "AnnotationTrace",
    "annotation_tracing",
    "describe_overlap",
]


@frozen(kw_only=True, auto_attribs=True)
class AnnotationTrace:
    """Record one overlap's derivation.

    Attributes:
        fact_events (tuple[FactEvent, ...]): Facts touched by the derivation, in materialisation
            order, followed by facts proved false through ``requires`` propagation.
        rule_events (tuple[RuleEvent, ...]): Every rule's verdict, in rank order.
        prune_events (tuple[PruneEvent, ...]): Every subsumption decision.
        emitted_terms (tuple[str, ...]): Terms remaining after subsumption.
        annotator_name (str): Registered annotator name.
        overlap_label (str): One-line description of the annotated overlap.
    """

    fact_events: tuple[FactEvent, ...] = field(default=())
    rule_events: tuple[RuleEvent, ...] = field(default=())
    prune_events: tuple[PruneEvent, ...] = field(default=())
    emitted_terms: tuple[str, ...] = field(default=())
    annotator_name: str = field(default="")
    overlap_label: str = field(default="")

    def events_by_origin(self, origin: FactOrigin) -> tuple[FactEvent, ...]:
        """Return the fact events with one origin, in trace order."""
        selected = []
        for event in self.fact_events:
            if event.origin is origin:
                selected.append(event)
        return tuple(selected)

    @property
    def computed_cost(self) -> int:
        """Return the declared cost charged by this derivation."""
        return sum(event.cost for event in self.fact_events if event.origin is FactOrigin.COMPUTED)

    @property
    def avoided_cost(self) -> int:
        """Return the declared cost avoided through ``requires`` propagation."""
        avoided = 0
        for event in self.fact_events:
            if event.origin is FactOrigin.FORCED_FALSE:
                avoided += event.cost
        return avoided

    def format_explanation(self) -> str:
        """Return this trace as a human-readable derivation.

        Returns:
            str: Facts, rule verdicts, pruning decisions, and emitted terms.
        """
        lines = [f"{self.overlap_label}          annotator: {self.annotator_name}"]

        computed = self.events_by_origin(FactOrigin.COMPUTED)
        prefilled = self.events_by_origin(FactOrigin.PREFILLED)
        cached = self.events_by_origin(FactOrigin.CACHED)
        forced_false = self.events_by_origin(FactOrigin.FORCED_FALSE)
        reads = sum(event.reads for event in self.fact_events)

        lines.append(
            f"  facts: {len(computed)} computed (cost {self.computed_cost}), {len(prefilled)} prefilled, "
            f"{len(forced_false)} forced False (cost {self.avoided_cost} avoided), "
            f"{len(cached)} cache hits, {reads} reads"
        )

        for event in self.fact_events:
            origin_text = event.origin.value
            if event.origin is FactOrigin.FORCED_FALSE and event.forced_by:
                origin_text = f"forced False via {', '.join(event.forced_by)}"
            cost_text = CostTier.label(event.cost) if event.origin is FactOrigin.COMPUTED else ""
            read_text = f"x{event.reads}" if event.reads else ""
            # A trailing space after every padded field (rather than relying on the padding alone)
            # keeps columns visually separated even when a value overflows its nominal width -- a long
            # fact name, or a tuple-valued fact's repr -- instead of the next column running into it.
            lines.append(
                f"    {event.name:<26} {str(event.value):<10} {origin_text:<30} "
                f"{cost_text:<9} {read_text}".rstrip()
            )

        lines.append("  rules:")
        for rule_event in self.rule_events:
            verdict = "FIRED" if rule_event.fired else "FALSE"
            lines.append(f"    [{rule_event.rank:>3}] {rule_event.term:<34} {verdict:<6} {rule_event.reason}")

        if self.prune_events:
            pruned_text = "; ".join(
                f"{event.dropped_term} <- subsumed by {event.subsumed_by}" for event in self.prune_events
            )
            lines.append(f"  pruned: {pruned_text}")
        else:
            lines.append("  pruned: (none)")

        emitted_text = ", ".join(self.emitted_terms) if self.emitted_terms else "(none)"
        lines.append(f"  emitted: {emitted_text}")

        return "\n".join(lines)


class _SuppressAnnotatorDebugFilter(logging.Filter):
    """Suppress annotator debug records when the caller renders traces."""

    def filter(self, record: logging.LogRecord) -> bool:
        """Return whether a record should be emitted."""
        return not (record.levelno == logging.DEBUG and record.name.startswith("ensembl.vep.annotators"))


@contextmanager
def annotation_tracing(*, quiet: bool) -> Iterator[None]:
    """Enable annotator tracing for the duration of a context.

    Tracing uses the annotator logger level instead of adding a second flag to fact evaluation.
    Quiet mode suppresses the annotator's debug rendering when the caller renders the same trace.

    Args:
        quiet (bool): Whether to suppress annotator debug records.
    """
    annotators_logger = logging.getLogger("ensembl.vep.annotators")
    previous_level = annotators_logger.level
    annotators_logger.setLevel(logging.DEBUG)

    # The filter goes on the handlers, not the logger: a logger's filters apply only to records
    # logged through that logger itself, never to its children's.
    handlers = logging.getLogger("ensembl.vep").handlers
    log_filter = _SuppressAnnotatorDebugFilter()
    if quiet:
        for handler in handlers:
            handler.addFilter(log_filter)
    try:
        yield
    finally:
        for handler in handlers:
            handler.removeFilter(log_filter)
        annotators_logger.setLevel(previous_level)


def describe_overlap(overlap: object) -> str:
    """Return a one-line label for a transcript and variant overlap.

    The object is inspected without importing ``TranscriptVariationAllele`` so the trace event
    layer remains independent of ``models``. The allele half is its ``display_name``, so a trace
    header names the variant the way every other command names it.

    Args:
        overlap (object): Object exposing transcript and overlapping-allele attributes.

    Returns:
        str: Label containing the feature and the allele's name when available.
    """
    feature = getattr(overlap, "feature", None)
    feature_label = getattr(feature, "stable_id", None) or repr(feature)
    allele = getattr(overlap, "overlapping_feature", None)
    return f"{feature_label} / {getattr(allele, 'display_name', None) or UNKNOWN_VALUE}"
