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
"""The finished derivation of one overlap, and how it is rendered.

An ``AnnotationTrace`` is immutable and self-describing: it names the overlap and the annotator it
came from, so a consumer needs nothing but the trace to render, diff or serialise it. It is never
built on the annotation path -- ``trace/reconstruct.py`` assembles one afterwards, from the memo
table the run left behind plus the little ``TraceRecorder`` gathered alongside it.

Fields a mode could not fill are empty tuples rather than None, so a renderer iterates and draws
nothing instead of branching.
"""

from __future__ import annotations

import logging
from collections.abc import Iterator
from contextlib import contextmanager

from attrs import field, frozen

from ensembl.vep.common import CostTier
from ensembl.vep.trace.events import FactEvent, FactOrigin, PruneEvent, RuleEvent, TraceRecorder

__all__ = [
    "AnnotationTrace",
    "describe_overlap",
    "format_explanation",
    "maybe_start_recorder",
    "tracing_enabled",
]


@frozen(kw_only=True, auto_attribs=True)
class AnnotationTrace:
    """One overlap's derivation: facts observed, rule verdicts, pruning, and the result.

    Attributes:
        fact_events (tuple[FactEvent, ...]): Every fact the derivation touched, memoised ones in
            first-materialisation order followed by the facts ``requires`` propagation proved
            False without computing.
        rule_events (tuple[RuleEvent, ...]): Every rule's verdict, in rank order.
        prune_events (tuple[PruneEvent, ...]): Every subsumption decision.
        emitted_terms (tuple[str, ...]): The final, post-pruning term list.
        annotator_name (str): The annotator's registry name, e.g. "consequence.linear".
        overlap_label (str): The one-line label of the overlap this explains (describe_overlap).
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
        """The declared cost this derivation actually charged."""
        return sum(event.cost for event in self.fact_events if event.origin is FactOrigin.COMPUTED)

    @property
    def avoided_cost(self) -> int:
        """The declared cost the caches and ``requires`` propagation kept it from charging."""
        avoided = 0
        for event in self.fact_events:
            if event.origin is FactOrigin.FORCED_FALSE:
                avoided += event.cost
        return avoided


def maybe_start_recorder(logger: logging.Logger) -> TraceRecorder | None:
    """Return a fresh TraceRecorder if ``logger`` is enabled for DEBUG, else None.

    The single gate every annotator's ``annotate_feature_overlap`` calls through, so the decision
    to trace is made once per overlap rather than per fact. Everything downstream of it is a plain
    ``is not None`` check.

    Args:
        logger (logging.Logger): The calling annotator module's own logger.

    Returns:
        TraceRecorder | None: A new recorder if DEBUG is enabled for ``logger``, else None.
    """
    return TraceRecorder() if logger.isEnabledFor(logging.DEBUG) else None


class _DropAnnotatorDebug(logging.Filter):
    """Drops the annotators' own DEBUG records while a caller renders the same traces itself."""

    def filter(self, record: logging.LogRecord) -> bool:
        """Return whether a record should be emitted."""
        return not (record.levelno == logging.DEBUG and record.name.startswith("ensembl.vep.annotators"))


@contextmanager
def tracing_enabled(*, quiet: bool) -> Iterator[None]:
    """Turn the annotators' tracing switch on for the duration of the block.

    The switch is the annotator module's own logger level, so this adds no second flag to the hot
    path. A caller that will render the traces itself passes ``quiet`` to keep the annotator's own
    per-overlap explanation line out of the log; printing the same derivation twice helps nobody.

    Args:
        quiet (bool): True to suppress the annotators' DEBUG records for the duration.
    """
    annotators_logger = logging.getLogger("ensembl.vep.annotators")
    previous_level = annotators_logger.level
    annotators_logger.setLevel(logging.DEBUG)

    # The filter goes on the handlers, not the logger: a logger's filters apply only to records
    # logged through that logger itself, never to its children's.
    handlers = logging.getLogger("ensembl.vep").handlers
    log_filter = _DropAnnotatorDebug()
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
    """Format a one-line label for a (transcript, variant allele) overlap.

    Duck-typed rather than importing ``TranscriptVariationAllele``, so this package stays free of a
    dependency on ``models`` and the fallbacks degrade for any object shaped closely enough.

    Args:
        overlap (object): Typically a TranscriptVariationAllele, exposing ``.feature`` (with
            ``.stable_id``) and ``.overlapping_feature`` (with ``.slice.seq_region_name``,
            ``.seq_region_start``, ``.reference_allele``, ``.alternate_allele``).

    Returns:
        str: ``"<feature stable_id> / <seq_region>:<position> <ref>><alt>"``, degrading field by
            field when an attribute is missing.
    """
    feature = getattr(overlap, "feature", None)
    feature_label = getattr(feature, "stable_id", None) or repr(feature)

    allele = getattr(overlap, "overlapping_feature", None)
    if allele is None:
        return f"{feature_label} / ?"

    allele_slice = getattr(allele, "slice", None)
    seq_region_name = getattr(allele_slice, "seq_region_name", None) if allele_slice is not None else None
    position = getattr(allele, "seq_region_start", None)
    if position is None:
        position = getattr(allele, "start", None)
    locus = f"{seq_region_name}:{position}" if seq_region_name is not None else f"pos {position}"

    reference_allele = getattr(allele, "reference_allele", "") or ""
    alternate_allele = getattr(allele, "alternate_allele", "") or ""
    allele_label = (
        f"{locus} {reference_allele}>{alternate_allele}" if (reference_allele or alternate_allele) else locus
    )

    return f"{feature_label} / {allele_label}"


def format_explanation(trace: AnnotationTrace) -> str:
    """Render a trace as a human-readable derivation.

    Args:
        trace (AnnotationTrace): The trace to render.

    Returns:
        str: A header, the facts observed (a summary line then one line each), every rule's
            verdict, every pruning decision, and the final emitted term list.
    """
    lines = [f"{trace.overlap_label}          annotator: {trace.annotator_name}"]

    computed = trace.events_by_origin(FactOrigin.COMPUTED)
    prefilled = trace.events_by_origin(FactOrigin.PREFILLED)
    cached = trace.events_by_origin(FactOrigin.CACHED)
    forced_false = trace.events_by_origin(FactOrigin.FORCED_FALSE)
    reads = sum(event.reads for event in trace.fact_events)

    lines.append(
        f"  facts: {len(computed)} computed (cost {trace.computed_cost}), {len(prefilled)} prefilled, "
        f"{len(forced_false)} forced False (cost {trace.avoided_cost} avoided), "
        f"{len(cached)} cache hits, {reads} reads"
    )

    for event in trace.fact_events:
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
    for rule_event in trace.rule_events:
        verdict = "FIRED" if rule_event.fired else "FALSE"
        lines.append(f"    [{rule_event.rank:>3}] {rule_event.term:<34} {verdict:<6} {rule_event.reason}")

    if trace.prune_events:
        pruned_text = "; ".join(
            f"{event.dropped_term} <- subsumed by {event.subsumed_by}" for event in trace.prune_events
        )
        lines.append(f"  pruned: {pruned_text}")
    else:
        lines.append("  pruned: (none)")

    emitted_text = ", ".join(trace.emitted_terms) if trace.emitted_terms else "(none)"
    lines.append(f"  emitted: {emitted_text}")

    return "\n".join(lines)
