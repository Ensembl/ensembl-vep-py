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
"""Measuring an annotator: the conformance gate that runs before any timing, and the timing itself.

The gate comes first on purpose. A strategy that is fast because it is wrong is not a result, so a
benchmark refuses to time anything until every candidate agrees with the reference on every
overlap. The gate returns a report rather than exiting: the runner decides what a mismatch means for
the process.
"""

from __future__ import annotations

import time
from collections.abc import Callable, Mapping, Sequence

from attrs import field, frozen

from ensembl.vep.annotators.base import Annotator
from ensembl.vep.annotators.registry import available_annotators
from ensembl.vep.models.overlap import OverlapConsequence, TranscriptVariationAllele
from ensembl.vep.stats.profile import count_fact_computations

__all__ = [
    "KERNELS_STRATEGY",
    "MISMATCH_LIMIT",
    "REFERENCE_STRATEGY",
    "GateReport",
    "Timing",
    "TimingPolicy",
    "annotator_extras",
    "conformance_gate",
    "consequence_strategies",
    "fact_cost",
    "native_kernels_available",
    "terms_of",
    "time_best",
]

REFERENCE_STRATEGY = "consequence.linear"
"""The strategy every other one must agree with."""
KERNELS_STRATEGY = "consequence.columnar-kernels"
"""The strategy that needs the compiled extension."""
MISMATCH_LIMIT = 20
"""Mismatches collected before the gate stops, so one broken strategy cannot flood the report."""

_ANNOTATOR_COUNTERS = (
    ("cache_hits", "candidate cache hits"),
    ("cache_misses", "candidate cache misses"),
    ("settled_early", "settled early"),
    ("needed_expensive", "needed expensive tier"),
    ("skipped_by_subsumption", "skipped by subsumption"),
)
_CACHE_COUNTERS = (
    (
        "zone_map_cache",
        "zone map",
        (
            ("maps_built", "maps built"),
            ("hit_count", "hits"),
            ("eviction_count", "evictions"),
            ("build_seconds", "build seconds"),
            ("facts_served", "facts served"),
        ),
    ),
    ("codon_cache", "codon cache", (("hit_count", "hits"), ("miss_count", "misses"))),
    ("mapper_index_cache", "mapper cache", (("hit_count", "hits"), ("miss_count", "misses"))),
)


def terms_of(consequences: Sequence[OverlapConsequence]) -> list[str]:
    """Return the terms of a consequence list, in order."""
    return [consequence.term for consequence in consequences]


def consequence_strategies() -> list[str]:
    """Return every registered ``consequence.*`` strategy, the reference first."""
    names = [name for name in available_annotators() if name.startswith("consequence.")]
    if REFERENCE_STRATEGY in names:
        names.remove(REFERENCE_STRATEGY)
        names.insert(0, REFERENCE_STRATEGY)
    return names


def native_kernels_available() -> bool:
    """Return whether the strategy that needs the compiled extension is registered."""
    return KERNELS_STRATEGY in available_annotators()


@frozen(kw_only=True, auto_attribs=True)
class GateReport:
    """What the conformance gate found.

    Attributes:
        reference (str): The strategy the candidates were checked against.
        candidates (tuple[str, ...]): The strategies checked.
        overlaps (int): Overlaps checked.
        mismatches (tuple[str, ...]): One line per disagreement found, at most ``MISMATCH_LIMIT``.
        truncated (bool): Whether the gate stopped early at the limit.
    """

    reference: str
    candidates: tuple[str, ...]
    overlaps: int
    mismatches: tuple[str, ...] = field(default=())
    truncated: bool = field(default=False)

    @property
    def ok(self) -> bool:
        """Return whether every candidate agreed on every overlap."""
        return not self.mismatches


def conformance_gate(
    overlaps: Sequence[TranscriptVariationAllele],
    *,
    reference: Annotator[OverlapConsequence],
    candidates: Mapping[str, Annotator[OverlapConsequence]],
    reference_name: str = REFERENCE_STRATEGY,
    limit: int = MISMATCH_LIMIT,
) -> GateReport:
    """Check every candidate against the reference, overlap by overlap.

    Args:
        overlaps (Sequence[TranscriptVariationAllele]): The overlaps to check.
        reference (Annotator[OverlapConsequence]): The annotator whose terms are right by definition.
        candidates (Mapping[str, Annotator[OverlapConsequence]]): Name -> annotator to check.
        reference_name (str): The reference's name, for the report.
        limit (int): Stop after this many mismatches.

    Returns:
        GateReport: What was found.
    """
    mismatches: list[str] = []
    for overlap in overlaps:
        expected = terms_of(reference.annotate_feature_overlap(overlap))
        for name, candidate in candidates.items():
            actual = terms_of(candidate.annotate_feature_overlap(overlap))
            if actual == expected:
                continue
            mismatches.append(
                f"{name} disagrees with {reference_name} on {overlap.feature.stable_id} at "
                f"{overlap.overlapping_feature.start}: {actual} != {expected}"
            )
            if len(mismatches) >= limit:
                return GateReport(
                    reference=reference_name,
                    candidates=tuple(candidates),
                    overlaps=len(overlaps),
                    mismatches=tuple(mismatches),
                    truncated=True,
                )
    return GateReport(
        reference=reference_name,
        candidates=tuple(candidates),
        overlaps=len(overlaps),
        mismatches=tuple(mismatches),
    )


@frozen(kw_only=True, auto_attribs=True)
class TimingPolicy:
    """How a measurement is taken.

    Attributes:
        repeats (int): Timed samples; the best is reported, as the least noisy.
        warmup (int): Untimed passes first, so the timed ones do not include warming caches.
    """

    repeats: int = field(default=3)
    warmup: int = field(default=1)


@frozen(kw_only=True, auto_attribs=True)
class Timing:
    """The samples one measurement took.

    Attributes:
        samples (tuple[float, ...]): Elapsed seconds per timed pass.
        warmup (int): Untimed passes that preceded them.
    """

    samples: tuple[float, ...]
    warmup: int = field(default=0)

    @property
    def best(self) -> float:
        """Return the fastest sample in seconds."""
        return min(self.samples)

    def ms_per_thousand(self, contexts: int) -> float:
        """Return the best sample as milliseconds per thousand contexts."""
        return self.best * 1_000_000 / max(1, contexts)


def time_best(run_once: Callable[[], object], policy: TimingPolicy) -> Timing:
    """Run a workload the policy's warmup times untimed, then time it the policy's repeats times."""
    for _ in range(policy.warmup):
        run_once()
    samples = []
    for _ in range(max(1, policy.repeats)):
        started = time.perf_counter()
        run_once()
        samples.append(time.perf_counter() - started)
    return Timing(samples=tuple(samples), warmup=policy.warmup)


def fact_cost(
    annotator: Annotator[OverlapConsequence], overlaps: Sequence[TranscriptVariationAllele]
) -> tuple[float, float]:
    """Return facts computed and declared cost charged, each per context, over one pass.

    Instruments every fact process-wide for the pass; see ``stats.profile.count_fact_computations``.
    """
    contexts = max(1, len(overlaps))
    with count_fact_computations() as counters:
        for overlap in overlaps:
            annotator.annotate_feature_overlap(overlap)
    return counters["computed"] / contexts, counters["cost"] / contexts


def annotator_extras(annotator: object) -> dict[str, object]:
    """Return whatever instrumentation an annotator and its caches expose, labelled.

    The counters live where the work is done: the staged strategies count on the annotator, while
    the zone map, codon and mapper caches count on themselves, reached through the annotator's
    context.
    """
    extras: dict[str, object] = {}
    for attribute, label in _ANNOTATOR_COUNTERS:
        if hasattr(annotator, attribute):
            extras[label] = getattr(annotator, attribute)
    context = getattr(annotator, "annotator_context", None)
    if context is None:
        return extras
    for cache_attribute, prefix, counters in _CACHE_COUNTERS:
        cache = getattr(context, cache_attribute, None)
        if cache is None:
            continue
        for attribute, label in counters:
            if hasattr(cache, attribute):
                extras[f"{prefix} {label}"] = getattr(cache, attribute)
    return extras
