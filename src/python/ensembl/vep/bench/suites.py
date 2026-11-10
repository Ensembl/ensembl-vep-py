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
"""The strategy x cache matrix: every configuration worth timing, and timing them.

A registered name is a *strategy*; the caches are orthogonal to it. Timing only registered names
would leave the configurations people actually run unmeasured, so a configuration is a strategy
plus a cache set, and the set is chosen on the command line rather than hardcoded.
"""

from __future__ import annotations

from collections.abc import Sequence
from functools import partial

from attrs import field, frozen

from ensembl.vep.annotators.base import Annotator
from ensembl.vep.annotators.factory import build_annotator
from ensembl.vep.bench.harness import TimingPolicy, annotator_extras, fact_cost, time_best
from ensembl.vep.models.overlap import OverlapConsequence, TranscriptVariationAllele

__all__ = [
    "CACHES",
    "NO_CACHES",
    "Configuration",
    "MatrixRow",
    "build_configurations",
    "parse_cache_set",
    "run_matrix",
]

CACHES = ("zonemap", "codon", "mapper")
"""The caches a configuration can switch on, in the order a label names them."""
NO_CACHES = "none"


def parse_cache_set(text: str) -> frozenset[str]:
    """Parse ``none`` or a ``+``-joined set of cache names.

    Raises:
        ValueError: If a name is not one of ``CACHES``.
    """
    names = [part.strip() for part in text.split("+") if part.strip()]
    if names == [NO_CACHES]:
        return frozenset()
    unknown = sorted(set(names) - set(CACHES))
    if unknown or not names:
        raise ValueError(
            f"unknown cache set {text!r}: expected {NO_CACHES!r} or a +-joined set of {', '.join(CACHES)}"
        )
    return frozenset(names)


@frozen(kw_only=True, auto_attribs=True)
class Configuration:
    """One strategy with one set of caches.

    Attributes:
        strategy (str): The registered strategy.
        caches (frozenset[str]): The caches switched on.
        batch_size (int): Chunk size, for the strategies that batch.
    """

    strategy: str
    caches: frozenset[str] = field(factory=frozenset)
    batch_size: int = field(default=256)

    @property
    def label(self) -> str:
        """Return the configuration's table label: the strategy, then its caches."""
        if not self.caches:
            return self.strategy
        return f"{self.strategy} +{'+'.join(name for name in CACHES if name in self.caches)}"

    def build(self) -> Annotator[OverlapConsequence]:
        """Return a fresh annotator with fresh caches, so no configuration is timed warm from another."""
        return build_annotator(
            annotator=self.strategy,
            batch_size=self.batch_size,
            zone_map="zonemap" in self.caches,
            codon_cache="codon" in self.caches,
            mapper_cache="mapper" in self.caches,
        )


def build_configurations(
    strategies: Sequence[str], cache_sets: Sequence[str], batch_size: int
) -> list[Configuration]:
    """Cross strategies with cache sets, strategies outermost, without repeats.

    Raises:
        ValueError: If a cache set does not parse.
    """
    parsed = [parse_cache_set(text) for text in cache_sets]
    configurations: list[Configuration] = []
    for strategy in strategies:
        for caches in parsed:
            configuration = Configuration(strategy=strategy, caches=caches, batch_size=batch_size)
            if configuration not in configurations:
                configurations.append(configuration)
    return configurations


@frozen(kw_only=True, auto_attribs=True)
class MatrixRow:
    """One configuration's measurements over one corpus.

    Attributes:
        corpus (str): The corpus measured over.
        configuration (str): The configuration's label.
        contexts (int): Overlaps in the corpus.
        ms_per_1k (float): Best wall clock per thousand overlaps, one overlap at a time.
        batch_ms_per_1k (float): The same, in chunks of the configuration's batch size.
        facts_per_context (float): Facts computed per overlap.
        cost_per_context (float): Declared cost charged per overlap -- the number that survives a
            port to another language, where wall clock does not.
        extras (dict[str, object]): Whatever instrumentation the annotator and caches expose.
    """

    corpus: str
    configuration: str
    contexts: int
    ms_per_1k: float
    batch_ms_per_1k: float
    facts_per_context: float
    cost_per_context: float
    extras: dict[str, object] = field(factory=dict)

    def to_json(self) -> dict[str, object]:
        """Return the row as plain JSON-ready data."""
        return {
            "corpus": self.corpus,
            "configuration": self.configuration,
            "contexts": self.contexts,
            "ms_per_1k": self.ms_per_1k,
            "batch_ms_per_1k": self.batch_ms_per_1k,
            "facts_per_context": self.facts_per_context,
            "cost_per_context": self.cost_per_context,
            "extras": dict(self.extras),
        }


def _run_single(
    annotator: Annotator[OverlapConsequence], overlaps: Sequence[TranscriptVariationAllele]
) -> None:
    """Annotate every overlap one at a time."""
    for overlap in overlaps:
        annotator.annotate_feature_overlap(overlap)


def _run_chunked(
    annotator: Annotator[OverlapConsequence], overlaps: Sequence[TranscriptVariationAllele], chunk_size: int
) -> None:
    """Annotate every overlap in chunks."""
    for start in range(0, len(overlaps), chunk_size):
        annotator.annotate_feature_overlaps(list(overlaps[start : start + chunk_size]))


def run_matrix(
    corpus: str,
    overlaps: Sequence[TranscriptVariationAllele],
    configurations: Sequence[Configuration],
    policy: TimingPolicy,
) -> list[MatrixRow]:
    """Measure every configuration over one corpus.

    Each measurement builds its own annotator: one for the single-overlap timing, a clean one for
    the fact count (instrumenting needs its own pass), and one for the batch timing.
    """
    rows = []
    for configuration in configurations:
        single = configuration.build()
        single_timing = time_best(partial(_run_single, single, overlaps), policy)

        counted = configuration.build()
        facts_per_context, cost_per_context = fact_cost(counted, overlaps)

        batched = configuration.build()
        size = max(1, configuration.batch_size)
        batch_timing = time_best(partial(_run_chunked, batched, overlaps, size), policy)

        rows.append(
            MatrixRow(
                corpus=corpus,
                configuration=configuration.label,
                contexts=len(overlaps),
                ms_per_1k=single_timing.ms_per_thousand(len(overlaps)),
                batch_ms_per_1k=batch_timing.ms_per_thousand(len(overlaps)),
                facts_per_context=facts_per_context,
                cost_per_context=cost_per_context,
                extras=annotator_extras(counted),
            )
        )
    return rows
