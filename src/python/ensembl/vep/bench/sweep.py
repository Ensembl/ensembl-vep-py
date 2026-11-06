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
"""Sweeping one numeric annotation setting and timing each value.

Every value is sanity-checked against the reference strategy before it is timed: a setting that
changes the answer is a bug, not a speed-up. The chunking follows the strategy. A Python batch
strategy chunks internally by its own ``batch_size``, so it is handed the whole corpus at once; the
native kernels chunk nowhere else, so for them the swept batch size is the chunk handed over.
"""

from __future__ import annotations

from collections.abc import Sequence

from attrs import field, frozen

from ensembl.vep.annotators.base import Annotator
from ensembl.vep.annotators.factory import build_annotator
from ensembl.vep.bench.harness import KERNELS_STRATEGY, REFERENCE_STRATEGY, TimingPolicy, terms_of, time_best
from ensembl.vep.caches.zone_map import DEFAULT_ZONE_MAP_BYTES
from ensembl.vep.config.sections.bench.sweep import SWEEPABLE_KNOBS
from ensembl.vep.models.overlap import OverlapConsequence, TranscriptVariationAllele

__all__ = [
    "DEFAULT_BATCH_STRATEGIES",
    "SANITY_CHECK_COUNT",
    "SWEEPABLE_KNOBS",
    "Sweep",
    "SweepPoint",
    "sweep_strategy",
]

SANITY_CHECK_COUNT = 64
"""Overlaps each swept value is checked on before it is timed."""
DEFAULT_BATCH_STRATEGIES = ("consequence.columnar-staged", KERNELS_STRATEGY)

_DEFAULT_BATCH_SIZE = 256
_DEFAULT_CODON_CACHE_ENTRIES = 65536
_DEFAULT_MAPPER_CACHE_BYTES = 64 * 1024 * 1024


@frozen(kw_only=True, auto_attribs=True)
class Sweep:
    """A setting and the values to try.

    Attributes:
        knob (str): One of ``SWEEPABLE_KNOBS``.
        values (tuple[int, ...]): The values, in order.
    """

    knob: str
    values: tuple[int, ...] = field(default=())

    def __attrs_post_init__(self) -> None:
        """Reject a setting that cannot be swept."""
        if self.knob not in SWEEPABLE_KNOBS:
            raise ValueError(f"cannot sweep {self.knob!r}; expected one of {', '.join(SWEEPABLE_KNOBS)}")


@frozen(kw_only=True, auto_attribs=True)
class SweepPoint:
    """One strategy at one value of the swept setting.

    Attributes:
        corpus (str): The corpus measured over.
        strategy (str): The strategy.
        knob (str): The swept setting.
        value (int): Its value for this point.
        ms_per_1k (float): Best wall clock per thousand overlaps.
    """

    corpus: str
    strategy: str
    knob: str
    value: int
    ms_per_1k: float

    @property
    def configuration(self) -> str:
        """Return the point's name in a results file."""
        return f"{self.strategy} {self.knob}={self.value}"

    def to_json(self) -> dict[str, object]:
        """Return the point as plain JSON-ready data."""
        return {
            "corpus": self.corpus,
            "configuration": self.configuration,
            "strategy": self.strategy,
            "knob": self.knob,
            "value": self.value,
            "ms_per_1k": self.ms_per_1k,
        }


def _annotator_at(strategy: str, knob: str, value: int) -> Annotator[OverlapConsequence]:
    """Build a strategy with one setting at the swept value, and the cache it needs switched on."""
    return build_annotator(
        annotator=strategy,
        batch_size=value if knob == "batch_size" else _DEFAULT_BATCH_SIZE,
        zone_map=knob == "zone_map_bytes",
        zone_map_bytes=value if knob == "zone_map_bytes" else DEFAULT_ZONE_MAP_BYTES,
        codon_cache=knob == "codon_cache_entries",
        codon_cache_entries=value if knob == "codon_cache_entries" else _DEFAULT_CODON_CACHE_ENTRIES,
        mapper_cache=knob == "mapper_cache_bytes",
        mapper_cache_bytes=value if knob == "mapper_cache_bytes" else _DEFAULT_MAPPER_CACHE_BYTES,
    )


def sweep_strategy(
    corpus: str,
    strategy: str,
    overlaps: Sequence[TranscriptVariationAllele],
    sweep: Sweep,
    policy: TimingPolicy,
) -> tuple[list[SweepPoint], list[str]]:
    """Time one strategy at every value of a setting, checking each value first.

    Args:
        corpus (str): The corpus's label.
        strategy (str): The strategy to sweep.
        overlaps (Sequence[TranscriptVariationAllele]): The corpus.
        sweep (Sweep): The setting and its values.
        policy (TimingPolicy): How to time.

    Returns:
        tuple[list[SweepPoint], list[str]]: The timed points, and one line per value that changed
            the answer -- which is not timed.
    """
    checked = list(overlaps[:SANITY_CHECK_COUNT])
    reference = build_annotator(annotator=REFERENCE_STRATEGY)
    expected = [terms_of(reference.annotate_feature_overlap(overlap)) for overlap in checked]

    points: list[SweepPoint] = []
    failures: list[str] = []
    for value in sweep.values:
        annotator = _annotator_at(strategy, sweep.knob, value)
        actual = [terms_of(consequences) for consequences in annotator.annotate_feature_overlaps(checked)]
        if actual != expected:
            failures.append(
                f"{strategy} at {sweep.knob}={value} disagrees with {REFERENCE_STRATEGY} "
                f"on the first {len(checked)} overlaps"
            )
            continue
        chunk_size = (
            value if sweep.knob == "batch_size" and strategy == KERNELS_STRATEGY else max(1, len(overlaps))
        )

        def run_once(annotator: Annotator[OverlapConsequence] = annotator, size: int = chunk_size) -> None:
            for start in range(0, len(overlaps), size):
                annotator.annotate_feature_overlaps(list(overlaps[start : start + size]))

        timing = time_best(run_once, policy)
        points.append(
            SweepPoint(
                corpus=corpus,
                strategy=strategy,
                knob=sweep.knob,
                value=value,
                ms_per_1k=timing.ms_per_thousand(len(overlaps)),
            )
        )
    return points, failures
