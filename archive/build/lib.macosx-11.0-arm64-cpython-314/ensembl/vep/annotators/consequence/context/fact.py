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
"""FactContext: the per-overlap memoisation unit the fact layer is evaluated over."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import TYPE_CHECKING, Generic, TypeVar

from attrs import define, field

from ensembl.vep.trace.events import FactOrigin, TraceRecorder

if TYPE_CHECKING:
    from ensembl.vep.annotators.consequence.context.annotator import AnnotatorContext
    from ensembl.vep.models.genic import Transcript
    from ensembl.vep.models.overlap import FeatureOverlap, TranscriptVariationAllele
    from ensembl.vep.models.variation import VariationAllele

__all__ = ["FactContext", "OverlapT"]

OverlapT = TypeVar("OverlapT", bound="FeatureOverlap", covariant=True)
"""The kind of overlap a context carries.

We need a covariant here, so that a ``FactContext[TranscriptVariationAllele]`` can be passed onto a fact
declared over ``FactContext[FeatureOverlap]``.
"""


@define(kw_only=True, auto_attribs=True)
class FactContext(Generic[OverlapT]):
    """One (variant, feature) overlap, the unit a fact rule set is evaluated over.

    Wraps a overlap object and holds relevant context such asthe memo table for computed facts
    and other memoised values.

    Attributes:
        overlap (OverlapT): The (variant, feature) overlap being annotated. OverlapT is a covariant type
            generic and therefore, it is recommended that overlap is not mutable here but not enforced here
            with ``attrs.setters.frozen`` because of cost (~5x).
        recorder (TraceRecorder | None): Opt-in trace gatherer, or None to disable tracing
            entirely. It collects only what a later reconstruction cannot recover -- read counts
            and which cache served which fact -- so the whole cost of tracing on this class is one
            ``is not None`` check per fact lookup. See ensembl.vep.trace.
        annotator_context (AnnotatorContext | None): The run-scoped state this overlap is being
            annotated under or None when the context is evaluated standalone (e.g., unit tests).
    """

    overlap: OverlapT
    recorder: TraceRecorder | None = field(default=None)
    annotator_context: AnnotatorContext | None = field(default=None)
    _memo: dict[str, object] = field(factory=dict, init=False, repr=False)
    _variant_span: tuple[int, int] | None = field(default=None, init=False, repr=False)
    _biotype: str | None = field(default=None, init=False, repr=False)

    @property
    def memo(self) -> dict[str, object]:
        """Public accessor for the fact memo table.

        Returns:
            dict[str, object]: The mutable memo table mapping fact name to computed value.
        """
        return self._memo

    @property
    def variant_span(self) -> tuple[int, int]:
        """The variant's genomic ``(min, max)`` bounds (memoised).

        An insertion is represented with ``end == start - 1``, so ``start`` alone is not reliably the
        lower bound, hence the min/max rather than a direct read.

        Returns:
            tuple[int, int]: variant genomic (min, max) span.
        """
        span = self._variant_span
        if span is None:
            variant = self.overlap.overlapping_feature
            span = (min(variant.start, variant.end), max(variant.start, variant.end))
            self._variant_span = span
        return span

    @property
    def biotype(self) -> str | None:
        """The overlapped feature's biotype (memoised).

        If transcript/feature carries no ``biotype`` of its own in this model then falls back
        to gene. Resolved once per overlap

        Returns:
            str | None: The biotype, or None if neither the feature nor its gene declares one.
        """
        if self._biotype is None:
            feature = self.overlap.feature
            biotype = getattr(feature, "biotype", None)
            if biotype is None:
                gene = getattr(feature, "gene", None)
                if gene is not None:
                    biotype = gene.biotype
            self._biotype = biotype
        return self._biotype

    @classmethod
    def for_overlaps(cls, overlaps: Iterable[OverlapT]) -> list[FactContext[OverlapT]]:
        """Return one context per overlap.

        What every caller that evaluates a whole corpus needs -- profiling, validation, benchmark
        gates -- rather than each rebuilding the same comprehension.
        """
        return [cls(overlap=overlap) for overlap in overlaps]

    def reset(self) -> None:
        """Clear the memo table and any other memoised state.

        Lets the context be reused for a new overlap.
        """
        self._memo.clear()
        self._variant_span = None
        self._biotype = None

    def prefill(self, values: Mapping[str, object], *, origin: FactOrigin = FactOrigin.PREFILLED) -> None:
        """Seed the memo table with already-known values, without charging their cost.

        Useful for an annotator that has derived some facts by cheaper means than calling their
        registered ``compute`` function, and wants later fact lookups to see those values as already memoised.

        Args:
            values (Mapping[str, object]): Fact name -> value pairs to seed the memo table with.
            origin (FactOrigin): Which uncharged source these values came from, for example
                FactOrigin.PREFILLED for a zone map or FactOrigin.CACHED for the shared-codon
                cache. Ignored entirely when self.recorder is None.
        """
        # Noted before the update, so only the names this call genuinely introduced are credited to
        # the cache. A codon-cache miss computes both facts through get_fact and then prefills the
        # same values back; those were computed, and the memo check is what keeps them so.
        if self.recorder is not None:
            self.recorder.note_served((name for name in values if name not in self._memo), origin)
        self._memo.update(values)
