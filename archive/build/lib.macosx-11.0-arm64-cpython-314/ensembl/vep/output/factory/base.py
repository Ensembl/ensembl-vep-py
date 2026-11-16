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
"""The format-neutral result of annotating one variant."""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from typing import Any

from attrs import define, field, frozen

from ensembl.vep.stats import RunSummary, RunSummaryCollector

__all__ = ["AnnotatedVariant", "AnnotatedVariants", "AnnotationRecord"]

AnnotationRecord = dict[str, object]
"""One VEP output hash, copied from a selected allele-feature annotation."""


def _as_iterator(variants: Iterable[AnnotatedVariant]) -> Iterator[AnnotatedVariant]:
    """Return a single-pass iterator over annotated variants."""
    return iter(variants)


@frozen(kw_only=True, auto_attribs=True)
class AnnotatedVariant:
    """Carry one annotated variant independently of its eventual output format.

    Attributes:
        id (str | None): The uploaded or stable variant identifier.
        seq_region_name (str): The sequence-region name.
        start (int): The inclusive genomic start.
        end (int): The inclusive genomic end.
        strand (int): The genomic strand.
        allele_string (str): The reference and alternate alleles separated by ``/``.
        input_columns (tuple[str, ...]): The original input columns, when present.
        annotations (tuple[AnnotationRecord, ...]): Selected VEP output hashes.
        most_severe_consequence (str | None): The most severe emitted consequence term.
        variant_class (str | None): The Sequence Ontology class term.
        assembly_name (str | None): The assembly name, when known.
        skipped (bool): Whether annotation was skipped for this input variant.
    """

    id: str | None
    seq_region_name: str
    start: int
    end: int
    strand: int
    allele_string: str
    input_columns: tuple[str, ...] = ()
    annotations: tuple[AnnotationRecord, ...] = ()
    most_severe_consequence: str | None = None
    variant_class: str | None = None
    assembly_name: str | None = None
    skipped: bool = False


@define(kw_only=True, auto_attribs=True)
class AnnotatedVariants(Iterator[AnnotatedVariant]):
    """Expose a closable annotation iterator and its completed run summary.

    Attributes:
        variants (Iterable[AnnotatedVariant]): Ordered annotated variants.
        collector (RunSummaryCollector): Summary state updated as variants become available.
    """

    _iterator: Iterator[AnnotatedVariant] = field(alias="variants", converter=_as_iterator)
    _collector: RunSummaryCollector = field(alias="collector", repr=False)
    _completed: bool = field(default=False, init=False)
    _closed: bool = field(default=False, init=False)
    _variant_count: int = field(default=0, init=False)

    def __iter__(self) -> AnnotatedVariants:
        """Return this single-pass annotation stream."""
        return self

    def __next__(self) -> AnnotatedVariant:
        """Return the next annotated variant in input order."""
        if self._closed:
            raise StopIteration
        try:
            variant = next(self._iterator)
        except StopIteration:
            self._completed = True
            self.close()
            raise
        self._variant_count += 1
        return variant

    def __enter__(self) -> AnnotatedVariants:
        """Return this stream for bounded consumption."""
        return self

    def __exit__(self, *_error: Any) -> None:
        """Close worker resources when consumption stops."""
        self.close()

    @property
    def completed(self) -> bool:
        """Return whether every input variant has been consumed."""
        return self._completed

    @property
    def variant_count(self) -> int:
        """Return the number of variants yielded to the consumer."""
        return self._variant_count

    @property
    def collector(self) -> RunSummaryCollector:
        """Return the collector after complete consumption.

        Raises:
            RuntimeError: If the stream has not been exhausted.
        """
        if not self._completed:
            raise RuntimeError("the run summary is available after the annotation stream is exhausted")
        return self._collector

    @property
    def summary(self) -> RunSummary:
        """Return the final summary after complete consumption.

        Raises:
            RuntimeError: If the stream has not been exhausted.
        """
        if not self._completed:
            raise RuntimeError("the run summary is available after the annotation stream is exhausted")
        return self._collector.finish()

    def close(self) -> None:
        """Release resources held by the underlying annotation iterator."""
        if self._closed:
            return
        close = getattr(self._iterator, "close", None)
        if close is not None:
            close()
        self._closed = True
