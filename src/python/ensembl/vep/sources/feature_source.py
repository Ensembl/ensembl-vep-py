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
"""The feature source contract.

Reading and model building are two halves. The record getters answer in whatever the source's own
records are -- what the format holds, before anything has been interpreted -- and
``records_to_features`` turns a set of them into models. The feature getters are those two
composed, which is why they are concrete here: a source implements the record half and gets the
feature half. A format that assembles one feature from several records (GFF3, GTF) overrides the
feature getters, because a record at a time is not enough to build anything complete.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Iterator
from typing import Generic, TypeVar
from attr import define, field

from ensembl.vep.caches import SliceCache, WeakSliceCache
from ensembl.vep.common import SourceType
from ensembl.vep.models import Feature
from ensembl.vep.sources.sequence_source import SequenceSource

F = TypeVar("F", bound=Feature)
R = TypeVar("R")
"""One of the source's own records, as its format holds it."""


@define(kw_only=True, auto_attribs=True)
class FeatureSource(ABC, Generic[R, F]):
    """Source of domain feature objects.

    Attributes:
        source_type (SourceType): Kind of backing source.
        sequence_source (SequenceSource | None): Optional reference sequence source.
        slice_cache (SliceCache): Shared genomic slice cache.
    """

    source_type: SourceType
    sequence_source: SequenceSource | None = field(default=None)
    slice_cache: SliceCache = field(factory=WeakSliceCache)

    @abstractmethod
    def get_all_records(self) -> list[R]:
        """Return every unconverted record held by this source."""

        raise NotImplementedError

    @abstractmethod
    def iter_all_records(self) -> Iterator[R]:
        """Yield every unconverted record held by this source."""

        raise NotImplementedError

    @abstractmethod
    def get_records_by_region(
        self,
        region: str | None = None,
        *,
        seq_region_name: str | None = None,
        start: int | None = None,
        end: int | None = None,
        **kwargs: object,
    ) -> list[R]:
        """Return unconverted records overlapping a genomic region.

        Args:
            region (str | None): The region (e.g., "chr1:1000-2000") to get records from.
            seq_region_name (str | None): The name of the region to get records from.
            start (int | None): The start position of the region to get records from.
            end (int | None): The end position of the region to get records from.
            **kwargs (object): What a format's own query accepts, such as the record types a GXF
                query may filter to before a line is parsed.

        Returns:
            list[R]: The records overlapping the specified region.
        """

        raise NotImplementedError

    @abstractmethod
    def iter_records_by_region(
        self,
        region: str | None = None,
        *,
        seq_region_name: str | None = None,
        start: int | None = None,
        end: int | None = None,
        **kwargs: object,
    ) -> Iterator[R]:
        """Yield unconverted records overlapping a genomic region.

        Args:
            region (str | None): The region (e.g., "chr1:1000-2000") to get records from.
            seq_region_name (str | None): The name of the region to get records from.
            start (int | None): The start position of the region to get records from.
            end (int | None): The end position of the region to get records from.
            **kwargs (object): What a format's own query accepts.

        Returns:
            Iterator[R]: The records overlapping the specified region.
        """

        raise NotImplementedError

    @abstractmethod
    def records_to_features(self, records: list[R]) -> list[F]:
        """Return the features described by records from this source.

        Args:
            records (list[R]): The records to convert, in the order the source read them.

        Returns:
            list[F]: The features they describe.
        """

        raise NotImplementedError

    def get_all_features(self) -> list[F]:
        """Return all features from this source."""

        return self.records_to_features(self.get_all_records())

    def iter_all_features(self) -> Iterator[F]:
        """Yield all features from this source.

        One record at a time, which only a format whose record carries a whole feature can do.
        A format that assembles a feature from several records overrides this.
        """

        for record in self.iter_all_records():
            yield from self.records_to_features([record])

    def get_features_by_region(
        self,
        region: str | None = None,
        *,
        seq_region_name: str | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> list[F]:
        """Return features overlapping a genomic region.

        Args:
            region (str | None): The region (e.g., "chr1:1000-2000") to get features from.
            seq_region_name (str | None): The name of the region to get features from.
            start (int | None): The start position of the region to get features from.
            end (int | None): The end position of the region to get features from.

        Returns:
            list[F]: A list of features overlapping the specified region.
        """

        return self.records_to_features(
            self.get_records_by_region(region, seq_region_name=seq_region_name, start=start, end=end)
        )

    def iter_features_by_region(
        self,
        region: str | None = None,
        *,
        seq_region_name: str | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> Iterator[F]:
        """Yield features overlapping a genomic region.

        Args:
            region (str | None): The region (e.g., "chr1:1000-2000") to get features from.
            seq_region_name (str | None): The name of the region to get features from.
            start (int | None): The start position of the region to get features from.
            end (int | None): The end position of the region to get features from.

        Returns:
            Iterator[F]: An iterator over features overlapping the specified region.
        """

        records = self.iter_records_by_region(region, seq_region_name=seq_region_name, start=start, end=end)
        for record in records:
            yield from self.records_to_features([record])

    def record_count(self, *, scan: bool = False) -> int | None:  # pylint: disable=unused-argument
        """Return the declared record count, if available.

        Concrete and None by default: a count is optional metadata, not something every source
        can answer. A source whose format carries one -- an indexed file, a table a query can
        count -- overrides this; the rest say they do not know, and a caller such as progress
        reporting does without rather than paying for a count.

        Args:
            scan (bool): Whether to count the hard way when nothing declares it, which for a
                    file-backed source means reading the whole input. Default is False.

        Returns:
            int | None: The record count, or None when this source cannot answer.
        """

        return None

    @abstractmethod
    def close(self) -> None:
        """Release resources held by this source (e.g., close file handles or DB connections)."""

        raise NotImplementedError

    def __enter__(self) -> FeatureSource[R, F]:
        return self

    def __exit__(self, *exc_info: object) -> None:
        self.close()
