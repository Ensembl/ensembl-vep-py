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
"""Abstract base for all feature sources."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Iterator, Generic, TypeVar
from attr import define

from ensembl.vep.common import SourceType
from ensembl.vep.models import Feature

F = TypeVar("F", bound="Feature")


@define(kw_only=True, auto_attribs=True)
class FeatureSource(ABC, Generic[F]):
    """Unified interface for all sources that produce Feature objects."""

    source_type: SourceType

    @abstractmethod
    def get_all_features(self) -> list[F]:
        """Return all features from this source."""

        raise NotImplementedError

    @abstractmethod
    def iter_all_features(self) -> Iterator[F]:
        """Lazily yield all features from this source."""

        raise NotImplementedError

    @abstractmethod
    def get_features_by_region(
        self,
        region: str | None = None,
        *,
        seq_region_name: str | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> list[F]:
        """Return features overlapping the given genomic region.

        Args:
            region (str | None): The region (e.g., "chr1:1000-2000") to get features from.
            seq_region_name (str | None): The name of the region to get features from.
            start (int | None): The start position of the region to get features from.
            end (int | None): The end position of the region to get features from.

        Returns:
            list[F]: A list of features overlapping the specified region.
        """

        raise NotImplementedError

    @abstractmethod
    def iter_features_by_region(
        self,
        region: str | None = None,
        *,
        seq_region_name: str | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> Iterator[F]:
        """Lazily yield features overlapping the given genomic region.

        Args:
            region (str | None): The region (e.g., "chr1:1000-2000") to get features from.
            seq_region_name (str | None): The name of the region to get features from.
            start (int | None): The start position of the region to get features from.
            end (int | None): The end position of the region to get features from.

        Returns:
            Iterator[F]: An iterator over features overlapping the specified region.
        """

        raise NotImplementedError

    @abstractmethod
    def close(self) -> None:
        """Release resources held by this source (e.g., close file handles or DB connections)."""

        raise NotImplementedError

    def __enter__(self) -> FeatureSource[F]:
        return self

    def __exit__(self, *args: object) -> None:
        self.close()
