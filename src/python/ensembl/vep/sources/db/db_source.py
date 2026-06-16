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
"""Abstract base for database-backed feature sources."""

from abc import abstractmethod
from typing import Generic, Iterator, TypeVar

from attrs import define, field

from ensembl.vep.common import SourceType
from ensembl.vep.models import Feature
from ensembl.vep.models.sequence.slice import Slice
from ensembl.vep.sources.feature_source import FeatureSource

D = TypeVar("D")
F = TypeVar("F", bound=Feature)


@define(kw_only=True, auto_attribs=True)
class DBSource(FeatureSource, Generic[D, F]):
    """Abstract base for database-backed feature sources.

    D: DBAdaptor type
    F: Output Feature type (e.g. Transcript)

    Concrete subclasses implement `_rows_to_features` to handle
    schema-specific query results and model construction.
    """

    adaptor: D
    source_type: SourceType = field(default=SourceType.DB, init=False)

    def get_all_features(self) -> list[F]:
        return self._rows_to_features(self.adaptor.fetch_all())

    def iter_all_features(self) -> Iterator[F]:
        raise NotImplementedError("DB sources do not yet support streaming iteration.")

    def get_features_by_region(
        self,
        seq_region_name: str | Slice | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> list[F]:
        if isinstance(seq_region_name, Slice):
            seq_region_name = seq_region_name.seq_region_name
        return self._rows_to_features(
            self.adaptor.fetch_by_region(seq_region_name=seq_region_name, start=start, end=end)
        )

    def iter_features_by_region(
        self,
        seq_region_name: str | Slice | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> Iterator[F]:
        raise NotImplementedError("DB sources do not yet support streaming iteration.")

    def close(self) -> None:
        pass

    @abstractmethod
    def _rows_to_features(self, rows: list) -> list[F]:
        """Convert database result rows to domain Feature objects."""
        ...
