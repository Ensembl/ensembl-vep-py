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
"""Abstract base for file-backed feature sources."""

import logging
from abc import abstractmethod
from collections.abc import Iterator
from typing import Generic, TypeVar
from attrs import define, field

from ensembl.vep.common import SourceType
from ensembl.vep.common.methods import format_region_string
from ensembl.vep.models import Feature
from ensembl.vep.parsers import Parser
from ensembl.vep.sources.feature_source import FeatureSource

P = TypeVar("P", bound=Parser)
F = TypeVar("F", bound=Feature)
logger = logging.getLogger(__name__)


@define(kw_only=True, auto_attribs=True)
class FileFeatureSource(FeatureSource, Generic[P, F]):
    """Abstract base for file-backed feature sources.

    P: Parser/adaptor type (e.g. GFF3Parser)
    F: Output Feature type (e.g. Transcript, Variation)

    Concrete subclasses implement `_records_to_features` to handle
    format-specific record parsing and model construction.
    """

    file_path: str
    adaptor: P
    source_type: SourceType = field(default=SourceType.FILE, init=False)

    def get_all_features(self) -> list[F]:
        """Returns all features from the provided file.

        Returns:
            list[F]: Every feature the file holds.

        Raises:
            NotImplementedError: If the adaptor cannot parse a whole file at once. Re-raised with
                the adaptor named, which the bare error from the parser does not say.
            FileNotFoundError: If the file does not exist. Allowed to propagate: its message
                already names the path.

        A library reports what went wrong and lets its caller decide what to do about it. These
        both used to call ``sys.exit(1)``, which ended the process from underneath any caller --
        including a test, and including a runner that had a better answer than dying.
        """

        try:
            records = self.adaptor.parse_file(self.file_path)
        except NotImplementedError as unsupported:
            raise NotImplementedError(
                f"{type(self.adaptor).__name__} does not support parse_file"
            ) from unsupported
        return self._records_to_features(records)

    def get_features_by_region(
        self,
        region: str | None = None,
        *,
        seq_region_name: str | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> list[F]:
        """Returns features overlapping the given genomic region.

        Args:
            region (str | None): The region (e.g., "chr1:1000-2000") to get features from.
            seq_region_name (str | None): The name of the sequence region (e.g., "chr1").
            start (int | None): The start position of the region.
            end (int | None): The end position of the region.

        Returns:
            list[F]: A list of features overlapping the specified region.
        """

        try:
            region = format_region_string(region, seq_region_name=seq_region_name, start=start, end=end)
        except ValueError:
            logger.debug(
                f"Cannot parse region: region={region}, "
                f"seq_region_name={seq_region_name}, start={start}, end={end}"
            )
            return []

        with self.adaptor.open(self.file_path, mode="query") as adaptor:
            if not adaptor.seek(region):
                return []
            records = adaptor.records(limit=None)

        return self._records_to_features(records)

    def iter_all_features(self) -> Iterator[F]:
        """Lazily yield all features from the provided file.

        Converts and yields each record as soon as it comes off the parser,
        rather than buffering the whole file. This default is only correct
        for formats where one record maps to one (or more) complete features
        (e.g. VCF); formats that need to assemble several records into one
        feature (e.g. GFF3) must override this.
        """

        with self.adaptor.open(self.file_path, mode="sequential") as adaptor:
            while (record := adaptor.next()) is not None:
                yield from self._records_to_features([record])

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
            seq_region_name (str | None): The name of the sequence region (e.g., "chr1").
            start (int | None): The start position of the region.
            end (int | None): The end position of the region.

        Returns:
            Iterator[F]: An iterator over features overlapping the specified region.
        """

        try:
            region = format_region_string(region, seq_region_name=seq_region_name, start=start, end=end)
        except ValueError:
            logger.debug(
                f"Cannot parse region: region={region}, "
                f"seq_region_name={seq_region_name}, start={start}, end={end}"
            )
            return

        with self.adaptor.open(self.file_path, mode="query") as adaptor:
            if not adaptor.seek(region):
                return
            while (record := adaptor.next()) is not None:
                yield from self._records_to_features([record])

    def close(self) -> None:
        """Release the attached sequence source, if any.

        ``FeatureSource`` owns the ``sequence_source`` field and ``__exit__`` calls ``close()``, but
        nothing in the chain ever released it -- so a source opened as a context manager leaked its
        sequence handle. Closing it here is what the field's ownership already implied.
        """

        if self.sequence_source is not None:
            self.sequence_source.close()

    @abstractmethod
    def _records_to_features(self, records: list) -> list[F]:
        """Convert format-specific raw records to domain Feature objects."""

        raise NotImplementedError
