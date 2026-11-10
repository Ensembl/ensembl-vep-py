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
"""The file-backed feature source contract."""

from __future__ import annotations

import logging
from collections.abc import Iterator
from typing import Generic, TypeVar
from attrs import define, field

from ensembl.vep.common import SourceType
from ensembl.vep.common.methods import format_region_string
from ensembl.vep.models import Feature
from ensembl.vep.parsers import Parser
from ensembl.vep.sources.feature_source import FeatureSource

P = TypeVar("P", bound=Parser)
R = TypeVar("R")
F = TypeVar("F", bound=Feature)
logger = logging.getLogger(__name__)


@define(kw_only=True, auto_attribs=True)
class FileFeatureSource(FeatureSource[R, F], Generic[P, R, F]):
    """Feature source backed by a parser-readable file.

    Attributes:
        file_path (str): Path to the feature file.
        adaptor (P): Parser for the file's record type.
        source_type (SourceType): File source marker.
    """

    file_path: str
    adaptor: P
    source_type: SourceType = field(default=SourceType.FILE, init=False)

    def get_all_records(self) -> list[R]:
        """Return every unconverted record in the source file.

        Returns:
            list[R]: Every record in the file.

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
            return self.adaptor.parse_file(self.file_path)
        except NotImplementedError as unsupported:
            raise NotImplementedError(
                f"{type(self.adaptor).__name__} does not support parse_file"
            ) from unsupported

    def iter_all_records(self) -> Iterator[R]:
        """Yield every unconverted record in the source file."""

        with self.adaptor.open(self.file_path, mode="sequential") as adaptor:
            while (record := adaptor.next()) is not None:
                yield record

    def get_records_by_region(
        self,
        region: str | None = None,
        *,
        seq_region_name: str | None = None,
        start: int | None = None,
        end: int | None = None,
        **kwargs: object,
    ) -> list[R]:
        """Return unconverted records from an indexed region query.

        Args:
            region (str | None): The region (e.g., "chr1:1000-2000") to read.
            seq_region_name (str | None): The name of the sequence region (e.g., "chr1").
            start (int | None): The start position of the region.
            end (int | None): The end position of the region.
            **kwargs (object): Unused here; a format whose query takes more accepts it.

        Returns:
            list[R]: The records the query reads, empty when the region is unreadable.
        """

        del kwargs
        try:
            query_region = format_region_string(region, seq_region_name=seq_region_name, start=start, end=end)
        except ValueError:
            logger.debug("Cannot parse region: %s %s:%s-%s", region, seq_region_name, start, end)
            return []

        with self.adaptor.open(self.file_path, mode="query") as adaptor:
            if not adaptor.seek(query_region):
                return []
            return adaptor.records(limit=None)

    def iter_records_by_region(
        self,
        region: str | None = None,
        *,
        seq_region_name: str | None = None,
        start: int | None = None,
        end: int | None = None,
        **kwargs: object,
    ) -> Iterator[R]:
        """Yield unconverted records from an indexed region query.

        Args:
            region (str | None): The region (e.g., "chr1:1000-2000") to read.
            seq_region_name (str | None): The name of the sequence region (e.g., "chr1").
            start (int | None): The start position of the region.
            end (int | None): The end position of the region.
            **kwargs (object): Unused here; a format whose query takes more accepts it.

        Returns:
            Iterator[R]: The records the query reads, in file order.
        """

        del kwargs
        try:
            query_region = format_region_string(region, seq_region_name=seq_region_name, start=start, end=end)
        except ValueError:
            logger.debug("Cannot parse region: %s %s:%s-%s", region, seq_region_name, start, end)
            return

        with self.adaptor.open(self.file_path, mode="query") as adaptor:
            if not adaptor.seek(query_region):
                return
            while (record := adaptor.next()) is not None:
                yield record

    def close(self) -> None:
        """Release the attached sequence source, if any.

        ``FeatureSource`` owns the ``sequence_source`` field and ``__exit__`` calls ``close()``, but
        nothing in the chain ever released it -- so a source opened as a context manager leaked its
        sequence handle. Closing it here is what the field's ownership already implied.
        """

        if self.sequence_source is not None:
            self.sequence_source.close()
