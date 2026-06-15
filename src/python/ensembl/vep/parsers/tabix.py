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
"""Tabix parser interface backed by Rust."""

from collections.abc import Iterator
from pathlib import Path
from functools import cached_property
from attrs import define, field

from ensembl.vep.common import format_region_string
from ensembl.vep.parsers.parser import Parser
from ensembl.vep.ensembl_vep_core import (  # noqa: F401,F403 pylint: disable=E0401,E0611
    TabixParser as _RustTabixParser,
)

__all__ = ["TabixParser"]


@define(kw_only=True, auto_attribs=True)
class TabixParser(Parser[str]):
    """Parser facade for tabix-indexed row streams.

    Attributes:
        format (str): The format of the file being parsed.
        file_path (str | Path | None): The path to the file to be parsed.
        tabix_path (str | Path | None): The path to the tabix executable.
        _parser (_RustTabixParser): The underlying Rust tabix parser.
    """

    format: str = field(init=False, default="tabix")
    file_path: str | Path | None = field(default=None)
    tabix_path: str | Path | None = field(default="tabix")
    _parser: _RustTabixParser = field(init=False, repr=False)

    def __attrs_post_init__(self) -> None:
        self._parser = _RustTabixParser(
            str(self.file_path) if self.file_path is not None else None,
            str(self.tabix_path) if self.tabix_path is not None else None,
        )

    @property
    def metadata(self) -> dict:
        """Public accessor for the _metadata.

        Returns:
            dict: The metadata dictionary.
        """

        return self._metadata

    @cached_property
    def _metadata(self) -> dict:
        """Metadata dictionary.

        Returns:
            dict: The metadata dictionary.
        """

        return dict(self._parser.metadata)

    def open(self, file_path: str | Path) -> "TabixParser":
        """Open a tabix-indexed BGZF file.

        Args:
            file_path (str | Path): The path to the file to be parsed.

        Returns:
            TabixParser: The parser instance.
        """

        self._parser.open(str(file_path))
        self._clear_metadata_cache()
        return self

    def close(self) -> None:
        """Close the parser."""

        self._parser.close()
        self._clear_metadata_cache()

    def seek(
        self,
        region: str | None = None,
        *,
        seq_region_name: str | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> bool:
        """Seek to a region and make next() stream matching raw rows.

        Args:
            region (str | None): The region to seek to, in the format "chr:start-end".
            seq_region_name (str | None): The sequence ID of the region to seek to.
            start (int | None): The start position of the region to seek to.
            end (int | None): The end position of the region to seek to.

        Returns:
            bool: True if the seek was successful, False otherwise.
        """

        return self._parser.seek(
            format_region_string(region, seq_region_name=seq_region_name, start=start, end=end)
        )

    def next(self) -> str | None:
        """Return the next raw row from the active region query.

        Returns:
            str | None: The next raw row from the active region query, or None if there are no more rows.
        """

        return self._parser.next()

    def records(self, limit: int | None = None) -> list[str]:
        """Return raw rows from the active region query.

        Args:
            limit (int | None): The maximum number of records to return. If None, return all records.

        Returns:
            list[str]: A list of raw rows from the active region query.
        """

        return self._parser.records(limit)

    def __iter__(self) -> Iterator[str]:
        return self

    def __next__(self) -> str:
        record = self.next()

        if record is None:
            raise StopIteration

        return record

    def __enter__(self) -> "TabixParser":
        return self

    def __exit__(self, *_exc_info: object) -> None:
        self.close()

    def _clear_metadata_cache(self) -> None:
        try:
            del self._metadata
        except AttributeError:
            pass
