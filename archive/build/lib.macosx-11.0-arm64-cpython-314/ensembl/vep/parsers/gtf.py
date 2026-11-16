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
"""GTF parser interface backed by the Rust noodles-gtf parser."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from attr import define, field, setters

from ensembl.vep.common import format_region_string
from ensembl.vep.parsers.parser import Parser
from ensembl.vep.ensembl_vep_core import (  # noqa: F401,F403 pylint: disable=E0401,E0611
    GTFParser as _RustGTFParser,
    GTFRecord,
)

__all__ = ["GTFParser", "GTFRecord"]


@define(kw_only=True, auto_attribs=True)
class GTFParser(Parser[GTFRecord]):
    """Parser facade for GTF files.

    Attributes:
        format (str): The format of the file being parsed.
        file_path (str | Path | None): The path to the file being parsed.
        tabix_path (str): The path to the tabix executable.
        mode (str): The reader mode, either "auto", "query" or "sequential".
                Default is "auto", which will attempt to use "query" mode
                and fall back to "sequential" if the file is not indexed.
    """

    format: str = field(init=False, default="gtf")
    file_path: str | Path | None = field(default=None)
    tabix_path: str | Path | None = field(default="tabix")
    mode: str = field(default="auto", on_setattr=setters.frozen)
    _parser: _RustGTFParser = field(init=False, repr=False)

    def __attrs_post_init__(self) -> None:
        self._parser = _RustGTFParser(
            str(self.file_path) if self.file_path is not None else None,
            str(self.tabix_path) if self.tabix_path is not None else None,
            self.mode,
        )

    @property
    def metadata(self) -> dict:
        """Public accessor for the _metadata.

        Returns:
            dict: The metadata dictionary.
        """

        return self._metadata

    @property
    def _metadata(self) -> dict:
        """Metadata dictionary containing GTF header lines.

        Returns:
            dict: A dictionary with keys "header".
        """

        return dict(self._parser.metadata)

    @property
    def reader_mode(self) -> str:
        """Returns the configured reader mode."""

        return self._parser.reader_mode

    @property
    def index_sequence_names(self) -> list[str]:
        """Reference sequence names declared by the tabix index, in index order.

        Empty when the file is unopened, read sequentially, or falls back to the system
        ``tabix`` executable.

        Returns:
            list[str]: The reference sequence names the index declares.
        """

        return self._parser.index_sequence_names

    def open(self, file_path: str | Path, mode: str = "auto") -> GTFParser:
        """Open a GTF file.

        Args:
            file_path (str | Path): The path to the GTF file to be parsed.

        Returns:
            GTFParser: The parser instance.
        """

        self.file_path = file_path
        self._clear_metadata_cache()
        self._parser.open(str(file_path), mode)
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
        """Seek to a region and make next() stream matching records.

        Args:
            region (str | None): The region to seek to, in the format "chr:start-end".
            seq_region_name (str | None): The name of the region to seek to..
            start (int | None): The start position of the region to seek to.
            end (int | None): The end position of the region to seek to.

        Returns:
            bool: True if the seek was successful, False otherwise.
        """

        return self._parser.seek(
            format_region_string(region, seq_region_name=seq_region_name, start=start, end=end)
        )

    def next(self) -> GTFRecord | None:
        """Returns the next GTF record, or None at EOF.

        Returns:
            GTFRecord | None: The next GTF record,
                    or None if the end of the file reached.
        """

        return self._parser.next()

    def records(
        self, limit: int | None = 10000, feature_types: Iterable[str] | None = None
    ) -> list[GTFRecord]:
        """Returns records from the current stream position.

        Args:
            limit (int | None): The maximum number of records to return. Default is 10000.
            feature_types (Iterable[str] | None): When given, only rows whose column 3 matches one
                of these values (case-insensitively) are returned. The Rust side applies this
                before a non-matching line is parsed into a record or converted into a Python
                object.

        Returns:
            list[GTFRecord]: A list of GTF records.
        """

        return self._parser.records(limit, None if feature_types is None else list(feature_types))

    def _clear_metadata_cache(self) -> None:
        try:
            del self._metadata
        except AttributeError:
            pass

    @staticmethod
    def parse_text(data: str, limit: int | None = None) -> list[GTFRecord]:
        """Parse GTF records from text.

        Args:
            data (str): The GTF text to parse.
            limit (int | None): The maximum number of records to return. If None, return all records.

        Returns:
            list[GTFRecord]: A list of GTF records.
        """

        return _RustGTFParser.parse_text(data, limit)

    @staticmethod
    def parse_file(file_path: str | Path, limit: int | None = None) -> list[GTFRecord]:
        """Parse GTF records from a file.

        Args:
            file_path (str | Path): The path to the GTF file to parse.
            limit (int | None): The maximum number of records to return. If None, return all records.

        Returns:
            list[GTFRecord]: A list of GTF records.
        """

        return _RustGTFParser.parse_file(str(file_path), limit)
