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
"""FASTA parser interface backed by the Rust parser."""

from __future__ import annotations

from pathlib import Path
from attr import define, field

from ensembl.vep.common import format_region_string
from ensembl.vep.parsers.parser import Parser
from ensembl.vep.ensembl_vep_core import (  # noqa: F401,F403 pylint: disable=E0401,E0611
    FastaParser as _RustFastaParser,
    FastaRecord,
)

__all__ = ["FastaParser", "FastaRecord"]


@define(kw_only=True, auto_attribs=True)
class FastaParser(Parser[FastaRecord]):
    """Parser facade for FASTA files.

    Attributes:
        format (str): The format of the file being parsed.
        file_path (str | Path | None): The path to the file being parsed.
        mode (str): The reader mode, either "auto", "query" or "sequential".
                Default is "auto", which will attempt to use "query" mode
                and fall back to "sequential" if the file is not indexed.
    """

    format: str = field(init=False, default="fasta")
    file_path: str | Path | None = field(default=None)
    mode: str = field(default="auto")
    _parser: _RustFastaParser = field(init=False, repr=False)

    def __attrs_post_init__(self) -> None:
        self._parser = _RustFastaParser(
            str(self.file_path) if self.file_path is not None else None,
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
        """Metadata dictionary containing FASTA sequence lengths.

        Returns:
            dict: A dictionary with sequence length metadata.
        """

        return dict(self._parser.metadata)

    @property
    def reader_mode(self) -> str:
        """Return the configured reader mode."""

        return self._parser.reader_mode

    def open(self, file_path: str | Path, mode: str = "auto") -> FastaParser:
        """Open a FASTA file.

        Args:
            file_path (str | Path): The path to the FASTA file to be parsed.
            mode (str): The mode in which to open the file.

        Returns:
            FastaParser: The parser instance.
        """

        self.file_path = file_path
        self._parser.open(str(file_path), mode)
        return self

    def close(self) -> None:
        """Close the parser."""

        self._parser.close()

    def seek(
        self,
        region: str | None = None,
        *,
        seq_region_name: str | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> bool:
        """Seek to a FASTA record or subsequence.

        Args:
            region (str | None): The region to seek to, in the format "chr:start-end".
            seq_region_name (str | None): The name of the region to seek to.
            start (int | None): The start position of the region to seek to.
            end (int | None): The end position of the region to seek to.

        Returns:
            bool: True if the seek was successful, False otherwise.
        """

        return self._parser.seek(
            format_region_string(region, seq_region_name=seq_region_name, start=start, end=end)
        )

    def next(self) -> FastaRecord | None:
        """Return the next FASTA record, or None at EOF.

        Returns:
            FastaRecord | None: The next FASTA record, or None if there are no more records.
        """

        return self._parser.next()

    def records(self, limit: int | None = 10000) -> list[FastaRecord]:
        """Return records from the current stream position.

        Args:
            limit (int | None): The maximum number of records to return. Default is 10000.

        Returns:
            list[FastaRecord]: A list of FASTA records.
        """

        return self._parser.records(limit)

    @staticmethod
    def parse_text(data: str, limit: int | None = None) -> list[FastaRecord]:
        """Parse FASTA records from text."""

        return _RustFastaParser.parse_text(data, limit)

    @staticmethod
    def parse_file(file_path: str | Path, limit: int | None = None) -> list[FastaRecord]:
        """Parse FASTA records from a file."""

        return _RustFastaParser.parse_file(str(file_path), limit)
