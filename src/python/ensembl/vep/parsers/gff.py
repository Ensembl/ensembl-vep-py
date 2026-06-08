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
"""GFF parser interface backed by the Rust noodles-gff parser."""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path
from functools import cached_property
from attr import define, field

from ensembl.vep.common import format_region_string
from ensembl.vep.parsers import Parser
from ensembl.vep.ensembl_vep_core import GFFParser as _RustGFFParser
from ensembl.vep.ensembl_vep_core import GFFRecord

__all__ = ["GFFParser", "GFFRecord"]


@define(kw_only=True, auto_attribs=True)
class GFFParser(Parser[GFFRecord]):
    """Parser facade for GFF files.
    
    Attributes: 
        format (str): The format of the file being parsed.
        file_path (str | Path | None): The path to the file being parsed.
        tabix_path (str): The path to the tabix executable.
        mode (str): The reader mode, either "auto", "query" or "sequential". Default is "auto", which will attempt to use "query" mode and fall back to "sequential" if the file is not indexed.
    """

    format: str = field(init=False, default="gff")
    file_path: str | Path | None = field(default=None)
    tabix_path: str | Path | None = field(default="tabix")
    mode: str = field(default="auto")
    _parser: _RustGFFParser = field(init=False, repr=False)

    def __attrs_post_init__(self) -> None:
        self._parser = _RustGFFParser(
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
    
    @cached_property
    def _metadata(self) -> dict:
        """Metadata dictionary containing GFF directives and header lines.

        Returns:
            dict: A dictionary with keys "directives" and "header", where "directives" is a list of GFF directives (lines starting with ##) and "header" is a list of header lines (lines starting with # but not ##).
        """

        return self._parser.metadata()

    @property
    def reader_mode(self) -> str:
        """Return the configured reader mode."""

        return self._parser.reader_mode

    def open(self, file_path: str | Path) -> GFFParser:
        """Open a GFF file.
        
        Args:
            file_path (str | Path): The path to the GFF file to be parsed.
            
        Returns:
            GFFParser: The parser instance.
        """

        self._parser.open(str(file_path))
        return self

    def close(self) -> None:
        """Close the parser."""

        self._parser.close()

    def seek(
        self,
        region: str | None = None,
        *,
        seqid: str | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> bool:
        """Seek to a region and make next() stream matching records.
        
        Args:
            region (str | None): The region to seek to, in the format "chr:start-end".
            seqid (str | None): The sequence ID of the region to seek to..
            start (int | None): The start position of the region to seek to.
            end (int | None): The end position of the region to seek to.

        Returns:
            bool: True if the seek was successful, False otherwise.
        """

        return self._parser.seek(format_region_string(region, seqid=seqid, start=start, end=end))

    def next(self) -> GFFRecord | None:
        """Return the next GFF record, or None at EOF/##FASTA.
        
        Returns:
            GFFRecord | None: The next GFF record, or None if the end of the file or a ##FASTA directive is reached.
        """

        return self._parser.next()

    def records(self, limit: int | None = 10000) -> list[GFFRecord]:
        """Return records from the current stream position.
        
        Args:
            limit (int | None): The maximum number of records to return. Default is 10000.

        Returns:
            list[GFFRecord]: A list of GFF records.
        """

        return self._parser.records(limit)

    def __iter__(self) -> Iterator[GFFRecord]:
        return self

    def __next__(self) -> GFFRecord:
        record = self.next()

        if record is None:
            raise StopIteration

        return record

    def __enter__(self) -> GFFParser:
        return self

    def __exit__(self, *_exc_info: object) -> None:
        self.close()

    @staticmethod
    def parse_text(data: str, limit: int | None = None) -> list[GFFRecord]:
        """Parse GFF records from text."""

        return _RustGFFParser.parse_text(data, limit)

    @staticmethod
    def parse_file(file_path: str | Path, limit: int | None = None) -> list[GFFRecord]:
        """Parse GFF records from a file."""

        return _RustGFFParser.parse_file(str(file_path), limit)
