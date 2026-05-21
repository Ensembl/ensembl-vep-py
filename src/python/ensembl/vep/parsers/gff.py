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

from collections.abc import Iterator
from pathlib import Path

from ensembl.vep.common import format_region_string
from ensembl.vep.ensembl_vep_core import GFFParser as _RustGFFParser
from ensembl.vep.ensembl_vep_core import GFFRecord


class GFFParser:
    """Parser facade matching the Python parser interface."""

    format = "gff"

    def __init__(
        self,
        file_path: str | Path | None = None,
        *,
        tabix_path: str = "tabix",
        reader: str = "auto",
    ) -> None:
        self.metadata: dict[str, str] = {}
        self._parser = _RustGFFParser(str(file_path) if file_path is not None else None, tabix_path, reader)

    @property
    def reader(self) -> str:
        """Return the configured reader mode."""

        return self._parser.reader

    def open(self, file_path: str | Path) -> "GFFParser":
        """Open a GFF file for streaming reads."""

        self._parser.open(str(file_path))
        return self

    def set_reader(self, reader: str) -> "GFFParser":
        """Set the reader mode."""

        self._parser.set_reader(reader)
        return self

    def close(self) -> None:
        """Close the parser."""

        self._parser.close()

    def next(self) -> GFFRecord | None:
        """Return the next GFF record, or None at EOF/##FASTA."""

        return self._parser.next()

    def records(self, limit: int | None = None) -> list[GFFRecord]:
        """Return records from the current stream position."""

        return self._parser.records(limit)

    def seek(
        self,
        region: str | None = None,
        *,
        seqid: str | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> bool:
        """Seek to a tabix region and make next() stream matching records."""

        return self._parser.seek(format_region_string(region, seqid=seqid, start=start, end=end))

    def __iter__(self) -> Iterator[GFFRecord]:
        return self

    def __next__(self) -> GFFRecord:
        record = self.next()

        if record is None:
            raise StopIteration

        return record

    def __enter__(self) -> "GFFParser":
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

__all__ = ["GFFParser", "GFFRecord"]
