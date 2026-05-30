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
from ensembl.vep.ensembl_vep_core import TabixParser as _RustTabixParser  # noqa: F401,F403 pylint: disable=E0401,E0611


@define(kw_only=True, auto_attribs=True)
class TabixParser(Parser[str]):
    """Parser facade for tabix-indexed row streams."""

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
        """Metadata dictionary."""

        return self._metadata

    @cached_property
    def _metadata(self) -> dict:
        """Metadata dictionary."""

        return dict(self._parser.metadata)

    def open(self, file_path: str | Path) -> "TabixParser":
        """Open a tabix-indexed BGZF file."""

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
        seqid: str | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> bool:
        """Seek to a region and make next() stream matching raw rows."""

        return self._parser.seek(format_region_string(region, seqid=seqid, start=start, end=end))

    def next(self) -> str | None:
        """Return the next raw row from the active region query."""

        return self._parser.next()

    def records(self, limit: int | None = None) -> list[str]:
        """Return raw rows from the active region query."""

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


__all__ = ["TabixParser"]
