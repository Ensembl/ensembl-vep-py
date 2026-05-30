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
"""Parser interface."""

from abc import ABC, abstractmethod
from typing import TypeVar, Generic
from attrs import define

T = TypeVar("T")


@define(kw_only=True, auto_attribs=True)
class Parser(ABC, Generic[T]):
    """Parser interface."""

    format: str

    @property
    def metadata(self) -> dict:
        """Metadata dictionary."""

        return self._metadata

    @property
    def _metadata(self) -> dict:
        """Metadata dictionary."""

        return {}

    @abstractmethod
    def open(self, file_path: str) -> "Parser":
        """Open a file for parsing."""

        raise NotImplementedError

    @abstractmethod
    def close(self) -> None:
        """Close the parser."""

        raise NotImplementedError

    @abstractmethod
    def seek(
        self,
        region: str | None = None,
        *,
        seqid: str | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> bool:
        """Seek to a region and make next() stream matching raw rows."""

        raise NotImplementedError

    @abstractmethod
    def next(self) -> T | None:
        """Return the next raw row from the active region query."""

        raise NotImplementedError

    def records(self, limit: int | None = None) -> list[T]:
        """Return all records from the active region query."""

        records: list[T] = []
        while limit is None or len(records) < limit:
            record = self.next()
            if record is None:
                break
            records.append(record)
        return records
