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
    """Parser interface.
    
    Attributes:
        format (str): The format of the file being parsed.
        metadata (dict): A dictionary containing metadata about the file being parsed.
    """

    format: str

    @property
    def metadata(self) -> dict:
        """Public accessor for the _metadata.
        
        Returns:
            dict: The metadata dictionary.
        """

        return self._metadata

    @property
    def _metadata(self) -> dict:
        """Metadata dictionary.
        
        Returns:
            dict: The metadata dictionary.
        """

        return {}

    @abstractmethod
    def open(self, file_path: str) -> "Parser":
        """Open a file for parsing.
        
        Args:
            file_path (str): The path to the file to be parsed.
            
        Returns:
            Parser: The parser instance.
        """

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
        """Seek to a region and make next() stream matching raw rows.
        
        Args:
            region (str | None): The region to seek to, in the format "chr:start-end".
            seqid (str | None): The sequence ID of the region to seek to.
            start (int | None): The start position of the region to seek to.
            end (int | None): The end position of the region to seek to.
            
        Returns:
            bool: True if the seek was successful, False otherwise.
        """

        raise NotImplementedError

    @abstractmethod
    def next(self) -> T | None:
        """Return the next raw row from the active region query.
        
        Returns:
            T | None: The next raw row, or None if there are no more rows.
        """

        raise NotImplementedError

    def records(self, limit: int | None = None) -> list[T]:
        """Return all records from the active region query.

        Args:
            limit (int | None): The maximum number of records to return. If None, return all records.

        Returns:
            list[T]: A list of records from the active region query.
        """

        records: list[T] = []
        while limit is None or len(records) < limit:
            record = self.next()
            if record is None:
                break
            records.append(record)
        return records
