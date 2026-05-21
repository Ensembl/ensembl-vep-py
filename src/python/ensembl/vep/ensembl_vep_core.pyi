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

"""Type stubs for the PyO3 extension module."""

from typing import Any


def sum_as_string(a: str, b: str) -> str: ...


class GFFRecord:
    seqid: str
    source: str
    feature_type: str
    start: int
    end: int
    score: float | None
    strand: str
    phase: int | None
    attributes: dict[str, tuple[str, ...]]
    raw_attributes: str


class GFFParser:
    format: str
    metadata: dict[str, Any]
    reader: str

    def __init__(
        self,
        file_path: str | None = None,
        tabix_path: str | None = None,
        reader: str | None = None,
    ) -> None: ...
    def open(self, file_path: str) -> None: ...
    def close(self) -> None: ...
    def set_reader(self, reader: str) -> None: ...
    def next(self) -> GFFRecord | None: ...
    def records(self, limit: int | None = None) -> list[GFFRecord]: ...
    def seek(self, region: str) -> bool: ...
    @staticmethod
    def parse_text(data: str, limit: int | None = None) -> list[GFFRecord]: ...
    @staticmethod
    def parse_file(file_path: str, limit: int | None = None) -> list[GFFRecord]: ...


class TabixParser:
    format: str
    metadata: dict[str, Any]

    def __init__(self, file_path: str | None = None, tabix_path: str | None = None) -> None: ...
    def open(self, file_path: str) -> None: ...
    def close(self) -> None: ...
    def seek(self, region: str) -> bool: ...
    def next(self) -> str | None: ...
    def records(self, limit: int | None = None) -> list[str]: ...
