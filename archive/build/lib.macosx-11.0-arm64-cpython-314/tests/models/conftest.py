# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Shared fixtures for model unit tests."""

from attrs import define, field
import pytest

from ensembl.vep.sources.sequence_source import SequenceSource


@define(kw_only=True, auto_attribs=True)
class DummySequenceSource(SequenceSource):
    """Minimal in-memory sequence source used by model tests."""

    data: dict[str, str] = field(factory=dict)
    calls: list[tuple[str, int, int]] = field(factory=list)
    closed: bool = field(default=False)

    def close(self) -> None:
        self.closed = True

    def _fetch_sequence_chunk(self, seq_region_name: str, start: int, end: int) -> str:
        self.calls.append((seq_region_name, start, end))
        sequence = self.data.get(seq_region_name)
        if sequence is None:
            return "N" * (end - start + 1)

        return sequence[start - 1 : end]

    def _get_seq_region_length(self, seq_region_name: str) -> int | None:
        sequence = self.data.get(seq_region_name)
        return len(sequence) if sequence is not None else None


@pytest.fixture
def dummy_sequence_source() -> DummySequenceSource:
    """Create an empty in-memory sequence source."""
    return DummySequenceSource()
