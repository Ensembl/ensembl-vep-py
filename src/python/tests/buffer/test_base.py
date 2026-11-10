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
"""Tests for the buffer contract."""

from __future__ import annotations

from collections.abc import Iterator

from attrs import define

from ensembl.vep.buffer import Buffer, FixedSizeInputBuffer, InputBuffer
from ensembl.vep.models.variation import Variation


@define(kw_only=True, auto_attribs=True)
class _EagerBuffer(InputBuffer[Variation]):
    """Buffer that takes everything left in one refill."""

    def next(self, items: Iterator[Variation]) -> None:
        """Take every remaining object."""

        self._fill(list(items), exhausted=True)


class TestBufferProtocol:
    """What a buffer has to provide."""

    def test_fixed_size_input_buffer_is_a_buffer(self) -> None:
        """Recognise the built-in buffer as a buffer."""
        assert isinstance(FixedSizeInputBuffer(), Buffer)

    def test_a_buffer_without_a_refill_policy_is_still_a_buffer(self) -> None:
        """Recognise the unspecialised input buffer too; its ``next`` just refuses."""
        assert isinstance(InputBuffer(), Buffer)

    def test_a_custom_policy_fills_the_batch(self, variation_list: list[Variation]) -> None:
        """Let a refill policy decide what one batch holds."""
        buffer: _EagerBuffer = _EagerBuffer()

        buffer.next(iter(variation_list))

        assert buffer.objects == tuple(variation_list)
        assert buffer.is_empty
