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
"""Fixed Sized Buffer."""

from __future__ import annotations

from attrs import define

# from ensembl.vep.config.base import Config
from ensembl.vep.buffers.buffer import InputBuffer
from ensembl.vep.buffers.registry import register_buffer
from ensembl.vep.models.variation.variation import Variation


@register_buffer("variation.source.fixed")
@define(kw_only=True, auto_attribs=True)
class FixedSizeInputBuffer(InputBuffer[Variation]):
    """Input buffer with a fixed size."""

    def next(self) -> None:
        """Reset and re-fill the buffer from the source.

        Returns:
            None
        """

        self._buffer = []
        seq_region_name = None

        while self._pre_buffer and len(self._buffer) < self.config.buffer_size:
            variant = self._pre_buffer.popleft()
            assert variant.slice is not None
            variant_seq_region_name = variant.slice.seq_region_name

            if seq_region_name is None:
                seq_region_name = variant_seq_region_name
            elif variant_seq_region_name != seq_region_name:
                self._pre_buffer.append(variant)
                return

            self._buffer.append(variant)

        if len(self._buffer) == self.config.buffer_size:
            return

        for variant in self.variation_iterator:
            assert variant.slice is not None
            variant_seq_region_name = variant.slice.seq_region_name

            if seq_region_name is None:
                seq_region_name = variant_seq_region_name
            elif variant_seq_region_name != seq_region_name:
                self._pre_buffer.append(variant)
                return

            if len(self._buffer) == self.config.buffer_size:
                self._pre_buffer.append(variant)
                return

            self._buffer.append(variant)

        self._is_empty = True

    def get_overlapping_variation(self, start: int, end: int) -> list[Variation]:
        """Return the variants in the buffer that overlap the supplied interval."""
        raise NotImplementedError
