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

from attrs import define, field

# from ensembl.vep.config.base import Config
from ensembl.vep.buffers.buffer import InputBuffer
from ensembl.vep.buffers.registry import register_buffer
from ensembl.vep.caches.interval_tree import IntervalTree
from ensembl.vep.models.variation.variation import Variation


@register_buffer("variation.source.fixed")
@define(kw_only=True, auto_attribs=True)
class FixedSizeInputBuffer(InputBuffer[Variation]):
    """Input buffer with a fixed size."""

    _interval_tree: IntervalTree | None = field(default=None, init=False)

    def next(self) -> None:
        """Reset and re-fill the buffer from the source.

        Returns:
            None
        """

        self._buffer = []
        self._interval_tree = None
        seq_region_name = None
        prev_start = None

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
            if prev_start is not None:
                self._check_for_unsorted_variants(variant, prev_start)
            prev_start = variant.start

        self._is_empty = True

    def get_overlapping_variation(self, start: int, end: int) -> list[Variation] | None:
        """Return the variants in the buffer that overlap the supplied interval.

        Args:
            start: start coordinate of the interval.
            end: end coordinate of the interval.

        Returns:
            A list of Variation objects from the current buffer that overlap the interval
            or None if the buffer is empty.

        Raises:
        """

        if not self._buffer:
            return None

        # get/create interval tree
        if not self._interval_tree:
            self._interval_tree = IntervalTree(features=self._buffer)

        # get the overlapping variation objects
        if start > end:
            start, end = (end, start)
        variations = self._interval_tree.get_overlapping_features(start, end)

        return variations

    def _check_for_unsorted_variants(self, var: Variation, prev_start: int) -> None:
        """Function to check if this variant appears before the previous one.

        There is some leeway to account for deletions or variants that have been minimised.

        Args:
            var: The current variant to test
            prev_start: The position of the previous variant

        Returns:
            None

        Raises:
            RuntimeError: If the threshold for an unsorted file is reached.
        """

        if var.start - prev_start < 0:
            if abs(var.start - prev_start) > self.config.max_not_ordered_variants_distance:
                self._count_not_ordered_variants += 1
                if self._count_not_ordered_variants > self.config.max_not_ordered_variants:
                    raise RuntimeError(
                        "Exiting the program. The input file appears to be unsorted. "
                        "Please sort by chromosome and by location and re-submit."
                    )
