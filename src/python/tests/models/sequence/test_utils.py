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

# pylint: disable=redefined-outer-name
"""Unit tests for sequence utils functions."""

from ensembl.vep.models.sequence.slice import Slice
from ensembl.vep.models.sequence.utils import get_or_build_slice
from tests.models.conftest import DummySequenceSource


class TestGetOrBuildSlice:
    """Tests for the shared, cached Slice builder."""

    def test_get_or_build_slice_uses_sequence_source_length(self) -> None:
        """Test the built slice uses the sequence source's region length."""

        source = DummySequenceSource()
        source.data = {"chr1": "ACGTACGTACGTACGTACGTACGTACGTACGT"}

        slice_cache: dict[str, Slice] = {}
        seq_slice = get_or_build_slice("chr1", slice_cache, source)  # pylint: disable=protected-access

        assert seq_slice.end == 32
        assert seq_slice.seq_region_length == 32
        assert seq_slice.seq == "ACGTACGTACGTACGTACGTACGTACGTACGT"

    def test_get_or_build_slice_without_sequence_source_uses_placeholder_end(self) -> None:
        """Test the built slice falls back to a placeholder end without a sequence source."""

        source = DummySequenceSource()

        seq_slice = get_or_build_slice("chr1", {}, source)  # pylint: disable=protected-access

        assert seq_slice.start == 1
        assert seq_slice.seq_region_length is None
        assert seq_slice.end == 2**31 - 1

    def test_get_or_build_slice_reuses_cached_instance(self) -> None:
        """Test the same Slice instance is returned for the same region."""

        source = DummySequenceSource()
        source.data = {"chr1": "ACGTACGTACGTACGTACGTACGTACGTACGT"}
        slice_cache: dict[str, Slice] = {}

        first = get_or_build_slice("chr1", slice_cache)  # pylint: disable=protected-access
        second = get_or_build_slice("chr1", slice_cache)  # pylint: disable=protected-access

        assert first is second

    def test_get_or_build_slice_different_regions_get_different_instances(self) -> None:
        """Test different regions produce distinct Slice instances."""

        source = DummySequenceSource()
        source.data = {
            "chr1": "ACGTACGTACGTACGTACGTACGTACGTACGT",
            "chr2": "GCTGGCTGGCTGGCTGGCTGGCTGGCTGGCTGGCTG",
        }
        slice_cache: dict[str, Slice] = {}

        chr1_slice = get_or_build_slice("chr1", slice_cache)  # pylint: disable=protected-access
        chr2_slice = get_or_build_slice("chr2", slice_cache)  # pylint: disable=protected-access

        assert chr1_slice is not chr2_slice
        assert chr1_slice.seq_region_name == "chr1"
        assert chr2_slice.seq_region_name == "chr2"
