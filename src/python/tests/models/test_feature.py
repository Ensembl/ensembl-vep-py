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
"""Unit tests for Feature class."""

import pytest

from ensembl.vep.common import Strand
from ensembl.vep.models import Feature
from ensembl.vep.models.sequence import Slice


class TestFeatureInitialization:
    """Test Feature initialization."""

    def test_initialization_without_keywords(self) -> None:
        """Test Feature initialization without keywords raises error."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = Feature(start=10, end=20, strand=Strand.FORWARD, seq=None)


class TestFeatureProperties:
    """Test property accessors."""

    @pytest.fixture
    def feature(self) -> Feature:
        """Create a Feature instance for testing."""
        return Feature(start=10, end=20, strand=Strand.FORWARD, slice=None)

    @pytest.fixture
    def slice(self) -> Slice:
        """Create a Slice instance for testing."""
        return Slice(
            coord_system="chromosome",
            seq_region_name="1",
            seq_region_length=1000,
            start=100,
            end=200,
            strand=Strand.FORWARD,
        )

    def test_core_properties(self, feature: Feature) -> None:
        """Test core properties are set correctly."""
        assert feature.start == 10
        assert feature.end == 20
        assert feature.strand == Strand.FORWARD
        assert feature.slice is None
        assert feature.seq is None

    def test_length(self, feature: Feature) -> None:
        """Test length setter calculates from start/end when None."""
        assert feature.length == 11  # 20 - 10 + 1

    def test_length_immutable(self, feature: Feature) -> None:
        """Test length cannot be directly set."""
        with pytest.raises(AttributeError, match="property 'length' of 'Feature' object has no setter"):
            feature.length = 100  # type: ignore[misc]

        with pytest.raises(
            AttributeError, match="'Feature' object has no attribute '_TestFeatureProperties__length'"
        ):
            feature.__length = 100

    def test_length_update(self, feature: Feature) -> None:
        """Test length updates when start or end changes."""
        feature.start = 15
        assert feature.length == 6  # 20 - 15 + 1

        feature.end = 25
        assert feature.length == 11  # 25 - 15 + 1

        feature.start = 10
        feature.end = 20
        assert feature.length == 11  # back to original

    def test_seq_region_start_no_slice(self, feature: Feature) -> None:
        """Test seq_region_start returns None when no slice attached."""
        feature.slice = None
        assert feature.seq_region_start is None

    def test_seq_region_start_with_slice(self, feature: Feature, slice: Slice) -> None:
        """Test seq_region_start property with slice."""
        feature.slice = slice
        assert feature.seq_region_start == 109  # slice.start + feature.start - 1, as positions are 1-based and inclusive

    def test_seq_region_start_immutable(self, feature: Feature) -> None:
        """Test seq_region_start cannot be directly set."""
        with pytest.raises(
            AttributeError, match="property 'seq_region_start' of 'Feature' object has no setter"
        ):
            feature.seq_region_start = None  # type: ignore[misc]

        with pytest.raises(
            AttributeError,
            match="'Feature' object has no attribute '_TestFeatureProperties__seq_region_start'",
        ):
            feature.__seq_region_start = None

    def test_seq_region_start_update(self, feature: Feature, slice: Slice) -> None:
        """Test seq_region_start updates when start changes."""
        feature.slice = slice
        assert feature.seq_region_start == 109  # slice.start + feature.start - 1, as positions are 1-based and inclusive

        feature.start = 15
        assert feature.seq_region_start == 114  # slice.start + feature.start - 1, as positions are 1-based and inclusive

        feature.start = 10
        assert feature.seq_region_start == 109  # back to original

    def test_seq_region_end_no_slice(self, feature: Feature) -> None:
        """Test seq_region_end returns None when no slice attached."""
        assert feature.seq_region_end is None

    def test_seq_region_end_with_slice(self, feature: Feature, slice: Slice) -> None:
        """Test seq_region_end property with slice."""
        feature.slice = slice

        # slice.start + feature.end - 1, as positions are 1-based and inclusive
        assert feature.seq_region_end == 119

    def test_seq_region_end_immutable(self, feature: Feature) -> None:
        """Test seq_region_end cannot be directly set."""
        with pytest.raises(
            AttributeError, match="property 'seq_region_end' of 'Feature' object has no setter"
        ):
            feature.seq_region_end = None  #type: ignore[misc]

        with pytest.raises(
            AttributeError, match="'Feature' object has no attribute '_TestFeatureProperties__seq_region_end'"
        ):
            feature.__seq_region_end = None

    def test_seq_region_end_update(self, feature: Feature, slice: Slice) -> None:
        """Test seq_region_end updates when end changes."""
        feature.slice = slice
        assert feature.seq_region_end == 119  # slice.start + feature.end - 1, as positions are 1-based and inclusive

        feature.end = 25
        assert feature.seq_region_end == 124  # slice.start + feature.end - 1, as positions are 1-based and inclusive

        feature.end = 20
        assert feature.seq_region_end == 119  # back to original

    def test_seq_region_strand_no_slice(self, feature: Feature) -> None:
        """Test seq_region_strand returns None when no slice attached."""
        assert feature.seq_region_strand is None

    def test_seq_region_strand_positive(self, feature: Feature, slice: Slice) -> None:
        """Test seq_region_strand property with slice."""
        feature.slice = slice
        assert feature.seq_region_strand == Strand.FORWARD  # slice.strand * feature.strand

    def test_seq_region_strand_negative(self, feature: Feature, slice: Slice) -> None:
        """Test seq_region_strand with negative strand feature."""
        feature.slice = slice
        feature.strand = Strand.REVERSE
        assert feature.seq_region_strand == Strand.REVERSE  # slice.strand * feature.strand

    def test_seq_region_strand_immutable(self, feature: Feature) -> None:
        """Test seq_region_strand cannot be directly set."""
        with pytest.raises(
            AttributeError, match="property 'seq_region_strand' of 'Feature' object has no setter"
        ):
            feature.seq_region_strand = None  #type: ignore[misc]

        with pytest.raises(
            AttributeError,
            match="'Feature' object has no attribute '_TestFeatureProperties__seq_region_strand'",
        ):
            feature.__seq_region_strand = None

    def test_seq_region_strand_update(self, feature: Feature, slice: Slice) -> None:
        """Test seq_region_strand updates when strand changes."""
        feature.slice = slice
        assert feature.seq_region_strand == Strand.FORWARD  # slice.strand * feature.strand

        feature.strand = Strand.REVERSE
        assert feature.seq_region_strand == Strand.REVERSE  # slice.strand * feature.strand

        feature.strand = Strand.FORWARD
        assert feature.seq_region_strand == Strand.FORWARD  # back to original


class TestFeatureEquality:
    """Test equality comparison."""

    @pytest.fixture
    def slice(self) -> Slice:
        """Create a Slice instance for testing."""
        return Slice(
            coord_system="chromosome",
            seq_region_name="1",
            seq_region_length=1000,
            start=100,
            end=200,
            strand=Strand.FORWARD,
        )

    def test_equal_features(self, slice: Slice) -> None:
        """Test two equal features are equal."""
        f1 = Feature(start=10, end=20, strand=Strand.FORWARD, slice=slice)
        f2 = Feature(start=10, end=20, strand=Strand.FORWARD, slice=slice)
        assert f1 == f2

    def test_different_start(self, slice: Slice) -> None:
        """Test features with different start are not equal."""
        f1 = Feature(start=10, end=20, strand=Strand.FORWARD, slice=slice)
        f2 = Feature(start=15, end=20, strand=Strand.FORWARD, slice=slice)
        assert f1 != f2

    def test_different_end(self, slice: Slice) -> None:
        """Test features with different end are not equal."""
        f1 = Feature(start=10, end=20, strand=Strand.FORWARD, slice=slice)
        f2 = Feature(start=10, end=25, strand=Strand.FORWARD, slice=slice)
        assert f1 != f2

    def test_different_strand(self, slice: Slice) -> None:
        """Test features with different strand are not equal."""
        f1 = Feature(start=10, end=20, strand=Strand.FORWARD, slice=slice)
        f2 = Feature(start=10, end=20, strand=Strand.REVERSE, slice=slice)
        assert f1 != f2

    def test_different_slice(self, slice: Slice) -> None:
        """Test features with different slices are not equal."""
        slice2 = Slice(
            coord_system="chromosome",
            seq_region_name="1",
            seq_region_length=1000,
            start=200,
            end=300,
            strand=Strand.FORWARD,
        )
        f1 = Feature(start=10, end=20, strand=Strand.FORWARD, slice=slice)
        f2 = Feature(start=10, end=20, strand=Strand.FORWARD, slice=slice2)
        assert f1.slice != f2.slice
        assert f1 != f2

    def test_different_type(self, slice: Slice) -> None:
        """Test features of different types are not equal."""
        f1 = Feature(start=10, end=20, strand=Strand.FORWARD, slice=slice)

        class AnotherFeature(Feature):
            pass

        f2 = AnotherFeature(start=10, end=20, strand=Strand.FORWARD, slice=slice)
        assert f1 != f2


class TestFeatureOverlaps:
    """Test overlaps method."""

    @pytest.fixture
    def slice(self) -> Slice:
        """Create a Slice instance for testing."""
        return Slice(
            coord_system="chromosome",
            seq_region_name="1",
            seq_region_length=1000,
            start=1,
            end=1000,
            strand=Strand.FORWARD,
        )

    def test_overlaps_same_slice(self, slice: Slice) -> None:
        """Test overlapping features on same slice."""
        f1 = Feature(start=10, end=20, strand=Strand.FORWARD, slice=slice)
        f2 = Feature(start=15, end=25, strand=Strand.FORWARD, slice=slice)
        assert f1.overlaps(f2) is True

    def test_overlaps_no_overlap(self, slice: Slice) -> None:
        """Test non-overlapping features."""
        f1 = Feature(start=10, end=20, strand=Strand.FORWARD, slice=slice)
        f2 = Feature(start=30, end=40, strand=Strand.FORWARD, slice=slice)
        assert f1.overlaps(f2) is False

    def test_overlaps_different_seq_regions(self, slice: Slice) -> None:
        """Test features on different sequence regions don't overlap."""
        slice1 = slice
        slice2 = Slice(
            coord_system="chromosome",
            seq_region_name="2",
            seq_region_length=1000,
            start=1,
            end=1000,
            strand=Strand.FORWARD,
        )

        f1 = Feature(start=10, end=20, strand=Strand.FORWARD, slice=slice1)
        f2 = Feature(start=10, end=20, strand=Strand.FORWARD, slice=slice2)
        assert f1.overlaps(f2) is False

    def test_overlaps_no_slice(self) -> None:
        """Test features without slices return None."""
        f1 = Feature(start=10, end=20, strand=Strand.FORWARD, slice=None)
        f2 = Feature(start=15, end=25, strand=Strand.FORWARD, slice=None)
        assert f1.overlaps(f2) is None


class TestFeatureFeatureSlice:
    """Test feature_slice method."""

    @pytest.fixture
    def slice(self) -> Slice:
        """Create a Slice instance for testing."""
        return Slice(
            coord_system="chromosome",
            seq_region_name="1",
            seq_region_length=1000,
            start=100,
            end=200,
            strand=Strand.FORWARD,
        )

    def test_feature_slice(self, slice: Slice) -> None:
        """Test feature_slice creates correct slice."""
        f = Feature(start=10, end=20, strand=Strand.FORWARD, slice=slice)
        feature_slice = f.feature_slice()

        assert feature_slice is not None
        assert feature_slice.seq_region_name == "1"
        assert feature_slice.start == 109  # seq_region_start
        assert feature_slice.end == 119  # seq_region_end
        assert feature_slice.strand == Strand.FORWARD  # seq_region_strand

    def test_feature_slice_no_slice(self) -> None:
        """Test feature_slice returns None when no slice attached."""
        f = Feature(start=10, end=20, strand=Strand.FORWARD, slice=None)
        assert f.feature_slice() is None


class TestFeatureTransform:
    """Test transform method."""

    ...


class TestFeatureTransfer:
    """Test transfer method."""

    ...
