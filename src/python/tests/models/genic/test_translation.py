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
"""Unit tests for Translation class."""

import pytest

from ensembl.vep.common import Strand
from ensembl.vep.models.genic.translation import Translation
from ensembl.vep.models.genic.transcript import Transcript


@pytest.fixture
def transcript() -> Transcript:
    """Transcript fixture for testing."""
    return Transcript(start=1, end=100, strand=Strand.FORWARD, stable_id="ENST00000000001")


@pytest.fixture
def translation() -> Translation:
    """Translation fixture for testing."""
    return Translation(start=10, end=90, stable_id="ENSP00000000001")


class TestTranslationInitialization:
    """Test Translation initialization."""

    def test_initialization_without_keywords(self) -> None:
        """Test that initializing without keywords raises error."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = Translation(10, end=90, stable_id="ENSP00000000001")  # type: ignore[misc]


class TestTranslationProperties:
    """Test Translation properties and cached properties."""

    def test_core_properties(self, translation: Translation) -> None:
        """Test that Translation core properties are set correctly."""
        assert translation.start == 10
        assert translation.end == 90
        assert translation.stable_id == "ENSP00000000001"

    def test_externally_settable_properties(self, translation: Translation, transcript: Transcript) -> None:
        """Test that externally settable properties can be set and accessed."""
        assert translation.version is None
        translation.version = 1
        assert translation.version == 1

        assert translation.seq is None
        translation.seq = "M" * 81
        assert translation.seq == "M" * 81

        assert translation.transcript is None
        translation.transcript = transcript
        assert translation.transcript == transcript

    def test_length_cached_property(self, translation: Translation) -> None:
        """Test that length cached property is calculated correctly and cached."""
        # Length should be calculated as end - start + 1
        assert translation.length == 81

        # Change start and end to invalidate cache
        translation.start = 20
        translation.end = 80

        # Length should be recalculated as 80 - 20 + 1 = 61
        assert translation.length == 61
