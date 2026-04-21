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
"""Unit tests for Exon class."""

import pytest

from ensembl.vep.common import ExonPhase, Strand
from ensembl.vep.models.genic.exon import Exon


@pytest.fixture
def exon() -> Exon:
    """Exon fixture for testing."""
    return Exon(start=10, end=50, strand=Strand.FORWARD, stable_id="ENSE00000000001")


class TestExonInitialization:
    """Test Exon initialization."""

    def test_initialization_without_keywords(self) -> None:
        """Test that initializing without keywords raises error."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = Exon(start=10, end=50, strand=Strand.FORWARD, stable_id="ENSE00000000001")


class TestExonProperties:
    """Test Exon properties and cached properties."""

    def test_core_properties(self, exon: Exon) -> None:
        """Test that Exon core properties are set correctly."""
        assert exon.start == 10
        assert exon.end == 50
        assert exon.strand == Strand.FORWARD
        assert exon.stable_id == "ENSE00000000001"

    def test_externally_settable_properties(self, exon: Exon) -> None:
        """Test that externally settable properties can be set and accessed."""
        assert exon.phase is None
        exon.phase = ExonPhase.ZERO
        assert exon.phase == ExonPhase.ZERO

        assert exon.end_phase is None
        exon.end_phase = ExonPhase.ONE
        assert exon.end_phase == ExonPhase.ONE

    def test_frame(self, exon: Exon) -> None:
        """Test that frame is calculated correctly based on phase and end_phase."""
        assert exon.frame is None

        exon.phase = ExonPhase.ZERO
        assert exon.frame == 1

        exon.phase = ExonPhase.ONE
        assert exon.frame == 0

        exon.phase = ExonPhase.TWO
        assert exon.frame == 2
