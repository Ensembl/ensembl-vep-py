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
"""Unit tests for UTR class."""

import pytest

from ensembl.vep.common import Strand
from ensembl.vep.models.genic.utr import UTR
from ensembl.vep.models.genic.transcript import Transcript


@pytest.fixture
def transcript() -> Transcript:
    """Transcript fixture for testing."""
    return Transcript(start=1, end=100, strand=Strand.FORWARD, stable_id="ENST00000000001")


@pytest.fixture
def utr() -> UTR:
    """UTR fixture for testing."""
    return UTR(start=10, end=20, strand=Strand.FORWARD)


class TestUTRInitialization:
    """Test UTR initialization."""

    def test_initialization_without_keywords(self) -> None:
        """Test that initializing without keywords raises error."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = UTR(start=10, end=20, strand=Strand.FORWARD)


class TestUTRProperties:
    """Test UTR properties."""

    def test_core_properties(self, utr: UTR) -> None:
        """Test that UTR core properties are set correctly."""
        assert utr.start == 10
        assert utr.end == 20
        assert utr.strand == 1

    def test_externally_settable_properties(self, utr: UTR, transcript: Transcript) -> None:
        """Test that externally settable properties can be set and accessed."""
        assert utr.transcript is None
        utr.transcript = transcript
        assert utr.transcript == transcript
