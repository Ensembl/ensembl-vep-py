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
"""Unit tests for cDNA class."""

import pytest

from ensembl.vep.common import Strand
from ensembl.vep.models.genic.cdna import cDNA
from ensembl.vep.models.genic.transcript import Transcript


@pytest.fixture
def transcript() -> Transcript:
    return Transcript(stable_id="ENST00000000001", start=1, end=100, strand=Strand.FORWARD)


@pytest.fixture
def cdna() -> cDNA:
    return cDNA(start=10, end=50, strand=Strand.FORWARD)


class TestcDNAInitialization:
    """Test cDNA initialization."""

    def test_initialization_without_keywords(self) -> None:
        """Test that initializing without keywords raises error."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = cDNA(start=10, end=50, strand=Strand.FORWARD)


class TestcDNAProperties:
    """Test cDNA properties."""

    def test_core_properties(self, cdna: cDNA) -> None:
        """Test that cDNA properties are set correctly."""
        assert cdna.start == 10
        assert cdna.end == 50
        assert cdna.strand == 1

    def test_externally_settable_properties(self, cdna: cDNA, transcript: Transcript) -> None:
        """Test that externally settable properties can be set and accessed."""
        assert cdna.transcript is None
        cdna.transcript = transcript
        assert cdna.transcript == transcript
