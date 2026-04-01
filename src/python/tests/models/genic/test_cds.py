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
"""Unit tests for CDS class."""

import pytest

from ensembl.vep.models.genic.cds import CDS
from ensembl.vep.models.genic.transcript import Transcript
from ensembl.vep.models.genic.translation import Translation

@pytest.fixture
def transcript():
    return Transcript(
        stable_id="ENST00000000001",
        start=1, 
        end=100, 
        strand=1
    )

@pytest.fixture
def translation(transcript):
    return Translation(
        start=1, 
        end=90,
        stable_id="ENSP00000000001"
    )

@pytest.fixture
def cds():
    return CDS(
        start=11, 
        end=40, 
        strand=1
    )

class TestCDSInitialization:
    """Test CDS initialization."""

    def test_initialization_without_keywords(self):
        """Test that initializing without keywords raises error."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = CDS(
                10, 
                end=50, 
                strand=1
            )


class TestCDSProperties:
    """Test CDS properties and cached properties."""

    def test_core_properties(self, cds, transcript, translation):
        """Test that CDS properties are set correctly."""
        assert cds.start == 11
        assert cds.end == 40
        assert cds.strand == 1

    def test_externally_settable_properties(self, cds, transcript, translation):
        """Test that externally settable properties can be set and accessed."""
        assert cds.transcript is None
        cds.transcript = transcript
        assert cds.transcript == transcript

        assert cds.translation is None
        cds.translation = translation
        assert cds.translation == translation

    def test_neucleotide_length(self, cds):
        """Test that CDS cached properties are accessed correctly."""
        assert cds.neucleotide_length is None

        cds.seq = "ATG" * 10
        assert cds.neucleotide_length == 30

    def test_protein_length(self, cds):
        """Test that CDS cached properties are accessed correctly."""
        # assert cds.protein_length is None

        # cds.seq = "ATG" * 10
        # assert cds.protein_length == 10

        # protein_length accesses neucleotide_length and gets into deadlock, workaround?
