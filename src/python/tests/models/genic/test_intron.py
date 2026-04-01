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
"""Unit tests for Intron class."""

import pytest

from ensembl.vep.models.genic.intron import Intron
from ensembl.vep.models.genic.exon import Exon

@pytest.fixture
def exon():
    return Exon(
        start=1, 
        end=10, 
        strand=1, 
        stable_id="ENSE00000000001"
    )

@pytest.fixture
def intron():
    return Intron(
        start=11, 
        end=19, 
        strand=1, 
    )


class TestIntronInitialization:
    """Test Intron initialization."""

    def test_initialization_without_keywords(self):
        """Test that initializing without keywords raises error."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = Intron(11, end=19, strand=1)

class TestIntronProperties:
    """Test Intron properties."""

    def test_core_properties(self, intron):
        """Test that Intron core properties are set correctly."""
        assert intron.start == 11
        assert intron.end == 19
        assert intron.strand == 1

    def test_externally_settable_properties(self, intron, exon):
        """Test that externally settable properties can be set and accessed."""
        assert intron.next_exon is None
        intron.next_exon = exon
        assert intron.next_exon == exon

        assert intron.prev_exon is None
        intron.prev_exon = exon
        assert intron.prev_exon == exon
