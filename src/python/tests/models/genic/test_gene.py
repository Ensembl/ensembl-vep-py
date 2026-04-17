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
"""Unit tests for Gene class."""

import pytest

from ensembl.vep.common import Strand
from ensembl.vep.models.genic.gene import Gene
from ensembl.vep.models.genic.transcript import Transcript


@pytest.fixture
def transcript() -> Transcript:
    return Transcript(stable_id="ENST00000000001", start=1, end=100, strand=Strand.FORWARD)


@pytest.fixture
def gene() -> Gene:
    return Gene(stable_id="ENSG00000000001", start=10, end=200, strand=Strand.FORWARD)


class TestGeneInitialization:
    """Test Gene initialization."""

    def test_initialization_without_keywords(self) -> None:
        """Test that initializing without keywords raises error."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = Gene("ENSG00000000001", start=10, end=200, strand=Strand.FORWARD)


class TestGeneProperties:
    """Test property accessors."""

    def test_core_properties(self, gene) -> None:
        """Test that core properties are set correctly."""
        assert gene.start == 10
        assert gene.end == 200
        assert gene.strand == 1
        assert gene.stable_id == "ENSG00000000001"

    def test_externally_settable_properties(self, gene, transcript) -> None:
        """Test that externally settable properties can be set and accessed."""
        assert gene.version is None
        gene.version = 2
        assert gene.version == 2

        assert gene.name is None
        gene.name = "Gene1"
        assert gene.name == "Gene1"

        assert gene.symbol is None
        gene.symbol = "G1"
        assert gene.symbol == "G1"

        assert gene.alternative_names == []
        gene.alternative_names.append("GeneOne")
        assert gene.alternative_names == ["GeneOne"]

        assert gene.transcripts == []
        gene.transcripts.append(transcript)
        assert gene.transcripts == [transcript]
