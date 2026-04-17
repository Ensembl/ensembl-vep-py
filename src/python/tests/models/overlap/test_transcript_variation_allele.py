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
"""Unit tests for TranscriptVariationAllele class."""

import pytest

from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.models.variation import VariationAllele
from ensembl.vep.models.genic import Transcript


@pytest.fixture
def transcript() -> Transcript:
    return Transcript(start=1, end=100, strand=1, stable_id="ENST00000000001")


@pytest.fixture
def variation_allele() -> VariationAllele:
    return VariationAllele(
        start=50,
        end=150,
        strand=1,
        reference_allele="A",
        alternate_allele="T",
    )


@pytest.fixture
def transcript_variation_allele(transcript, variation_allele) -> TranscriptVariationAllele:
    return TranscriptVariationAllele(feature=transcript, overlapping_feature=variation_allele)


class TestTranscriptVariationAlleleInitialization:
    """Test TranscriptVariationAllele initialization."""

    def test_initialization_without_keywords(self, transcript, variation_allele) -> None:
        """Test that initializing without keywords raises error."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = TranscriptVariationAllele(transcript, variation_allele)


class TestTranscriptVariationAlleleProperties:
    """Test TranscriptVariationAllele properties."""

    def test_core_properties(self, transcript_variation_allele, transcript, variation_allele) -> None:
        """Test that feature and overlapping_feature properties are set correctly."""
        assert transcript_variation_allele.feature == transcript
        assert transcript_variation_allele.overlapping_feature == variation_allele

    def test_default_properties(self, transcript_variation_allele) -> None:
        """Test that default properties are set correctly."""
        assert transcript_variation_allele.is_reference is False

    def test_cdna_locus(self, transcript_variation_allele) -> None:
        """Test cdna_locus property."""
        pass
