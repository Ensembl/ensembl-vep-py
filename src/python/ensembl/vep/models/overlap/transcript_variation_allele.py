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
"""TranscriptVariationAllele class."""

from __future__ import annotations

from typing import TYPE_CHECKING
from functools import cached_property
from threading import RLock
from attr import define, field

from ensembl.vep.common import invalidates
from .feature_overlap import FeatureOverlap

if TYPE_CHECKING:
    from ensembl.vep.models.genic import Transcript
    from ensembl.vep.models.variation import VariationAllele
    from ensembl.vep.common import Locus


@define(kw_only=True, auto_attribs=True)
class TranscriptVariationAllele(FeatureOverlap):
    """
    Represents the overlap between a transcript and a variation allele.

    Attributes:
        feature (Transcript):
            The transcript that is getting overlapped.
        overlapping_feature (VariationAllele):
            The variation allele that is overlapping.
        is_reference (bool):
            Whether the overlapping variation allele is the reference allele.
        cdna_locus (Locus):
            The locus of the variation allele on the cDNA sequence of the transcript.
        cds_locus (Locus | None):
            The locus of the variation allele on the CDS sequence of the transcript.
        protein_locus (Locus | None):
            The locus of the variation allele on the protein sequence of the transcript.
        distance_to_transcript (int | None):
            The distance from the variation allele to the transcript,
            or None if the variation allele overlaps the transcript.
        exon_number (int | None):
            The exon number of the transcript that is overlapped by the variation allele,
            or None if the variation allele does not overlap any exon.
        intron_number (int | None):
            The intron number of the transcript that is overlapped by the variation allele,
            or None if the variation allele does not overlap any intron.
    """

    feature: Transcript
    overlapping_feature: VariationAllele
    is_reference: bool = field(default=False)
    cdna_locus: Locus = field(init=False, on_setattr=invalidates("__hgvs_transcript"))
    cds_locus: Locus | None = field(
        default=None, init=False, on_setattr=invalidates("__hgvs_transcript", "__hgvs_protein")
    )
    protein_locus: Locus | None = field(default=None, init=False, on_setattr=invalidates("__hgvs_protein"))
    distance_to_transcript: int | None = field(default=None, init=False)
    exon_number: int | None = field(default=None, init=False)
    intron_number: int | None = field(default=None, init=False)
    _lock: RLock = field(factory=RLock, init=False, repr=False)

    @property
    def shifted(self) -> bool:
        """
        Public accessor for __shifted.

        Returns:
            bool: True if the allele has been shifted, False otherwise.
        """
        return self.__shifted

    @cached_property
    def __shifted(self) -> bool:
        """
        Returns whether the allele has been shifted from its original position.

        Returns:
            bool: True if the allele has been shifted, False otherwise.
        """
        return False

    @property
    def shifted_cdna_locus(self) -> Locus | None:
        """
        Public accessor for __shifted_cdna_locus.

        Returns:
            The cDNA locus of the variation allele after shifting,
            or None if not applicable.
        """
        return self.__shifted_cdna_locus

    @cached_property
    def __shifted_cdna_locus(self) -> Locus | None:
        """
        Returns the cDNA locus of the variation allele after shifting, or None if not applicable.

        Returns:
            The cDNA locus of the variation allele after shifting,
            or None if not applicable.
        """
        return None

    @property
    def shifted_cds_locus(self) -> Locus | None:
        """
        Public accessor for __shifted_cds_locus.

        Returns:
            The CDS locus of the variation allele after shifting,
            or None if not applicable.
        """
        return self.__shifted_cds_locus

    @cached_property
    def __shifted_cds_locus(self) -> Locus | None:
        """
        Returns the CDS locus of the variation allele after shifting, or None if not applicable.

        Returns:
            LThe CDS locus of the variation allele after shifting,
            or None if not applicable.
        """
        return None

    @property
    def shifted_protein_locus(self) -> Locus | None:
        """
        Public accessor for __shifted_protein_locus.

        Returns:
            The protein locus of the variation allele after shifting,
            or None if not applicable.
        """
        return self.__shifted_protein_locus

    @cached_property
    def __shifted_protein_locus(self) -> Locus | None:
        """
        Returns the protein locus of the variation allele after shifting, or None if not applicable.

        Returns:
            The protein locus of the variation allele after shifting,
            or None if not applicable.
        """
        return self.__shifted_protein_locus

    @property
    def original_cdna_locus(self) -> Locus | None:
        """Public accessor for __original_cdna_locus.

        Returns:
            The original cDNA locus of the variation allele before shifting,
            or None if not applicable.
        """
        return self.__original_cdna_locus

    @cached_property
    def __original_cdna_locus(self) -> Locus | None:
        """Returns the original cDNA locus of the variation allele before shifting, or None if not applicable.

        Returns:
            The original cDNA locus of the variation allele before shifting,
            or None if not applicable.
        """
        return None

    @property
    def original_cds_locus(self) -> Locus | None:
        """Public accessor for __original_cds_locus.

        Returns:
            The original CDS locus of the variation allele before shifting,
            or None if not applicable.
        """
        return self.__original_cds_locus

    @cached_property
    def __original_cds_locus(self) -> Locus | None:
        """Returns the original CDS locus of the variation allele before shifting, or None if not applicable.

        Returns:
            The original CDS locus of the variation allele before shifting,
            or None if not applicable.
        """
        return None

    @property
    def original_protein_locus(self) -> Locus | None:
        """Public accessor for __original_protein_locus.

        Returns:
            The original protein locus of the variation allele before shifting,
            or None if not applicable.
        """
        return self.__original_protein_locus

    @cached_property
    def __original_protein_locus(self) -> Locus | None:
        """The original protein locus of the variation allele before shifting,
        or None if not applicable.

        Returns:
            The original protein locus of the variation allele before shifting,
            or None if not applicable.
        """
        raise NotImplementedError("Feature is not implemented yet.")

    def shift(self) -> None:
        """Shifts the variation allele to its most 3' position on the transcript according to HGVS rules."""
        raise NotImplementedError("Feature is not implemented yet.")

    def hgvs_genomic(self) -> str | None:
        """Returns the HGVS genomic notation for this variation allele."""
        raise NotImplementedError("Feature is not implemented yet.")

    def hgvs_transcript(self) -> str | None:
        """Returns the HGVS transcript notation for this variation allele."""
        raise NotImplementedError("Feature is not implemented yet.")

    def hgvs_protein(self) -> str | None:
        """Returns the HGVS protein notation for this variation allele."""
        raise NotImplementedError("Feature is not implemented yet.")
