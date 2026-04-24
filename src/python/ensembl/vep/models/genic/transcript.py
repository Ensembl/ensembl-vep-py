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
"""Transcript class."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING
from functools import cached_property
from threading import Lock
from attrs import define, field

from ensembl.vep.models import Feature
from ensembl.vep.models.genic import CDS, UTR
from ensembl.vep.common import Attribute, invalidates

if TYPE_CHECKING:
    from ensembl.vep.models.genic import Gene, Exon, Intron, cDNA, Translation

logger = logging.getLogger(__name__)


@define(kw_only=True, auto_attribs=True)
class Transcript(Feature):
    """
    Represents a transcript.

    Attributes:
        stable_id (str):
          Ensembl stable transcript identifier (e.g. "ENST00000366667").
        version (int | None):
          Version of the stable identifier.
        gene (Gene | None):
          The gene to which this transcript belongs.
        spliced_exons (list[Exon]):
          List of exons that are spliced together to form the transcript.
        introns (list[Intron]):
          List of introns in the transcript.
        translation (Translation | None):
          The translation associated with this transcript,
          if applicable.
        alternative_translations (list[Translation] | None):
          List of alternative translations associated with this transcript,
          if applicable.
        attributes (list[Attribute] | None):
          List of attributes associated with this transcript.
    """

    stable_id: str
    version: int | None = field(default=None, init=False)
    gene: Gene | None = field(default=None, init=False)
    spliced_exons: list[Exon] = field(
        factory=list,
        init=False,
        on_setattr=invalidates(
            "__end_exon_index",
            "__spliced_seq",
            "__cds",
            "__cdna_coding_start",
            "__cdna_coding_end",
            "__cds_start",
            "__cds_end",
            "__five_prime_utr",
            "__three_prime_utr",
            "__translated_seq",
        ),
    )
    introns: list[Intron] = field(factory=list, init=False)
    translation: Translation | None = field(
        default=None,
        init=False,
        on_setattr=invalidates(
            "__cds",
            "__cdna_coding_start",
            "__cdna_coding_end",
            "__cds_start",
            "__cds_end",
            "__five_prime_utr",
            "__three_prime_utr",
            "__translated_seq",
        ),
    )
    alternative_translations: list[Translation] | None = field(factory=list, init=False)
    attributes: list[Attribute] | None = field(factory=list, init=False)
    _lock: Lock = field(factory=Lock, init=False, repr=False)

    @property
    def start_exon_index(self) -> int:
        """Returns the index of the exon containing the start codon.

        Returns:
            int: Index of the exon containing the start codon.
        """
        return 0

    @property
    def end_exon_index(self) -> int:
        """Public accessor for __end_exon_index.

        Returns:
            int: Index of the exon containing the stop codon.
        """
        return self.__end_exon_index

    @cached_property
    def __end_exon_index(self) -> int:
        """Returns the index of the exon containing the stop codon.

        Returns:
            int: Index of the exon containing the stop codon.
        """
        return len(self.spliced_exons) - 1

    @property
    def spliced_seq(self) -> str | None:
        """Public accessor for __spliced_seq.

        Returns:
            The spliced sequence of the transcript,
            or None if sequence information is not available for all exons.
        """
        return self.__spliced_seq

    @cached_property
    def __spliced_seq(self) -> str | None:
        """
        The spliced sequence of the transcript,
        or None if sequence information is not available for all exons.

        Returns:
            The spliced sequence of the transcript,
            or None if sequence information is not available for all exons.
        """
        with self._lock:
            seq = ""
            warning_template = (
                "Exon %s has no sequence information,"
                " adding N instead of the exact sequence for this exon"
                " when constructing spliced sequence for transcript %s."
                " Spliced sequence may not be accurate."
            )

            for exon in self.spliced_exons:
                if exon.seq is None:
                    logger.warning(warning_template, exon.stable_id, self.stable_id)
                    seq += "N" * exon.length
                else:
                    seq += exon.seq
            return seq if seq else None

    @property
    def cds(self) -> CDS | None:
        """Public accessor for __cds.

        Returns:
            CDS | None: The CDS of the transcript, or None if not applicable.
        """
        return self.__cds

    @cached_property
    def __cds(self) -> CDS | None:
        """
        Returns the CDS associated with this transcript, or None if not applicable.

        Returns:
            CDS | None: The CDS associated with this transcript, or None if not applicable.
        """
        with self._lock:
            if self.translation is None:
                logger.warning(
                    "Cannot compute CDS for transcript %s because translation information is not available.",
                    self.stable_id,
                )
                return None

            seq = ""
            cds_start: int | None = None
            cds_end: int | None = None
            warning_template = (
                "Exon %s has no sequence information,"
                " adding N instead of the exact sequence for this exon"
                " when constructing CDS sequence for transcript %s."
                " Transcript sequence may not be accurate."
            )

            for idx, exon in enumerate(self.spliced_exons[self.start_exon_index : self.end_exon_index + 1]):
                if exon.seq is None:
                    logger.warning(warning_template, exon.stable_id, self.stable_id)
                    exon_coding_seq = "N" * exon.length
                else:
                    exon_coding_seq = exon.seq

                if idx == self.start_exon_index:
                    exon_coding_seq = exon_coding_seq[self.translation.exon_coding_start - 1 :]
                    if self.translation.exon_coding_start is not None:
                        cds_start = exon.start + self.translation.exon_coding_start - 1
                    else:
                        cds_start = exon.start
                if idx == self.end_exon_index:
                    exon_coding_seq = exon_coding_seq[: self.translation.exon_coding_end]
                    if self.translation.exon_coding_end is not None:
                        cds_end = exon.start + self.translation.exon_coding_end - 1
                    else:
                        cds_end = exon.end

                seq += exon_coding_seq

            if seq == "" or cds_start is None or cds_end is None:
                return None

            cds = CDS(start=cds_start, end=cds_end, strand=self.strand, slice=self.slice)
            cds.seq = seq
            cds.transcript = self
            cds.translation = self.translation

            return cds

    @cached_property
    def cdna(self) -> cDNA | None:
        """Returns the cDNA sequence for the transcript, or None if not applicable."""
        raise NotImplementedError("Feature is not implemented yet.")

    @property
    def cdna_coding_start(self) -> int | None:
        """Public accessor for __cdna_coding_start.

        Returns:
            int | None: cDNA coordinate of the start of the coding sequence, or None if not applicable.
        """
        return self.__cdna_coding_start

    @cached_property
    def __cdna_coding_start(self) -> int | None:
        """
        Returns the start of the coding sequence with respect to the cDNA, or None if not applicable.

        Returns:
            int | None: cDNA coordinate of the start of the coding sequence, or None if not applicable.
        """
        with self._lock:
            if self.translation is None:
                return None

            start = 1
            for idx, exon in enumerate(self.spliced_exons):
                if idx < self.start_exon_index:
                    # if non-translateable then add the whole exon length
                    start += exon.length
                else:
                    start += self.translation.exon_coding_start - 1
                    break

            return start

    @property
    def cdna_coding_end(self) -> int | None:
        """Public accessor for __cdna_coding_end.

        Returns:
            int | None: cDNA coordinate of the end of the coding sequence, or None if not applicable.
        """
        return self.__cdna_coding_end

    @cached_property
    def __cdna_coding_end(self) -> int | None:
        """
        Returns the end of the coding sequence with respect to the cDNA, or None if not applicable.

        Returns:
            int | None: cDNA coordinate of the end of the coding sequence, or None if not applicable.
        """
        with self._lock:
            if self.translation is None:
                return None

            end = 0
            for idx, exon in enumerate(self.spliced_exons):
                if idx < self.end_exon_index:
                    # if non-translateable then add the whole exon length
                    end += exon.length
                else:
                    if self.translation.exon_coding_end is None:
                        end += exon.length
                    else:
                        end += self.translation.exon_coding_end
                    break

            return end

    @property
    def cds_start(self) -> int | None:
        """Public accessor for __cds_start.

        Returns:
            Coordinate of the start of the CDS with respect to the transcript,
            or None if not applicable.
        """
        return self.__cds_start

    @cached_property
    def __cds_start(self) -> int | None:
        """
        Returns the start of the CDS with respect to the transcript, or None if not applicable.

        Returns:
            Coordinate of the start of the CDS with respect to the transcript,
            or None if not applicable.
        """
        with self._lock:
            if self.spliced_exons is None or len(self.spliced_exons) == 0 or self.translation is None:
                return None

            start_exon = self.spliced_exons[self.start_exon_index]
            # cds start relative to the slice
            cds_start = start_exon.start + self.translation.exon_coding_start - 1
            # cds start relative to the transcript
            return cds_start - self.start + 1

    @property
    def cds_end(self) -> int | None:
        """Public accessor for __cds_end.

        Returns:
            Coordinate of the end of the CDS with respect to the transcript,
            or None if not applicable.
        """
        return self.__cds_end

    @cached_property
    def __cds_end(self) -> int | None:
        """
        Returns the end of the CDS with respect to the transcript, or None if not applicable.

        Returns:
            Coordinate of the end of the CDS with respect to the transcript,
            or None if not applicable.
        """
        with self._lock:
            if self.spliced_exons is None or len(self.spliced_exons) == 0 or self.translation is None:
                return None

            end_exon = self.spliced_exons[self.end_exon_index]
            # cds end relative to the slice
            if self.translation.exon_coding_end is None:
                cds_end = end_exon.end
            else:
                cds_end = end_exon.start + self.translation.exon_coding_end - 1
            # cds end relative to the transcript
            return cds_end - self.start + 1

    @property
    def five_prime_utr(self) -> UTR | None:
        """Public accessor for __five_prime_utr.

        Returns:
            The 5' UTR of the transcript, or None if not applicable.
        """
        return self.__five_prime_utr

    @cached_property
    def __five_prime_utr(self) -> UTR | None:
        """Returns the 5' UTR of the transcript, or None if not applicable.

        Returns:
            The 5' UTR of the transcript, or None if not applicable.
        """
        with self._lock:
            if self.cds_start is None or self.cdna_coding_start is None:
                warning_template = (
                    "Cannot compute five_prime_utr for transcript %s"
                    " because cds_start or cdna_coding_start is None."
                )

                logger.warning(warning_template, self.stable_id)
                return None

            if self.spliced_seq is None or len(self.spliced_seq) < self.cdna_coding_start:
                warning_template = (
                    "Cannot compute five_prime_utr for transcript %s"
                    " because spliced_seq is None or too short."
                )

                logger.warning(warning_template, self.stable_id)
                return None

            end = self.cds_start - 1
            # spliced sequence is the same as cDNA sequence
            seq = self.spliced_seq[: self.cdna_coding_start - 1]
            start = end - len(seq) + 1

            utr = UTR(start=start, end=end, strand=self.strand, slice=self.slice)
            utr.transcript = self
            return utr

    @property
    def three_prime_utr(self) -> UTR | None:
        """Public accessor for __three_prime_utr.

        Returns:
            The 3' UTR of the transcript, or None if not applicable.
        """
        return self.__three_prime_utr

    @cached_property
    def __three_prime_utr(self) -> UTR | None:
        """Returns the 3' UTR of the transcript, or None if not applicable.

        Returns:
            The 3' UTR of the transcript, or None if not applicable.
        """
        with self._lock:
            if self.cds_end is None or self.cdna_coding_end is None:
                warning_template = (
                    "Cannot compute three_prime_utr for transcript %s"
                    " because cds_end or cdna_coding_end is None."
                )

                logger.warning(
                    warning_template,
                    self.stable_id,
                )
                return None

            if self.spliced_seq is None or len(self.spliced_seq) <= self.cdna_coding_end:
                warning_template = (
                    "Cannot compute three_prime_utr for transcript %s"
                    " because spliced_seq is None or too short."
                )

                logger.warning(
                    warning_template,
                    self.stable_id,
                )
                return None

            start = self.cds_end + 1
            # spliced sequence is the same as cDNA sequence
            seq = self.spliced_seq[self.cdna_coding_end :]
            end = start + len(seq) - 1

            utr = UTR(start=start, end=end, strand=self.strand, slice=self.slice)
            utr.transcript = self
            return utr

    @property
    def translated_seq(self) -> str | None:
        """Public accessor for __translated_seq.

        Returns:
            str | None: The translated sequence of the transcript, or None if not applicable.
        """
        return self.__translated_seq

    @cached_property
    def __translated_seq(self) -> str | None:
        """Returns the translated sequence of the transcript, or None if not applicable.

        Returns:
            str | None: The translated sequence of the transcript, or None if not applicable.
        """
        with self._lock:
            seq: str | None

            if self.translation is not None:
                seq = self.translation.seq
            elif self.cds is None:
                # to avoid pylint warning, # TBD: need codon table
                seq = "dummy"
            else:
                warning_template = (
                    "Cannot compute translated_seq for transcript %s"
                    " because translation and CDS sequence are both not available."
                )

                logger.warning(warning_template, self.stable_id)
                seq = None

            return seq

    def __eq__(self, other: object) -> bool:
        """Compare two Transcripts using various criteria.

        Args:
            other (object): Another object to compare against.


        Returns:
            bool: True if equals, False otherwise.
        """
        if not isinstance(other, Transcript):
            return False
        return self.stable_id == other.stable_id and self.version == other.version
