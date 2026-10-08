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
from collections.abc import Mapping
from typing import TYPE_CHECKING, ClassVar
from functools import cached_property
from threading import RLock
from attrs import define, field

from ensembl.vep.models import Feature
from ensembl.vep.models.genic.cds import CDS
from ensembl.vep.models.genic.intron import Intron
from ensembl.vep.models.genic.translation import Translation
from ensembl.vep.models.genic.utr import UTR
from ensembl.vep.models.sequence import Mapper, MapperIndex
from ensembl.vep.common import Attribute, ExonPhase, Strand, invalidates, is_resolved
from ensembl.vep.common.codons import translate_dna

if TYPE_CHECKING:
    from ensembl.vep.models.genic import Gene, Exon, CDNA

logger = logging.getLogger(__name__)


_CODING_DERIVATIONS: tuple[str, ...] = (
    "_coding_start",
    "_coding_end",
    "_first_coding_exon_index",
    "_last_coding_exon_index",
    "_exon_coding_start",
    "_exon_coding_end",
    "_cds_start_phase",
    "_cdna_coding_start",
    "_cdna_coding_end",
    "_cds",
    "_five_prime_utr",
    "_three_prime_utr",
    "_coding_sequence",
    "_translated_seq",
    "_cds_genomic_spans",
    "_five_prime_utr_genomic_spans",
    "_three_prime_utr_genomic_spans",
    "_cds_mapper",
    "_cdna_mapper",
)
"""Everything downstream of a coding boundary. Seeding any coding-signal source makes all of
this stale. 

Here _coding_start/_coding_end are the first things to go, they are what every other entry below
is computed from, and leaving them cached is what made a seed applied after a read silently have no
effect.

And _cds and _translation are deliberately absent as an explicitly seeded translation must survive
a later cds seed and vice versa, the same exception ``spliced_exons``' own invalidates() list makes
for it.
"""


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
        biotype (str | None):
          Biotype of the transcript (e.g. "protein_coding").
        spliced_exons (list[Exon]):
          List of exons that are spliced together to form the transcript, in transcript (5' to
          3') order, therefore, index 0 is always the first exon.
        introns (list[Intron]):
          Introns between consecutive spliced exons, in transcript order. Derived from
          ``spliced_exons`` on first access rather than loaded; seed with ``seed_attrs``.
        translation (Translation | None):
          The translation associated with this transcript, explicitly assigned or derived from
          whatever coding-boundary information is available.
        alternative_translations (list[Translation] | None):
          List of alternative translations associated with this transcript,
          if applicable.
        attributes (list[Attribute] | None):
          List of attributes associated with this transcript.
    """

    stable_id: str
    version: int | None = field(default=None, init=False)
    gene: Gene | None = field(default=None, init=False)
    biotype: str | None = field(default=None, init=False)
    spliced_exons: list[Exon] = field(
        factory=list,
        init=False,
        on_setattr=invalidates(
            *_CODING_DERIVATIONS,
            "_spliced_seq",
            "_translated_seq",
            "_introns",
        ),
    )
    alternative_translations: list[Translation] | None = field(factory=list, init=False)
    attributes: list[Attribute] | None = field(factory=list, init=False)
    _lock: RLock = field(factory=RLock, init=False, repr=False)

    SEED_INVALIDATES: ClassVar[Mapping[str, tuple[str, ...]]] = {
        "cds": _CODING_DERIVATIONS,
        "translation": _CODING_DERIVATIONS,
        "five_prime_utr": _CODING_DERIVATIONS,
        "three_prime_utr": _CODING_DERIVATIONS,
    }
    """Seeding any coding-signal source invalidates the whole coding derivation closure."""

    @property
    def introns(self) -> list[Intron]:
        """Public accessor for _introns.

        Returns:
            list[Intron]: The introns between consecutive spliced exons, in transcript (5' to 3')
                order.
        """
        return self._introns

    @cached_property
    def _introns(self) -> list[Intron]:
        """Derive one intron per pair of consecutive spliced exons.

        Every adjacent pair gets an intron, including a pair with no genomic gap between them --
        such an intron has ``start > end`` and records adjacency only; callers that care about gap
        size check that themselves.

        Returns:
            list[Intron]: The introns, in transcript (5' to 3') order.
        """
        with self._lock:
            introns: list[Intron] = []
            for prev_exon, next_exon in zip(self.spliced_exons, self.spliced_exons[1:]):
                first, second = (
                    (next_exon, prev_exon) if self.strand == Strand.REVERSE else (prev_exon, next_exon)
                )
                intron = Intron(
                    start=first.end + 1, end=second.start - 1, strand=self.strand, slice=self.slice
                )
                intron.assign_without_hooks(prev_exon=prev_exon, next_exon=next_exon)
                introns.append(intron)
            return introns

    # --------------------------------------------------------------------------
    # Coding boundaries derivation:
    # If a transcript is protein coding, then we should have coding boundaries. The coding boundary is
    # captured by coding_start and coding_end and can be derived from either of the following:
    #       - CDS (attribute: cds)
    #       - UTRs (attributes: five_prime_utr and three_prime_utr)
    #       - Translation (attribute: translation)
    #       - first and last coding exon and coding start and end relative to the respective exon (attributes:
    # first_coding_exon_index, last_coding_exon_index, exon_coding_start and exon_coding_end)
    #
    # Each of the above set of attributes can be derived from any of the other set but one of them must be
    # seeded. Otherwise the transcript would not have enough information to determine its coding boundaries
    # and be treated as non coding.
    #
    # The mutual dependency here (exon_coding_start can derive from cds, cds derives from
    # exon_coding_start) never actually recurses: because the point of contention is in coding_start and
    # coding_end function body and there we only peek at the underlying slot using is_resolved() and never
    # instigates a computation which can recurse back to the same property.
    #
    # If any of the attributes directly reference the other attribute in their function body, then
    # the other attribute should not also directly reference the first attribute in its function body. In
    # such case is_resolved() should be used to break recursive calls.
    # --------------------------------------------------------------------------

    @property
    def coding_start(self) -> int | None:
        """Public accessor for _coding_start.

        Returns:
            int | None: Start of the coding region or None if not applicable.
        """
        return self._coding_start

    @cached_property
    def _coding_start(self) -> int | None:
        """
        Returns the start of the coding region.

        Returns:
            int | None: Start of the coding region or None if not applicable.
        """
        with self._lock:
            if is_resolved(self, "_cds") and self.cds is not None:
                cds = self.cds
                return cds.start if self.strand == Strand.FORWARD else cds.end

            if is_resolved(self, "_five_prime_utr") and self.five_prime_utr is not None:
                five_prime_utr = self.five_prime_utr
                return five_prime_utr.end + 1 if self.strand == Strand.FORWARD else five_prime_utr.start - 1

            if is_resolved(self, "_translation") and self.translation is not None:
                translation = self.translation
                return (
                    self.start + translation.start - 1
                    if self.strand == Strand.FORWARD
                    else self.end - translation.start + 1
                )

            if is_resolved(self, "_first_coding_exon_index") and is_resolved(self, "_exon_coding_start"):
                first_coding_exon_index = self.first_coding_exon_index
                exon_coding_start = self.exon_coding_start
                if not self.spliced_exons or first_coding_exon_index is None or exon_coding_start is None:
                    return None
                start_exon = self.spliced_exons[first_coding_exon_index]

                return (
                    start_exon.start + exon_coding_start - 1
                    if self.strand == Strand.FORWARD
                    else start_exon.end - exon_coding_start + 1
                )

            # No coding signal from any source.
            return None

    @property
    def coding_end(self) -> int | None:
        """Public accessor for _coding_end.

        Returns:
            int | None: End of the coding region, or None if not applicable.
        """
        return self._coding_end

    @cached_property
    def _coding_end(self) -> int | None:
        """
        Returns the end of the coding region.

        Returns:
            int | None: End of the coding region or None if not applicable.
        """
        with self._lock:
            if is_resolved(self, "_cds") and self.cds is not None:
                cds = self.cds
                return cds.end if self.strand == Strand.FORWARD else cds.start

            if is_resolved(self, "_three_prime_utr") and self.three_prime_utr is not None:
                three_prime_utr = self.three_prime_utr
                return three_prime_utr.start - 1 if self.strand == Strand.FORWARD else three_prime_utr.end + 1

            if is_resolved(self, "_translation") and self.translation is not None:
                translation = self.translation
                return (
                    self.start + translation.end - 1
                    if self.strand == Strand.FORWARD
                    else self.end - translation.end + 1
                )

            if is_resolved(self, "_last_coding_exon_index") and is_resolved(self, "_exon_coding_end"):
                last_coding_exon_index = self.last_coding_exon_index
                exon_coding_end = self.exon_coding_end
                if not self.spliced_exons or last_coding_exon_index is None or exon_coding_end is None:
                    return None
                end_exon = self.spliced_exons[last_coding_exon_index]

                return (
                    end_exon.start + exon_coding_end - 1
                    if self.strand == Strand.FORWARD
                    else end_exon.end - exon_coding_end + 1
                )

            # No coding signal from any source.
            return None

    @property
    def first_coding_exon_index(self) -> int | None:
        """Public accessor for _first_coding_exon_index.

        Returns:
            int | None: Index into spliced_exons of the exon where coding starts,
                or None if the transcript has no coding signal at all.
        """
        return self._first_coding_exon_index

    @cached_property
    def _first_coding_exon_index(self) -> int | None:
        """Derives the first coding exon's index.

        Returns:
            int | None: Index into spliced_exons, or None if not coding.
        """
        with self._lock:
            if self.coding_start is not None:
                target_position = self.coding_start
                for exon_index, exon in enumerate(self.spliced_exons):
                    if exon.start <= target_position <= exon.end:
                        return exon_index
            return None

    @property
    def last_coding_exon_index(self) -> int | None:
        """Public accessor for _last_coding_exon_index.

        Returns:
            int | None: Index into spliced_exons of the exon where coding ends, or None
                if the transcript has no coding signal at all.
        """
        return self._last_coding_exon_index

    @cached_property
    def _last_coding_exon_index(self) -> int | None:
        """Derives the last coding exon's index.

        Returns:
            int | None: Index into spliced_exons, or None if not coding.
        """
        with self._lock:
            if self.coding_end is not None:
                target_position = self.coding_end
                for exon_index, exon in enumerate(self.spliced_exons):
                    if exon.start <= target_position <= exon.end:
                        return exon_index
            return None

    @property
    def exon_coding_start(self) -> int | None:
        """Public accessor for _exon_coding_start.

        Returns:
            int | None: 1-based offset, within the first coding exon, of the first coding base,
                or None if the transcript has no coding signal at all.
        """
        return self._exon_coding_start

    @cached_property
    def _exon_coding_start(self) -> int | None:
        """Derives the coding start offset from the first coding exon.

        Returns:
            int | None: 1-based exon-relative offset, or None if not coding.
        """
        with self._lock:
            if self.coding_start is not None and self.first_coding_exon_index is not None:
                target_position = self.coding_start
                start_exon = self.spliced_exons[self.first_coding_exon_index]
                return (
                    (target_position - start_exon.start + 1)
                    if self.strand == Strand.FORWARD
                    else (start_exon.end - target_position + 1)
                )
            return None

    @property
    def exon_coding_end(self) -> int | None:
        """Public accessor for _exon_coding_end.

        Returns:
            int | None: 1-based offset, within the last coding exon, of the last coding base, or
                None if the transcript has no coding signal at all.
        """
        return self._exon_coding_end

    @cached_property
    def _exon_coding_end(self) -> int | None:
        """Derives the coding end offset from the last coding exon.

        Returns:
            int | None: 1-based exon-relative offset, or None if not coding.
        """
        with self._lock:
            if self.coding_end is not None and self.last_coding_exon_index is not None:
                target_position = self.coding_end
                end_exon = self.spliced_exons[self.last_coding_exon_index]
                return (
                    (target_position - end_exon.start + 1)
                    if self.strand == Strand.FORWARD
                    else (end_exon.end - target_position + 1)
                )
            return None

    @property
    def cds_start_phase(self) -> int:
        """Public accessor for _cds_start_phase.

        Returns:
            int: Bases of the first codon missing from the CDS's 5' end, 0 when it is complete.
        """
        return self._cds_start_phase

    @cached_property
    def _cds_start_phase(self) -> int:
        """Derive the start phase from the first coding exon, as VEP does.

        An exon's phase counts the bases of its first codon that lie before it. On the first coding
        exon that can only mean the CDS begins part-way through a codon, and VEP pads the
        translateable sequence with that many ``N`` (``Bio::EnsEMBL::Transcript::translateable_seq``).
        No phase, or an unknown one, means a complete 5' end.

        Returns:
            int: 0, 1 or 2.
        """
        with self._lock:
            first_coding_exon_index = self.first_coding_exon_index
            if first_coding_exon_index is None:
                return 0
            phase = self.spliced_exons[first_coding_exon_index].phase
            if phase in (ExonPhase.ONE, ExonPhase.TWO):
                return phase.value
            return 0

    @property
    def cdna_coding_start(self) -> int | None:
        """Public accessor for _cdna_coding_start.

        Returns:
            int | None: cDNA coordinate of the start of the coding sequence, or None if not applicable.
        """
        return self._cdna_coding_start

    @cached_property
    def _cdna_coding_start(self) -> int | None:
        """
        Returns the start of the coding sequence with respect to the cDNA, or None if not applicable.

        Returns:
            int | None: cDNA coordinate of the start of the coding sequence, or None if not applicable.
        """
        with self._lock:
            first_coding_exon_index = self.first_coding_exon_index
            exon_coding_start = self.exon_coding_start
            if first_coding_exon_index is None or exon_coding_start is None:
                return None

            start = 1
            for exon_index, exon in enumerate(self.spliced_exons):
                if exon_index < first_coding_exon_index:
                    # if non-translateable then add the whole exon length
                    start += exon.length
                else:
                    start += exon_coding_start - 1
                    break

            return start

    @property
    def cdna_coding_end(self) -> int | None:
        """Public accessor for _cdna_coding_end.

        Returns:
            int | None: cDNA coordinate of the end of the coding sequence, or None if not applicable.
        """
        return self._cdna_coding_end

    @cached_property
    def _cdna_coding_end(self) -> int | None:
        """
        Returns the end of the coding sequence with respect to the cDNA, or None if not applicable.

        Returns:
            int | None: cDNA coordinate of the end of the coding sequence, or None if not applicable.
        """
        with self._lock:
            last_coding_exon_index = self.last_coding_exon_index
            if last_coding_exon_index is None:
                return None
            exon_coding_end = self.exon_coding_end

            end = 0
            for exon_index, exon in enumerate(self.spliced_exons):
                if exon_index < last_coding_exon_index:
                    # if non-translateable then add the whole exon length
                    end += exon.length
                else:
                    if exon_coding_end is None:
                        end += exon.length
                    else:
                        end += exon_coding_end
                    break

            return end

    @property
    def cds(self) -> CDS | None:
        """Public accessor for _cds.

        Returns:
            CDS | None: The CDS of the transcript, or None if not applicable.
        """
        return self._cds

    @cached_property
    def _cds(self) -> CDS | None:
        """
        Returns the CDS associated with this transcript, or None if not applicable.

        Returns:
            CDS | None: The CDS associated with this transcript, or None if not applicable.
        """
        with self._lock:
            genomic_spans = self.cds_genomic_spans
            first_coding_exon_index = self.first_coding_exon_index
            last_coding_exon_index = self.last_coding_exon_index
            if genomic_spans is None or first_coding_exon_index is None or last_coding_exon_index is None:
                logger.warning(
                    f"Cannot compute CDS for transcript {self.stable_id} because no coding"
                    " boundary information is available."
                )
                return None

            coding_exons = self.spliced_exons[first_coding_exon_index : last_coding_exon_index + 1]
            if len(coding_exons) != len(genomic_spans):
                # An exon clipped to nothing means the coding offsets are inconsistent with the exons.
                logger.warning(
                    f"Cannot compute CDS for transcript {self.stable_id} because its coding boundaries"
                    " do not fit its exons."
                )
                return None

            blocks: list[str] = []
            for exon, (block_start, block_end) in zip(coding_exons, genomic_spans):
                if exon.seq is None:
                    logger.warning(
                        f"Exon {exon.stable_id} has no sequence information,"
                        f" adding N instead of the exact sequence for this exon"
                        f" when constructing CDS sequence for transcript {self.stable_id}."
                        f" Transcript sequence may not be accurate."
                    )
                    exon_seq = "N" * exon.length
                else:
                    exon_seq = exon.seq

                # exon sequence is already reverse complemented
                if self.strand == Strand.REVERSE:
                    blocks.append(exon_seq[exon.end - block_end : exon.end - block_start + 1])
                else:
                    blocks.append(exon_seq[block_start - exon.start : block_end - exon.start + 1])

            # Blocks are in transcript order and the sequence has already been reverse complemented
            # so we can just join the sequences together
            seq = "".join(blocks)

            cds_start = min(block_start for block_start, _ in genomic_spans)
            cds_end = max(block_end for _, block_end in genomic_spans)
            cds = CDS(start=cds_start, end=cds_end, strand=self.strand, slice=self.slice)
            cds.seq = seq
            cds.transcript = self
            if is_resolved(self, "_translation"):
                cds.translation = self.translation

            return cds

    @cached_property
    def cdna(self) -> CDNA | None:
        """Returns the cDNA sequence for the transcript, or None if not applicable."""
        raise NotImplementedError("Feature is not implemented yet.")

    @property
    def five_prime_utr(self) -> UTR | None:
        """Public accessor for _five_prime_utr.

        Returns:
            The 5' UTR of the transcript, or None if not applicable.
        """
        return self._five_prime_utr

    @cached_property
    def _five_prime_utr(self) -> UTR | None:
        """Returns the 5' UTR of the transcript, or None if not applicable.

        Returns:
            UTR | None: The 5' UTR, or None if the transcript is not coding or has no 5' UTR.
        """
        with self._lock:
            return self._utr_envelope(five_prime=True)

    @property
    def three_prime_utr(self) -> UTR | None:
        """Public accessor for _three_prime_utr.

        Returns:
            The 3' UTR of the transcript, or None if not applicable.
        """
        return self._three_prime_utr

    @cached_property
    def _three_prime_utr(self) -> UTR | None:
        """Returns the 3' UTR of the transcript, or None if not applicable.

        Returns:
            UTR | None: The 3' UTR, or None if the transcript is not coding or has no 3' UTR.
        """
        with self._lock:
            return self._utr_envelope(five_prime=False)

    def _utr_envelope(self, *, five_prime: bool) -> UTR | None:
        """Build the genomic envelope of one UTR.

        We are calling it envelope because it is just the genomic span of the UTR and does not capture
        full information. The actual UTR may span multiple exons, and therefore there are segments in the
        between the envelope's range which is not part of the UTR but an intron.

        Args:
            five_prime (bool): True if we want the 5' UTR, False if we want the 3' UTR.

        Returns:
            UTR | None: The UTR feature, or None when the transcript is not coding or the
                requested side is empty (the CDS runs to that end of the transcript).
        """
        if not self.spliced_exons:
            return None

        if five_prime:
            if self.coding_start is None:
                return None
            if self.strand == Strand.FORWARD:
                start, end = self.spliced_exons[0].start, self.coding_start - 1
            else:
                start, end = self.coding_start + 1, self.spliced_exons[0].end
        else:
            if self.coding_end is None:
                return None
            if self.strand == Strand.FORWARD:
                start, end = self.coding_end + 1, self.spliced_exons[-1].end
            else:
                start, end = self.spliced_exons[-1].start, self.coding_end - 1

        if start > end:
            # The coding boundary runs all the way to this side of the transcript -- no UTR here.
            logger.debug(
                f"No {'five' if five_prime else 'three'}_prime_utr for transcript"
                f" {self.stable_id}: the coding region runs to that end of the transcript."
            )
            return None

        utr = UTR(start=start, end=end, strand=self.strand, slice=self.slice)
        utr.transcript = self
        return utr

    @property
    def translation(self) -> Translation | None:
        """Public accessor for _translation.

        Returns:
            Translation | None: The transcript's translation, seeded or derived, or None if the
                transcript has no coding signal at all.
        """
        return self._translation

    @cached_property
    def _translation(self) -> Translation | None:
        """Returns the translation associated with this transcript, or None if not applicable.

        Returns:
            Translation | None: The derived translation, or None if not coding.
        """
        if self.coding_start is None or self.coding_end is None:
            return None
        coding_start, coding_end = self.coding_start, self.coding_end

        # 1-based, measured from this transcript's own 5' end -- the exact inverse of
        # _coding_start/_coding_end's translation-fallback branches above.
        if self.strand == Strand.FORWARD:
            translation_start = coding_start - self.start + 1
            translation_end = coding_end - self.start + 1
        else:
            translation_start = self.end - coding_start + 1
            translation_end = self.end - coding_end + 1

        exon_coding_start = self.exon_coding_start if is_resolved(self, "_exon_coding_start") else None
        exon_coding_end = self.exon_coding_end if is_resolved(self, "_exon_coding_end") else None

        translation = Translation(
            start=translation_start,
            end=translation_end,
            stable_id=f"{self.stable_id}_derived_translation",
            exon_coding_start=exon_coding_start,
            exon_coding_end=exon_coding_end,
        )
        translation.transcript = self
        coding_sequence = self.coding_sequence
        if coding_sequence:
            translation.seq = translate_dna(coding_sequence)
            translation.derived_from_coding_sequence = True
        return translation

    @property
    def spliced_seq(self) -> str | None:
        """Public accessor for _spliced_seq.

        Returns:
            The spliced sequence of the transcript,
            or None if sequence information is not available for all exons.
        """
        return self._spliced_seq

    @cached_property
    def _spliced_seq(self) -> str | None:
        """
        The spliced sequence of the transcript,
        or None if sequence information is not available for all exons.

        Returns:
            The spliced sequence of the transcript,
            or None if sequence information is not available for all exons.
        """
        with self._lock:
            seq = ""

            for exon in self.spliced_exons:
                if exon.seq is None:
                    logger.warning(
                        f"Exon {exon.stable_id} has no sequence information,"
                        f" adding N instead of the exact sequence for this exon"
                        f" when constructing spliced sequence for transcript {self.stable_id}."
                        f" Spliced sequence may not be accurate."
                    )
                    seq += "N" * exon.length
                else:
                    seq += exon.seq
            return seq if seq else None

    @property
    def coding_sequence(self) -> str | None:
        """Public accessor for _coding_sequence.

        Returns:
            str | None: The full CDS nucleotide sequence in transcript-sense, or
                None if the transcript is not coding or sequence is unavailable.
        """
        return self._coding_sequence

    @cached_property
    def _coding_sequence(self) -> str | None:
        """Returns the full CDS nucleotide sequence in transcript-sense.

        Prefixed with ``cds_start_phase`` ``N``, so that it lines up with ``cds_mapper`` offsets.

        Returns:
            str | None: CDS nucleotide sequence in transcript-sense orientation, or None if the
                transcript is not coding or no sequence is available from either source.
        """
        with self._lock:
            padding = "N" * self.cds_start_phase
            cds = self.cds
            if cds is not None and cds.seq is not None:
                return padding + cds.seq

            genomic_spans = self.cds_genomic_spans
            if genomic_spans is None or self.slice is None:
                return None

            blocks: list[str] = [padding]
            for block_start, block_end in genomic_spans:
                block = self.slice.subseq(block_start, block_end, self.strand)
                if block is None:
                    return None
                blocks.append(block)
            return "".join(blocks)

    @property
    def translated_seq(self) -> str | None:
        """Public accessor for _translated_seq.

        Returns:
            str | None: The translated sequence of the transcript, or None if not applicable.
        """
        return self._translated_seq

    @cached_property
    def _translated_seq(self) -> str | None:
        """Returns the translated sequence of the transcript, or None if not applicable.

        Returns:
            str | None: The translated sequence of the transcript, or None if not applicable.
        """
        with self._lock:
            if self.translation is not None and self.translation.seq is not None:
                return self.translation.seq

            coding_sequence = self.coding_sequence
            if coding_sequence is None:
                logger.warning(
                    f"Cannot compute translated_seq for transcript {self.stable_id} because "
                    f"translation and CDS sequence are both not available."
                )
                return None

            return translate_dna(coding_sequence)

    # ----------------------------------------------------------------------
    # Genomic <-> CDS and Genomic <-> cDNA mapping.
    # ----------------------------------------------------------------------

    @property
    def cds_genomic_spans(self) -> list[tuple[int, int]] | None:
        """Public accessor for _cds_genomic_spans.

        Returns:
            list[tuple[int, int]] | None: Transcript (5' to 3') order exon blocks clipped to the
                CDS, or None if the transcript is not coding.
        """
        return self._cds_genomic_spans

    @cached_property
    def _cds_genomic_spans(self) -> list[tuple[int, int]] | None:
        """Returns the exon blocks clipped to the CDS, in transcript (5' to 3') order.

        Both endpoints of every block are always given low-to-high regardless of strand, but the
        list itself is in transcript order, matching ``spliced_exons``' own convention.

        Returns:
            list[tuple[int, int]] | None: Transcript-order CDS-clipped exon blocks, or None if
                the transcript has no coding signal or no exons.
        """
        with self._lock:
            first_coding_exon_index, exon_coding_start = self.first_coding_exon_index, self.exon_coding_start
            last_coding_exon_index, exon_coding_end = self.last_coding_exon_index, self.exon_coding_end
            if (
                first_coding_exon_index is None
                or last_coding_exon_index is None
                or exon_coding_start is None
                or exon_coding_end is None
                or not self.spliced_exons
            ):
                return None

            genomic_blocks: list[tuple[int, int]] = []
            for exon_index, exon in enumerate(self.spliced_exons):
                if exon_index < first_coding_exon_index or exon_index > last_coding_exon_index:
                    continue

                block_start, block_end = exon.start, exon.end

                if exon_index == first_coding_exon_index:
                    if self.strand == Strand.REVERSE:
                        block_end = exon.end - (exon_coding_start - 1)
                    else:
                        block_start = exon.start + exon_coding_start - 1

                if exon_index == last_coding_exon_index:
                    if self.strand == Strand.REVERSE:
                        block_start = exon.end - exon_coding_end + 1
                    else:
                        block_end = exon.start + exon_coding_end - 1

                if block_start <= block_end:
                    genomic_blocks.append((block_start, block_end))

            if not genomic_blocks:
                return None

            return genomic_blocks

    @property
    def five_prime_utr_genomic_spans(self) -> list[tuple[int, int]] | None:
        """Public accessor for _five_prime_utr_genomic_spans.

        Returns:
            list[tuple[int, int]] | None: Transcript (5' to 3') order exon blocks clipped to the
                5' UTR, or None if the transcript is not coding or has no 5' UTR.
        """
        return self._five_prime_utr_genomic_spans

    @cached_property
    def _five_prime_utr_genomic_spans(self) -> list[tuple[int, int]] | None:
        """Returns the exon blocks clipped to the 5' UTR, in transcript (5' to 3') order.

        The mirror of cds_genomic_spans for the region upstream of the CDS in transcript sense:
        every exon strictly before the first coding exon, plus the portion of the first coding
        exon itself that precedes the CDS. Both endpoints of every block are always given
        low-to-high regardless of strand, matching cds_genomic_spans' own convention.

        Returns:
            list[tuple[int, int]] | None: Transcript-order 5' UTR exon blocks, or None if the
                transcript has no coding signal or no 5' UTR.
        """
        with self._lock:
            if not self.cds_genomic_spans:
                return None
            first_coding_exon_index, exon_coding_start = self.first_coding_exon_index, self.exon_coding_start
            if first_coding_exon_index is None or exon_coding_start is None:
                return None

            blocks: list[tuple[int, int]] = [
                (exon.start, exon.end) for exon in self.spliced_exons[:first_coding_exon_index]
            ]

            if exon_coding_start > 1:
                first_coding_exon = self.spliced_exons[first_coding_exon_index]
                if self.strand == Strand.REVERSE:
                    blocks.append((first_coding_exon.end - exon_coding_start + 2, first_coding_exon.end))
                else:
                    blocks.append((first_coding_exon.start, first_coding_exon.start + exon_coding_start - 2))

            return blocks if blocks else None

    @property
    def three_prime_utr_genomic_spans(self) -> list[tuple[int, int]] | None:
        """Public accessor for _three_prime_utr_genomic_spans.

        Returns:
            list[tuple[int, int]] | None: Transcript (5' to 3') order exon blocks clipped to the
                3' UTR, or None if the transcript is not coding or has no 3' UTR.
        """
        return self._three_prime_utr_genomic_spans

    @cached_property
    def _three_prime_utr_genomic_spans(self) -> list[tuple[int, int]] | None:
        """Returns the exon blocks clipped to the 3' UTR, in transcript (5' to 3') order.

        The mirror of _five_prime_utr_genomic_spans for the region downstream of the CDS in
        transcript sense: the portion of the last coding exon itself that follows the CDS, plus
        every exon strictly after it.

        Returns:
            list[tuple[int, int]] | None: Transcript-order 3' UTR exon blocks, or None if the
                transcript has no coding signal or no 3' UTR.
        """
        with self._lock:
            if not self.cds_genomic_spans:
                return None
            last_coding_exon_index, exon_coding_end = self.last_coding_exon_index, self.exon_coding_end
            if last_coding_exon_index is None or exon_coding_end is None:
                return None

            blocks: list[tuple[int, int]] = []
            if exon_coding_end is not None:
                last_coding_exon = self.spliced_exons[last_coding_exon_index]
                if self.strand == Strand.REVERSE:
                    block_end = last_coding_exon.end - exon_coding_end
                    if last_coding_exon.start <= block_end:
                        blocks.append((last_coding_exon.start, block_end))
                else:
                    block_start = last_coding_exon.start + exon_coding_end
                    if block_start <= last_coding_exon.end:
                        blocks.append((block_start, last_coding_exon.end))

            blocks.extend((exon.start, exon.end) for exon in self.spliced_exons[last_coding_exon_index + 1 :])

            return blocks if blocks else None

    @property
    def cds_mapper(self) -> Mapper | None:
        """Public accessor for _cds_mapper.

        Offsets are **CDS-relative**: offset 0 is the first base of the start codon, and UTR bases
        map to nothing at all. A CDS with an incomplete 5' end starts at offset ``cds_start_phase``
        instead, the offsets before it being VEP's ``N`` padding. Use this for anything expressed in
        coding coordinates, e.g., codon windows, protein offsets. Its sibling ``cdna_mapper`` covers
        the whole spliced transcript instead; the two disagree by ``cdna_coding_start - 1`` over the
        CDS, and only ``cdna_mapper`` can place a UTR position.

        Returns:
            Mapper | None: The genomic<->CDS<->protein coordinate mapper for this transcript's
                CDS, or None if the transcript is not coding.
        """
        return self._cds_mapper

    @cached_property
    def _cds_mapper(self) -> Mapper | None:
        """Builds the Mapper for this transcript's CDS, or None if it has none.

        Returns:
            Mapper | None: The coordinate mapper, or None if the transcript is not coding.
        """
        genomic_spans = self.cds_genomic_spans
        if genomic_spans is None or self.cds_length is None:
            return None
        return Mapper(
            genomic_spans=genomic_spans,
            strand=self.strand,
            coding_start_offset=0,
            # Offsets in this mapper start after the padding, so the last coding one is that far
            # along from cds_length's own unpadded count.
            coding_end_offset=self.cds_start_phase + self.cds_length - 1,
            padding=self.cds_start_phase,
        )

    @property
    def cdna_mapper(self) -> Mapper | None:
        """Public accessor for _cdna_mapper.

        Offsets are **cDNA-relative**: offset 0 is the transcript's first base, so 5' UTR, CDS and
        3' UTR positions all map. Use this to place a variant anywhere on the mature transcript.

        Returns:
            Mapper | None: The genomic<->cDNA coordinate mapper for this transcript's whole
                spliced sequence (5' UTR + CDS + 3' UTR), or None if it has no exons.
        """
        return self._cdna_mapper

    @cached_property
    def _cdna_mapper(self) -> Mapper | None:
        """Builds the Mapper for this transcript's whole spliced sequence, or None if it has no
        exons.

        The whole-transcript sibling of ``cds_mapper``, same class and mechanism, just built from
        every spliced exon rather than the CDS-clipped spans.

        A non-coding transcript still gets a mapper, since a cDNA position is meaningful without a
        CDS -- VEP reports one for a non-coding transcript variant. Its coding range is left empty
        (``coding_end_offset`` below ``coding_start_offset``), so the codon and protein methods
        answer None rather than reading a frame that does not exist.

        Returns:
            Mapper | None: The coordinate mapper, or None if the transcript has no exons.
        """
        if not self.spliced_exons:
            return None
        genomic_spans = [(exon.start, exon.end) for exon in self.spliced_exons]
        cdna_coding_start, cdna_coding_end = self.cdna_coding_start, self.cdna_coding_end
        if cdna_coding_start is None or cdna_coding_end is None:
            # No coding region: an empty range, which the codon and protein methods answer None to.
            cdna_coding_start, cdna_coding_end = 0, -1
        return Mapper(
            genomic_spans=genomic_spans,
            strand=self.strand,
            coding_start_offset=cdna_coding_start,
            coding_end_offset=cdna_coding_end,
        )

    @property
    def cds_length(self) -> int | None:
        """Number of nucleotides in the CDS.

        Returns:
            int | None: Number of nucleotides in the CDS, or None if the transcript is not coding.
        """
        genomic_spans = self.cds_genomic_spans
        if genomic_spans is None:
            return None
        return sum(end - start + 1 for start, end in genomic_spans)

    def genomic_to_cds_offset(
        self, genomic_position: int, *, accelerator: MapperIndex | None = None
    ) -> int | None:
        """Map a genomic position to its 0-based CDS offset.

        Args:
            genomic_position (int): A genomic position, 1-based.
            accelerator (MapperIndex | None): Optional O(1)-lookup accelerator for the mapping.

        Returns:
            int | None: The 0-based offset into the CDS, or None if the position is not coding.
        """
        mapper = self.cds_mapper
        return None if mapper is None else mapper.genomic_to_offset(genomic_position, accelerator=accelerator)

    def cds_offset_to_genomic(self, cds_offset: int, *, accelerator: MapperIndex | None = None) -> int | None:
        """Map a 0-based CDS offset to its genomic position.

        Args:
            cds_offset (int): 0-based offset into the CDS.
            accelerator (MapperIndex | None): Optional O(1)-lookup accelerator for the mapping.

        Returns:
            int | None: The genomic position, or None if the offset is out of range or the transcript
                is not coding.
        """
        mapper = self.cds_mapper
        return None if mapper is None else mapper.offset_to_genomic(cds_offset, accelerator=accelerator)

    def start_codon_genomic_spans(self) -> list[tuple[int, int]] | None:
        """Genomic spans of the three start-codon bases.

        A list rather than a single ``(min, max)`` envelope because a start codon whose three bases
        split across an exon-exon junction can have more than one genomic span here.

        Returns:
            list[tuple[int, int]] | None: Ascending genomic spans, or None if the transcript is
                not coding.
        """
        mapper = self.cds_mapper
        return None if mapper is None else mapper.codon_genomic_spans(0)

    def stop_codon_genomic_spans(self) -> list[tuple[int, int]] | None:
        """Genomic spans of the three stop-codon bases.

        A list rather than a single ``(min, max)`` envelope.

        Returns:
            list[tuple[int, int]] | None: Ascending genomic spans, or None if not coding or the
                CDS is shorter than 3 bases.
        """
        mapper = self.cds_mapper
        if mapper is None or mapper.length < 3:
            return None
        return mapper.codon_genomic_spans(mapper.length - 1)

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
