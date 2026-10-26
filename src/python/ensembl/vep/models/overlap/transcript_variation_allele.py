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

from ensembl.vep.common import Readthrough, SequenceEdit, Strand, invalidates, reverse_complement
from ensembl.vep.common.constants import EMPTY_ALLELE
from ensembl.vep.common.codons import translate_dna
from ensembl.vep.models.overlap.feature_overlap import FeatureOverlap

if TYPE_CHECKING:
    from ensembl.vep.models.genic import Transcript
    from ensembl.vep.models.sequence import Mapper
    from ensembl.vep.models.variation import VariationAllele
    from ensembl.vep.common import Locus

MAX_READTHROUGH_BASES = 1000
"""Cap on how far ``Readthrough.BEYOND_TRANSCRIPT`` reads past a transcript's own boundary."""

Position = tuple[int | None, int | None]
"""1-based positions. Primarily used for a variant's 5' and 3' ends on a transcript sequence, None for an
end that lies off it"""


def _nucleotide_sequence_up_to_stop(nucleotide_sequence: str) -> tuple[str, str]:
    """Truncate ``nucleotide_sequence`` to end at (and include) its first stop codon, if any.

    Args:
        nucleotide_sequence (str): Nucleotide sequence, transcript-sense.

    Returns:
        tuple[str, str]: The input truncated at (and including) its first stop codon, or
            unchanged if no stop codon is found within it, paired with its translation
            (truncated the same way). ``translate_dna`` translates strictly in fixed 3-base
            chunks from position 0 with no state carried between codons, so slicing its result
            at the stop residue is exactly ``translate_dna`` of the truncated nucleotides --
            re-translating the truncated sequence would repeat work for the same answer.
    """
    peptide = translate_dna(nucleotide_sequence)
    stop_index = peptide.find("*")
    if stop_index == -1:
        return nucleotide_sequence, peptide
    return nucleotide_sequence[: (stop_index + 1) * 3], peptide[: stop_index + 1]


def _edited_boundary(original_boundary: int, edit_start: int, edit_end: int, length_change: int) -> int:
    """Returns the offset of a boundary after an edit.

    The new offset can be used to clip the edited sequence to only have the coding sequence. If the edit
        1. starts after the boundary, the boundary is unaffected.
        2. ends before the boundary, the boundary is shifted by the length change of the edit.
        3. overlaps the boundary, the boundary is anchored to the start of the edit.

    Args:
        original_boundary (int): The boundary's own offset before the edit (0-based, into spliced_seq).
        edit_start (int): Start of the edited half-open range (0-based, into spliced_seq).
        edit_end (int): End of the edited half-open range (exclusive).
        length_change (int): Net length change of the edit (negative for deletion).

    Returns:
        int: The boundary's offset in the edited sequence.
    """
    if original_boundary < edit_start:
        return original_boundary
    if original_boundary >= edit_end:
        return original_boundary + length_change
    return edit_start


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
        cdna_position (Position | None):
            The variant's position on the spliced transcript, or None if it touches no exon.
        cds_position (Position | None):
            The variant's position on the CDS, counted from a whole first codon, or None if
            it does not touch the CDS.
        protein_position (Position | None):
            The residues the variant's two ends fall in, or None if it does not touch the CDS.
        distance_to_transcript (int | None):
            The distance from the variation allele to the transcript,
            or None if the variation allele overlaps the transcript.
        exon_number (int | None):
            The exon number of the transcript that is overlapped by the variation allele,
            or None if the variation allele does not overlap any exon.
        intron_number (int | None):
            The intron number of the transcript that is overlapped by the variation allele,
            or None if the variation allele does not overlap any intron.
        reference_codon (str | None):
            The reference allele's whole-codon window in the CDS (transcript-sense), covering
            every codon touched by the variant, or None if the variant does not overlap the CDS.
            Cased for display as VEP writes it: lower case but for the allele's own bases.
        alternate_codon (str | None):
            The alternate allele's whole-codon window in the CDS (transcript-sense), i.e.
            reference_codon with the variant allele spliced in, or None if not applicable. Cased
            the same way, so a deletion's is all lower case.
        reference_peptide (str | None):
            The amino acid translation of reference_codon, or None if not applicable.
        alternate_peptide (str | None):
            The amino acid translation of alternate_codon, or None if not applicable.
        reference_cds (str | None):
            The whole reference CDS (transcript-sense, includes the start and stop codon), or
            None if the transcript has no CDS at all.
        alternate_cds (str | None):
            The whole reference CDS (transcript-sense) with the variant allele spliced in, or None
            if not applicable.
        reference_cds_peptide (str | None):
            The amino acid translation of reference_cds, or None if reference_cds is None.
        alternate_cds_peptide (str | None):
            The amino acid translation of alternate_cds, or None if alternate_cds is None.
        readthrough (Readthrough):
            How far translation of alternate_cds may run when the variant destroys the stop codon.
    """

    feature: Transcript
    overlapping_feature: VariationAllele
    is_reference: bool = field(default=False)
    readthrough: Readthrough = field(
        default=Readthrough.NONE, on_setattr=invalidates("_full_cds_alleles", "_full_peptide_alleles")
    )
    distance_to_transcript: int | None = field(default=None, init=False)
    exon_number: int | None = field(default=None, init=False)
    intron_number: int | None = field(default=None, init=False)
    _lock: RLock = field(factory=RLock, init=False, repr=False)

    def _position_on(self, mapper: Mapper | None) -> Position | None:
        """Place the variant's 5' and 3' ends on ``mapper``'s sequence, 1-based.

        An insertion's ends are its two flanking bases. As in VEP, where the Ensembl mapper answers
        a gap for an end off the sequence, that end is None; a variant covering none of ``mapper``'s
        spans has no position at all.

        Args:
            mapper (Mapper | None): The mapper whose sequence to place the variant on.

        Returns:
            Position | None: The two ends, or None if the variant misses every span.
        """
        if mapper is None:
            return None
        variant = self.overlapping_feature
        variant_min, variant_max = min(variant.start, variant.end), max(variant.start, variant.end)
        if mapper.genomic_span_to_offsets(variant_min, variant_max) is None:
            return None

        if self.feature.strand == Strand.REVERSE:
            five_prime_end, three_prime_end = variant_max, variant_min
        else:
            five_prime_end, three_prime_end = variant_min, variant_max
        five_prime_offset = mapper.genomic_to_offset(five_prime_end)
        three_prime_offset = mapper.genomic_to_offset(three_prime_end)
        return (
            None if five_prime_offset is None else five_prime_offset + 1,
            None if three_prime_offset is None else three_prime_offset + 1,
        )

    @property
    def cdna_position(self) -> Position | None:
        """Public accessor for _cdna_position.

        Returns:
            Position | None: The variant's position on the spliced transcript, or None if
                it touches no exon.
        """
        return self._cdna_position

    @cached_property
    def _cdna_position(self) -> Position | None:
        """Place the variant on the spliced transcript.

        Returns:
            Position | None: The position, or None if the variant touches no exon.
        """
        return self._position_on(self.feature.cdna_mapper)

    @property
    def cds_position(self) -> Position | None:
        """Public accessor for _cds_position.

        Returns:
            Position | None: The variant's position on the CDS, or None if it does not
                touch the CDS.
        """
        return self._cds_position

    @cached_property
    def _cds_position(self) -> Position | None:
        """Place the variant on the CDS.

        CDS offsets already start after any ``cds_start_phase`` padding, so a CDS with an
        incomplete 5' end is counted from a whole first codon, as VEP adds the start exon's phase.

        Returns:
            Position | None: The position, or None if the variant does not touch the CDS.
        """
        return self._position_on(self.feature.cds_mapper)

    @property
    def protein_position(self) -> Position | None:
        """Public accessor for _protein_position.

        Returns:
            Position | None: The residues the variant's two ends fall in, or None if it
                does not touch the CDS.
        """
        return self._protein_position

    @cached_property
    def _protein_position(self) -> Position | None:
        """Derive the residue positions from the CDS positions, three bases a residue.

        Returns:
            Position | None: The residue positions, or None if the variant does not touch
                the CDS.
        """
        cds_position = self.cds_position
        if cds_position is None:
            return None
        residues: list[int | None] = []
        for position in cds_position:
            residues.append(None if position is None else (position - 1) // 3 + 1)
        return residues[0], residues[1]

    def _edit_on(self, mapper: Mapper | None) -> SequenceEdit | None:
        """Express the edit caused by this variant in ``mapper``'s offset space.

        The geometry shared by every downstream reconstruction, parameterised by which mapper does
        the genomic->offset conversion. Uses pure coordinate and allele-string arithmetic with no
        sequence involved (which makes ``cds_length_change`` cheap).

        Args:
            mapper (Mapper | None): The mapper whose offset space the edit is expressed in.

        Returns:
            SequenceEdit | None: The edit, or None when the variant covers none of ``mapper``'s
                spans (so there is nothing to replace).
        """
        with self._lock:
            if mapper is None:
                return None

            variant = self.overlapping_feature
            variant_min, variant_max = min(variant.start, variant.end), max(variant.start, variant.end)
            footprint = mapper.genomic_span_to_offsets(variant_min, variant_max)
            if footprint is None:
                return None
            low_offset, high_offset = footprint

            alternate_allele = (variant.alternate_allele or "").replace(EMPTY_ALLELE, "")
            reference_allele_length = len((variant.reference_allele or "").replace(EMPTY_ALLELE, ""))
            transcript_sense_alternate = (
                alternate_allele
                if self.feature.strand == Strand.FORWARD
                else reverse_complement(alternate_allele)
            )

            if reference_allele_length == 0:
                # Pure insertion: splice in immediately after the anchor offset.
                edit_start = edit_end = low_offset + 1
            else:
                edit_start, edit_end = low_offset, high_offset + 1

            return SequenceEdit(
                footprint=footprint,
                start=edit_start,
                end=edit_end,
                replacement=transcript_sense_alternate,
            )

    @cached_property
    def _cds_edit(self) -> SequenceEdit | None:
        """This variant as an edit in CDS offset space (clipped to the coding sequence).

        Returns:
            SequenceEdit | None: The edit, or None if the variant does not overlap the CDS.
        """
        return self._edit_on(self.feature.cds_mapper)

    @cached_property
    def _cdna_edit(self) -> SequenceEdit | None:
        """This variant as an edit in cDNA offset space (clipped to the exons).

        Returns:
            SequenceEdit | None: The edit, or None if the variant touches no exon at all.
        """
        return self._edit_on(self.feature.cdna_mapper)

    @property
    def reference_codon(self) -> str | None:
        """Returns the codons variant touches in the reference sequence.

        Returns:
            str | None: The reference allele's whole-codon window (transcript-sense), or None if
                the variant does not overlap the CDS.
        """
        codons = self._codon_window_alleles
        return None if codons is None else codons[0]

    @property
    def alternate_codon(self) -> str | None:
        """Returns the codons variant touches with the variant allele spliced in.

        Returns:
            str | None: reference_codon with the variant allele spliced in, or None if not
                applicable.
        """
        codons = self._codon_window_alleles
        return None if codons is None else codons[1]

    @cached_property
    def _codon_window_alleles(self) -> tuple[str, str] | None:
        """Returns (reference_codon, alternate_codon), the whole-codon window around the variant,
        in transcript-sense, with the variant allele spliced in for the alternate half.

        An insertion that falls between two codons has an empty reference window and the inserted
        bases as its alternate one, as VEP's does.

        Each is lower case but for its own allele's bases, which are upper case: VEP's
        ``display_codon``. An empty allele (the reference side of an insertion, the alternate side
        of a deletion) leaves its codon all lower case. Anything computing on a codon upper-cases it
        first.
        Returns:
            tuple[str, str] | None: (reference_codon, alternate_codon), or None if the variant
                does not overlap the CDS.
        """
        with self._lock:
            transcript = self.feature
            mapper = transcript.cds_mapper
            edit = self._cds_edit
            if mapper is None or edit is None:
                return None

            window = mapper.codon_window(edit.footprint)
            if window is None:
                return None

            coding_sequence = transcript.coding_sequence
            if coding_sequence is None:
                return None

            if edit.start == edit.end and edit.start % 3 == 0:
                # An insertion between two codons: VEP's window runs from the residue after the
                # insertion point back to the one before it, which is empty, leaving the inserted
                # bases alone on the alternate side (``TranscriptVariationAllele::codon``).
                return "", edit.replacement.upper()

            codon_start, codon_end = window
            leading_bases = coding_sequence[codon_start : edit.start].lower()
            trailing_bases = coding_sequence[edit.end : codon_end].lower()
            reference_bases = coding_sequence[edit.start : edit.end].upper()
            reference_codon = leading_bases + reference_bases + trailing_bases
            alternate_codon = leading_bases + edit.replacement.upper() + trailing_bases
            return reference_codon, alternate_codon

    @property
    def reference_peptide(self) -> str | None:
        """Returns the amino acid translation of the reference_codon.

        Returns:
            str | None: The amino acid translation of reference_codon, or None if not applicable.
        """
        peptides = self._peptide_window_alleles
        return None if peptides is None else peptides[0]

    @property
    def alternate_peptide(self) -> str | None:
        """Returns the amino acid translation of the alternate_codon.

        Returns:
            str | None: The amino acid translation of alternate_codon, or None if not applicable.
        """
        peptides = self._peptide_window_alleles
        return None if peptides is None else peptides[1]

    @cached_property
    def _peptide_window_alleles(self) -> tuple[str, str] | None:
        """Returns (reference_peptide, alternate_peptide), translated from _codon_window_alleles.

        Returns:
            tuple[str, str] | None: (reference_peptide, alternate_peptide), or None if the
                variant does not overlap the CDS.
        """
        codons = self._codon_window_alleles
        if codons is None:
            return None
        reference_codon, alternate_codon = codons
        return translate_dna(reference_codon.upper()), translate_dna(alternate_codon.upper())

    @property
    def cds_length_change(self) -> int | None:
        """Net change in **coding** length the variant makes, in bases.

        Positive for an insertion, negative for a deletion, zero for a substitution. Derived from the CDS
        edit geometry rather than the cDNA. Any intronic portion of the variant's genomic span is already
        discounted, since the offsets come from the mapper's block-overlap clipping.

        Returns:
            int | None: The net coding length change, or None when the variant does not overlap
                the CDS and so changes no coding length at all.
        """
        edit = self._cds_edit
        return None if edit is None else edit.length_change

    @property
    def reference_cds(self) -> str | None:
        """Returns the whole reference CDS.

        Returns:
            str | None: The whole reference CDS (transcript-sense), or None if the transcript has
                no CDS at all.
        """
        alleles = self._full_cds_alleles
        return None if alleles is None else alleles[0]

    @property
    def alternate_cds(self) -> str | None:
        """Returns the whole CDS with the variant applied.

        Returns:
            str | None: The whole CDS with the variant applied, or None if the transcript has no CDS
                at all.
        """
        alleles = self._full_cds_alleles
        return None if alleles is None else alleles[1]

    @cached_property
    def _full_cds_alleles(self) -> tuple[str, str, str | None] | None:
        """Returns (reference_cds, alternate_cds, alternate_peptide), the whole CDS, in
        transcript-sense, with the variant applied for the alternate half.

        For alternate CDS edit, replaces the half-open transcript-wide cDNA offset that the variant overlapped
        with the strand-corrected alternate allele in the transcript.spliced_seq (5' UTR + CDS + 3' UTR). By
        default, we assume there is no intron retention, exon skipping or readthrough. Only the
        coding sequence
        is considered.

        Rather than just the CDS it uses cDNA sequence to optionally allow readthrough into the UTR. If not
        opted-for, the CDS is truncated to include only the coding sequence.

        ``reference_cds`` reuses ``transcript.coding_sequence`` in place of re-slicing ``spliced_seq``
        when its length matches -- both are assembled from the same coding exons, so the content
        agrees (proved corpus-wide by tests, not re-checked per overlap). Both sides carry the
        transcript's ``cds_start_phase`` ``N`` padding, as its CDS offsets do.
        ``alternate_peptide`` is the translation ``_nucleotide_sequence_up_to_stop`` already produced
        while truncating ``alternate_cds`` at its first stop codon, so ``_full_peptide_alleles`` can
        reuse it instead of re-translating. It is None only when the variant touches no exon at all,
        the one path that leaves ``alternate_cds`` untranslated here.

        Returns:
            tuple[str, str, str | None] | None: (reference_cds, alternate_cds, alternate_peptide),
                or None if the transcript has no CDS at all (non-coding, or missing sequence
                data). A variant that changes nothing about the CDS (entirely within a UTR,
                entirely intronic, or entirely outside the transcript) yields two identical
                sequences and no alternate peptide.
        """
        with self._lock:
            transcript = self.feature
            cdna_mapper = transcript.cdna_mapper
            spliced_seq = transcript.spliced_seq
            cds_length = transcript.cds_length
            cdna_coding_start = transcript.cdna_coding_start
            if cdna_mapper is None or spliced_seq is None or not cds_length or cdna_coding_start is None:
                return None

            # 0-based start/end position for slicing. cds_length counts only real bases, so the
            # padding is prefixed separately; the spliced sequence does not carry it either.
            padding = "N" * transcript.cds_start_phase
            cdna_coding_start0 = cdna_coding_start - 1
            cdna_coding_end0 = cdna_coding_start0 + cds_length
            coding_sequence = transcript.coding_sequence
            if coding_sequence is not None and len(coding_sequence) == cds_length + len(padding):
                reference_cds = coding_sequence
            else:
                reference_cds = padding + spliced_seq[cdna_coding_start0:cdna_coding_end0]

            edit = self._cdna_edit
            if edit is None:
                # The variant touches no exon at all (purely intronic, or entirely outside the
                # transcript). Splicing is assumed to proceed normally, so the CDS is unchanged.
                return reference_cds, reference_cds, None
            edit_start, edit_end, transcript_sense_alternate = edit.start, edit.end, edit.replacement

            alternate_spliced_seq = (
                spliced_seq[:edit_start] + transcript_sense_alternate + spliced_seq[edit_end:]
            )

            length_change = edit.length_change

            # Clip the edited sequence boundary to contain only the coding sequence.
            anchor = _edited_boundary(cdna_coding_start0, edit_start, edit_end, length_change)
            if self.readthrough >= Readthrough.UTR:
                # Translate on through whatever follows the coding sequence in the edited cDNA.
                end_bound = len(alternate_spliced_seq)
            else:
                end_bound = _edited_boundary(cdna_coding_end0, edit_start, edit_end, length_change)

            alternate_tail = padding + alternate_spliced_seq[anchor:end_bound]
            alternate_cds, alternate_peptide = _nucleotide_sequence_up_to_stop(alternate_tail)

            if self.readthrough >= Readthrough.BEYOND_TRANSCRIPT and "*" not in alternate_peptide:
                # No in-frame stop anywhere in the edited cDNA, so carry on into raw genomic
                # sequence past the transcript's own boundary.
                flank = self._downstream_flanking_sequence(MAX_READTHROUGH_BASES)
                if flank:
                    alternate_cds, alternate_peptide = _nucleotide_sequence_up_to_stop(alternate_tail + flank)

            return reference_cds, alternate_cds, alternate_peptide

    def _downstream_flanking_sequence(self, length: int) -> str:
        """Up to ``length`` bases of raw genomic sequence past this transcript's own 3' boundary
        (transcript sense), read directly from the attached Slice.

        Only reached at ``Readthrough.BEYOND_TRANSCRIPT``. It trusts the attached slice already covers
        what is needed rather than issuing a fresh sequence-source fetch.

        Args:
            length (int): How many bases to try to read.

        Returns:
            str: Up to ``length`` bases of transcript-sense sequence past the transcript's own
                boundary, or ``""`` if none is available.
        """
        transcript = self.feature
        slice_ = transcript.slice
        if length <= 0 or slice_ is None or slice_.seq is None:
            return ""
        if transcript.strand == Strand.FORWARD:
            search_end = min(slice_.end, transcript.end + length)
            if search_end <= transcript.end:
                return ""
            return slice_.seq[transcript.end : search_end]
        search_start = max(slice_.start, transcript.start - length)
        if search_start >= transcript.start:
            return ""
        return reverse_complement(slice_.seq[search_start - 1 : transcript.start - 1])

    @property
    def reference_cds_peptide(self) -> str | None:
        """Returns the amino acid translation of the reference CDS.

        Returns:
            str | None: The amino acid translation of reference_cds, or None if not applicable.
        """
        peptides = self._full_peptide_alleles
        return None if peptides is None else peptides[0]

    @property
    def alternate_cds_peptide(self) -> str | None:
        """Returns the amino acid translation of the alternate CDS.

        Returns:
            str | None: The amino acid translation of alternate_cds, or None if not applicable.
        """
        peptides = self._full_peptide_alleles
        return None if peptides is None else peptides[1]

    @cached_property
    def _full_peptide_alleles(self) -> tuple[str, str] | None:
        """Returns (reference_cds_peptide, alternate_cds_peptide), translated from
        _full_cds_alleles.

        ``alternate_cds_peptide`` reuses the peptide ``_full_cds_alleles`` already produced while
        truncating ``alternate_cds`` at its first stop codon, translating it here only for the one
        path that leaves it untranslated (the variant touches no exon at all). See
        ``_reference_cds_peptide`` for how the reference half is sourced.

        Returns:
            tuple[str, str] | None: (reference_cds_peptide, alternate_cds_peptide), or None if
                the variant does not overlap the CDS.
        """
        alleles = self._full_cds_alleles
        if alleles is None:
            return None
        reference_cds, alternate_cds, alternate_peptide = alleles
        if alternate_peptide is None:
            alternate_peptide = translate_dna(alternate_cds)
        return self._reference_cds_peptide(reference_cds), alternate_peptide

    def _reference_cds_peptide(self, reference_cds: str) -> str:
        """Returns the translation of the reference CDS, reusing ``transcript.translated_seq``
        when it is provably ``translate_dna(reference_cds)``.

        ``translated_seq`` returns ``translation.seq`` when set. For every source in this codebase
        today that is itself ``translate_dna(cds.seq)`` (GFF/GTF-built translations seed
        ``seq=None`` and fall back the same way), and ``Transcript._translation`` marks exactly
        that case via ``Translation.derived_from_coding_sequence``. A hypothetical future source
        could instead seed a curated peptide (for example, one edited for selenocysteine or without
        a trailing stop), which the stop facts compare positionally against this CDS -- so reuse is
        guarded to the cases where that cannot have happened: the transcript's coding sequence
        matches ``reference_cds`` in length, and its translation is either unset, has no seq of its
        own, or is marked as derived from this same coding sequence.

        Checking the marker rather than *when* ``_translation`` resolved (the previous guard) means
        reuse no longer depends on evaluation order: a transcript with several overlapping variants
        benefits from the very first one resolving ``_translation`` the plain way, instead of only
        the first overlap getting to reuse it and every later one paying to retranslate.

        Args:
            reference_cds (str): The reference CDS this TVA computed (possibly itself already
                reused from ``transcript.cds.seq``).

        Returns:
            str: The reference CDS peptide.
        """
        transcript = self.feature
        coding_sequence = transcript.coding_sequence
        if coding_sequence is None or len(coding_sequence) != len(reference_cds):
            return translate_dna(reference_cds)

        translation = transcript.translation
        if translation is None or translation.seq is None or translation.derived_from_coding_sequence:
            translated_seq = transcript.translated_seq
            if translated_seq is not None:
                return translated_seq

        return translate_dna(reference_cds)

    @property
    def shifted(self) -> bool:
        """
        Public accessor for _shifted.

        Returns:
            bool: True if the allele has been shifted, False otherwise.
        """
        return self._shifted

    @cached_property
    def _shifted(self) -> bool:
        """
        Returns whether the allele has been shifted from its original position.

        Returns:
            bool: True if the allele has been shifted, False otherwise.
        """
        return False

    @property
    def shifted_cdna_locus(self) -> Locus | None:
        """
        Public accessor for _shifted_cdna_locus.

        Returns:
            Locus | None: The cDNA locus of the variation allele after shifting,
                or None if not applicable.
        """
        return self._shifted_cdna_locus

    @cached_property
    def _shifted_cdna_locus(self) -> Locus | None:
        """
        Returns the cDNA locus of the variation allele after shifting, or None if not applicable.

        Returns:
            Locus | None: The cDNA locus of the variation allele after shifting,
                or None if not applicable.
        """
        return None

    @property
    def shifted_cds_locus(self) -> Locus | None:
        """
        Public accessor for _shifted_cds_locus.

        Returns:
            Locus | None: The CDS locus of the variation allele after shifting,
                or None if not applicable.
        """
        return self._shifted_cds_locus

    @cached_property
    def _shifted_cds_locus(self) -> Locus | None:
        """
        Returns the CDS locus of the variation allele after shifting, or None if not applicable.

        Returns:
            Locus | None: The CDS locus of the variation allele after shifting,
                or None if not applicable.
        """
        return None

    @property
    def shifted_protein_locus(self) -> Locus | None:
        """
        Public accessor for _shifted_protein_locus.

        Returns:
            Locus | None: The protein locus of the variation allele after shifting,
                or None if not applicable.
        """
        return self._shifted_protein_locus

    @cached_property
    def _shifted_protein_locus(self) -> Locus | None:
        """
        Returns the protein locus of the variation allele after shifting, or None if not applicable.

        Returns:
            Locus | None: The protein locus of the variation allele after shifting,
                or None if not applicable.
        """
        return None

    @property
    def original_cdna_locus(self) -> Locus | None:
        """Public accessor for _original_cdna_locus.

        Returns:
            Locus | None: The original cDNA locus of the variation allele before shifting,
                or None if not applicable.
        """
        return self._original_cdna_locus

    @cached_property
    def _original_cdna_locus(self) -> Locus | None:
        """Returns the original cDNA locus of the variation allele before shifting, or None if not applicable.

        Returns:
            Locus | None: The original cDNA locus of the variation allele before shifting,
                or None if not applicable.
        """
        return None

    @property
    def original_cds_locus(self) -> Locus | None:
        """Public accessor for _original_cds_locus.

        Returns:
            Locus | None: The original CDS locus of the variation allele before shifting,
                or None if not applicable.
        """
        return self._original_cds_locus

    @cached_property
    def _original_cds_locus(self) -> Locus | None:
        """Returns the original CDS locus of the variation allele before shifting, or None if not applicable.

        Returns:
            Locus | None: The original CDS locus of the variation allele before shifting,
                or None if not applicable.
        """
        return None

    @property
    def original_protein_locus(self) -> Locus | None:
        """Public accessor for _original_protein_locus.

        Returns:
            Locus | None: The original protein locus of the variation allele before shifting,
                or None if not applicable.
        """
        return self._original_protein_locus

    @cached_property
    def _original_protein_locus(self) -> Locus | None:
        """The original protein locus of the variation allele before shifting,
        or None if not applicable.

        Returns:
            Locus | None: The original protein locus of the variation allele before shifting,
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
