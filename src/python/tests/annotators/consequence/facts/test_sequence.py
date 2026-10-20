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
"""Tests for consequence sequence facts."""

from __future__ import annotations

from typing import cast

import pytest

from ensembl.vep.annotators.consequence.facts import sequence
from ensembl.vep.annotators.consequence.facts.base import REGISTRY, get_fact
from ensembl.vep.annotators.consequence.facts.helper import get_full_cds_peptide_alleles
from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.common import CostTier, Readthrough, Strand
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.models.variation import VariationAllele


def _genomic_base(transcript: Transcript, position: int) -> str:
    """Return the raw (forward-strand) base at a genomic (slice-relative) position."""
    assert transcript.slice is not None
    assert transcript.slice.seq is not None
    return transcript.slice.seq[position - 1]


def _create_transcript_variation_allele(
    transcript: Transcript,
    *,
    start: int,
    end: int,
    reference: str,
    alternate: str,
    readthrough: Readthrough = Readthrough.NONE,
) -> TranscriptVariationAllele:
    """Create a transcript variation allele at explicit genomic coordinates."""
    variation_allele = VariationAllele(
        start=start,
        end=end,
        strand=Strand.FORWARD,
        slice=transcript.slice,
        reference_allele=reference,
        alternate_allele=alternate,
    )
    return TranscriptVariationAllele(
        feature=transcript,
        overlapping_feature=variation_allele,
        readthrough=readthrough,
    )


def _substitution_context(transcript: Transcript, position: int) -> FactContext:
    """A single-base substitution context at one genomic position."""
    reference_base = _genomic_base(transcript, position)
    alternate_base = "C" if reference_base != "C" else "G"
    overlap = _create_transcript_variation_allele(
        transcript, start=position, end=position, reference=reference_base, alternate=alternate_base
    )
    return FactContext(overlap=overlap)


def _seeded_context(
    transcript: Transcript, position: int, fact_name: str, peptides: tuple[str, str]
) -> FactContext:
    """A substitution context at ``position`` whose peptide-pair fact is seeded with ``peptides``.

    The guard atoms are left to compute against the real transcript, so the position must be one
    that satisfies the fact's own ``requires``; only the peptide pair is supplied, because
    reaching an arbitrary pair from genomic coordinates alone would need a bespoke transcript per
    case. Seeding is the same mechanism the shared-codon cache uses.

    Args:
        transcript (Transcript): The transcript to overlap.
        position (int): Genomic position of the substitution.
        fact_name (str): ``_peptide_alleles`` or ``_full_cds_peptide_alleles``.
        peptides (tuple[str, str]): The (reference, alternate) peptide pair to seed.

    Returns:
        FactContext: The seeded context.
    """
    context = _substitution_context(transcript, position)
    context.prefill({fact_name: peptides})
    return context


def _mid_cds_position(transcript: Transcript) -> int:
    """A genomic position well inside the CDS, clear of both the start and the stop codon."""
    position = transcript.cds_offset_to_genomic(300)
    assert position is not None
    return position


def _start_codon_position(transcript: Transcript) -> int:
    """The lowest genomic position of the start codon."""
    spans = transcript.start_codon_genomic_spans()
    assert spans is not None
    return spans[0][0]


def _stop_codon_position(transcript: Transcript) -> int:
    """The lowest genomic position of the stop codon."""
    spans = transcript.stop_codon_genomic_spans()
    assert spans is not None
    return spans[0][0]


class TestCodonAndPeptideWindowFacts:
    """Facts comparing only the reference and alternate peptide window around the variant."""

    @pytest.mark.parametrize(
        ("peptides", "expected"),
        [
            pytest.param(("M", "M"), False, id="identical"),
            pytest.param(("M", "L"), True, id="one_residue_differs"),
            pytest.param(("MK", "MK"), False, id="identical_pair"),
            pytest.param(("MK", "MKA"), True, id="alternate_is_longer"),
        ],
    )
    def test_peptide_changed(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        peptides: tuple[str, str],
        expected: bool,
    ) -> None:
        """The fact is True whenever the two peptide windows are not the same string."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        position = _mid_cds_position(transcript)
        context = _seeded_context(transcript, position, "_peptide_alleles", peptides)
        assert get_fact(context, "peptide_changed") is expected

    @pytest.mark.parametrize(
        ("peptides", "expected"),
        [
            pytest.param(("M", "L"), True, id="both_one_residue"),
            pytest.param(("MKA", "LRG"), True, id="both_three_residues"),
            pytest.param(("M", "ML"), False, id="alternate_is_longer"),
            pytest.param(("MKA", "MK"), False, id="alternate_is_shorter"),
        ],
    )
    def test_peptide_same_length(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        peptides: tuple[str, str],
        expected: bool,
    ) -> None:
        """The fact compares lengths only, not content."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        position = _mid_cds_position(transcript)
        context = _seeded_context(transcript, position, "_peptide_alleles", peptides)
        assert get_fact(context, "peptide_same_length") is expected

    @pytest.mark.parametrize(
        ("peptides", "expected"),
        [
            pytest.param(("M", "*"), True, id="alternate_is_a_stop"),
            pytest.param(("M", "L*"), True, id="alternate_ends_with_a_stop"),
            pytest.param(("M", "*L"), True, id="alternate_starts_with_a_stop"),
            pytest.param(("M", "L"), False, id="alternate_has_no_stop"),
            pytest.param(("*", "M"), False, id="only_the_reference_has_a_stop"),
        ],
    )
    def test_alt_peptide_has_ter(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        peptides: tuple[str, str],
        expected: bool,
    ) -> None:
        """The fact reads the alternate window only, wherever the stop sits in it."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        position = _mid_cds_position(transcript)
        context = _seeded_context(transcript, position, "_peptide_alleles", peptides)
        assert get_fact(context, "alt_peptide_has_ter") is expected

    @pytest.mark.parametrize(
        ("peptides", "expected"),
        [
            pytest.param(("*", "M"), True, id="reference_is_a_stop"),
            pytest.param(("L*", "M"), True, id="reference_ends_with_a_stop"),
            pytest.param(("*L", "M"), True, id="reference_starts_with_a_stop"),
            pytest.param(("L", "M"), False, id="reference_has_no_stop"),
            pytest.param(("M", "*"), False, id="only_the_alternate_has_a_stop"),
        ],
    )
    def test_ref_peptide_has_ter(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        peptides: tuple[str, str],
        expected: bool,
    ) -> None:
        """The fact reads the reference window only, wherever the stop sits in it."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        position = _mid_cds_position(transcript)
        context = _seeded_context(transcript, position, "_peptide_alleles", peptides)
        assert get_fact(context, "ref_peptide_has_ter") is expected

    @pytest.mark.parametrize(
        ("peptides", "expected"),
        [
            pytest.param(("X", "M"), True, id="reference_is_ambiguous"),
            pytest.param(("M", "X"), True, id="alternate_is_ambiguous"),
            pytest.param(("MX", "ML"), True, id="ambiguous_residue_mid_window"),
            pytest.param(("M", "L"), False, id="neither_is_ambiguous"),
        ],
    )
    def test_peptide_has_x(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        peptides: tuple[str, str],
        expected: bool,
    ) -> None:
        """An ambiguous residue in either window is enough."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        position = _mid_cds_position(transcript)
        context = _seeded_context(transcript, position, "_peptide_alleles", peptides)
        assert get_fact(context, "peptide_has_X") is expected

    @pytest.mark.parametrize(
        ("peptides", "expected"),
        [
            pytest.param(("*", "M"), True, id="reference_is_exactly_a_stop"),
            pytest.param(("*L", "M"), False, id="reference_starts_with_a_stop"),
            pytest.param(("L*", "M"), False, id="reference_ends_with_a_stop"),
            pytest.param(("", "M"), False, id="reference_is_empty"),
            pytest.param(("M", "*"), False, id="only_the_alternate_is_a_stop"),
        ],
    )
    def test_ref_peptide_is_ter(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        peptides: tuple[str, str],
        expected: bool,
    ) -> None:
        """ "Is" means the whole window, not merely containing a stop."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        position = _mid_cds_position(transcript)
        context = _seeded_context(transcript, position, "_peptide_alleles", peptides)
        assert get_fact(context, "ref_peptide_is_ter") is expected

    @pytest.mark.parametrize(
        ("peptides", "expected"),
        [
            pytest.param(("M", "*"), True, id="alternate_is_exactly_a_stop"),
            pytest.param(("M", "*L"), False, id="alternate_starts_with_a_stop"),
            pytest.param(("M", "L*"), False, id="alternate_ends_with_a_stop"),
            pytest.param(("M", ""), False, id="alternate_is_empty"),
            pytest.param(("*", "M"), False, id="only_the_reference_is_a_stop"),
        ],
    )
    def test_alt_peptide_is_ter(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        peptides: tuple[str, str],
        expected: bool,
    ) -> None:
        """ "Is" means the whole window, not merely containing a stop."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        position = _mid_cds_position(transcript)
        context = _seeded_context(transcript, position, "_peptide_alleles", peptides)
        assert get_fact(context, "alt_peptide_is_ter") is expected

    @pytest.mark.parametrize(
        ("peptides", "expected"),
        [
            pytest.param(("*", "M"), True, id="reference_is_a_stop"),
            pytest.param(("*LK", "M"), True, id="reference_starts_with_a_stop"),
            pytest.param(("L*", "M"), False, id="stop_is_not_first"),
            pytest.param(("", "M"), False, id="reference_is_empty"),
            pytest.param(("M", "*L"), False, id="only_the_alternate_starts_with_a_stop"),
        ],
    )
    def test_ref_peptide_startswith_ter(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        peptides: tuple[str, str],
        expected: bool,
    ) -> None:
        """Only a stop in the reference window's first position counts."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        position = _mid_cds_position(transcript)
        context = _seeded_context(transcript, position, "_peptide_alleles", peptides)
        assert get_fact(context, "ref_peptide_startswith_ter") is expected

    @pytest.mark.parametrize(
        ("peptides", "expected"),
        [
            pytest.param(("M", "AM"), True, id="residue_added_before"),
            pytest.param(("MK", "AGMK"), True, id="residues_added_before"),
            pytest.param(("M", "MA"), False, id="residue_added_after"),
            pytest.param(("M", "L"), False, id="reference_not_kept"),
            pytest.param(("", "A"), False, id="empty_reference_is_excluded"),
        ],
    )
    def test_insertion_before_start(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        peptides: tuple[str, str],
        expected: bool,
    ) -> None:
        """The reference window survives at the end of the alternate one."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        position = _start_codon_position(transcript)
        context = _seeded_context(transcript, position, "_peptide_alleles", peptides)
        assert get_fact(context, "insertion_before_start") is expected

    @pytest.mark.parametrize(
        ("peptides", "expected"),
        [
            pytest.param(("M", "MA"), True, id="residue_added_after"),
            pytest.param(("MK", "MKAG"), True, id="residues_added_after"),
            pytest.param(("M", "AM"), False, id="residue_added_before"),
            pytest.param(("M", "L"), False, id="reference_not_kept"),
            pytest.param(("", "A"), False, id="empty_reference_is_excluded"),
        ],
    )
    def test_insertion_after_start(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        peptides: tuple[str, str],
        expected: bool,
    ) -> None:
        """The reference window survives at the start of the alternate one."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        position = _start_codon_position(transcript)
        context = _seeded_context(transcript, position, "_peptide_alleles", peptides)
        assert get_fact(context, "insertion_after_start") is expected

    def test_an_unchanged_window_reads_as_an_insertion_on_both_sides(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """Neither insertion fact separates "inserted nothing" from an insertion on its own side.

        ``str.startswith`` and ``str.endswith`` are both true of an identical string, so the two
        facts are only opposites once the alternate window is genuinely longer.
        """
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        position = _start_codon_position(transcript)
        context = _seeded_context(transcript, position, "_peptide_alleles", ("M", "M"))

        assert get_fact(context, "insertion_before_start") is True
        assert get_fact(context, "insertion_after_start") is True

    @pytest.mark.parametrize(
        "fact_name",
        [
            "peptide_changed",
            "peptide_same_length",
            "alt_peptide_has_ter",
            "ref_peptide_has_ter",
            "peptide_has_X",
            "ref_peptide_is_ter",
            "alt_peptide_is_ter",
            "ref_peptide_startswith_ter",
            "insertion_before_start",
            "insertion_after_start",
        ],
    )
    def test_window_facts_are_false_for_an_intronic_variant(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript], fact_name: str
    ) -> None:
        """The ``overlap_cds`` guard settles every window fact without translating anything."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        intron = transcript.introns[0]
        context = _substitution_context(transcript, (intron.start + intron.end) // 2)

        assert get_fact(context, "overlap_cds") is False
        assert get_fact(context, fact_name) is False

    def test_peptide_alleles_translate_the_real_codon_window(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """Without seeding for test purpose, the window comes from the transcript's own sequence."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        context = _substitution_context(transcript, _start_codon_position(transcript))

        codons = cast("tuple[str, str] | None", get_fact(context, "_codon_alleles"))
        peptides = cast("tuple[str, str] | None", get_fact(context, "_peptide_alleles"))
        assert codons is not None and peptides is not None
        assert codons[0] == "ATG"
        assert peptides[0] == "M"

    def test_stop_gained_sets_peptide_facts(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """Stop gained sets peptide facts."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        # Try every base of the second and third codons against every alternate base: with 6
        # positions x 3 alternates, at least one is virtually certain to introduce a stop codon
        # (~50% of single-base substitutions to a random sense codon produce a stop or not,
        # per position; trying this many combinations makes flakiness negligible).
        found_stop_gain = False
        for cds_offset in range(3, 9):
            target_position = transcript.cds_offset_to_genomic(cds_offset)
            assert target_position is not None
            reference_base = _genomic_base(transcript, target_position)
            for alternate_base in "ACGT":
                if alternate_base == reference_base:
                    continue
                overlap = _create_transcript_variation_allele(
                    transcript,
                    start=target_position,
                    end=target_position,
                    reference=reference_base,
                    alternate=alternate_base,
                )
                context = FactContext(overlap=overlap)
                if get_fact(context, "alt_peptide_has_ter") and not get_fact(context, "ref_peptide_has_ter"):
                    found_stop_gain = True
                    assert overlap.reference_peptide is not None
                    assert overlap.alternate_peptide is not None
                    assert "*" in overlap.alternate_peptide
                    break
            if found_stop_gain:
                break
        assert found_stop_gain, "expected at least one base substitution to introduce a stop codon"


class TestInframeEdgeFacts:
    """What VEP asks before calling an indel inframe: did it only add or remove whole codons."""

    @pytest.mark.parametrize(
        ("codons", "expected"),
        [
            pytest.param(("ATGCTC", "ATG"), True, id="alternate_is_the_start"),
            pytest.param(("ATGCTC", "CTC"), True, id="alternate_is_the_end"),
            pytest.param(("ATGCTC", ""), True, id="whole_window_deleted"),
            pytest.param(("ATCACG", "ATG"), True, id="whole_codon_from_inside"),
            pytest.param(("ATCACGG", "ATG"), False, id="inside_but_not_whole_codons"),
            pytest.param(("TGGGATAAGTTC", "GAACTT"), False, id="neither_edge_nor_inside"),
        ],
    )
    def test_alt_matches_ref_codon_edge(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        codons: tuple[str, str],
        expected: bool,
    ) -> None:
        """A deletion is inframe when what is left sits at an edge, or whole codons came out."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        position = _mid_cds_position(transcript)
        context = _seeded_context(transcript, position, "_codon_alleles", codons)

        assert get_fact(context, "alt_matches_ref_codon_edge") is expected

    @pytest.mark.parametrize(
        ("peptides", "expected"),
        [
            pytest.param(("I", "IR"), True, id="reference_opens_the_alternate"),
            pytest.param(("I", "RI"), True, id="reference_closes_the_alternate"),
            pytest.param(("", "T"), True, id="no_reference_residues"),
            pytest.param(("T", "NA"), False, id="the_window_residue_changed"),
            pytest.param(("TL", "XL*TL"), False, id="only_what_precedes_a_stop_counts"),
            pytest.param(("T", "TA*"), True, id="an_added_stop_keeps_the_reference"),
        ],
    )
    def test_alt_peptide_matches_ref_edge(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        peptides: tuple[str, str],
        expected: bool,
    ) -> None:
        """An insertion is inframe when the reference residues survive at one end of the window."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        position = _mid_cds_position(transcript)
        context = _seeded_context(transcript, position, "_peptide_alleles", peptides)

        assert get_fact(context, "alt_peptide_matches_ref_edge") is expected


class TestFullCdsPeptideAllelesCostTier:
    """``_full_cds_peptide_alleles`` is tagged with the dedicated, more expensive tier."""

    def test_registered_at_transcript_sequence_tier(self) -> None:
        """The full-CDS rebuild/translate costs more than the narrow codon-window tier, so it
        sorts after ``peptide_has_X`` in a rule's compiled plan (see the columnar strategy's
        ``TestStopGainedPlanOrder``)."""
        registered_fact = REGISTRY["_full_cds_peptide_alleles"]
        assert registered_fact.cost == CostTier.TRANSCRIPT_SEQUENCE
        assert CostTier.TRANSCRIPT_SEQUENCE > CostTier.SEQUENCE


class TestFullCdsPeptideFacts:
    """Facts comparing the whole reference and alternate CDS translations."""

    @pytest.mark.parametrize(
        ("peptides", "expected"),
        [
            pytest.param(("MAKG*", "MAKG*"), False, id="unchanged"),
            pytest.param(("MAKG*", "LAKG*"), True, id="initiator_differs"),
            pytest.param(("MAKG*", "MTKG*"), False, id="second_residue_differs"),
            pytest.param(("MAKG*", "MATG*"), False, id="third_residue_differs"),
            pytest.param(("MAKG*", "MAKT*"), False, id="fourth_residue_differs"),
        ],
    )
    def test_cds_start_lost(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        peptides: tuple[str, str],
        expected: bool,
    ) -> None:
        """The comparison window is the first three residues, not the initiator alone."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        position = _start_codon_position(transcript)
        context = _seeded_context(transcript, position, "_full_cds_peptide_alleles", peptides)
        assert get_fact(context, "cds_start_lost") is expected

    @pytest.mark.parametrize(
        ("peptides", "expected"),
        [
            pytest.param(("MAKG*", "MAKG*"), True, id="unchanged"),
            pytest.param(("MAKG*", "LAKG*"), False, id="initiator_differs"),
            pytest.param(("MAKG*", "MTKG*"), True, id="second_residue_differs"),
            pytest.param(("MAKG*", "MATG*"), True, id="third_residue_differs"),
            pytest.param(("MAKG*", "MAKT*"), True, id="fourth_residue_differs"),
        ],
    )
    def test_cds_start_retained(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        peptides: tuple[str, str],
        expected: bool,
    ) -> None:
        """Retention is computed from the peptides directly, not as ``not cds_start_lost``."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        position = _start_codon_position(transcript)
        context = _seeded_context(transcript, position, "_full_cds_peptide_alleles", peptides)
        assert get_fact(context, "cds_start_retained") is expected

    @pytest.mark.parametrize(
        ("peptides", "expected"),
        [
            pytest.param(("MAK*", "MAKQ"), True, id="stop_replaced_by_a_residue"),
            pytest.param(("MAK*", "MAK"), True, id="alternate_ends_before_the_stop"),
            pytest.param(("MAK*", "MAK*"), False, id="stop_unchanged"),
            pytest.param(("MAK*", "MTK*"), False, id="prefix_differs_so_not_local"),
            pytest.param(("MAKQ", "MAKL"), False, id="reference_has_no_stop"),
        ],
    )
    def test_cds_stop_lost(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        peptides: tuple[str, str],
        expected: bool,
    ) -> None:
        """Tested at the reference's own stop position, with the prefix required to match."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        position = _stop_codon_position(transcript)
        context = _seeded_context(transcript, position, "_full_cds_peptide_alleles", peptides)
        assert get_fact(context, "cds_stop_lost") is expected

    @pytest.mark.parametrize(
        ("peptides", "expected"),
        [
            pytest.param(("MAK*", "MAK*"), True, id="stop_unchanged"),
            pytest.param(("MAK*", "MAK*G"), True, id="stop_kept_with_trailing_residues"),
            pytest.param(("MAK*", "MAKQ"), False, id="stop_replaced_by_a_residue"),
            pytest.param(("MAK*", "MAK"), False, id="alternate_ends_before_the_stop"),
            pytest.param(("MAK*", "MTK*"), False, id="prefix_differs_so_not_local"),
            pytest.param(("MAKQ", "MAKQ"), False, id="reference_has_no_stop"),
        ],
    )
    def test_cds_stop_retained(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        peptides: tuple[str, str],
        expected: bool,
    ) -> None:
        """The stop must be at the very same reading-frame position, behind the same prefix."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        position = _stop_codon_position(transcript)
        context = _seeded_context(transcript, position, "_full_cds_peptide_alleles", peptides)
        assert get_fact(context, "cds_stop_retained") is expected

    @pytest.mark.parametrize(
        ("peptides", "expected"),
        [
            pytest.param(("MAKQ*", "MA*Q*"), True, id="new_stop_before_the_reference_stop"),
            pytest.param(("MAKQ*", "*AKQ*"), True, id="new_stop_at_the_first_residue"),
            pytest.param(("MAKQ*", "MAKQ*"), False, id="unchanged"),
            pytest.param(("MAKQ*", "MAKQL"), False, id="alternate_has_no_stop"),
            pytest.param(("MAKQL", "MA*QL"), False, id="reference_has_no_stop"),
        ],
    )
    def test_cds_stop_gain(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        peptides: tuple[str, str],
        expected: bool,
    ) -> None:
        """A stop gain is a stop strictly earlier than the reference's own."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        position = _mid_cds_position(transcript)
        context = _seeded_context(transcript, position, "_full_cds_peptide_alleles", peptides)
        assert get_fact(context, "cds_stop_gain") is expected

    def test_start_codon_lost_and_retained_is_mutually_exclusive(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """Start codon lost and retained is mutually exclusive."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        start_spans = transcript.start_codon_genomic_spans()
        assert start_spans is not None
        first_base_position = start_spans[0][0]
        reference_base = _genomic_base(transcript, first_base_position)
        assert reference_base == "A"  # first base of ATG

        lost_context = FactContext(
            overlap=_create_transcript_variation_allele(
                transcript, start=first_base_position, end=first_base_position, reference="A", alternate="C"
            )
        )
        assert get_fact(lost_context, "cds_start_lost") is True
        assert get_fact(lost_context, "cds_start_retained") is False

    def test_start_codon_lost_reverse_strand(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """The forward-only test above pins the common case; this pins the same fact on the
        reverse strand, where the transcript-sense first base is the *highest* genomic position
        in the start codon span, not the lowest."""
        transcript = synthetic_transcripts_by_stable_id["ENST0002"]
        assert transcript.strand == Strand.REVERSE
        start_spans = transcript.start_codon_genomic_spans()
        assert start_spans is not None
        first_base_position = start_spans[-1][1]
        reference_base = _genomic_base(transcript, first_base_position)
        assert reference_base == "T"  # complement of 'A', the first base of ATG in transcript sense

        lost_context = FactContext(
            overlap=_create_transcript_variation_allele(
                transcript, start=first_base_position, end=first_base_position, reference="T", alternate="C"
            )
        )
        assert get_fact(lost_context, "cds_start_lost") is True
        assert get_fact(lost_context, "cds_start_retained") is False

    def test_stop_codon_synonymous_change_retained(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """Changing the middle base of TAA to G gives TGA -- still a stop codon, same position."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        stop_spans = transcript.stop_codon_genomic_spans()
        assert stop_spans is not None
        middle_base_position = stop_spans[0][0] + 1
        reference_base = _genomic_base(transcript, middle_base_position)
        assert reference_base == "A"  # middle base of TAA

        context = FactContext(
            overlap=_create_transcript_variation_allele(
                transcript, start=middle_base_position, end=middle_base_position, reference="A", alternate="G"
            )
        )
        assert get_fact(context, "cds_stop_retained") is True
        assert get_fact(context, "cds_stop_lost") is False

    def test_synonymous_mid_cds_snv_does_not_touch_start_or_stop(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """A synonymous SNV mid-CDS (GCG -> GCA, both Ala) touches neither the start nor the
        stop codon."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        genomic_position = transcript.cds_offset_to_genomic(278)
        assert genomic_position is not None
        assert _genomic_base(transcript, genomic_position) == "G"

        overlap = _create_transcript_variation_allele(
            transcript, start=genomic_position, end=genomic_position, reference="G", alternate="A"
        )
        context = FactContext(overlap=overlap)
        for name in ("cds_start_lost", "cds_start_retained", "cds_stop_lost", "cds_stop_retained"):
            assert get_fact(context, name) is False, name
            assert sequence.__dict__[name](FactContext(overlap=overlap)) is False, f"{name} (unguarded)"

    def test_inframe_mid_cds_deletion_does_not_touch_start_or_stop(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """A whole-codon (3bp, in-frame) deletion mid-CDS touches neither the start nor the stop
        codon."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        genomic_start = transcript.cds_offset_to_genomic(306)
        genomic_end = transcript.cds_offset_to_genomic(308)
        assert genomic_start is not None and genomic_end is not None

        overlap = _create_transcript_variation_allele(
            transcript, start=genomic_start, end=genomic_end, reference="TGG", alternate="-"
        )
        context = FactContext(overlap=overlap)
        for name in ("cds_start_lost", "cds_start_retained", "cds_stop_lost", "cds_stop_retained"):
            assert get_fact(context, name) is False, name
            assert sequence.__dict__[name](FactContext(overlap=overlap)) is False, f"{name} (unguarded)"

    def test_frameshift_before_stop_does_not_touch_start_or_stop_codon_directly(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """A 1bp deletion 10 codons upstream of the stop shifts everything downstream of it, but
        doesn't overlap the start or stop codon's own genomic span -- these facts must stay False
        regardless of what the frameshift does to the peptide, since overlap_start_codon/
        overlap_stop_codon (what they're actually gated on) are False here."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        genomic_position = transcript.cds_offset_to_genomic(519)
        assert genomic_position is not None

        overlap = _create_transcript_variation_allele(
            transcript, start=genomic_position, end=genomic_position, reference="G", alternate="-"
        )
        context = FactContext(overlap=overlap)
        assert get_fact(context, "overlap_start_codon") is False
        assert get_fact(context, "overlap_stop_codon") is False
        for name in ("cds_start_lost", "cds_start_retained", "cds_stop_lost", "cds_stop_retained"):
            assert get_fact(context, name) is False, name
            assert sequence.__dict__[name](FactContext(overlap=overlap)) is False, f"{name} (unguarded)"


class TestSequenceLengthFacts:
    """Facts that are arithmetic over the variant's coding-length change."""

    @pytest.mark.parametrize(
        ("reference_length", "alternate_length", "expected"),
        [
            pytest.param(1, 1, True, id="substitution_changes_nothing"),
            pytest.param(3, 3, True, id="three_base_substitution"),
            pytest.param(1, 0, False, id="one_base_deletion"),
            pytest.param(3, 0, True, id="three_base_deletion"),
            pytest.param(0, 1, False, id="one_base_insertion"),
            pytest.param(0, 3, True, id="three_base_insertion"),
            pytest.param(0, 4, False, id="four_base_insertion"),
        ],
    )
    def test_cds_length_change_inframe(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        reference_length: int,
        alternate_length: int,
        expected: bool,
    ) -> None:
        """The net coding-length change is in frame when it is a multiple of three."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        low_position = _mid_cds_position(transcript)
        # An insertion carries no reference base, and is spelled with end == start - 1.
        if reference_length == 0:
            start, end = low_position + 1, low_position
            reference = ""
        else:
            start, end = low_position, low_position + reference_length - 1
            reference = "".join(
                _genomic_base(transcript, start + offset) for offset in range(reference_length)
            )

        overlap = _create_transcript_variation_allele(
            transcript, start=start, end=end, reference=reference, alternate="G" * alternate_length
        )
        context = FactContext(overlap=overlap)
        assert get_fact(context, "overlap_cds") is True
        assert get_fact(context, "cds_length_change_inframe") is expected

    def test_cds_length_change_inframe_is_false_outside_the_cds(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """An intronic variant has no coding-length change to be in frame."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        intron = transcript.introns[0]
        context = _substitution_context(transcript, (intron.start + intron.end) // 2)

        assert get_fact(context, "overlap_cds") is False
        assert get_fact(context, "cds_length_change_inframe") is False


class TestReadthroughFacts:
    """The full-CDS facts read whatever readthrough level the overlap carries.

    Raising the level must not disturb a variant that leaves the stop codon alone, and must not
    turn a genuine stop loss into a retention just because translation later runs into some other
    stop downstream.
    """

    @staticmethod
    def _four_facts(context: FactContext) -> dict[str, bool]:
        return {
            name: cast("bool", get_fact(context, name))
            for name in ("cds_start_lost", "cds_start_retained", "cds_stop_lost", "cds_stop_retained")
        }

    def test_the_default_readthrough_level_is_none(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """An overlap built without a readthrough level does not read through its stop codon."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        position = _stop_codon_position(transcript)
        overlap = _create_transcript_variation_allele(
            transcript, start=position, end=position, reference="T", alternate="C"
        )
        assert overlap.readthrough is Readthrough.NONE

    def test_raising_readthrough_level_never_shortens_the_alternate_peptide(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """Translation of a stop-losing variant runs at least as far at each higher level."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        position = _stop_codon_position(transcript)

        lengths = []
        for level in Readthrough:
            overlap = _create_transcript_variation_allele(
                transcript, start=position, end=position, reference="T", alternate="C", readthrough=level
            )
            peptides = get_full_cds_peptide_alleles(FactContext(overlap=overlap))
            assert peptides is not None
            lengths.append(len(peptides[1]))

        assert lengths == sorted(lengths)

    def test_the_reference_peptide_is_readthrough_level_invariant(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """Readthrough follows the alternate translation; the reference keeps its own stop."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        position = _stop_codon_position(transcript)

        reference_peptides = set()
        for level in Readthrough:
            overlap = _create_transcript_variation_allele(
                transcript, start=position, end=position, reference="T", alternate="C", readthrough=level
            )
            peptides = get_full_cds_peptide_alleles(FactContext(overlap=overlap))
            assert peptides is not None
            reference_peptides.add(peptides[0])

        assert len(reference_peptides) == 1

    def test_localisation_negatives_are_readthrough_level_invariant(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """A synonymous mid-CDS SNV touches neither the start nor the stop codon, so it leaves the
        stop intact -- every level must agree, and all four facts stay False."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        genomic_position = transcript.cds_offset_to_genomic(278)
        assert genomic_position is not None

        results = []
        for level in Readthrough:
            overlap = _create_transcript_variation_allele(
                transcript,
                start=genomic_position,
                end=genomic_position,
                reference="G",
                alternate="A",
                readthrough=level,
            )
            results.append(self._four_facts(FactContext(overlap=overlap)))

        assert results[0] == results[1] == results[2]
        assert results[0] == dict.fromkeys(results[0], False)

    def test_stop_loss_survives_readthrough_into_a_downstream_stop(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """The stop codon *was* destroyed; that translation then runs on and meets some other stop
        in the 3' UTR does not make it retained.

        ``cds_stop_lost`` tests the reference's own stop *position*, not "no stop anywhere in the
        alternate", so a coincidental in-frame stop further downstream cannot flip it. Only the
        readthrough levels put such a stop there at all: under ``NONE`` translation halts at the
        CDS end and the alternate carries no stop, which is why the case is asserted per level.
        """
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        stop_spans = transcript.stop_codon_genomic_spans()
        assert stop_spans is not None
        first_base_position = stop_spans[0][0]
        assert _genomic_base(transcript, first_base_position) == "T"  # first base of TAA

        for level in Readthrough:
            context = FactContext(
                overlap=_create_transcript_variation_allele(
                    transcript,
                    start=first_base_position,
                    end=first_base_position,
                    reference="T",
                    alternate="C",
                    readthrough=level,
                )
            )
            peptides = get_full_cds_peptide_alleles(context)
            assert peptides is not None
            reference_peptide, alternate_peptide = peptides
            reference_stop = reference_peptide.find("*")
            assert reference_stop != -1

            alternate_stop = alternate_peptide.find("*")
            if level is Readthrough.NONE:
                assert alternate_stop == -1, level

            if alternate_stop != -1:
                assert alternate_stop > reference_stop, level

            assert get_fact(context, "cds_stop_lost") is True, level
            assert get_fact(context, "cds_stop_retained") is False, level

    def test_readthrough_actually_extends_the_peptide_the_facts_see(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """Guards against the option being silently inert at the fact layer: the same stop-losing
        variant must yield a strictly longer alternate peptide once readthrough is on."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        stop_spans = transcript.stop_codon_genomic_spans()
        assert stop_spans is not None
        position = stop_spans[0][0]

        def alternate_peptide(level: Readthrough) -> str:
            overlap = _create_transcript_variation_allele(
                transcript, start=position, end=position, reference="T", alternate="C", readthrough=level
            )
            peptides = get_full_cds_peptide_alleles(FactContext(overlap=overlap))
            assert peptides is not None
            return peptides[1]

        assert len(alternate_peptide(Readthrough.UTR)) > len(alternate_peptide(Readthrough.NONE))
