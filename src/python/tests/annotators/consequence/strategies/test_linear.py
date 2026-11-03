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
"""Tests for the linear consequence annotator."""

from __future__ import annotations

import pytest

from ensembl.vep.annotators.consequence.strategies.linear import LinearFactsAnnotator
from ensembl.vep.caches.codon_cache import CodonEditLRUCache
from ensembl.vep.caches.zone_map import ZoneMapCache
from ensembl.vep.common import Strand
from ensembl.vep.common.codons import translate_dna
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.models.variation import VariationAllele

# (overlap key, expected term list) pins for consequence.linear, derived from the biology of the
# synthetic corpus (see ensembl.vep.corpus.synthetic.default.build_genome):
#   ENST0001 forward protein_coding, ENST0002 reverse protein_coding,
#   ENST0003 forward lincRNA-like (non-coding), ENST0004 forward NMD.
# Keys match boundary_overlaps' "<transcript stable_id>::<variant name>" convention.
EXPECTED_TERMS_BY_OVERLAP_KEY: dict[str, list[str]] = {
    # -- forward-strand protein-coding transcript (ENST0001) --
    "ENST0001::ENST0001_start_codon_pos1050": ["start_lost"],
    "ENST0001::ENST0001_stop_codon_pos3197": ["stop_lost"],
    "ENST0001::ENST0001_cds3_del1": ["frameshift_variant"],
    "ENST0001::ENST0001_cds3_del3": ["inframe_deletion"],
    "ENST0001::ENST0001_intron0_donor+1": ["splice_donor_variant"],
    "ENST0001::ENST0001_intron0_acceptor-1": ["splice_acceptor_variant"],
    "ENST0001::ENST0001_intron0_acceptor-17": ["splice_polypyrimidine_tract_variant", "intron_variant"],
    "ENST0001::ENST0001_utr5_edge": ["5_prime_UTR_variant"],
    "ENST0001::ENST0001_utr3_edge": ["3_prime_UTR_variant"],
    "ENST0001::ENST0001_insertion_at_start": ["5_prime_UTR_variant"],
    "ENST0001::ENST0001_upstream5000": ["upstream_gene_variant"],
    "ENST0001::ENST0001_downstream5000": ["downstream_gene_variant"],
    # -- reverse-strand protein-coding transcript (ENST0002): 5'/3' mirrored vs. ENST0001 --
    "ENST0002::ENST0002_start_codon_pos12100": ["start_lost"],
    "ENST0002::ENST0002_stop_codon_pos10102": ["stop_lost"],
    "ENST0002::ENST0002_intron0_donor+1": ["splice_donor_variant"],
    "ENST0002::ENST0002_intron0_acceptor-1": ["splice_acceptor_variant"],
    "ENST0002::ENST0002_utr5_edge": ["5_prime_UTR_variant"],
    "ENST0002::ENST0002_utr3_edge": ["3_prime_UTR_variant"],
    "ENST0002::ENST0002_insertion_at_start": ["3_prime_UTR_variant"],
    # -- non-coding transcript (ENST0003) --
    "ENST0003::ENST0003_noncoding_exon": ["non_coding_transcript_exon_variant"],
    "ENST0003::ENST0003_intron0_donor+1": ["splice_donor_variant", "non_coding_transcript_variant"],
    # -- NMD-biotype transcript (ENST0004): coding terms plus NMD_transcript_variant alongside --
    "ENST0004::ENST0004_start_codon_pos30050": ["start_lost", "NMD_transcript_variant"],
    "ENST0004::ENST0004_stop_codon_pos31197": ["stop_lost", "NMD_transcript_variant"],
    "ENST0004::ENST0004_cds3_del3": ["inframe_deletion", "NMD_transcript_variant"],
    "ENST0004::ENST0004_intron0_donor+1": ["splice_donor_variant", "NMD_transcript_variant"],
}


class _NotATranscriptOverlap:
    """A stand-in object that is not a TranscriptVariationAllele, for the type-check test."""


@pytest.fixture(scope="session")
def linear_annotator() -> LinearFactsAnnotator:
    """Create a linear consequence annotator fixture."""
    return LinearFactsAnnotator()


class TestLinearAnnotatorBoundaryCorpus:
    """Pin exact expected terms for a curated set of boundary variants (plan Section 9)."""

    @pytest.mark.parametrize(
        ("overlap_key", "expected_terms"),
        list(EXPECTED_TERMS_BY_OVERLAP_KEY.items()),
        ids=list(EXPECTED_TERMS_BY_OVERLAP_KEY.keys()),
    )
    def test_expected_terms_for_named_boundary_overlap(
        self,
        linear_annotator: LinearFactsAnnotator,
        boundary_overlaps: dict[str, TranscriptVariationAllele],
        overlap_key: str,
        expected_terms: list[str],
    ) -> None:
        """Test that each named boundary overlap produces the expected terms."""
        overlap = boundary_overlaps[overlap_key]
        consequences = linear_annotator.annotate_feature_overlap(overlap)
        assert [consequence.term for consequence in consequences] == expected_terms

    def test_every_boundary_overlap_annotates_without_error(
        self,
        linear_annotator: LinearFactsAnnotator,
        boundary_overlaps: dict[str, TranscriptVariationAllele],
    ) -> None:
        """Every boundary overlap annotates without error."""
        for overlap in boundary_overlaps.values():
            consequences = linear_annotator.annotate_feature_overlap(overlap)
            assert isinstance(consequences, list)

    def test_no_overlap_ever_produces_the_default_term_itself(
        self,
        linear_annotator: LinearFactsAnnotator,
        boundary_overlaps: dict[str, TranscriptVariationAllele],
    ) -> None:
        """The annotator returns [] rather than [intergenic_variant]; the caller applies the default."""
        for overlap in boundary_overlaps.values():
            consequences = linear_annotator.annotate_feature_overlap(overlap)
            assert "intergenic_variant" not in {consequence.term for consequence in consequences}


class TestFiveThousandBaseFlankBoundary:
    """The 4999/5000/5001 upstream/downstream boundary, exact to the base.

    Built directly (rather than through the auto-generated boundary corpus) against ENST0002,
    which starts at genomic 10000 -- far enough from the chromosome's own start that clamping
    the flank position to stay within the chromosome (as the corpus builder does for variants
    that would otherwise fall off the left edge) cannot interfere with the boundary itself.
    """

    @pytest.mark.parametrize(
        ("distance", "expect_upstream"),
        [(4999, True), (5000, True), (5001, False)],
        ids=["4999bp_is_upstream", "5000bp_is_upstream", "5001bp_is_not_upstream"],
    )
    def test_reverse_strand_upstream_boundary(
        self,
        linear_annotator: LinearFactsAnnotator,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        distance: int,
        expect_upstream: bool,
    ) -> None:
        """Test the inclusive upstream boundary on a reverse-strand transcript."""
        # ENST0002 is reverse strand, so transcript.end is its 5' (upstream) side -- the higher
        # genomic coordinate. Place the variant beyond transcript.end for "upstream" there.
        transcript = synthetic_transcripts_by_stable_id["ENST0002"]
        position = transcript.end + distance
        variant = VariationAllele(
            start=position,
            end=position,
            strand=Strand.FORWARD,
            slice=transcript.slice,
            reference_allele="A",
            alternate_allele="T",
        )
        overlap = TranscriptVariationAllele(feature=transcript, overlapping_feature=variant)
        consequences = linear_annotator.annotate_feature_overlap(overlap)
        terms = {consequence.term for consequence in consequences}
        if expect_upstream:
            assert terms == {"upstream_gene_variant"}
        else:
            assert terms == set()


class TestMultiBaseSubstitution:
    """A substitution longer than one base, which VEP calls missense like any other."""

    def test_an_amino_acid_change_is_missense(
        self,
        linear_annotator: LinearFactsAnnotator,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
    ) -> None:
        """Two bases replaced inside one codon, changing its residue: missense, not protein altering."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        coding_sequence = transcript.coding_sequence
        assert coding_sequence is not None
        codon_offset = 30
        reference_codon = coding_sequence[codon_offset : codon_offset + 3]
        alternate_bases = next(
            first + second
            for first in "ACGT"
            for second in "ACGT"
            if translate_dna(first + second + reference_codon[2]) not in ("*", translate_dna(reference_codon))
        )
        start = transcript.cds_offset_to_genomic(codon_offset)
        assert start is not None and transcript.cds_offset_to_genomic(codon_offset + 1) == start + 1
        variation_allele = VariationAllele(
            start=start,
            end=start + 1,
            strand=Strand.FORWARD,
            slice=transcript.slice,
            reference_allele=reference_codon[:2],
            alternate_allele=alternate_bases,
        )
        overlap = TranscriptVariationAllele(feature=transcript, overlapping_feature=variation_allele)

        consequences = linear_annotator.annotate_feature_overlap(overlap)

        assert [consequence.term for consequence in consequences] == ["missense_variant"]


class TestTypeChecking:
    """Test overlap type validation."""

    def test_rejects_non_transcript_variation_allele_overlaps(
        self, linear_annotator: LinearFactsAnnotator
    ) -> None:
        """Test that non-transcript variation allele overlaps raise TypeError."""
        with pytest.raises(TypeError):
            linear_annotator.annotate_feature_overlap(_NotATranscriptOverlap())  # type: ignore[arg-type]


class TestGuardDigestCache:
    """Test guard digest caching."""

    def test_repeated_identical_guard_flags_hit_the_cache(
        self, boundary_overlaps: dict[str, TranscriptVariationAllele]
    ) -> None:
        """Repeated identical guard flags hit the cache."""
        annotator = LinearFactsAnnotator()
        first_key = next(iter(boundary_overlaps))
        overlap = boundary_overlaps[first_key]

        annotator.annotate_feature_overlap(overlap)
        assert annotator.cache_misses == 1
        assert annotator.cache_hits == 0

        annotator.annotate_feature_overlap(overlap)
        assert annotator.cache_hits == 1


class TestOptionalCaches:
    """Tests for the optional strategy caches."""

    def test_zone_map_is_used_without_changing_terms(
        self, boundary_overlaps: dict[str, TranscriptVariationAllele]
    ) -> None:
        """Test the enabled zone map serves facts without changing the result."""
        overlap = boundary_overlaps["ENST0001::ENST0001_cds3_del1"]
        expected_terms = [
            consequence.term for consequence in LinearFactsAnnotator().annotate_feature_overlap(overlap)
        ]
        zone_map = ZoneMapCache()
        annotator = LinearFactsAnnotator(zone_map=zone_map)

        consequences = annotator.annotate_feature_overlap(overlap)

        assert [consequence.term for consequence in consequences] == expected_terms
        assert annotator.annotator_context.zone_map_cache is zone_map
        assert zone_map.maps_built == 1
        assert zone_map.facts_served > 0

    def test_codon_cache_is_reused_without_changing_terms(
        self, boundary_overlaps: dict[str, TranscriptVariationAllele]
    ) -> None:
        """Test the enabled codon cache serves a repeated coding edit."""
        overlap = boundary_overlaps["ENST0001::ENST0001_cds3_del1"]
        codon_cache = CodonEditLRUCache()
        annotator = LinearFactsAnnotator(codon_cache=codon_cache)

        consequences = annotator.annotate_feature_overlaps([overlap, overlap])

        assert consequences[1] == consequences[0]
        assert annotator.annotator_context.codon_cache is codon_cache
        assert codon_cache.miss_count == 1
        assert codon_cache.hit_count == 1
