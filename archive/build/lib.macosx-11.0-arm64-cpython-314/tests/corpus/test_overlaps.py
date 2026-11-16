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
"""Tests for the generated overlap source."""

from __future__ import annotations

import hashlib
import json

import pytest

import ensembl.vep.annotators  # noqa: F401  pylint: disable=unused-import

# ^ keeps the annotators-before-caches import order.
from ensembl.vep.corpus.base import as_corpus
from ensembl.vep.corpus.combine import combine_corpora
from ensembl.vep.corpus.overlaps import KINDS, SyntheticOverlapSource
from ensembl.vep.corpus.spec import CorpusSpec
from ensembl.vep.corpus.synthetic import (
    SyntheticCorpus,
    boundary_variant_alleles,
    build_synthetic_genome,
    build_transcript_variation_alleles,
    fuzz_variant_alleles,
    heavy_locus_variant_alleles,
)
from ensembl.vep.corpus.validation import validate_corpus
from ensembl.vep.models.sequence import Slice


class TestCorpora:
    """Which corpora a kind draws from."""

    def test_legacy_source_adapts_to_corpus_and_rejects_other_objects(self) -> None:
        """The retired facade has one explicit bridge and unrelated values do not duck-type in."""
        source = SyntheticOverlapSource(kind="default", fuzz=0)
        assert as_corpus(source) is source.corpus
        with pytest.raises(TypeError, match="expected a Corpus"):
            as_corpus(object())

    def test_combining_no_corpora_is_rejected(self) -> None:
        """A combined corpus must still contain a concrete dataset."""
        with pytest.raises(ValueError, match="at least one"):
            combine_corpora(CorpusSpec(kind="synthetic.all"), ())

    @pytest.mark.parametrize("kind,expected", [("default", 1), ("heavy", 1), ("both", 2)])
    def test_the_kind_decides_how_many_are_built(self, kind: str, expected: int) -> None:
        """``both`` is the only one that builds more than one."""
        source = SyntheticOverlapSource(kind=kind, fuzz=0, step=600)
        assert len(source.corpora) == expected

    def test_every_declared_kind_works(self) -> None:
        """A kind offered on the command line must build something."""
        for kind in KINDS:
            assert SyntheticOverlapSource(kind=kind, fuzz=0, step=600).corpora

    def test_an_unknown_kind_is_rejected(self) -> None:
        """Named in the message, since it reaches the user as a usage error."""
        with pytest.raises(ValueError, match="nope"):
            _ = SyntheticOverlapSource(kind="nope").corpora


class TestStability:
    """Walking a source twice must give the same overlaps."""

    def test_the_corpora_are_built_once(self) -> None:
        """A conformance gate then a timing pass must see the same overlaps."""
        source = SyntheticOverlapSource(kind="default", fuzz=0)
        assert source.corpora is source.corpora

    def test_two_walks_agree(self) -> None:
        """named_overlaps is a generator, so it must not consume the corpus."""
        source = SyntheticOverlapSource(kind="default", fuzz=0)
        first = [name for name, _overlap in source.named_overlaps()]
        second = [name for name, _overlap in source.named_overlaps()]
        assert first == second

    def test_the_fuzz_count_is_honoured(self) -> None:
        """The knob a caller uses to keep a check fast."""
        small = SyntheticOverlapSource(kind="default", fuzz=0)
        larger = SyntheticOverlapSource(kind="default", fuzz=50)
        assert len(larger.corpora[0].named_variants) == len(small.corpora[0].named_variants) + 50
        assert len(list(larger.named_overlaps())) > len(list(small.named_overlaps()))


class TestGeneratorValidation:
    """Public low-level generators reject invalid sizes directly."""

    def test_negative_fuzz_count_is_rejected(self) -> None:
        """Calling the generator directly has the same invariant as ``CorpusSpec``."""
        synthetic_slice, _transcripts = build_synthetic_genome()
        with pytest.raises(ValueError, match="variant_count"):
            fuzz_variant_alleles(synthetic_slice, variant_count=-1)

    def test_non_positive_heavy_step_is_rejected(self) -> None:
        """The error is domain-specific rather than leaking ``range`` behavior."""
        synthetic_slice, _transcripts = build_synthetic_genome()
        with pytest.raises(ValueError, match="step"):
            heavy_locus_variant_alleles(synthetic_slice, step=0)

    def test_negative_overlap_window_is_rejected(self) -> None:
        """A negative feature window is never a meaningful matching policy."""
        synthetic_slice, transcripts = build_synthetic_genome()
        variants = fuzz_variant_alleles(synthetic_slice, variant_count=1)
        with pytest.raises(ValueError, match="window"):
            build_transcript_variation_alleles(variants, transcripts, window=-1)


class TestCorpusIntegrity:
    """The two views of a corpus cannot silently disagree."""

    def test_named_variants_are_required(self) -> None:
        """Legacy overlap-only construction fails loudly instead of exporting an empty VCF."""
        source = SyntheticOverlapSource(kind="default", fuzz=0)
        corpus = source.corpora[0]
        with pytest.raises(TypeError, match="named_variants"):
            SyntheticCorpus(  # type: ignore[call-arg]
                synthetic_slice=corpus.synthetic_slice,
                transcripts=corpus.transcripts,
                named_overlaps=corpus.named_overlaps,
            )

    def test_sequence_region_names_are_unique_across_objects(self) -> None:
        """A flat-file corpus cannot serialise two distinct regions under the same name."""
        corpus = SyntheticOverlapSource(kind="default", fuzz=0).corpora[0]
        sequence_region = corpus.synthetic_slice
        duplicate_name = Slice(
            coord_system=sequence_region.coord_system,
            seq_region_name=sequence_region.seq_region_name,
            start=sequence_region.start,
            end=sequence_region.end,
            strand=sequence_region.strand,
        )
        with pytest.raises(ValueError, match="duplicate corpus sequence region name"):
            validate_corpus((sequence_region, duplicate_name), (), (), ())

    def test_duplicate_variant_names_are_rejected(self) -> None:
        """A name is the stable join key and therefore must be unique."""
        corpus = SyntheticOverlapSource(kind="default", fuzz=0).corpora[0]
        first = corpus.named_variants[0]
        with pytest.raises(ValueError, match="duplicate"):
            SyntheticCorpus(
                synthetic_slice=corpus.synthetic_slice,
                named_variants=(first, first),
            )

    def test_mutable_input_collections_are_frozen(self) -> None:
        """A caller passing lists cannot mutate cached corpus membership afterwards."""
        corpus = SyntheticOverlapSource(kind="default", fuzz=0).corpora[0]
        rebuilt = SyntheticCorpus(
            synthetic_slice=corpus.synthetic_slice,
            named_variants=list(corpus.named_variants),
            transcripts=list(corpus.transcripts),
            named_overlaps=list(corpus.named_overlaps),
        )
        assert isinstance(rebuilt.named_variants, tuple)
        assert isinstance(rebuilt.transcripts, tuple)
        assert isinstance(rebuilt.named_overlaps, tuple)

    def test_an_overlap_must_reference_its_named_variant(self) -> None:
        """Features cannot be attached to another same-named allele by accident."""
        corpus = SyntheticOverlapSource(kind="default", fuzz=0).corpora[0]
        overlap_name, overlap = corpus.named_overlaps[0]
        other_allele = corpus.named_variants[1][1]
        with pytest.raises(ValueError, match="does not refer"):
            SyntheticCorpus(
                synthetic_slice=corpus.synthetic_slice,
                named_variants=((overlap_name, other_allele),),
                transcripts=corpus.transcripts,
                named_overlaps=((overlap_name, overlap),),
            )

    def test_an_overlap_name_must_resolve_to_a_variant(self) -> None:
        """Dangling overlap join keys are rejected at construction."""
        corpus = SyntheticOverlapSource(kind="default", fuzz=0).corpora[0]
        _overlap_name, overlap = corpus.named_overlaps[0]
        with pytest.raises(ValueError, match="has no corpus variant"):
            SyntheticCorpus(
                synthetic_slice=corpus.synthetic_slice,
                named_variants=corpus.named_variants,
                transcripts=corpus.transcripts,
                named_overlaps=(("missing", overlap),),
            )

    def test_an_overlap_feature_must_belong_to_the_corpus(self) -> None:
        """The overlap view cannot introduce a transcript absent from the corpus."""
        corpus = SyntheticOverlapSource(kind="default", fuzz=0).corpora[0]
        overlap_name, overlap = corpus.named_overlaps[0]
        with pytest.raises(ValueError, match="transcript outside"):
            SyntheticCorpus(
                synthetic_slice=corpus.synthetic_slice,
                named_variants=corpus.named_variants,
                named_overlaps=((overlap_name, overlap),),
            )

    def test_duplicate_transcript_membership_is_rejected(self) -> None:
        """A transcript appears at most once in a corpus."""
        corpus = SyntheticOverlapSource(kind="default", fuzz=0).corpora[0]
        transcript = corpus.transcripts[0]
        with pytest.raises(ValueError, match="duplicate corpus transcript"):
            SyntheticCorpus(
                synthetic_slice=corpus.synthetic_slice,
                named_variants=corpus.named_variants,
                transcripts=(transcript, transcript),
            )

    def test_duplicate_variant_transcript_overlap_is_rejected(self) -> None:
        """One variant-feature pair has one overlap record."""
        corpus = SyntheticOverlapSource(kind="default", fuzz=0).corpora[0]
        overlap = corpus.named_overlaps[0]
        with pytest.raises(ValueError, match="duplicate overlap"):
            SyntheticCorpus(
                synthetic_slice=corpus.synthetic_slice,
                named_variants=corpus.named_variants,
                transcripts=corpus.transcripts,
                named_overlaps=(overlap, overlap),
            )


class TestBoundaryVariantStability:
    """The boundary suite is an ordered seeded fixture, not merely a set of cases."""

    def test_the_default_boundary_variant_stream_is_unchanged(self) -> None:
        """Names, coordinates, alleles, order, and seeded alternate bases remain stable."""
        synthetic_slice, transcripts = build_synthetic_genome()
        signature = [
            (name, allele.start, allele.end, allele.reference_allele, allele.alternate_allele)
            for name, allele in boundary_variant_alleles(synthetic_slice, transcripts)
        ]
        encoded_signature = json.dumps(signature, separators=(",", ":")).encode()

        assert len(signature) == 655
        assert hashlib.sha256(encoded_signature).hexdigest() == (
            "e2df39a1494dad381324a15d2e58eb94a5048efeeb7116e7601b4b6e37dd40b0"
        )
