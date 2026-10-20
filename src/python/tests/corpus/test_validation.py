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
"""Tests for corpus membership validation."""

from __future__ import annotations

import pytest

from ensembl.vep.corpus.combine import combine_corpora
from ensembl.vep.corpus.spec import CorpusSpec
from ensembl.vep.corpus.synthetic import SyntheticCorpus
from ensembl.vep.corpus.synthetic.default import default_corpus
from ensembl.vep.corpus.synthetic.heavy import heavy_corpus
from ensembl.vep.corpus.validation import validate_corpus
from ensembl.vep.models.sequence import Slice


class TestCorpusValidation:
    """Cross-collection relationship checks."""

    def test_rejects_duplicate_sequence_region_names(self) -> None:
        """A flat-file corpus cannot serialise separate objects under one name."""
        corpus = default_corpus(fuzz_count=0)
        sequence_region = corpus.synthetic_slice
        duplicate = Slice(
            coord_system=sequence_region.coord_system,
            seq_region_name=sequence_region.seq_region_name,
            start=sequence_region.start,
            end=sequence_region.end,
            strand=sequence_region.strand,
        )
        with pytest.raises(ValueError, match="duplicate corpus sequence region name"):
            validate_corpus((sequence_region, duplicate), (), (), ())

    def test_rejects_an_overlap_without_its_named_variant(self) -> None:
        """Overlap names are stable joins back to the variant collection."""
        corpus = default_corpus(fuzz_count=0)
        _name, overlap = corpus.named_overlaps[0]
        with pytest.raises(ValueError, match="has no corpus variant"):
            validate_corpus(
                corpus.sequence_regions, corpus.named_variants, corpus.transcripts, (("missing", overlap),)
            )

    def test_rejects_duplicate_variant_names(self) -> None:
        """Variant names remain unique stable join keys."""
        corpus = default_corpus(fuzz_count=0)
        first = corpus.named_variants[0]
        with pytest.raises(ValueError, match="duplicate"):
            SyntheticCorpus(
                synthetic_slice=corpus.synthetic_slice,
                named_variants=(first, first),
            )

    def test_rejects_an_overlap_attached_to_another_allele(self) -> None:
        """An overlap must reference the exact allele registered under its name."""
        corpus = default_corpus(fuzz_count=0)
        overlap_name, overlap = corpus.named_overlaps[0]
        other_allele = corpus.named_variants[1][1]
        with pytest.raises(ValueError, match="does not refer"):
            SyntheticCorpus(
                synthetic_slice=corpus.synthetic_slice,
                named_variants=((overlap_name, other_allele),),
                transcripts=corpus.transcripts,
                named_overlaps=((overlap_name, overlap),),
            )

    def test_rejects_an_overlap_feature_outside_the_corpus(self) -> None:
        """The overlap view cannot introduce a transcript absent from membership."""
        corpus = default_corpus(fuzz_count=0)
        overlap_name, overlap = corpus.named_overlaps[0]
        with pytest.raises(ValueError, match="transcript outside"):
            SyntheticCorpus(
                synthetic_slice=corpus.synthetic_slice,
                named_variants=corpus.named_variants,
                named_overlaps=((overlap_name, overlap),),
            )

    def test_rejects_an_overlap_feature_with_different_slice_to_feature(self) -> None:
        """The overlapped allele and transcript can't be on different slices."""
        corpus_default = default_corpus(fuzz_count=0)
        corpus_heavy = heavy_corpus()
        corpus = combine_corpora(CorpusSpec(), (corpus_default, corpus_heavy))
        sequence_regions = corpus.sequence_regions
        _overlap_name, overlap = corpus.named_overlaps[0]
        overlap.feature.slice = (
            sequence_regions[1] if overlap.feature.slice == sequence_regions[0] else sequence_regions[0]
        )
        with pytest.raises(ValueError, match="transcript and variant on different sequence regions"):
            validate_corpus(
                sequence_regions, corpus.named_variants, corpus.transcripts, corpus.named_overlaps
            )

    def test_rejects_duplicate_transcript_membership(self) -> None:
        """One transcript object appears at most once in a corpus."""
        corpus = default_corpus(fuzz_count=0)
        transcript = corpus.transcripts[0]
        with pytest.raises(ValueError, match="duplicate corpus transcript"):
            SyntheticCorpus(
                synthetic_slice=corpus.synthetic_slice,
                named_variants=corpus.named_variants,
                transcripts=(transcript, transcript),
            )

    def test_rejects_duplicate_variant_transcript_overlaps(self) -> None:
        """One variant and transcript pair has one overlap record."""
        corpus = default_corpus(fuzz_count=0)
        overlap = corpus.named_overlaps[0]
        with pytest.raises(ValueError, match="duplicate overlap"):
            SyntheticCorpus(
                synthetic_slice=corpus.synthetic_slice,
                named_variants=corpus.named_variants,
                transcripts=corpus.transcripts,
                named_overlaps=(overlap, overlap),
            )
