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

"""Unit tests for codon caches."""

from __future__ import annotations

from collections.abc import Callable

from ensembl.vep.annotators.consequence.context import AnnotatorContext
from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.annotators.consequence.strategies.columnar_staged import ColumnarStagedAnnotator
from ensembl.vep.annotators.registry import get_annotator
from ensembl.vep.corpus.synthetic import SyntheticCorpus
from ensembl.vep.caches.codon_cache import (
    CodonCache,
    CodonEdit,
    CodonEditKey,
    CodonEditLRUCache,
)
from ensembl.vep.caches.mapper_cache import MapperIndexLRUCache
from ensembl.vep.caches.zone_map import ZoneMapCache
from ensembl.vep.common import Strand
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.models.variation import VariationAllele


def _build_transcript_variation_allele(
    transcript: Transcript, *, start: int, end: int, reference: str, alternate: str
) -> TranscriptVariationAllele:
    """Create a transcript variation allele at the supplied genomic coordinates."""
    variation_allele = VariationAllele(
        start=start,
        end=end,
        strand=Strand.FORWARD,
        slice=transcript.slice,
        reference_allele=reference,
        alternate_allele=alternate,
    )
    return TranscriptVariationAllele(feature=transcript, overlapping_feature=variation_allele)


def _build_codon_edit() -> CodonEdit:
    """Create a codon edit for cache testing."""
    return CodonEdit(
        reference_codon="ATG", alternate_codon="ACG", reference_peptide="M", alternate_peptide="T"
    )


class TestCodonEditLRUCache:
    """Tests for the codon/peptide-allele LRU cache implementation."""

    def test_cache_implements_the_codon_cache_interface(self) -> None:
        """Test the concrete cache satisfies the protocol at runtime."""
        cache = CodonEditLRUCache()
        assert isinstance(cache, CodonCache)

    def test_hit_avoids_reload(self) -> None:
        """Test a cached key does not invoke the loader again."""
        cache = CodonEditLRUCache()
        key: CodonEditKey = CodonEditKey("chr1", (100, 101, 102), 0, 1, 100, 100, "A", "C")
        load_calls = 0

        def loader() -> CodonEdit:
            nonlocal load_calls
            load_calls += 1
            return _build_codon_edit()

        first = cache.get_codon_edit(key, loader)
        second = cache.get_codon_edit(key, loader)

        assert first == second
        assert load_calls == 1
        assert cache.hit_count == 1
        assert cache.miss_count == 1

    def test_lru_eviction_discards_oldest_entry(self) -> None:
        """Test the cache evicts the least-recently-used entry once over capacity."""
        cache = CodonEditLRUCache(max_entries=2)
        load_calls: list[CodonEditKey] = []

        def make_loader(key: CodonEditKey) -> Callable[[], CodonEdit]:
            def loader() -> CodonEdit:
                load_calls.append(key)
                return _build_codon_edit()

            return loader

        key_a: CodonEditKey = CodonEditKey("chr1", (1, 2, 3), 0, 1, 1, 1, "A", "C")
        key_b: CodonEditKey = CodonEditKey("chr1", (4, 5, 6), 0, 1, 4, 4, "A", "C")
        key_c: CodonEditKey = CodonEditKey("chr1", (7, 8, 9), 0, 1, 7, 7, "A", "C")

        cache.get_codon_edit(key_a, make_loader(key_a))
        cache.get_codon_edit(key_b, make_loader(key_b))
        # key_a is now the least-recently-used entry; inserting key_c evicts it.
        cache.get_codon_edit(key_c, make_loader(key_c))
        cache.get_codon_edit(key_a, make_loader(key_a))

        assert load_calls == [key_a, key_b, key_c, key_a]

    def test_disabled_cache_always_loads(self) -> None:
        """Test a disabled cache bypasses storage entirely."""
        cache = CodonEditLRUCache(enabled=False)
        load_calls = 0
        key: CodonEditKey = CodonEditKey("chr1", (100, 101, 102), 0, 1, 100, 100, "A", "C")

        def loader() -> CodonEdit:
            nonlocal load_calls
            load_calls += 1
            return _build_codon_edit()

        cache.get_codon_edit(key, loader)
        cache.get_codon_edit(key, loader)

        assert load_calls == 2
        assert cache.hit_count == 0
        assert cache.miss_count == 2

    def test_clear_removes_all_cached_entries(self) -> None:
        """Test clear() forces every subsequent lookup to reload."""
        cache = CodonEditLRUCache()
        load_calls = 0
        key: CodonEditKey = CodonEditKey("chr1", (100, 101, 102), 0, 1, 100, 100, "A", "C")

        def loader() -> CodonEdit:
            nonlocal load_calls
            load_calls += 1
            return _build_codon_edit()

        cache.get_codon_edit(key, loader)
        cache.clear()
        cache.get_codon_edit(key, loader)

        assert load_calls == 2


class TestCodonEditKey:
    """Tests for the codon edit cache key."""

    def test_none_when_variant_touches_no_codon(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """Test a deep-intronic variant (no CDS footprint) has no cache key."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        overlap = _build_transcript_variation_allele(
            transcript, start=1600, end=1600, reference="A", alternate="C"
        )
        fact_context = FactContext(overlap=overlap)
        assert CodonEditKey.build(fact_context) is None

    def test_none_when_transcript_has_no_cds_mapper(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """Test a non-coding transcript without a CDS mapper has no cache key."""
        transcript = synthetic_transcripts_by_stable_id["ENST0003"]
        overlap = _build_transcript_variation_allele(
            transcript,
            start=transcript.start,
            end=transcript.start,
            reference="A",
            alternate="C",
        )

        assert transcript.cds_mapper is None
        assert CodonEditKey.build(FactContext(overlap=overlap)) is None

    def test_same_edit_produces_the_same_key(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """Test identical edits produce equal cache keys."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        first = FactContext(
            overlap=_build_transcript_variation_allele(
                transcript, start=1051, end=1051, reference="A", alternate="C"
            )
        )
        second = FactContext(
            overlap=_build_transcript_variation_allele(
                transcript, start=1051, end=1051, reference="A", alternate="C"
            )
        )

        assert CodonEditKey.build(first) == CodonEditKey.build(second)

    def test_mapper_cache_accelerates_genomic_position_lookup(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """Test cache-key construction reuses the transcript's mapper index."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        overlap = _build_transcript_variation_allele(
            transcript, start=1051, end=1051, reference="A", alternate="C"
        )
        mapper_cache = MapperIndexLRUCache()
        annotator_context = AnnotatorContext(mapper_index_cache=mapper_cache)

        first_key = CodonEditKey.build(FactContext(overlap=overlap, annotator_context=annotator_context))
        second_key = CodonEditKey.build(FactContext(overlap=overlap, annotator_context=annotator_context))

        assert first_key == second_key
        assert mapper_cache.miss_count == 1
        assert mapper_cache.hit_count == 1

    def test_same_footprint_different_alternate_allele_produces_distinct_keys(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """Test alternate alleles at the same position produce distinct keys."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        alternate_c_context = FactContext(
            overlap=_build_transcript_variation_allele(
                transcript, start=1051, end=1051, reference="A", alternate="C"
            )
        )
        alternate_g_context = FactContext(
            overlap=_build_transcript_variation_allele(
                transcript, start=1051, end=1051, reference="A", alternate="G"
            )
        )

        key_c = CodonEditKey.build(alternate_c_context)
        key_g = CodonEditKey.build(alternate_g_context)

        assert key_c is not None
        assert key_g is not None
        assert key_c != key_g

    def test_same_alternate_allele_different_offset_in_window_produces_distinct_keys(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """Test equal alleles at different codon offsets produce distinct keys."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        # 1050, 1051, 1052 are the three bases of ENST0001's first codon (CDS starts at 1050).
        first_base = FactContext(
            overlap=_build_transcript_variation_allele(
                transcript, start=1050, end=1050, reference="A", alternate="T"
            )
        )
        second_base = FactContext(
            overlap=_build_transcript_variation_allele(
                transcript, start=1051, end=1051, reference="A", alternate="T"
            )
        )

        key_first = CodonEditKey.build(first_base)
        key_second = CodonEditKey.build(second_base)

        assert key_first is not None
        assert key_second is not None
        assert key_first != key_second
        assert key_first.codon_window_genomic_positions == key_second.codon_window_genomic_positions
        assert key_first.seq_region_name == key_second.seq_region_name
        assert key_first.variant_min != key_second.variant_min


class TestCodonCacheIntegration:
    """Tests for codon cache integration with consequence annotation."""

    def test_cached_and_uncached_annotations_match(self, full_heavy_locus_corpus: SyntheticCorpus) -> None:
        """Test cached and uncached annotations match linear annotation."""
        reference = get_annotator("consequence.linear")
        annotator_without_cache = ColumnarStagedAnnotator(zone_map=ZoneMapCache())
        cache = CodonEditLRUCache()
        annotator_with_cache = ColumnarStagedAnnotator(zone_map=ZoneMapCache(), codon_cache=cache)

        mismatches: list[str] = []
        for overlap_name, overlap in full_heavy_locus_corpus.named_overlaps:
            reference_terms = [c.term for c in reference.annotate_feature_overlap(overlap)]
            without_cache_terms = [c.term for c in annotator_without_cache.annotate_feature_overlap(overlap)]
            with_cache_terms = [c.term for c in annotator_with_cache.annotate_feature_overlap(overlap)]
            if without_cache_terms != reference_terms or with_cache_terms != reference_terms:
                mismatches.append(
                    f"{overlap_name}: reference={reference_terms} without_cache={without_cache_terms} "
                    f"with_cache={with_cache_terms}"
                )

        assert not mismatches, "\n".join(mismatches[:20])
        assert cache.hit_count > 0

    def test_cache_disabled_behaves_like_no_cache(self, small_default_corpus: SyntheticCorpus) -> None:
        """Test an enabled=False cache produces the same output as no cache at all."""
        annotator_without_cache = ColumnarStagedAnnotator(zone_map=ZoneMapCache())
        annotator_with_disabled_cache = ColumnarStagedAnnotator(
            zone_map=ZoneMapCache(), codon_cache=CodonEditLRUCache(enabled=False)
        )

        for _overlap_name, overlap in small_default_corpus.named_overlaps:
            without_terms = [c.term for c in annotator_without_cache.annotate_feature_overlap(overlap)]
            disabled_terms = [c.term for c in annotator_with_disabled_cache.annotate_feature_overlap(overlap)]
            assert without_terms == disabled_terms


class TestCodonCachePrefill:
    """Tests for prefilling consequence facts from the codon cache."""

    def test_non_coding_transcript_never_queries_the_cache(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """Test a non-coding transcript's overlap never reaches the cache at all."""
        transcript = synthetic_transcripts_by_stable_id["ENST0003"]
        assert transcript.translation is None
        overlap = _build_transcript_variation_allele(
            transcript, start=20010, end=20010, reference="A", alternate="C"
        )
        cache = CodonEditLRUCache()
        annotator = ColumnarStagedAnnotator(zone_map=ZoneMapCache(), codon_cache=cache)

        annotator.annotate_feature_overlap(overlap)

        assert cache.hit_count == 0
        assert cache.miss_count == 0
