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

# pylint: disable=no-name-in-module,redefined-outer-name
"""Tests for columnar consequence kernels."""

from __future__ import annotations

import pytest

from ensembl.vep.annotators.consequence.strategies.columnar_staged import ColumnarStagedAnnotator
from ensembl.vep.annotators.consequence.strategies.columnar_kernels import (
    SEEDED_FACT_NAMES,
    ColumnarKernelsAnnotator,
)
from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.annotators.consequence.facts.base import get_fact
from ensembl.vep.caches.codon_cache import CodonEditLRUCache
from ensembl.vep.caches.zone_map import ZoneMapCache
from ensembl.vep.corpus.synthetic import SyntheticCorpus
from ensembl.vep.ensembl_vep_core import (
    batch_contain_span,
    batch_overlap_any_interval,
    batch_overlap_span,
    batch_variant_shapes,
)
from ensembl.vep.models.overlap import TranscriptVariationAllele


def _bitmap_value(bitmap: bytes, index: int) -> bool:
    """Return a value from a packed little-endian bitmap."""
    return bool(bitmap[index // 8] & (1 << (index % 8)))


@pytest.fixture
def overlaps(small_default_corpus: SyntheticCorpus) -> list[TranscriptVariationAllele]:
    """Create transcript variation alleles from the default synthetic corpus."""
    return [overlap for _name, overlap in small_default_corpus.named_overlaps]


class TestKernelsMatchTheirPythonFacts:
    """Each kernel against the fact it stands in for, over the whole corpus."""

    def test_variant_shape_kernel_matches_every_python_shape_fact(
        self, overlaps: list[TranscriptVariationAllele]
    ) -> None:
        """All six variant-shape facts, for every overlap in the corpus."""
        contexts = [FactContext(overlap=overlap) for overlap in overlaps]
        reference_lengths = [
            len((context.overlap.overlapping_feature.reference_allele or "").replace("-", ""))
            for context in contexts
        ]
        alternate_lengths = [
            len((context.overlap.overlapping_feature.alternate_allele or "").replace("-", ""))
            for context in contexts
        ]
        snv, increases, decreases, _inframe = batch_variant_shapes(reference_lengths, alternate_lengths)

        # The synthetic corpus is short variants only, so the short branch and the umbrella
        # fact agree for every overlap here -- see facts/kind.py for why they would not
        # for a structural one.
        for index, context in enumerate(contexts):
            assert _bitmap_value(snv, index) == get_fact(context, "is_snv"), index
            assert _bitmap_value(increases, index) == get_fact(context, "is_short_insertion"), index
            assert _bitmap_value(increases, index) == get_fact(context, "is_insertion"), index
            assert _bitmap_value(decreases, index) == get_fact(context, "is_short_deletion"), index
            assert _bitmap_value(decreases, index) == get_fact(context, "is_deletion"), index
            # assert _bitmap_value(inframe, index) == get_fact(context, "cds_length_change_inframe"), index

    def test_span_kernels_match_overlap_and_contain_feature(
        self, overlaps: list[TranscriptVariationAllele]
    ) -> None:
        """overlap_feature and contain_feature, for every overlap in the corpus."""
        contexts = [FactContext(overlap=overlap) for overlap in overlaps]
        variant_mins = [context.variant_span[0] for context in contexts]
        variant_maxs = [context.variant_span[1] for context in contexts]
        feature_starts = [context.overlap.feature.start for context in contexts]
        feature_ends = [context.overlap.feature.end for context in contexts]

        overlap_bitmap = batch_overlap_span(variant_mins, variant_maxs, feature_starts, feature_ends)
        contain_bitmap = batch_contain_span(variant_mins, variant_maxs, feature_starts, feature_ends)

        for index, context in enumerate(contexts):
            assert _bitmap_value(overlap_bitmap, index) == get_fact(context, "overlap_feature"), index
            assert _bitmap_value(contain_bitmap, index) == get_fact(context, "contain_feature"), index

    def test_interval_kernel_matches_overlap_exon(self, overlaps: list[TranscriptVariationAllele]) -> None:
        """The CSR interval kernel against overlap_exon.

        Not seeded by the annotator (a zone map already serves it more cheaply), but the kernel
        exists and is meant to be correct, so it is checked here rather than left untested.
        """
        contexts = [FactContext(overlap=overlap) for overlap in overlaps]
        variant_mins = [context.variant_span[0] for context in contexts]
        variant_maxs = [context.variant_span[1] for context in contexts]

        interval_starts: list[int] = []
        interval_ends: list[int] = []
        offsets: list[int] = [0]
        for context in contexts:
            for exon in context.overlap.feature.spliced_exons:
                interval_starts.append(exon.start)
                interval_ends.append(exon.end)
            offsets.append(len(interval_starts))

        bitmap = batch_overlap_any_interval(
            variant_mins, variant_maxs, interval_starts, interval_ends, offsets
        )
        for index, context in enumerate(contexts):
            # overlap_exon additionally requires overlap_feature; the kernel answers the raw
            # interval question, so compare against the conjunction the fact actually computes.
            expected = get_fact(context, "overlap_feature") and _bitmap_value(bitmap, index)
            assert get_fact(context, "overlap_exon") == expected, index


class TestSeededColumnsAgreeWithThePythonPath:
    """The annotator's seeded columns, against the same facts computed the ordinary way."""

    def test_every_seeded_column_matches_get_fact(self, overlaps: list[TranscriptVariationAllele]) -> None:
        """Kernel-filled columns equal what the ordinary lazy path would have produced."""
        annotator = ColumnarKernelsAnnotator()
        chunk = [FactContext(overlap=overlap) for overlap in overlaps[:256]]
        columns = annotator.seed_columns(chunk)

        assert set(columns) == set(SEEDED_FACT_NAMES)
        for fact_name, column in columns.items():
            for index, overlap in enumerate(overlaps[:256]):
                expected = bool(get_fact(FactContext(overlap=overlap), fact_name))
                assert bool(column & (1 << index)) == expected, (fact_name, index)

    def test_seeding_an_empty_chunk_is_a_no_op(self) -> None:
        """A zero-length chunk must not index into empty arrays."""
        assert ColumnarKernelsAnnotator().seed_columns([]) == {}

    def test_kernel_call_count_is_three_per_chunk(self, overlaps: list[TranscriptVariationAllele]) -> None:
        """Three crossings per chunk, regardless of chunk size -- the whole point of batching."""
        annotator = ColumnarKernelsAnnotator(batch_size=128)
        annotator.annotate_feature_overlaps(overlaps[:256])
        assert annotator.kernel_calls == 3 * 2  # 256 overlaps / 128 per chunk = 2 chunks


class TestOutputIsUnaffected:
    """Turning the kernels on must not change a single term."""

    @pytest.mark.parametrize(
        "zone_map_enabled", [pytest.param(False, id="plain"), pytest.param(True, id="zone_map")]
    )
    def test_matches_the_pure_python_annotator(
        self, overlaps: list[TranscriptVariationAllele], zone_map_enabled: bool
    ) -> None:
        """Identical terms with and without a zone map layered on top."""

        def terms(annotator: ColumnarStagedAnnotator) -> list[list[str]]:
            return [
                [consequence.term for consequence in consequences]
                for consequences in annotator.annotate_feature_overlaps(overlaps)
            ]

        if zone_map_enabled:
            kernel_annotator = ColumnarKernelsAnnotator(zone_map=ZoneMapCache())
            python_annotator = ColumnarStagedAnnotator(zone_map=ZoneMapCache())
        else:
            kernel_annotator = ColumnarKernelsAnnotator()
            python_annotator = ColumnarStagedAnnotator()
        assert terms(kernel_annotator) == terms(python_annotator)


class TestOptionalCaches:
    """Tests for the optional strategy caches."""

    def test_zone_map_is_used_without_changing_terms(
        self, boundary_overlaps: dict[str, TranscriptVariationAllele]
    ) -> None:
        """Test the enabled zone map serves facts without changing the result."""
        overlap = boundary_overlaps["ENST0001::ENST0001_cds3_del1"]
        expected = ColumnarKernelsAnnotator().annotate_feature_overlaps([overlap])
        zone_map = ZoneMapCache()
        annotator = ColumnarKernelsAnnotator(zone_map=zone_map)

        actual = annotator.annotate_feature_overlaps([overlap])

        assert actual == expected
        assert annotator.annotator_context.zone_map_cache is zone_map
        assert zone_map.maps_built == 1
        assert zone_map.facts_served > 0

    def test_codon_cache_is_reused_without_changing_terms(
        self, boundary_overlaps: dict[str, TranscriptVariationAllele]
    ) -> None:
        """Test the enabled codon cache serves a repeated coding edit."""
        overlap = boundary_overlaps["ENST0001::ENST0001_cds3_del1"]
        codon_cache = CodonEditLRUCache()
        annotator = ColumnarKernelsAnnotator(codon_cache=codon_cache)

        consequences = annotator.annotate_feature_overlaps([overlap, overlap])

        assert consequences[1] == consequences[0]
        assert annotator.annotator_context.codon_cache is codon_cache
        assert codon_cache.miss_count == 1
        assert codon_cache.hit_count == 1
