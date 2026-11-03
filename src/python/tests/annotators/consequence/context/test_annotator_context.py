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
"""Tests for the ``AnnotatorContext`` class."""

from __future__ import annotations

import pytest

from ensembl.vep.annotators.base import Annotator
from ensembl.vep.annotators.consequence.strategies.linear import LinearFactsAnnotator
from ensembl.vep.annotators.consequence.context import AnnotatorContext
from ensembl.vep.caches.codon_cache import CodonEditLRUCache
from ensembl.vep.caches.zone_map import ZoneMapCache
from ensembl.vep.corpus.synthetic import SyntheticCorpus
from ensembl.vep.models.overlap import TranscriptVariationAllele


@pytest.fixture
def overlaps(small_default_corpus: SyntheticCorpus) -> list[TranscriptVariationAllele]:
    """Create a sample of transcript variation alleles for cache testing."""
    return [overlap for _name, overlap in small_default_corpus.named_overlaps][:200]


def _annotated_terms(annotator: Annotator, overlaps: list[TranscriptVariationAllele]) -> list[list[str]]:
    """Return the ordered consequence terms for each transcript variation allele."""
    return [
        [consequence.term for consequence in annotator.annotate_feature_overlap(overlap)]
        for overlap in overlaps
    ]


class TestDefaultIsolation:
    """With no context passed, each annotator owns a private one."""

    def test_two_annotators_do_not_share_by_default(self) -> None:
        """Separate constructions get separate contexts."""
        first, second = LinearFactsAnnotator(), LinearFactsAnnotator()
        assert first.annotator_context is not second.annotator_context

    def test_a_private_context_starts_empty(self) -> None:
        """A fresh context carries no caches and no counts."""
        context = LinearFactsAnnotator().annotator_context
        assert context.codon_cache is None
        assert context.digest_candidate_cache == {}
        assert context.cache_hits == 0
        assert context.cache_misses == 0


class TestSharing:
    """Passing one context to two annotators shares its caches, without changing what they produce."""

    def test_candidate_cache_entries_are_visible_to_the_second_annotator(
        self, overlaps: list[TranscriptVariationAllele]
    ) -> None:
        """The second annotator inherits a warm candidate cache and so records only hits."""
        context = AnnotatorContext()
        _annotated_terms(LinearFactsAnnotator(annotator_context=context), overlaps)
        entries_after_first = len(context.digest_candidate_cache)
        assert entries_after_first > 0

        misses_after_first = context.cache_misses
        _annotated_terms(LinearFactsAnnotator(annotator_context=context), overlaps)

        # The second pass sees exactly the same guard digests, so it adds no new entries and
        # misses nothing -- which is only possible because the cache is on the context.
        assert len(context.digest_candidate_cache) == entries_after_first
        assert context.cache_misses == misses_after_first

    def test_zone_maps_are_built_once_across_two_annotators(
        self, overlaps: list[TranscriptVariationAllele]
    ) -> None:
        """Two zone-map annotators on one context build each transcript's map once between them."""
        context = AnnotatorContext()
        first = LinearFactsAnnotator(annotator_context=context, zone_map=ZoneMapCache())
        _annotated_terms(first, overlaps)
        zone_map_cache = context.zone_map_cache
        assert zone_map_cache is not None
        maps_after_first = zone_map_cache.maps_built
        assert maps_after_first > 0

        second = LinearFactsAnnotator(annotator_context=context)
        _annotated_terms(second, overlaps)
        assert zone_map_cache.maps_built == maps_after_first

    def test_codon_edits_are_computed_once_across_two_annotators(
        self, overlaps: list[TranscriptVariationAllele]
    ) -> None:
        """Two annotators on one context translate each distinct codon edit once between them."""
        context = AnnotatorContext()
        codon_cache = CodonEditLRUCache()
        _annotated_terms(LinearFactsAnnotator(annotator_context=context, codon_cache=codon_cache), overlaps)
        assert context.codon_cache is codon_cache
        misses_after_first = codon_cache.miss_count
        hits_after_first = codon_cache.hit_count
        assert misses_after_first > 0

        _annotated_terms(LinearFactsAnnotator(annotator_context=context), overlaps)

        # The second pass asks for exactly the same codon edits, so it translates none of them
        # again -- which is only possible because the cache is on the context.
        assert codon_cache.miss_count == misses_after_first
        assert codon_cache.hit_count > hits_after_first

    def test_counters_are_shared_not_summed_per_annotator(
        self, overlaps: list[TranscriptVariationAllele]
    ) -> None:
        """Both annotators report the context's counters, so they agree."""
        context = AnnotatorContext()
        first = LinearFactsAnnotator(annotator_context=context)
        second = LinearFactsAnnotator(annotator_context=context)
        _annotated_terms(first, overlaps)
        _annotated_terms(second, overlaps)
        assert first.cache_hits == second.cache_hits == context.cache_hits
        assert first.cache_misses == second.cache_misses == context.cache_misses

    def test_shared_and_isolated_contexts_produce_identical_terms(
        self, overlaps: list[TranscriptVariationAllele]
    ) -> None:
        """A warm shared context yields exactly what a cold private one does."""
        shared = AnnotatorContext()
        _annotated_terms(LinearFactsAnnotator(annotator_context=shared), overlaps)
        # The second annotator sees a warm candidate cache and a warm zone map cache and
        # the third does not see any caches.
        warm_context_terms = _annotated_terms(LinearFactsAnnotator(annotator_context=shared), overlaps)
        isolated_context_terms = _annotated_terms(LinearFactsAnnotator(), overlaps)
        assert warm_context_terms == isolated_context_terms

    @pytest.mark.parametrize(
        ("use_zone_map", "use_codon_cache"),
        [
            pytest.param(True, False, id="zone_map"),
            pytest.param(False, True, id="codon_cache"),
            pytest.param(True, True, id="both"),
        ],
    )
    def test_a_cached_annotator_matches_a_plain_one_under_a_shared_context(
        self,
        overlaps: list[TranscriptVariationAllele],
        use_zone_map: bool,
        use_codon_cache: bool,
    ) -> None:
        """Sharing a context between a cached and an uncached annotator keeps both correct.

        The candidate caches are keyed on cheap-guard values, which are a property of the overlap
        and the rule set rather than of which caches are switched on -- so a cached scan and a
        plain one can safely populate each other's. This asserts that, since getting it wrong
        would show up as one serving the other the wrong candidate rules.

        The caches are built inside the test rather than in the table, so no two cases share one.
        """
        context = AnnotatorContext()
        cached = LinearFactsAnnotator(
            annotator_context=context,
            zone_map=ZoneMapCache() if use_zone_map else None,
            codon_cache=CodonEditLRUCache() if use_codon_cache else None,
        )
        cached_terms = _annotated_terms(cached, overlaps)
        plain_terms = _annotated_terms(LinearFactsAnnotator(annotator_context=context), overlaps)
        assert cached_terms == plain_terms


class TestResetCandidateCaches:
    """Clearing is the caller's escape hatch when the rule set changes underneath a context."""

    def test_reset_drops_entries_and_counters(self, overlaps: list[TranscriptVariationAllele]) -> None:
        """Both candidate caches and both counters go back to empty."""
        context = AnnotatorContext()
        _annotated_terms(LinearFactsAnnotator(annotator_context=context), overlaps)
        assert context.digest_candidate_cache

        context.reset_candidate_caches()
        assert context.digest_candidate_cache == {}
        assert context.bitmask_candidate_cache == {}
        assert context.cache_hits == 0
        assert context.cache_misses == 0

    def test_reset_leaves_the_zone_map_cache_alone(self, overlaps: list[TranscriptVariationAllele]) -> None:
        """Zone maps depend on the transcript, not the rule set, so they survive a reset."""
        context = AnnotatorContext()
        annotator = LinearFactsAnnotator(annotator_context=context, zone_map=ZoneMapCache())
        _annotated_terms(annotator, overlaps)
        zone_map_cache = context.zone_map_cache
        assert zone_map_cache is not None
        maps_built = zone_map_cache.maps_built

        context.reset_candidate_caches()
        _annotated_terms(annotator, overlaps)
        assert zone_map_cache.maps_built == maps_built

    def test_reset_leaves_the_codon_cache_alone(self, overlaps: list[TranscriptVariationAllele]) -> None:
        """Codon edits depend on the sequence, not the rule set, so they survive a reset."""
        context = AnnotatorContext()
        codon_cache = CodonEditLRUCache()
        annotator = LinearFactsAnnotator(annotator_context=context, codon_cache=codon_cache)
        _annotated_terms(annotator, overlaps)
        misses = codon_cache.miss_count
        assert misses > 0

        context.reset_candidate_caches()
        _annotated_terms(annotator, overlaps)
        assert codon_cache.miss_count == misses
