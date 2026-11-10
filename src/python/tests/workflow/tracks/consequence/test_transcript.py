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
"""Tests for the transcript consequence track."""

from __future__ import annotations

import random
from collections.abc import Iterator
from pathlib import Path
from types import SimpleNamespace
from typing import Any, cast

import pytest

from ensembl.vep.buffer import FixedSizeInputBuffer
from ensembl.vep.caches import Interval, TranscriptPool
from ensembl.vep.common import Readthrough
from ensembl.vep.config import Config
from ensembl.vep.config.sections.common.annotation import AnnotationConfig
from ensembl.vep.config.sections.common.execution import ExecutionConfig
from ensembl.vep.config.sections.common.source import SourceConfig
from ensembl.vep.models.overlap import FeatureOverlap
from ensembl.vep.models.variation import VariationAllele
from ensembl.vep.workflow.tracks.base import TrackResources, TrackResult
from ensembl.vep.workflow.tracks.consequence import transcript as transcript_module
from ensembl.vep.workflow.tracks.consequence.transcript import (
    TranscriptCandidateLookup,
    TranscriptConsequenceTrack,
    _CoverageIndex,
    build_transcript_overlap_finder,
    build_transcript_tracks,
    transcript_candidate_window,
)


def _file_config(**overrides: Any) -> Config:
    """Return a configuration naming the files a run needs."""
    source = overrides.pop("source", None) or SourceConfig(
        input_file=Path("variants.vcf"), gff=Path("transcripts.gff3"), fasta=Path("reference.fa")
    )
    return Config(source=source, **overrides)


def _resources(config: Config) -> TrackResources:
    """Return track resources over no reference sequence."""
    return TrackResources(
        source_config=config.source,
        sequence_source=None,
        synonyms=cast(Any, SimpleNamespace()),
    )


def _allele(position: int, region: str = "chr1") -> VariationAllele:
    """Return an allele-like object at one position."""
    return cast(
        VariationAllele,
        SimpleNamespace(
            start=position,
            end=position,
            slice=SimpleNamespace(seq_region_name=region),
        ),
    )


class _TranscriptSource:
    """Record regional source requests."""

    def __init__(self) -> None:
        """Initialise the query record."""
        self.calls: list[tuple[str, int, int]] = []
        # Assigned by TranscriptCandidateLookup, mirroring GXFTranscriptSource.transcript_pool --
        # declared here so static analysis sees the attribute this fake stands in for.
        self.transcript_pool: TranscriptPool | None = None
        self.transcript_pool_poisoned = False
        self.closed = False

    def get_features_by_region(
        self,
        *,
        seq_region_name: str,
        start: int,
        end: int,
    ) -> list[SimpleNamespace]:
        """Return one complete feature for each requested region."""
        self.calls.append((seq_region_name, start, end))
        return [
            SimpleNamespace(
                stable_id=f"{seq_region_name}:{start}-{end}",
                start=start,
                end=end,
                slice=SimpleNamespace(seq_region_name=seq_region_name),
            )
        ]

    def get_all_features(
        self,
    ) -> list[SimpleNamespace]:  # pragma: no cover - region index is present in this test
        """Fail if the indexed loading path regresses to a whole source read."""
        raise AssertionError("unexpected whole-source fallback")

    def close(self) -> None:
        """Record that the track released the source."""
        self.closed = True


class _StreamingTranscriptSource:
    """Expose only streaming whole-source transcript loading."""

    def __init__(self) -> None:
        """Initialise the iteration count and the pool attributes a lookup assigns."""
        self.iterations = 0
        self.transcript_pool: TranscriptPool | None = None
        self.transcript_pool_poisoned = False

    @staticmethod
    def get_all_features() -> list[SimpleNamespace]:
        """Fail if fallback materialises the complete source."""
        raise AssertionError("unexpected materialised whole-source load")

    def iter_features_by_regions(
        self,
        intervals_by_region: dict[str, list[tuple[int, int]]],
    ) -> Iterator[SimpleNamespace]:
        """Yield only complete transcripts intersecting the requested intervals."""
        self.iterations += 1
        for sequence_region in ("chr1", "chr2"):
            intervals = intervals_by_region.get(sequence_region, [])
            if any(100 <= end and 200 >= start for start, end in intervals):
                yield SimpleNamespace(
                    stable_id=f"transcript:{sequence_region}",
                    start=100,
                    end=200,
                    slice=SimpleNamespace(seq_region_name=sequence_region),
                )


class _MultiTranscriptSource:
    """Return two fixed, ordered transcript-like objects for every region query."""

    def __init__(self) -> None:
        """Initialise the pool attributes a lookup assigns."""
        self.transcript_pool: TranscriptPool | None = None
        self.transcript_pool_poisoned = False

    def get_features_by_region(self, *, seq_region_name: str, start: int, end: int) -> list[SimpleNamespace]:
        """Return the region's two features in a fixed order, ignoring the query span."""
        del start, end
        return [
            SimpleNamespace(
                stable_id=f"{seq_region_name}:A",
                start=100,
                end=150,
                slice=SimpleNamespace(seq_region_name=seq_region_name),
            ),
            SimpleNamespace(
                stable_id=f"{seq_region_name}:B",
                start=160,
                end=200,
                slice=SimpleNamespace(seq_region_name=seq_region_name),
            ),
        ]


def _lookup(
    source: object, *, window: int = 0, pool: TranscriptPool | None = None
) -> TranscriptCandidateLookup:
    """Return a lookup over a stub source."""
    return TranscriptCandidateLookup(
        source=cast(Any, source),
        window=window,
        pool=pool if pool is not None else TranscriptPool(),
    )


class TestSharedPairingParts:
    """The window and the pairing rule a run and ``vep explain`` must agree on."""

    @pytest.mark.parametrize(
        ("distance", "expected"),
        [
            pytest.param((0,), 0, id="no-flank"),
            pytest.param((5000,), 5000, id="symmetric"),
            pytest.param((100, 4000), 4000, id="downstream-wider"),
            pytest.param((4000, 100), 4000, id="upstream-wider"),
        ],
    )
    def test_candidate_window_is_the_wider_distance(self, distance: tuple[int, ...], expected: int) -> None:
        """the fetch window covers whichever direction a transcript's strand makes apply."""
        config = Config(source=SourceConfig(distance=distance))

        assert transcript_candidate_window(config) == expected

    def test_finder_carries_the_configured_distances_and_readthrough(self) -> None:
        """the pairing rule is built from the configuration, not from defaults."""
        config = Config(
            source=SourceConfig(distance=(100, 4000)),
            annotation=AnnotationConfig(readthrough="utr"),
        )

        finder = build_transcript_overlap_finder(config, lambda _buffer: [])

        assert finder.upstream_distance == 100
        assert finder.downstream_distance == 4000
        assert finder.readthrough is Readthrough.UTR

    def test_finder_pairs_the_candidates_it_is_given(self) -> None:
        """a caller's own candidate lookup is what the finder pairs, with no source of its own."""
        config = Config(source=SourceConfig(distance=(0,)))
        asked: list[object] = []

        def candidates(input_buffer: object) -> list[Any]:
            asked.append(input_buffer)
            return []

        finder = build_transcript_overlap_finder(config, candidates)
        buffer: FixedSizeInputBuffer[VariationAllele] = FixedSizeInputBuffer(objects=())

        assert not list(finder.find(buffer))
        assert asked == [buffer]


class TestTrackBuilding:
    """What ``build_transcript_tracks`` wires together."""

    @staticmethod
    def _stub_source(monkeypatch: pytest.MonkeyPatch, source: object) -> None:
        """Answer the track's source builder with a stub."""
        monkeypatch.setattr(transcript_module, "build_transcript_source", lambda _config, **_shared: source)
        monkeypatch.setattr(transcript_module, "build_annotator", lambda **_options: SimpleNamespace())

    def test_one_track_is_built_with_the_configured_batch_size(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """a run always gets exactly one transcript track."""
        self._stub_source(monkeypatch, _TranscriptSource())
        config = _file_config(annotation=AnnotationConfig(batch_size=64))

        (track,) = build_transcript_tracks(config, _resources(config))

        assert isinstance(track, TranscriptConsequenceTrack)
        assert track.identifier == "transcript"
        assert track.batch_size == 64

    def test_pool_uses_the_configured_byte_budget(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """the pool the track attaches uses ``execution.transcript_pool_bytes``."""
        source = _TranscriptSource()
        self._stub_source(monkeypatch, source)
        config = _file_config(execution=ExecutionConfig(transcript_pool_bytes=12345))

        build_transcript_tracks(config, _resources(config))

        assert isinstance(source.transcript_pool, TranscriptPool)
        assert source.transcript_pool.max_bytes == 12345

    def test_configured_flank_reaches_the_candidate_fetch(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """the track's lookup queries the source over the configured flank."""
        source = _TranscriptSource()
        self._stub_source(monkeypatch, source)
        config = _file_config(source=SourceConfig(gff=Path("transcripts.gff3"), distance=(6000,)))

        (track,) = build_transcript_tracks(config, _resources(config))
        cast(TranscriptConsequenceTrack, track).lookup.candidates(
            FixedSizeInputBuffer(objects=(_allele(10000),))
        )

        assert source.calls == [("chr1", 4000, 16000)]

    @pytest.mark.parametrize("enabled", [False, True])
    def test_fact_cost_override_reaches_the_annotator(
        self, enabled: bool, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """``annotation.fact_cost_override`` is read from the configuration unchanged."""
        received: list[dict[str, Any]] = []

        monkeypatch.setattr(
            transcript_module, "build_transcript_source", lambda _config, **_shared: _TranscriptSource()
        )

        def capture_options(**options: Any) -> SimpleNamespace:
            received.append(options)
            return SimpleNamespace()

        monkeypatch.setattr(transcript_module, "build_annotator", capture_options)
        config = _file_config(annotation=AnnotationConfig(fact_cost_override=enabled))

        build_transcript_tracks(config, _resources(config))

        assert received[0]["fact_cost_override"] is enabled


class TestCandidateLoading:
    """Worker-local complete-feature cache loading."""

    def test_contained_and_neighbouring_queries_reuse_and_extend_cached_regions(self) -> None:
        """A contained buffer is a cache hit while an adjacent buffer loads its own range."""
        source = _TranscriptSource()
        lookup = _lookup(source)

        lookup._load_region("chr1", [Interval(start=100, end=200)])  # pylint: disable=protected-access
        lookup._load_region("chr1", [Interval(start=120, end=180)])  # pylint: disable=protected-access
        lookup._load_region("chr1", [Interval(start=201, end=300)])  # pylint: disable=protected-access

        assert source.calls == [("chr1", 100, 200), ("chr1", 201, 300)]
        assert len(lookup._cache) == 2  # pylint: disable=protected-access

    def test_many_alleles_share_one_bounded_source_query(self) -> None:
        """A same-region input buffer must not issue one source query per allele."""
        source = _TranscriptSource()
        lookup = _lookup(source)
        alleles = tuple(_allele(position) for position in range(100, 200))

        lookup.candidates(FixedSizeInputBuffer(objects=alleles))

        assert source.calls == [("chr1", 100, 199)]

    def test_sparse_alleles_keep_disjoint_source_queries(self) -> None:
        """Distant same-region alleles must not fetch the unneeded interval between them."""
        source = _TranscriptSource()
        lookup = _lookup(source)
        alleles = tuple(_allele(position) for position in (1_000_000, 100))

        lookup.candidates(FixedSizeInputBuffer(objects=alleles))

        assert source.calls == [("chr1", 100, 100), ("chr1", 1_000_000, 1_000_000)]

    def test_indexed_payload_releases_previous_transcript_state(self) -> None:
        """A new indexed payload starts with an empty feature cache and coverage map."""
        lookup = _lookup(_TranscriptSource())
        # pylint: disable=protected-access
        lookup._cache.add("chr1", cast(Any, object()), Interval(start=100, end=200))
        coverage = _CoverageIndex()
        coverage.add(Interval(start=50, end=250))
        lookup._queried_regions["chr1"] = coverage  # pylint: disable=protected-access

        lookup.begin_payload()

        assert len(lookup._cache) == 0  # pylint: disable=protected-access
        assert not lookup._queried_regions  # pylint: disable=protected-access

    def test_whole_source_fallback_retains_only_matching_transcripts(self) -> None:
        """Stream an unindexed source without retaining unrelated transcript graphs."""
        source = _StreamingTranscriptSource()
        lookup = _lookup(source)

        # pylint: disable=protected-access
        lookup._load_complete_source({"chr1": [Interval(start=150, end=150)]})

        assert source.iterations == 1
        assert len(lookup._cache) == 1
        assert len(lookup._cache.get_objects("chr2", [Interval(start=150, end=150)])) == 0
        assert lookup._uses_complete_source_fallback

        lookup._load_complete_source({"chr1": [Interval(start=150, end=150)]})

        assert source.iterations == 1

    def test_unindexed_source_falls_back_once_and_keeps_serving(self) -> None:
        """A source with no region index is streamed instead of failing the payload."""

        class _UnindexedSource(_StreamingTranscriptSource):
            """Refuse region queries the way an unindexed GXF does."""

            file_path = "transcripts.gff3"

            @staticmethod
            def get_features_by_region(**_query: Any) -> list[SimpleNamespace]:
                """Fail the way a source with no region index fails."""
                raise OSError("no region index")

        source = _UnindexedSource()
        lookup = _lookup(source)

        candidates = lookup.candidates(FixedSizeInputBuffer(objects=(_allele(150),)))

        assert source.iterations == 1
        assert [transcript.stable_id for transcript in candidates] == ["transcript:chr1"]
        assert lookup._uses_complete_source_fallback  # pylint: disable=protected-access

    def test_closing_the_lookup_closes_its_source(self) -> None:
        """the track's source is released through the lookup that owns it."""
        source = _TranscriptSource()

        _lookup(source).close()

        assert source.closed


class TestOrderAcrossPayloads:
    """The hard CSQ-order invariant: insertion order into ``BinnedFeatureCache`` must not depend
    on how many payloads a lookup has already served."""

    def test_two_payloads_over_one_gene_yield_the_same_order_as_a_fresh_lookup(self) -> None:
        """a lookup that has already served one payload orders a later, repeat payload
        over the same gene exactly as a brand-new lookup would."""
        source = _MultiTranscriptSource()
        # Overlaps both fake transcripts (A: 100-150, B: 160-200) -- the fake ignores the query
        # bounds anyway, so this interval is reused both to trigger the load and to retrieve.
        interval = [Interval(start=100, end=200)]

        fresh = _lookup(source)
        fresh._load_region("chr1", interval)  # pylint: disable=protected-access
        fresh_order = [
            transcript.stable_id
            for transcript in fresh._cache.get_objects("chr1", interval)  # pylint: disable=protected-access
        ]

        reused = _lookup(source)
        reused._load_region("chr1", interval)  # pylint: disable=protected-access  # payload 1
        reused.begin_payload()
        reused._load_region("chr1", interval)  # pylint: disable=protected-access  # payload 2
        reused_order = [
            transcript.stable_id
            for transcript in reused._cache.get_objects("chr1", interval)  # pylint: disable=protected-access
        ]

        assert reused_order == fresh_order == ["chr1:A", "chr1:B"]


class TestPoolLifetime:
    """One ``TranscriptPool`` per worker, kept across payloads."""

    def test_transcript_pool_is_not_cleared_by_begin_payload(self) -> None:
        """``begin_payload`` clears the per-payload cache but keeps the pool intact."""
        source = _TranscriptSource()
        pool = TranscriptPool()
        lookup = _lookup(source, pool=pool)
        pool.put(
            ("chr1", "ENST0001", 100, 200), cast(Any, SimpleNamespace(start=100, end=200, spliced_exons=[]))
        )

        lookup.begin_payload()

        assert source.transcript_pool is pool
        assert len(pool) == 1


class TestPoolObservability:
    """Debug-level per-payload logging of the pool's running counters."""

    def test_logs_pool_counters_at_debug_level(self, caplog: pytest.LogCaptureFixture) -> None:
        """the log line carries hits/misses/evictions and the poisoned flag."""
        pool = TranscriptPool()
        pool.hit_count = 3
        pool.miss_count = 2
        lookup = _lookup(_TranscriptSource(), pool=pool)

        with caplog.at_level("DEBUG", logger=transcript_module.__name__):
            lookup.log_statistics()

        (record,) = caplog.records
        assert "transcript pool" in record.message
        assert "hits=3" in record.message
        assert "misses=2" in record.message
        assert "poisoned=False" in record.message

    def test_logs_zero_counters_when_pooling_is_disabled(self, caplog: pytest.LogCaptureFixture) -> None:
        """a disabled pool still reports, with nothing served."""
        lookup = _lookup(_TranscriptSource(), pool=TranscriptPool(max_bytes=0))

        with caplog.at_level("DEBUG", logger=transcript_module.__name__):
            lookup.log_statistics()

        (record,) = caplog.records
        assert "hits=0" in record.message
        assert "misses=0" in record.message


class TestAttachingResults:
    """Where the track puts what its annotator produced."""

    @staticmethod
    def _track() -> TranscriptConsequenceTrack:
        """Return a track whose finder and annotator are never called."""
        return TranscriptConsequenceTrack(
            identifier="transcript",
            lookup=_lookup(_TranscriptSource()),
            overlap_finder=cast(Any, SimpleNamespace()),
            annotator=cast(Any, SimpleNamespace()),
        )

    def test_consequences_reach_the_allele_through_its_overlap(self) -> None:
        """an overlap is attached to the allele it was found for, carrying its annotations."""
        allele = SimpleNamespace(transcript_variation_alleles=[])
        variation = SimpleNamespace(variation_alleles=[allele])
        overlap = transcript_module.TranscriptVariationAllele(
            feature=cast(Any, SimpleNamespace()), overlapping_feature=cast(Any, allele)
        )
        result: TrackResult[Any] = TrackResult(
            track_identifier="transcript", overlap=overlap, annotations=("missense",)
        )

        self._track().attach(cast(Any, [variation]), [result])

        assert allele.transcript_variation_alleles == [overlap]
        assert overlap.overlap_consequences == ["missense"]

    def test_rendered_window_annotations_are_released(self) -> None:
        """Attaching no results detaches the overlaps a rendered window retained."""
        allele = SimpleNamespace(transcript_variation_alleles=[object()])
        variation = SimpleNamespace(variation_alleles=[allele])

        self._track().attach(cast(Any, [variation]), ())

        assert allele.transcript_variation_alleles == []

    def test_an_overlap_of_another_kind_is_refused(self) -> None:
        """a result this track cannot have produced is reported, not dropped."""
        overlap = FeatureOverlap(
            feature=cast(Any, SimpleNamespace()), overlapping_feature=cast(Any, SimpleNamespace())
        )
        result: TrackResult[Any] = TrackResult(track_identifier="transcript", overlap=overlap, annotations=())

        with pytest.raises(TypeError, match="cannot attach"):
            self._track().attach((), [result])

    def test_closing_the_track_closes_its_source(self) -> None:
        """the workflow closing a track reaches the source it owns."""
        track = self._track()

        track.close()

        assert cast(Any, track.lookup.source).closed


class TestCoverageIndex:
    """Merged queried-range coverage index."""

    def test_empty_index_contains_nothing(self) -> None:
        """a fresh index has no covered range."""
        coverage = _CoverageIndex()

        assert not coverage.contains(Interval(start=1, end=10))

    def test_added_interval_is_contained(self) -> None:
        """an added interval, and any sub-range of it, is reported as covered."""
        coverage = _CoverageIndex()
        coverage.add(Interval(start=10, end=20))

        assert coverage.contains(Interval(start=10, end=20))
        assert coverage.contains(Interval(start=12, end=18))

    def test_query_outside_added_interval_is_not_contained(self) -> None:
        """a disjoint query is not reported as covered."""
        coverage = _CoverageIndex()
        coverage.add(Interval(start=10, end=20))

        assert not coverage.contains(Interval(start=21, end=30))
        assert not coverage.contains(Interval(start=1, end=9))

    def test_touching_intervals_are_not_merged(self) -> None:
        """two intervals with no coordinate gap stay separate entries.

        This is the conservative choice the class docstring documents: merging touching
        intervals was flagged as a risk to re-derive carefully rather than assume, so a query
        that straddles the touching boundary is answered exactly as the old ``any()`` scan over
        individual (unmerged) intervals would have answered it -- not covered, even though the
        two queries together did cover every position in the combined span.
        """
        coverage = _CoverageIndex()
        coverage.add(Interval(start=1, end=10))
        coverage.add(Interval(start=11, end=20))

        assert coverage.contains(Interval(start=1, end=10))
        assert coverage.contains(Interval(start=11, end=20))
        assert not coverage.contains(Interval(start=5, end=15))

    def test_strictly_overlapping_intervals_are_merged(self) -> None:
        """intervals sharing a coordinate merge into one covered span."""
        coverage = _CoverageIndex()
        coverage.add(Interval(start=1, end=10))
        coverage.add(Interval(start=5, end=20))

        assert coverage.contains(Interval(start=2, end=18))

    def test_gap_of_one_between_intervals_is_not_merged(self) -> None:
        """a one-position gap keeps two intervals from being treated as covering it."""
        coverage = _CoverageIndex()
        coverage.add(Interval(start=1, end=10))
        coverage.add(Interval(start=12, end=20))

        assert not coverage.contains(Interval(start=9, end=13))
        assert not coverage.contains(Interval(start=11, end=11))

    def test_nested_interval_is_absorbed_without_growing_the_index(self) -> None:
        """an interval fully inside an existing one does not add a new entry."""
        coverage = _CoverageIndex()
        coverage.add(Interval(start=1, end=100))
        coverage.add(Interval(start=10, end=20))

        assert coverage.contains(Interval(start=10, end=20))
        assert len(coverage._starts) == 1  # pylint: disable=protected-access

    def test_insertion_interval_covers_only_its_own_point(self) -> None:
        """a zero-width insertion interval is tracked by its flanking lookup point."""
        coverage = _CoverageIndex()
        insertion = Interval(start=11, end=10)  # start == end + 1: an insertion before position 11
        coverage.add(insertion)

        assert coverage.contains(insertion)
        assert not coverage.contains(Interval(start=1, end=9))

    def test_adding_a_disjoint_earlier_interval_does_not_drop_a_later_one(self) -> None:
        """inserting before an existing entry keeps that entry rather than overwriting it."""
        coverage = _CoverageIndex()
        coverage.add(Interval(start=200, end=300))
        coverage.add(Interval(start=10, end=30))

        assert coverage.contains(Interval(start=10, end=30))
        assert coverage.contains(Interval(start=200, end=300))

    def test_widening_an_interval_forward_does_not_shrink_it(self) -> None:
        """merging with a wider, already-seen interval keeps the wider bound."""
        coverage = _CoverageIndex()
        coverage.add(Interval(start=10, end=100))
        coverage.add(Interval(start=10, end=50))

        assert coverage.contains(Interval(start=10, end=100))

    def test_new_interval_merges_forward_into_an_existing_wider_one(self) -> None:
        """a new interval that overlaps a following, wider entry keeps that entry's reach."""
        coverage = _CoverageIndex()
        coverage.add(Interval(start=20, end=100))
        coverage.add(Interval(start=10, end=30))

        assert coverage.contains(Interval(start=10, end=100))

    def test_matches_brute_force_reference_over_random_sequences(self) -> None:
        """contains/add against a naive list-of-intervals reference, many random sequences."""
        random_generator = random.Random(20240914)
        for _ in range(200):
            coverage = _CoverageIndex()
            reference = _BruteForceCoverage()
            for _ in range(30):
                start = random_generator.randint(1, 200)
                end = start + random_generator.randint(0, 40)
                if random_generator.random() < 0.5:
                    coverage.add(Interval(start=start, end=end))
                    reference.add(start, end)
                else:
                    assert coverage.contains(Interval(start=start, end=end)) == reference.contains(start, end)


class _BruteForceCoverage:
    """Naive list-of-intervals reference for ``_CoverageIndex``, used only by the random test.

    Merges an added interval with every existing one it strictly overlaps (shares a coordinate),
    mirroring ``_CoverageIndex``'s merge policy without the sorted/bisect machinery.
    """

    def __init__(self) -> None:
        self._intervals: list[tuple[int, int]] = []

    def add(self, start: int, end: int) -> None:
        """Merge a new interval into the tracked set."""
        for interval_start, interval_end in self._intervals:
            if interval_start <= end and interval_end >= start:
                start = min(start, interval_start)
                end = max(end, interval_end)
        self._intervals = [
            (interval_start, interval_end)
            for interval_start, interval_end in self._intervals
            if not (interval_start <= end and interval_end >= start)
        ]
        self._intervals.append((start, end))

    def contains(self, start: int, end: int) -> bool:
        """Return whether some tracked interval fully contains ``[start, end]``."""
        return any(
            interval_start <= start and interval_end >= end
            for interval_start, interval_end in self._intervals
        )
