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
"""Tests for the default runner's buffered worker."""

from __future__ import annotations

import pickle
import random
from collections.abc import Iterator
from functools import partial
from pathlib import Path
from types import SimpleNamespace
from typing import Any, cast

import pytest

from ensembl.vep.buffer import InputBuffer
from ensembl.vep.caches import BinnedFeatureCache, Interval, TranscriptPool
from ensembl.vep.config import Config
from ensembl.vep.config.sections.common.annotation import AnnotationConfig
from ensembl.vep.config.sections.common.execution import ExecutionConfig
from ensembl.vep.config.sections.common.source import SourceConfig
from ensembl.vep.output.factory import AnnotatedVariant
from ensembl.vep.runner.run import worker as worker_module
from ensembl.vep.runner.run.worker import (
    AnnotatedBufferResult,
    BufferedWorker,
    _CoverageIndex,
    annotate_buffer,
    build_buffered_worker,
)
from ensembl.vep.stats import RunSummaryDelta


class _TranscriptSource:
    """Record regional source requests."""

    def __init__(self) -> None:
        """Initialise the query record."""
        self.calls: list[tuple[str, int, int]] = []
        # Assigned by BufferedWorker.__init__, mirroring GXFTranscriptSource.transcript_pool --
        # declared here so static analysis sees the attribute this fake stands in for.
        self.transcript_pool: TranscriptPool | None = None
        self.transcript_pool_poisoned = False

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


class _StreamingTranscriptSource:
    """Expose only streaming whole-source transcript loading."""

    iterations: int = 0

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


def _stub_sources(monkeypatch: pytest.MonkeyPatch, source: object) -> None:
    """Build the worker's sources from a stub transcript source instead of real files.

    The worker asks ``sources.build`` for its sources; these tests are about what it does with
    them, so each builder is answered with the stub the test wants.
    """

    monkeypatch.setattr(
        worker_module, "build_sequence_source", lambda _config: SimpleNamespace(unknown_regions=set())
    )
    monkeypatch.setattr(worker_module, "build_annotation_sources", lambda _config, **_shared: [source])
    monkeypatch.setattr(worker_module, "build_input_source", lambda _config, **_shared: SimpleNamespace())
    monkeypatch.setattr(worker_module, "transcript_source_of", lambda sources: sources[0])


class TestBufferedTranscriptLoading:
    """Worker-local complete-feature cache loading."""

    def test_worker_entry_points_are_process_safe(self) -> None:
        """Keep worker construction and execution serialisable for spawn processes."""
        worker_factory = partial(build_buffered_worker, Config())

        assert pickle.loads(pickle.dumps(worker_factory)).func is build_buffered_worker
        assert pickle.loads(pickle.dumps(annotate_buffer)) is annotate_buffer

    def test_worker_results_are_process_safe(self) -> None:
        """Serialise detached annotations without carrying their live model graph."""
        variant = AnnotatedVariant(
            id="rs1",
            seq_region_name="1",
            start=10,
            end=10,
            strand=1,
            allele_string="A/G",
            annotations=({"Consequence": ["missense_variant"]},),
        )
        result = AnnotatedBufferResult(
            index=3,
            variants=(variant,),
            summary_delta=RunSummaryDelta(variants=1),
        )

        restored = pickle.loads(pickle.dumps(result))

        assert restored == result

    @pytest.mark.parametrize(
        "flank", [pytest.param(4000, id="below-default"), pytest.param(6000, id="above-default")]
    )
    def test_configured_flank_reaches_compiler_and_candidate_fetch(
        self,
        flank: int,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Give workflow compilation and source retrieval the same non-default flank."""
        source = _TranscriptSource()
        compiled_configs: list[Config] = []
        compiled_workflow = SimpleNamespace()

        def capture_workflow(config: Config, context: Any) -> SimpleNamespace:
            del context
            compiled_configs.append(config)
            return compiled_workflow

        _stub_sources(monkeypatch, source)
        monkeypatch.setattr(worker_module, "build_annotator", lambda **_options: SimpleNamespace())
        monkeypatch.setattr(worker_module, "compile_workflow", capture_workflow)
        config = Config(
            source=SourceConfig(
                input_file=Path("variants.vcf"),
                gff=Path("transcripts.gff3"),
                fasta=Path("reference.fa"),
                distance=flank,
            )
        )
        worker = BufferedWorker(config)
        allele = SimpleNamespace(
            start=10000,
            end=10000,
            slice=SimpleNamespace(seq_region_name="chr1"),
        )

        worker._transcript_candidates(InputBuffer(objects=(allele,)))  # pylint: disable=protected-access

        assert compiled_configs == [config]
        assert worker.workflow is compiled_workflow
        assert source.calls == [("chr1", 10000 - flank, 10000 + flank)]

    def test_contained_and_neighbouring_queries_reuse_and_extend_cached_regions(self) -> None:
        """A contained buffer is a cache hit while an adjacent buffer loads its own range."""
        source = _TranscriptSource()
        worker = object.__new__(BufferedWorker)
        worker.config = Config(source=SourceConfig(distance=0))
        worker.transcript_source = cast(Any, source)
        worker.transcript_cache = BinnedFeatureCache()
        worker._queried_regions = {}  # pylint: disable=protected-access
        worker._uses_complete_source_fallback = False  # pylint: disable=protected-access

        worker._load_transcripts("chr1", [Interval(start=100, end=200)])  # pylint: disable=protected-access
        worker._load_transcripts("chr1", [Interval(start=120, end=180)])  # pylint: disable=protected-access
        worker._load_transcripts("chr1", [Interval(start=201, end=300)])  # pylint: disable=protected-access

        assert source.calls == [("chr1", 100, 200), ("chr1", 201, 300)]
        assert len(worker.transcript_cache) == 2

    def test_many_alleles_share_one_bounded_source_query(self) -> None:
        """A same-region input buffer must not issue one source query per allele."""
        source = _TranscriptSource()
        worker = object.__new__(BufferedWorker)
        worker.config = Config(source=SourceConfig(distance=0))
        worker.transcript_source = cast(Any, source)
        worker.transcript_cache = BinnedFeatureCache()
        worker._queried_regions = {}  # pylint: disable=protected-access
        worker._uses_complete_source_fallback = False  # pylint: disable=protected-access

        alleles = tuple(
            SimpleNamespace(
                start=position,
                end=position,
                slice=SimpleNamespace(seq_region_name="chr1"),
            )
            for position in range(100, 200)
        )
        worker._transcript_candidates(InputBuffer(objects=alleles))  # pylint: disable=protected-access

        assert source.calls == [("chr1", 100, 199)]

    def test_sparse_alleles_keep_disjoint_source_queries(self) -> None:
        """Distant same-region alleles must not fetch the unneeded interval between them."""
        source = _TranscriptSource()
        worker = object.__new__(BufferedWorker)
        worker.config = Config(source=SourceConfig(distance=0))
        worker.transcript_source = cast(Any, source)
        worker.transcript_cache = BinnedFeatureCache()
        worker._queried_regions = {}  # pylint: disable=protected-access
        worker._uses_complete_source_fallback = False  # pylint: disable=protected-access
        alleles = tuple(
            SimpleNamespace(
                start=position,
                end=position,
                slice=SimpleNamespace(seq_region_name="chr1"),
            )
            for position in (1_000_000, 100)
        )

        worker._transcript_candidates(InputBuffer(objects=alleles))  # pylint: disable=protected-access

        assert source.calls == [("chr1", 100, 100), ("chr1", 1_000_000, 1_000_000)]

    def test_rendered_window_annotations_are_released(self) -> None:
        """Clearing a completed window detaches overlaps from each allele."""
        allele = SimpleNamespace(transcript_variation_alleles=[object()])
        variation = SimpleNamespace(variation_alleles=[allele])

        BufferedWorker._clear_transient_annotations(  # pylint: disable=protected-access
            cast(list[Any], [variation])
        )

        assert allele.transcript_variation_alleles == []

    def test_indexed_payload_releases_previous_transcript_state(self) -> None:
        """A new indexed payload starts with an empty feature cache and coverage map."""
        worker = object.__new__(BufferedWorker)
        worker.transcript_cache = BinnedFeatureCache()
        worker.transcript_cache.add(
            "chr1",
            object(),
            Interval(start=100, end=200),
        )
        coverage = _CoverageIndex()
        coverage.add(Interval(start=50, end=250))
        worker._queried_regions = {"chr1": coverage}  # pylint: disable=protected-access
        worker._uses_complete_source_fallback = False  # pylint: disable=protected-access

        worker._begin_payload()  # pylint: disable=protected-access

        assert len(worker.transcript_cache) == 0
        assert not worker._queried_regions  # pylint: disable=protected-access

    def test_unindexed_payload_releases_its_filtered_source_fallback(self) -> None:
        """An unindexed source does not retain one payload's transcripts for the next."""
        worker = object.__new__(BufferedWorker)
        worker.transcript_cache = BinnedFeatureCache()
        worker.transcript_cache.add(
            "chr1",
            object(),
            Interval(start=100, end=200),
        )
        region_coverage = _CoverageIndex()
        region_coverage.add(Interval(start=50, end=250))
        coverage = {"chr1": region_coverage}
        worker._queried_regions = coverage  # pylint: disable=protected-access
        worker._uses_complete_source_fallback = True  # pylint: disable=protected-access

        worker._begin_payload()  # pylint: disable=protected-access

        assert len(worker.transcript_cache) == 0
        assert not coverage

    def test_whole_source_fallback_retains_only_matching_transcripts(self) -> None:
        """Stream an unindexed source without retaining unrelated transcript graphs."""
        source = _StreamingTranscriptSource()
        worker = object.__new__(BufferedWorker)
        worker.transcript_source = cast(Any, source)
        worker.transcript_cache = BinnedFeatureCache()
        worker._queried_regions = {}  # pylint: disable=protected-access
        worker._uses_complete_source_fallback = False  # pylint: disable=protected-access

        worker._load_complete_transcript_source(  # pylint: disable=protected-access
            {"chr1": [Interval(start=150, end=150)]}
        )

        assert source.iterations == 1
        assert len(worker.transcript_cache) == 1
        assert len(worker.transcript_cache.get_objects("chr2", [Interval(start=150, end=150)])) == 0
        assert worker._uses_complete_source_fallback  # pylint: disable=protected-access

        worker._load_complete_transcript_source(  # pylint: disable=protected-access
            {"chr1": [Interval(start=150, end=150)]}
        )

        assert source.iterations == 1


class _MultiTranscriptSource:
    """Return two fixed, ordered transcript-like objects for every region query."""

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


class TestOrderAcrossPayloads:
    """The hard CSQ-order invariant: insertion order into ``BinnedFeatureCache`` must not depend
    on how many payloads a worker has already processed."""

    @staticmethod
    def _new_worker(source: _MultiTranscriptSource) -> BufferedWorker:
        worker = object.__new__(BufferedWorker)
        worker.config = Config(source=SourceConfig(distance=0))
        worker.transcript_source = cast(Any, source)
        worker.transcript_cache = BinnedFeatureCache()
        worker._queried_regions = {}  # pylint: disable=protected-access
        worker._uses_complete_source_fallback = False  # pylint: disable=protected-access
        return worker

    def test_two_payloads_over_one_gene_yield_the_same_order_as_a_fresh_worker(self) -> None:
        """Test a worker that has already processed one payload orders a later, repeat payload
        over the same gene exactly as a brand-new worker would."""
        source = _MultiTranscriptSource()
        # Overlaps both fake transcripts (A: 100-150, B: 160-200) -- the fake ignores the query
        # bounds anyway, so this interval is reused both to trigger the load and to retrieve.
        interval = [Interval(start=100, end=200)]

        fresh_worker = self._new_worker(source)
        fresh_worker._load_transcripts("chr1", interval)  # pylint: disable=protected-access
        fresh_order = [
            transcript.stable_id for transcript in fresh_worker.transcript_cache.get_objects("chr1", interval)
        ]

        reused_worker = self._new_worker(source)
        reused_worker._load_transcripts("chr1", interval)  # pylint: disable=protected-access  # payload 1
        reused_worker._begin_payload()  # pylint: disable=protected-access
        reused_worker._load_transcripts("chr1", interval)  # pylint: disable=protected-access  # payload 2
        reused_order = [
            transcript.stable_id
            for transcript in reused_worker.transcript_cache.get_objects("chr1", interval)
        ]

        assert reused_order == fresh_order == ["chr1:A", "chr1:B"]


class TestCoverageIndex:
    """Tests for the merged, sorted queried-range tracker behind ``_queried_regions``."""

    def test_empty_index_contains_nothing(self) -> None:
        """Test a fresh index has no covered range."""
        coverage = _CoverageIndex()

        assert not coverage.contains(Interval(start=1, end=10))

    def test_added_interval_is_contained(self) -> None:
        """Test an added interval, and any sub-range of it, is reported as covered."""
        coverage = _CoverageIndex()
        coverage.add(Interval(start=10, end=20))

        assert coverage.contains(Interval(start=10, end=20))
        assert coverage.contains(Interval(start=12, end=18))

    def test_query_outside_added_interval_is_not_contained(self) -> None:
        """Test a disjoint query is not reported as covered."""
        coverage = _CoverageIndex()
        coverage.add(Interval(start=10, end=20))

        assert not coverage.contains(Interval(start=21, end=30))
        assert not coverage.contains(Interval(start=1, end=9))

    def test_touching_intervals_are_not_merged(self) -> None:
        """Test two intervals with no coordinate gap stay separate entries.

        This is the conservative choice the module docstring documents: merging touching
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
        """Test intervals sharing a coordinate merge into one covered span."""
        coverage = _CoverageIndex()
        coverage.add(Interval(start=1, end=10))
        coverage.add(Interval(start=5, end=20))

        assert coverage.contains(Interval(start=2, end=18))

    def test_gap_of_one_between_intervals_is_not_merged(self) -> None:
        """Test a one-position gap keeps two intervals from being treated as covering it."""
        coverage = _CoverageIndex()
        coverage.add(Interval(start=1, end=10))
        coverage.add(Interval(start=12, end=20))

        assert not coverage.contains(Interval(start=9, end=13))
        assert not coverage.contains(Interval(start=11, end=11))

    def test_nested_interval_is_absorbed_without_growing_the_index(self) -> None:
        """Test an interval fully inside an existing one does not add a new entry."""
        coverage = _CoverageIndex()
        coverage.add(Interval(start=1, end=100))
        coverage.add(Interval(start=10, end=20))

        assert coverage.contains(Interval(start=10, end=20))
        assert len(coverage._starts) == 1  # pylint: disable=protected-access

    def test_insertion_interval_covers_only_its_own_point(self) -> None:
        """Test a zero-width insertion interval is tracked by its flanking lookup point."""
        coverage = _CoverageIndex()
        insertion = Interval(start=11, end=10)  # start == end + 1: an insertion before position 11
        coverage.add(insertion)

        assert coverage.contains(insertion)
        assert not coverage.contains(Interval(start=1, end=9))

    def test_adding_a_disjoint_earlier_interval_does_not_drop_a_later_one(self) -> None:
        """Test inserting before an existing entry keeps that entry rather than overwriting it."""
        coverage = _CoverageIndex()
        coverage.add(Interval(start=200, end=300))
        coverage.add(Interval(start=10, end=30))

        assert coverage.contains(Interval(start=10, end=30))
        assert coverage.contains(Interval(start=200, end=300))

    def test_widening_an_interval_forward_does_not_shrink_it(self) -> None:
        """Test merging with a wider, already-seen interval keeps the wider bound."""
        coverage = _CoverageIndex()
        coverage.add(Interval(start=10, end=100))
        coverage.add(Interval(start=10, end=50))

        assert coverage.contains(Interval(start=10, end=100))

    def test_new_interval_merges_forward_into_an_existing_wider_one(self) -> None:
        """Test a new interval that overlaps a following, wider entry keeps that entry's reach."""
        coverage = _CoverageIndex()
        coverage.add(Interval(start=20, end=100))
        coverage.add(Interval(start=10, end=30))

        assert coverage.contains(Interval(start=10, end=100))

    def test_matches_brute_force_reference_over_random_sequences(self) -> None:
        """Test contains/add against a naive list-of-intervals reference, many random sequences."""
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


class TestTranscriptPoolWiring:
    """The worker builds one ``TranscriptPool`` per process and keeps it across payloads."""

    @staticmethod
    def _build_worker(
        monkeypatch: pytest.MonkeyPatch, source: _TranscriptSource, config: Config
    ) -> BufferedWorker:
        _stub_sources(monkeypatch, source)
        monkeypatch.setattr(worker_module, "build_annotator", lambda **_options: SimpleNamespace())
        monkeypatch.setattr(worker_module, "compile_workflow", lambda _config, _context: SimpleNamespace())
        return BufferedWorker(config)

    def test_worker_wires_a_pool_with_the_configured_byte_budget(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Test the pool the worker attaches uses ``execution.transcript_pool_bytes``."""
        source = _TranscriptSource()
        config = Config(
            source=SourceConfig(
                input_file=Path("variants.vcf"), gff=Path("transcripts.gff3"), fasta=Path("reference.fa")
            ),
            execution=ExecutionConfig(transcript_pool_bytes=12345),
        )

        self._build_worker(monkeypatch, source, config)

        assert isinstance(source.transcript_pool, TranscriptPool)
        assert source.transcript_pool.max_bytes == 12345

    def test_transcript_pool_is_not_cleared_by_begin_payload(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Test ``_begin_payload`` clears the per-payload cache but keeps the pool intact."""
        source = _TranscriptSource()
        config = Config(
            source=SourceConfig(
                input_file=Path("variants.vcf"), gff=Path("transcripts.gff3"), fasta=Path("reference.fa")
            )
        )
        worker = self._build_worker(monkeypatch, source, config)
        pool = source.transcript_pool
        assert pool is not None
        pool.put(
            ("chr1", "ENST0001", 100, 200), cast(Any, SimpleNamespace(start=100, end=200, spliced_exons=[]))
        )

        worker._begin_payload()  # pylint: disable=protected-access

        assert worker.transcript_source.transcript_pool is pool
        assert len(pool) == 1


class TestFactCostOverrideWiring:
    """The worker hands the configured fact-cost choice to the annotator builder."""

    @pytest.mark.parametrize("enabled", [False, True])
    def test_worker_forwards_fact_cost_override(self, enabled: bool, monkeypatch: pytest.MonkeyPatch) -> None:
        """``annotation.fact_cost_override`` reaches ``build_annotator`` unchanged."""
        received: list[dict[str, Any]] = []

        def capture_options(**options: Any) -> SimpleNamespace:
            received.append(options)
            return SimpleNamespace()

        _stub_sources(monkeypatch, _TranscriptSource())
        monkeypatch.setattr(worker_module, "build_annotator", capture_options)
        monkeypatch.setattr(worker_module, "compile_workflow", lambda _config, _context: SimpleNamespace())
        config = Config(
            source=SourceConfig(
                input_file=Path("variants.vcf"), gff=Path("transcripts.gff3"), fasta=Path("reference.fa")
            ),
            annotation=AnnotationConfig(fact_cost_override=enabled),
        )

        BufferedWorker(config)

        assert received[0]["fact_cost_override"] is enabled


class TestTranscriptPoolObservability:
    """Debug-level per-payload logging of the pool's running counters."""

    def test_logs_pool_counters_at_debug_level(self, caplog: pytest.LogCaptureFixture) -> None:
        """Test the log line carries hits/misses/evictions and the poisoned flag."""
        worker = object.__new__(BufferedWorker)
        pool = TranscriptPool()
        pool.hit_count = 3
        pool.miss_count = 2
        worker.transcript_pool = pool
        worker.transcript_source = cast(Any, SimpleNamespace(transcript_pool_poisoned=False))

        with caplog.at_level("DEBUG", logger=worker_module.__name__):
            worker._log_transcript_pool_stats()  # pylint: disable=protected-access

        (record,) = caplog.records
        assert "transcript pool" in record.message
        assert "hits=3" in record.message
        assert "misses=2" in record.message
        assert "poisoned=False" in record.message

    def test_logs_zero_counters_when_pooling_is_disabled(self, caplog: pytest.LogCaptureFixture) -> None:
        """Test a disabled pool still reports, with nothing served."""
        worker = object.__new__(BufferedWorker)
        worker.transcript_pool = TranscriptPool(max_bytes=0)
        worker.transcript_source = cast(Any, SimpleNamespace(transcript_pool_poisoned=False))

        with caplog.at_level("DEBUG", logger=worker_module.__name__):
            worker._log_transcript_pool_stats()  # pylint: disable=protected-access

        (record,) = caplog.records
        assert "hits=0" in record.message
        assert "misses=0" in record.message
