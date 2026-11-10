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
"""Tests for payload execution."""

from __future__ import annotations

import multiprocessing
from collections.abc import Callable, Generator, Iterable, Iterator
from concurrent.futures import Future
from typing import cast

import pytest

from ensembl.vep.common.enums import StartMethod
from ensembl.vep.config import load
from ensembl.vep.config.sections.common.execution import ExecutionConfig
from ensembl.vep.workflow import executor as executor_module
from ensembl.vep.workflow.executor import ExecutionSpec, Executor, PayloadIterator


def _region_of(record: str) -> str:
    """Return the sequence region a ``<region>:<name>`` test record sits on."""
    return record.partition(":")[0]


class _Multiplier:
    """Worker multiplying primitive payloads by a configured value."""

    built = 0

    def __init__(self, factor: int) -> None:
        """Store the configured factor and count this construction."""
        self.factor = factor
        type(self).built += 1

    def run(self, payload: int) -> int:
        """Return the payload multiplied by this worker's factor."""
        return self.factor * payload


def _build_multiplier() -> _Multiplier:
    """Build a process-safe worker."""
    return _Multiplier(factor=2)


def _executor(payloads: Iterable[int], specification: ExecutionSpec) -> Executor[int, int]:
    """Return an executor multiplying its payloads under one execution specification."""
    return Executor(
        specification=specification,
        worker_factory=_build_multiplier,
        payload_iterator=payloads,
    )


def _execute(payloads: Iterable[int], specification: ExecutionSpec) -> Iterator[int]:
    """Return the ordered results of one execution specification."""
    return _executor(payloads, specification).execute()


class _FakeFuture(Future[object]):
    """Record cancellation requests for completed test futures."""

    cancel_requested: bool = False

    def cancel(self) -> bool:
        """Record cancellation before delegating to the standard future."""
        self.cancel_requested = True
        return super().cancel()


class _FakeProcessPoolExecutor:
    """Record submitted futures and shutdown requests."""

    instances: list[_FakeProcessPoolExecutor] = []

    def __init__(self, **_kwargs: object) -> None:
        """Record this executor instance."""
        self.submitted_futures: list[_FakeFuture] = []
        self.shutdown_calls: list[tuple[bool, bool]] = []
        self.terminate_calls = 0
        self.instances.append(self)

    def submit(
        self,
        _function: Callable[[object], object],
        payload: object,
    ) -> Future[object]:
        """Return a completed future for a payload."""
        future = _FakeFuture()
        future.set_result(payload)
        self.submitted_futures.append(future)
        return future

    def shutdown(self, wait: bool = True, *, cancel_futures: bool = False) -> None:
        """Record the requested shutdown behaviour."""
        self.shutdown_calls.append((wait, cancel_futures))

    def terminate_workers(self) -> None:
        """Record immediate termination of active worker processes."""
        self.terminate_calls += 1


@pytest.fixture(name="process_pool")
def fixture_process_pool(
    monkeypatch: pytest.MonkeyPatch,
) -> type[_FakeProcessPoolExecutor]:
    """Replace process execution with a deterministic executor."""
    _FakeProcessPoolExecutor.instances.clear()
    monkeypatch.setattr(executor_module, "ProcessPoolExecutor", _FakeProcessPoolExecutor)
    return _FakeProcessPoolExecutor


class TestExecutionSpec:
    """Execution option translation."""

    def test_fork_spelling_sets_the_process_count(self) -> None:
        """``--fork`` on the command line reaches the specification as ``processes``."""
        config = load(argv=["--fork", "3"], section_names=("execution",), environ={})

        assert ExecutionSpec.from_config(config.execution).processes == 3

    def test_a_parallel_run_forks_where_forking_is_safe(self) -> None:
        """A forked worker inherits the loaded package, which is what keeps the parallel path fast."""
        specification = ExecutionSpec(processes=2)

        context = specification.multiprocessing_context()

        if "fork" in multiprocessing.get_all_start_methods():
            assert context is not None
            assert context.get_start_method() == "fork"
        else:
            assert context is None

    def test_a_single_worker_keeps_the_platform_default(self) -> None:
        """Nothing is started, so the context does not matter."""
        assert ExecutionSpec(processes=1).multiprocessing_context() is None

    def test_a_threaded_run_does_not_fork(self) -> None:
        """Forking a process that already has threads is unsafe, so the platform decides."""
        assert ExecutionSpec(processes=2, threads=2).multiprocessing_context() is None

    def test_an_explicit_start_method_wins(self) -> None:
        """A start-method option is more specific than the automatic choice."""
        specification = ExecutionSpec(processes=2, start_method=StartMethod.SPAWN)

        context = specification.multiprocessing_context()

        assert context is not None
        assert context.get_start_method() == "spawn"

    def test_explicit_context_is_available(self) -> None:
        """A requested start method builds its matching multiprocessing context."""
        specification = ExecutionSpec(start_method=StartMethod.SPAWN)

        assert specification.multiprocessing_context() is not None

    def test_default_queue_matches_the_worker_count(self) -> None:
        """An omitted limit occupies each worker without declaring a hard memory bound."""
        specification = ExecutionSpec(processes=3)

        assert specification.max_in_flight is None
        assert specification.effective_max_in_flight == 3
        assert ExecutionConfig().max_in_flight is None

    def test_explicit_queue_limit_remains_a_hard_cap(self) -> None:
        """A caller may trade worker utilisation for a stricter memory bound."""
        specification = ExecutionSpec(processes=3, max_in_flight=1)

        assert specification.effective_max_in_flight == 1

    def test_non_positive_explicit_queue_limit_is_rejected(self) -> None:
        """An explicit cap must still leave room for one submitted payload."""
        with pytest.raises(ValueError, match="max_in_flight"):
            ExecutionConfig(max_in_flight=0)

    @pytest.mark.parametrize(
        "build_specification,message",
        [
            (lambda: ExecutionSpec(processes=0), "processes"),
            (lambda: ExecutionSpec(threads=0), "threads"),
            (lambda: ExecutionSpec(buffer_size=0), "buffer_size"),
            (lambda: ExecutionSpec(annotation_window=0), "annotation_window"),
            (lambda: ExecutionSpec(max_in_flight=0), "max_in_flight"),
        ],
    )
    def test_direct_specification_rejects_unusable_limits(
        self,
        build_specification: Callable[[], ExecutionSpec],
        message: str,
    ) -> None:
        """The exported execution API cannot silently discard queued work."""

        with pytest.raises(ValueError, match=message):
            build_specification()


class TestExecutor:
    """Ordered execution across runtime backends."""

    def test_an_uninitialised_process_refuses_a_payload(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """A worker process that never ran its initialiser has no worker to execute against."""
        monkeypatch.setattr(executor_module, "_WORKER", None)

        with pytest.raises(RuntimeError, match="not initialised"):
            executor_module._submit_payload_to_worker(1)  # pylint: disable=protected-access

    def test_serial_execution_builds_one_worker_for_every_payload(self) -> None:
        """Serial execution runs the same worker protocol, against one worker it builds itself."""
        _Multiplier.built = 0

        results = list(_execute([1, 2, 3], ExecutionSpec()))

        assert results == [2, 4, 6]
        assert _Multiplier.built == 1

    def test_serial_execution_closes_input_after_early_close(self) -> None:
        """Closing results also closes an input iterator that owns resources."""
        input_closed = False

        def payloads() -> Generator[int, None, None]:
            """Yield payloads while exposing finalisation to the assertion."""
            nonlocal input_closed
            try:
                yield from [1, 2, 3]
            finally:
                input_closed = True

        results = _execute(payloads(), ExecutionSpec())

        assert next(results) == 2
        closable_results = cast(Generator[int, None, None], results)
        closable_results.close()

        assert input_closed

    def test_process_execution_preserves_input_order(self) -> None:
        """A parent consumes completed work in payload order."""
        results = _execute(
            [1, 2, 3], ExecutionSpec(processes=2, start_method=StartMethod.SPAWN, max_in_flight=2)
        )

        assert list(results) == [2, 4, 6]

    def test_process_execution_accepts_empty_input(self) -> None:
        """Discovering exhaustion before submission completes without looking up a future."""
        results = _execute([], ExecutionSpec(processes=2, start_method=StartMethod.SPAWN))

        assert not list(results)

    def test_process_execution_waits_after_normal_exhaustion(
        self,
        process_pool: type[_FakeProcessPoolExecutor],
    ) -> None:
        """Normal consumption waits for every submitted payload to finish."""
        results = _execute([1, 2, 3], ExecutionSpec(processes=2, max_in_flight=2))

        assert list(results) == [1, 2, 3]

        pool = process_pool.instances[0]
        assert pool.shutdown_calls == [(True, False)]
        assert not any(future.cancelled() for future in pool.submitted_futures)

    def test_process_execution_cancels_queued_work_after_early_close(
        self,
        process_pool: type[_FakeProcessPoolExecutor],
    ) -> None:
        """Closing the result iterator does not wait for its queued payloads."""
        results = _execute([1, 2, 3], ExecutionSpec(processes=2, max_in_flight=2))

        assert next(results) == 1
        closable_results = cast(Generator[int, None, None], results)
        closable_results.close()

        pool = process_pool.instances[0]
        assert len(pool.submitted_futures) == 2
        assert not pool.shutdown_calls
        assert pool.terminate_calls == 1
        assert pool.submitted_futures[1].cancel_requested

    def test_process_execution_terminates_workers_when_input_close_fails(
        self,
        process_pool: type[_FakeProcessPoolExecutor],
    ) -> None:
        """Input cleanup errors cannot leave the worker executor running."""

        def payloads() -> Generator[int, None, None]:
            """Raise while releasing the simulated input resource."""
            try:
                yield from [1, 2, 3]
            finally:
                raise RuntimeError("input close failed")

        results = _execute(payloads(), ExecutionSpec(processes=2, max_in_flight=2))

        assert next(results) == 1
        closable_results = cast(Generator[int, None, None], results)
        with pytest.raises(RuntimeError, match="input close failed"):
            closable_results.close()

        pool = process_pool.instances[0]
        assert pool.terminate_calls == 1

    @pytest.mark.skipif(
        "fork" not in multiprocessing.get_all_start_methods(),
        reason="the platform does not provide fork",
    )
    def test_fork_execution_preserves_input_order(self) -> None:
        """The Perl-compatible process backend follows the ordered worker contract."""
        results = _execute(
            [1, 2, 3], ExecutionSpec(processes=2, start_method=StartMethod.FORK, max_in_flight=1)
        )

        assert list(results) == [2, 4, 6]


class TestPayloadIterator:
    """Grouping a source's records into the payloads an executor runs."""

    def test_a_payload_names_one_sequence_region(self) -> None:
        """Close a payload at a region change, however much room it had left."""
        records = ["1:a", "1:b", "2:c"]

        payloads = PayloadIterator(records=records, batch_size=10, region_of=_region_of)

        assert list(payloads) == [("1:a", "1:b"), ("2:c",)]

    def test_a_payload_is_closed_once_it_is_full(self) -> None:
        """Close a payload at its size, even while the sequence region continues."""
        records = ["1:a", "1:b", "1:c"]

        payloads = PayloadIterator(records=records, batch_size=2, region_of=_region_of)

        assert list(payloads.payloads()) == [("1:a", "1:b"), ("1:c",)]

    def test_recurrent_regions_are_kept_apart(self) -> None:
        """Group region order 1, 2, 1 as it arrives, without sorting or combining it."""
        records = ["1:a", "2:b", "1:c"]

        payloads = PayloadIterator(records=records, batch_size=10, region_of=_region_of)

        assert list(payloads) == [("1:a",), ("2:b",), ("1:c",)]

    def test_records_that_name_no_region_are_grouped_by_size_alone(self) -> None:
        """Without a region to group by, a payload is closed only by its size."""
        payloads = PayloadIterator(records=["1:a", "2:b", "1:c"], batch_size=2)

        assert list(payloads) == [("1:a", "2:b"), ("1:c",)]

    def test_no_records_make_no_payload(self) -> None:
        """An empty input yields nothing rather than one empty payload."""
        assert not list(PayloadIterator(records=[], batch_size=10, region_of=_region_of))

    def test_a_non_positive_batch_size_is_refused(self) -> None:
        """Require each payload to have capacity for a record, before anything is read."""
        with pytest.raises(ValueError, match="at least one"):
            PayloadIterator(records=["1:a"], batch_size=0, region_of=_region_of)

    def test_leaving_the_payloads_early_closes_the_records(self) -> None:
        """Reading the records is what holds the input open, so this iterator closes them."""
        records_closed = False

        def records() -> Generator[str, None, None]:
            """Yield records while exposing finalisation to the assertion."""
            nonlocal records_closed
            try:
                yield from ["1:a", "1:b"]
            finally:
                records_closed = True

        payloads = cast(
            Generator[tuple[str, ...], None, None],
            PayloadIterator(records=records(), batch_size=1).payloads(),
        )

        assert next(payloads) == ("1:a",)
        payloads.close()

        assert records_closed
