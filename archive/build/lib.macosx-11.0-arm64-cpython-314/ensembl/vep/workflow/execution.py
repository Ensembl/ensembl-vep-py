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
"""Execute ordered workflow payloads in serial or worker processes."""

from __future__ import annotations

import multiprocessing
from collections.abc import Callable, Iterable, Iterator
from concurrent.futures import Future, ProcessPoolExecutor
from typing import Generic, TypeVar, cast

from attrs import define

from ensembl.vep.common.enums import StartMethod
from ensembl.vep.config.sections.common.execution import ExecutionConfig

__all__ = ["ExecutionSpec", "ordered_map"]

PayloadT = TypeVar("PayloadT")
ResultT = TypeVar("ResultT")
WorkerT = TypeVar("WorkerT")


@define(frozen=True, kw_only=True, auto_attribs=True)
class ExecutionSpec:
    """Describe how worker-local annotation runtimes are executed.

    The payload and result types are deliberately absent from this value object. They must be
    serialisable transport values whenever ``processes`` exceeds one; live model objects remain
    inside a worker.

    Attributes:
        processes (int): Number of independent worker processes.
        start_method (StartMethod): Multiprocessing context, or ``auto`` for the platform default.
        threads (int): Threads available to independent tracks in each worker.
        buffer_size (int): Maximum records in one sequence-region input buffer.
        annotation_window (int): Maximum variations retaining overlap graphs at one time.
        max_in_flight (int | None): Explicit submitted-payload cap, or None to match workers.
    """

    processes: int = 1
    start_method: StartMethod = StartMethod.AUTO
    threads: int = 1
    buffer_size: int = 5000
    annotation_window: int = 1024
    max_in_flight: int | None = None

    def __attrs_post_init__(self) -> None:
        """Reject values that could drop work or create an unusable executor."""

        if self.processes < 1:
            raise ValueError("processes must be at least one")
        if self.threads < 1:
            raise ValueError("threads must be at least one")
        if self.buffer_size < 1:
            raise ValueError("buffer_size must be at least one")
        if self.annotation_window < 1:
            raise ValueError("annotation_window must be at least one")
        if self.max_in_flight is not None and self.max_in_flight < 1:
            raise ValueError("max_in_flight must be at least one")

    @classmethod
    def from_config(cls, config: ExecutionConfig) -> ExecutionSpec:
        """Build a specification from the declared execution options."""
        return cls(
            processes=config.processes,
            start_method=config.start_method,
            threads=config.threads,
            buffer_size=config.buffer_size,
            annotation_window=config.annotation_window,
            max_in_flight=config.max_in_flight,
        )

    @property
    def effective_max_in_flight(self) -> int:
        """Return the explicit queue cap, defaulting to one payload per worker."""
        if self.max_in_flight is not None:
            return self.max_in_flight
        return self.processes

    def multiprocessing_context(self) -> multiprocessing.context.BaseContext | None:
        """Return the multiprocessing context to run workers in.

        ``auto`` means fork wherever a forked worker is safe, which is what keeps the parallel path
        fast: a forked worker inherits the loaded package, its catalogue and its caches, while a
        spawned one re-imports all of it. Measured on 10,000 ClinVar variants across four workers,
        spawning cost about 30% of wall time. Threads are the exception -- forking a process that
        already has them is unsafe -- so a threaded run takes the platform default instead.
        """
        if self.start_method is not StartMethod.AUTO:
            return multiprocessing.get_context(self.start_method.value)
        if self.processes == 1 or self.threads > 1:
            return None
        if StartMethod.FORK.value not in multiprocessing.get_all_start_methods():
            return None
        return multiprocessing.get_context(StartMethod.FORK.value)


@define(kw_only=True, auto_attribs=True)
class _WorkerState(Generic[WorkerT, PayloadT, ResultT]):
    """Hold a worker runtime constructed after a process starts."""

    factory: Callable[[], WorkerT]
    execute: Callable[[WorkerT, PayloadT], ResultT]
    runtime: WorkerT | None = None

    def initialise(self) -> None:
        """Construct the process-local runtime once."""
        self.runtime = self.factory()

    def run(self, payload: PayloadT) -> ResultT:
        """Run one primitive payload through the local runtime."""
        if self.runtime is None:
            self.initialise()
        assert self.runtime is not None
        return self.execute(self.runtime, payload)


_WORKER_STATE: _WorkerState[object, object, object] | None = None


def _initialise_process_worker(state: _WorkerState[object, object, object]) -> None:
    """Install and initialise a process-local worker state."""
    global _WORKER_STATE  # pylint: disable=global-statement
    _WORKER_STATE = state
    _WORKER_STATE.initialise()


def _run_process_payload(payload: object) -> object:
    """Execute one payload in the process-local worker runtime."""
    if _WORKER_STATE is None:
        raise RuntimeError("worker process was not initialised")
    return _WORKER_STATE.run(payload)


def ordered_map(
    payloads: Iterable[PayloadT],
    *,
    specification: ExecutionSpec,
    worker_factory: Callable[[], WorkerT],
    execute: Callable[[WorkerT, PayloadT], ResultT],
) -> Iterator[ResultT]:
    """Yield results in payload order while limiting queued work.

    Args:
        payloads (Iterable[PayloadT]): Pickle-safe work descriptions in input order.
        specification (ExecutionSpec): Process and queue configuration.
        worker_factory (Callable[[], WorkerT]): Construct one worker-local runtime.
        execute (Callable[[WorkerT, PayloadT], ResultT]): Run one payload against that runtime.

    Yields:
        Iterator[ResultT]: Results in the same order as ``payloads``.
    """
    state: _WorkerState[WorkerT, PayloadT, ResultT] = _WorkerState(
        factory=worker_factory,
        execute=execute,
    )
    payload_iterator = iter(payloads)
    if specification.processes == 1:
        state.initialise()
        try:
            for payload in payload_iterator:
                yield state.run(payload)
        finally:
            _close_iterator(payload_iterator)
        return

    context = specification.multiprocessing_context()
    process_state = cast(_WorkerState[object, object, object], state)
    executor = ProcessPoolExecutor(
        max_workers=specification.processes,
        mp_context=context,
        initializer=_initialise_process_worker,
        initargs=(process_state,),
    )
    pending: dict[int, Future[object]] = {}
    next_submit = 0
    next_yield = 0
    completed = False

    try:
        exhausted = False
        while not exhausted or pending:
            while not exhausted and len(pending) < specification.effective_max_in_flight:
                try:
                    payload = next(payload_iterator)
                except StopIteration:
                    exhausted = True
                    break
                pending[next_submit] = executor.submit(_run_process_payload, payload)
                next_submit += 1

            if not pending:
                break
            future = pending.pop(next_yield)
            yield cast(ResultT, future.result())
            next_yield += 1
        completed = True
    finally:
        try:
            _close_iterator(payload_iterator)
        finally:
            if completed:
                executor.shutdown(wait=True)
            else:
                for future in pending.values():
                    future.cancel()
                terminate_workers = getattr(executor, "terminate_workers", None)
                if terminate_workers is not None:
                    terminate_workers()
                else:
                    executor.shutdown(wait=False, cancel_futures=True)


def _close_iterator(iterator: Iterator[object]) -> None:
    """Close a payload iterator that owns input resources."""
    close = getattr(iterator, "close", None)
    if close is not None:
        close()
