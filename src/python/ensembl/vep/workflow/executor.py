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
"""Payload execution: who runs the work, and how the work is grouped into payloads.

Three things, one per step of a run's spine. ``ExecutionSpec`` is what the configuration asked
for. ``PayloadIterator`` groups a source's records into the payloads one unit of work holds.
``Executor`` runs them against a ``Worker``, serially or across processes, and yields the results
in payload order whatever order they finish in.
"""

from __future__ import annotations

import multiprocessing
from collections.abc import Callable, Iterable, Iterator
from concurrent.futures import Future, ProcessPoolExecutor
from typing import Any, Generic, TypeVar, cast

from attrs import define, field

from ensembl.vep.common.enums import StartMethod
from ensembl.vep.config.sections.common.execution import ExecutionConfig
from ensembl.vep.workflow.worker.base import Worker

__all__ = ["ExecutionSpec", "Executor", "PayloadIterator"]

PayloadT = TypeVar("PayloadT")
RecordT = TypeVar("RecordT")
ResultT = TypeVar("ResultT")


@define(frozen=True, kw_only=True, auto_attribs=True)
class ExecutionSpec:
    """Worker-local annotation execution settings.

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
class PayloadIterator(Generic[RecordT]):
    """One source's records, grouped into the payloads an executor runs.

    Here rather than in the source being read: a source says what a record is and which region it
    sits on, and how many records make one unit of work is an execution setting (``buffer_size``).

    Attributes:
        records (Iterable[RecordT]): The records to group, in input order. Closed with this
            iterator, since reading them is what holds the input open.
        batch_size (int): The most records one payload may hold.
        region_of (Callable[[RecordT], str] | None): Name the sequence region a record sits on.
            A payload then never spans two regions, which is what a worker's annotation tracks
            prime from -- one region's features, fetched once for the whole payload. None groups
            by size alone, for records that name no region.
    """

    records: Iterable[RecordT]
    batch_size: int
    region_of: Callable[[RecordT], str] | None = field(default=None)

    def __attrs_post_init__(self) -> None:
        """Reject a size no payload could hold a record at.

        Raises:
            ValueError: If ``batch_size`` is below one.
        """

        if self.batch_size < 1:
            raise ValueError("batch_size must be at least one")

    def __iter__(self) -> Iterator[tuple[RecordT, ...]]:
        """Return this iterator's payloads, so an executor can be given it as it is."""

        return self.payloads()

    def payloads(self) -> Iterator[tuple[RecordT, ...]]:
        """Yield payloads of at most ``batch_size`` records, in input order.

        Yields:
            Iterator[tuple[RecordT, ...]]: Each payload, closed by its size or by a region change.
        """

        region: str | None = None
        payload: list[RecordT] = []
        try:
            for record in self.records:
                record_region = self.region_of(record) if self.region_of is not None else None
                if payload and (record_region != region or len(payload) == self.batch_size):
                    yield tuple(payload)
                    payload = []
                region = record_region
                payload.append(record)
            if payload:
                yield tuple(payload)
        finally:
            _stop_payload_iterator(self.records)


@define(kw_only=True, auto_attribs=True)
class Executor(Generic[PayloadT, ResultT]):
    """Ordered payload executor with one worker per process.

    A worker is built after its process starts and then reused, so what the executor transports is
    the factory rather than the worker: a worker holds open sources and warm caches that neither
    cross a process boundary nor survive a rebuild. The serial path runs the same way, against a
    worker of its own in the calling process, so there is one execution contract rather than two.

    Attributes:
        specification (ExecutionSpec): Process and queue configuration.
        worker_factory (Callable[[], Worker[PayloadT, ResultT]]): Construct one worker. It is what
            crosses to a worker process, so it has to pickle; the worker it builds does not.
        payload_iterator (Iterable[PayloadT]): The payloads to run, in input order, each
            serialisable whenever ``processes`` exceeds one. A run passes a ``PayloadIterator``.
    """

    specification: ExecutionSpec
    worker_factory: Callable[[], Worker[PayloadT, ResultT]]
    payload_iterator: Iterable[PayloadT]

    def execute(self) -> Iterator[ResultT]:
        """Yield every payload's result in payload order, while limiting queued work.

        Yields:
            Iterator[ResultT]: Results in the order ``payload_iterator`` produced the payloads.
        """
        payloads = iter(self.payload_iterator)
        if self.specification.processes == 1:
            yield from self._execute_serial(payloads)
        else:
            yield from self._execute_parallel(payloads)

    def _execute_serial(self, payloads: Iterator[PayloadT]) -> Iterator[ResultT]:
        """Yield every result from one worker in the calling process."""
        _start_worker(self.worker_factory)
        try:
            for payload in payloads:
                yield cast(ResultT, _submit_payload_to_worker(payload))
        finally:
            _stop_payload_iterator(payloads)

    def _execute_parallel(self, payloads: Iterator[PayloadT]) -> Iterator[ResultT]:
        """Yield results in payload order from a pool of worker processes."""
        pool = ProcessPoolExecutor(
            max_workers=self.specification.processes,
            mp_context=self.specification.multiprocessing_context(),
            initializer=_start_worker,
            initargs=(self.worker_factory,),
        )
        pending: dict[int, Future[object]] = {}
        next_submit = 0
        next_yield = 0
        completed = False

        try:
            exhausted = False
            while not exhausted or pending:
                while not exhausted and len(pending) < self.specification.effective_max_in_flight:
                    try:
                        payload = next(payloads)
                    except StopIteration:
                        exhausted = True
                        break
                    pending[next_submit] = pool.submit(_submit_payload_to_worker, payload)
                    next_submit += 1

                if not pending:
                    break
                future = pending.pop(next_yield)
                yield cast(ResultT, future.result())
                next_yield += 1
            completed = True
        finally:
            try:
                _stop_payload_iterator(payloads)
            finally:
                _shut_down(pool, pending, completed=completed)


_WORKER: Worker[Any, Any] | None = None
"""This process's worker, built once before it is given any payload.

Per process rather than per executor, because that is what it is: a worker process is initialised
once and then runs whatever arrives. The serial path builds it the same way, in the calling
process, which is why one run at a time per process is the contract."""


def _start_worker(worker_factory: Callable[[], Worker[Any, Any]]) -> None:
    """Build this process's worker, replacing one an earlier run left behind."""
    global _WORKER  # pylint: disable=global-statement
    _WORKER = worker_factory()


def _submit_payload_to_worker(payload: object) -> object:
    """Run one payload against this process's worker.

    Raises:
        RuntimeError: If the process was never given a worker to run against.
    """
    if _WORKER is None:
        raise RuntimeError("worker process was not initialised")
    return _WORKER.run(payload)


def _shut_down(
    pool: ProcessPoolExecutor,
    pending: dict[int, Future[object]],
    *,
    completed: bool,
) -> None:
    """Wait for a finished pool, or drop the queued work of an abandoned one."""
    if completed:
        pool.shutdown(wait=True)
        return
    for future in pending.values():
        future.cancel()
    terminate_workers = getattr(pool, "terminate_workers", None)
    if terminate_workers is not None:
        terminate_workers()
    else:
        pool.shutdown(wait=False, cancel_futures=True)


def _stop_payload_iterator(iterator: Iterable[object]) -> None:
    """Close a payload or record iterator that owns input resources."""
    close = getattr(iterator, "close", None)
    if close is not None:
        close()
