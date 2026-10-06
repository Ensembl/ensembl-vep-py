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
"""Options controlling annotation-worker execution."""

from __future__ import annotations

from typing import ClassVar

from attrs import define, field

from ensembl.vep.common.enums import StartMethod
from ensembl.vep.config.option_spec import option
from ensembl.vep.config.sections.base import ConfigSection

__all__ = ["ExecutionConfig"]

# Mirrors ensembl.vep.caches.feature_pool.DEFAULT_TRANSCRIPT_POOL_BYTES (the source of truth
# every other call site imports). Not imported here so this config layer stays free of the
# domain-layer dependencies caches.feature_pool pulls in transitively (models.genic); kept in
# sync by test_transcript_pool_bytes_default_matches_the_caches_module.
_DEFAULT_TRANSCRIPT_POOL_BYTES = 128 * 1024 * 1024


@define(frozen=True, kw_only=True, auto_attribs=True)
class ExecutionConfig(ConfigSection):
    """Configure worker processes and optional independent-track threads.

    Attributes:
        processes (int): Number of worker processes.
        fork (int | None): Compatibility spelling for ``processes``.
        start_method (StartMethod): Multiprocessing start method, or ``auto``.
        threads (int): Threads available to independent annotation tracks in one worker.
        buffer_size (int): Maximum records in one sequence-region input buffer.
        annotation_window (int): Maximum variations retaining overlap graphs at one time.
        max_in_flight (int | None): Explicit submitted-buffer cap, or None to match worker count.
        transcript_pool_bytes (int): Byte budget for the per-worker cross-payload transcript pool.
            ``0`` disables pooling; every payload then rebuilds its own transcripts, matching
            behaviour before Stage 4.
    """

    section_name: ClassVar[str] = "execution"

    processes: int = field(default=1, metadata=option("worker process count"))
    fork: int | None = field(
        default=None,
        metadata=option("worker process count (Perl VEP compatibility)"),
    )
    start_method: StartMethod = field(
        default=StartMethod.AUTO,
        metadata=option("multiprocessing start method"),
    )
    threads: int = field(default=1, metadata=option("independent annotation-track thread count"))
    buffer_size: int = field(default=5000, metadata=option("records per input buffer"))
    annotation_window: int = field(
        default=1024,
        metadata=option("variations annotated together within an input buffer"),
    )
    max_in_flight: int | None = field(
        default=None,
        metadata=option("maximum submitted input buffers"),
    )
    transcript_pool_bytes: int = field(
        default=_DEFAULT_TRANSCRIPT_POOL_BYTES,
        metadata=option("transcript pool byte budget (0 disables)", metavar="BYTES"),
    )

    def __attrs_post_init__(self) -> None:
        """Reject execution settings that cannot produce a useful worker pool."""
        if self.processes < 1:
            raise ValueError("processes must be at least one")
        if self.fork is not None and self.fork < 1:
            raise ValueError("fork must be at least one")
        if self.fork is not None and self.processes != 1 and self.processes != self.fork:
            raise ValueError("processes and fork must agree when both are given")
        if self.threads < 1:
            raise ValueError("threads must be at least one")
        if self.buffer_size < 1:
            raise ValueError("buffer_size must be at least one")
        if self.annotation_window < 1:
            raise ValueError("annotation_window must be at least one")
        if self.max_in_flight is not None and self.max_in_flight < 1:
            raise ValueError("max_in_flight must be at least one")
        if self.transcript_pool_bytes < 0:
            raise ValueError("transcript_pool_bytes must be at least zero")

    @property
    def worker_processes(self) -> int:
        """Return the effective process count, applying the compatibility option."""
        return self.fork if self.fork is not None else self.processes
