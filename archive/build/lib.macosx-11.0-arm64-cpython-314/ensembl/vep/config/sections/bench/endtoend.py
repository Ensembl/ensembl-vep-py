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
"""Options for timing whole runs of the command line against each other."""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

from attrs import define, field

from ensembl.vep.config.option_spec import option
from ensembl.vep.config.sections.base import ConfigSection

__all__ = ["BenchEndToEndConfig"]


@define(frozen=True, kw_only=True, auto_attribs=True)
class BenchEndToEndConfig(ConfigSection):
    """Which configurations to run end to end, where to run them, and what to leave behind.

    Attributes:
        strategies (tuple[str, ...]): Strategies to time, or empty for every registered
            ``consequence.*`` strategy. ``consequence.linear`` is run first, as the reference the
            speed-up column and the output gate are against.
        caches (tuple[str, ...]): Cache sets to cross every strategy with: ``none``, or a
            ``+``-joined set of ``zonemap``, ``codon`` and ``mapper``. The default matches
            ``bench matrix``, so a label means the same run in both tables.
        processes (tuple[int, ...]): Worker process counts to cross every configuration with. An
            axis only an end-to-end benchmark has: the in-process verbs never fork.
        batch_size (int): Overlaps per batch, passed to each run.
        work_dir (Path | None): Where the corpus, its indexes and the run outputs go. Defaults to
            ``corpus.corpus_dir`` when that is set, and to a temporary directory otherwise.
        clean (bool): Remove the files this run created when it finishes. Files that were already
            there are left alone either way.
    """

    section_name: ClassVar[str] = "bench_endtoend"

    strategies: tuple[str, ...] = field(
        default=(),
        metadata=option("strategies to time (default: every consequence strategy)", metavar="NAME,..."),
    )
    caches: tuple[str, ...] = field(
        default=("none", "zonemap", "zonemap+codon"),
        metadata=option("cache sets to cross each strategy with", metavar="SET,..."),
    )
    processes: tuple[int, ...] = field(
        default=(1,),
        metadata=option("worker process counts to time each configuration with", metavar="N,..."),
    )
    batch_size: int = field(default=256, metadata=option("overlaps per batch", metavar="N"))
    work_dir: Path | None = field(
        default=None, metadata=option("directory for the corpus, its indexes and the run outputs")
    )
    clean: bool = field(default=True, metadata=option("remove the files this run created"))

    def __attrs_post_init__(self) -> None:
        """Reject counts that cannot be run.

        Raises:
            ValueError: If no process count was given, or any of them is below one, or the batch
                size is below one.
        """
        if not self.processes:
            raise ValueError("processes must name at least one worker count")
        if any(count < 1 for count in self.processes):
            raise ValueError("every process count must be at least one")
        if self.batch_size < 1:
            raise ValueError("batch_size must be at least one")
