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
"""Options for timing strategies and cache configurations against each other."""

from __future__ import annotations

from typing import ClassVar

from attrs import define, field

from ensembl.vep.config.option_spec import option
from ensembl.vep.config.sections.base import ConfigSection

__all__ = ["BenchMatrixConfig"]


@define(frozen=True, kw_only=True, auto_attribs=True)
class BenchMatrixConfig(ConfigSection):
    """Which strategies, crossed with which cache sets.

    Attributes:
        strategies (tuple[str, ...]): Strategies to time, or empty for every registered
            ``consequence.*`` strategy. ``consequence.linear`` is always timed first, as the
            reference the speed-up column is against.
        caches (tuple[str, ...]): Cache sets to cross every strategy with: ``none``, or a
            ``+``-joined set of ``zonemap``, ``codon`` and ``mapper``. Each configuration gets
            fresh caches, so none is timed against another's warmed ones.
        batch_size (int): Chunk size for the batch timing.
    """

    section_name: ClassVar[str] = "bench_matrix"

    strategies: tuple[str, ...] = field(
        default=(),
        metadata=option("strategies to time (default: every consequence strategy)", metavar="NAME,..."),
    )
    caches: tuple[str, ...] = field(
        default=("none", "zonemap", "zonemap+codon"),
        metadata=option("cache sets to cross each strategy with", metavar="SET,..."),
    )
    batch_size: int = field(default=256, metadata=option("chunk size for the batch timing", metavar="N"))
