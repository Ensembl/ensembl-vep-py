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
"""Options for sweeping one numeric annotation setting."""

from __future__ import annotations

from typing import ClassVar

from attrs import define, field

from ensembl.vep.config.option_spec import option
from ensembl.vep.config.sections.base import ConfigSection

__all__ = ["SWEEPABLE_KNOBS", "BenchSweepConfig"]

SWEEPABLE_KNOBS = ("batch_size", "codon_cache_entries", "mapper_cache_bytes", "zone_map_bytes")
"""The numeric annotation settings a sweep can vary."""


@define(frozen=True, kw_only=True, auto_attribs=True)
class BenchSweepConfig(ConfigSection):
    """Which setting to sweep, over which values, for which strategies.

    Attributes:
        strategies (tuple[str, ...]): Strategies to sweep, or empty for the batch strategies
            available here -- ``consequence.columnar-staged``, and ``consequence.columnar-kernels``
            when the native extension is built.
        sweep_knob (str): The setting to vary.
        sweep_values (tuple[int, ...]): The values to try, in order.
    """

    section_name: ClassVar[str] = "bench_sweep"

    strategies: tuple[str, ...] = field(
        default=(), metadata=option("strategies to sweep (default: the batch strategies)", metavar="NAME,...")
    )
    sweep_knob: str = field(default="batch_size", metadata=option("setting to vary", choices=SWEEPABLE_KNOBS))
    sweep_values: tuple[int, ...] = field(
        default=(32, 64, 128, 256, 512, 1024), metadata=option("values to try", metavar="N,N,...")
    )
