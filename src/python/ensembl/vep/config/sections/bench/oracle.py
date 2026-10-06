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
"""Options for comparing against the Perl ``VariationEffect.pm`` oracle."""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

from attrs import define, field

from ensembl.vep.config.option_spec import option
from ensembl.vep.config.sections.base import ConfigSection

__all__ = ["BenchOracleConfig"]


@define(frozen=True, kw_only=True, auto_attribs=True)
class BenchOracleConfig(ConfigSection):
    """Where the Perl oracle lives, and what to compare against it.

    Attributes:
        perl_modules (Path | None): ensembl-variation's ``modules`` directory. Without it the
            Python strategies are still compared against each other, and the Perl side is skipped
            with the reason.
        perl_oracle_script (Path): The oracle script.
        strategies (tuple[str, ...]): Strategies to compare, or empty for ``consequence.linear``
            and, when built, ``consequence.columnar-kernels``.
        oracle_number (int): Workload executions per timed repeat.
    """

    section_name: ClassVar[str] = "bench_oracle"

    perl_modules: Path | None = field(default=None, metadata=option("ensembl-variation modules directory"))
    perl_oracle_script: Path = field(
        default=Path("bench/variation_effect.pl"), metadata=option("the Perl oracle script")
    )
    strategies: tuple[str, ...] = field(
        default=(),
        metadata=option("strategies to compare (default: linear and the kernels)", metavar="NAME,..."),
    )
    oracle_number: int = field(
        default=200, metadata=option("workload executions per timed repeat", metavar="N")
    )
