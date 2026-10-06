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
"""Options every benchmark verb shares."""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

from attrs import define, field

from ensembl.vep.config.option_spec import option
from ensembl.vep.config.report import REPORT_FORMATS
from ensembl.vep.config.sections.base import ConfigSection

__all__ = ["BenchConfig"]


@define(frozen=True, kw_only=True, auto_attribs=True)
class BenchConfig(ConfigSection):
    """How every benchmark is measured and reported. Read by every ``bench`` verb.

    Attributes:
        repeats (int): Timed samples per configuration; the best is reported.
        warmup (int): Untimed passes before the first sample.
        no_gate (bool): Skip the conformance gate, which otherwise refuses to time anything until
            every strategy agrees with the reference. For re-timing a tree already known to conform.
        results_file (Path | None): Write the measured results here as JSON.
        baseline_file (Path | None): Compare against results previously written here.
        regression_threshold (float): Fractional slowdown against the baseline that counts as a
            regression.
        format (str): Table format.
    """

    section_name: ClassVar[str] = "bench"

    repeats: int = field(default=3, metadata=option("timed samples per configuration", metavar="N"))
    warmup: int = field(default=1, metadata=option("untimed passes before timing", metavar="N"))
    no_gate: bool = field(default=False, metadata=option("skip the conformance gate before timing"))
    results_file: Path | None = field(default=None, metadata=option("write results JSON here"))
    baseline_file: Path | None = field(default=None, metadata=option("compare against this results JSON"))
    regression_threshold: float = field(
        default=0.05, metadata=option("fractional slowdown counted as a regression")
    )
    output_format: str | None = field(
        default=None,
        metadata=option("report format (default: text)", choices=REPORT_FORMATS),
    )
