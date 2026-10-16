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
"""Statistics and profiling options."""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

from attrs import define, field

from ensembl.vep.config.option_spec import option
from ensembl.vep.config.sections.base import ConfigSection

__all__ = ["StatsConfig"]


@define(frozen=True, kw_only=True, auto_attribs=True)
class StatsConfig(ConfigSection):
    """What ``vep stats`` reads and what it writes.

    Attributes:
        output_file (Path | None): Report path, or None for stdout.
        output_format (str): Report format. Distinct from ``--format``, which is the input
            format of the file it summarises.
        profile (bool): Also collect the fact-level profile. This one genuinely costs -- it wraps
            every fact's compute function -- so it is never on by default.
    """

    section_name: ClassVar[str] = "stats"

    output_file: Path | None = field(default=None, metadata=option("report output file"))
    output_format: str = field(
        default="html", metadata=option("report format", choices=("html", "text", "json"))
    )
    profile: bool = field(default=False, metadata=option("collect fact-level profiling (slower)"))
