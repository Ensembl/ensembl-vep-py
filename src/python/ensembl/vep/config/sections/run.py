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
"""Options for the default ``run`` path and its summary file."""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

from attrs import define, field

from ensembl.vep.config.option_spec import option
from ensembl.vep.config.sections.base import ConfigSection

__all__ = ["RunConfig"]


@define(frozen=True, kw_only=True, auto_attribs=True)
class RunConfig(ConfigSection):
    """Where annotation output and the run summary go.

    The summary is produced by the run itself and written unless ``no_stats`` is set, matching
    Perl VEP. It is cheap enough to default on: aggregating it costs roughly a microsecond per
    overlap, against a fact layer that is called about ninety times more often.

    Attributes:
        output_file (Path): Annotation output path.
        output_format (str): Annotation output format.
        force_overwrite (bool): Overwrite ``output_file`` if it already exists.
        no_progress (bool): Suppress interactive variant progress.
        no_stats (bool): Suppress the run summary entirely.
        stats_file (Path | None): Summary path, or None to derive ``<output_file>_summary.html``.
        stats_format (str): Summary format.
    """

    section_name: ClassVar[str] = "run"

    output_file: Path = field(
        factory=lambda: Path("variant_effect_output.txt"),
        metadata=option("annotation output file"),
    )
    output_format: str = field(default="vcf", metadata=option("registered annotation output format"))
    force_overwrite: bool = field(default=False, metadata=option("overwrite an existing output file"))
    no_progress: bool = field(default=False, metadata=option("do not show variant progress"))
    no_stats: bool = field(default=False, metadata=option("do not write a run summary"))
    stats_file: Path | None = field(default=None, metadata=option("run summary output file"))
    stats_format: str = field(
        default="html", metadata=option("run summary format", choices=("html", "text", "json"))
    )
