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

_STDOUT_SPELLINGS = frozenset({"stdout", "-"})


def _as_output_file(value: Path | str | None) -> Path | None:
    """Map the stdout spellings, and an empty value, to ``None``; read anything else as a path."""
    if value is None:
        return None
    text = str(value)
    if text.lower() in _STDOUT_SPELLINGS:
        return None
    return Path(text)


@define(frozen=True, kw_only=True, auto_attribs=True)
class RunConfig(ConfigSection):
    """Where annotation output and the run summary go.

    The summary is produced by the run itself and written unless ``no_stats`` is set, matching
    Perl VEP. It is cheap enough to default on: aggregating it costs roughly a microsecond per
    overlap, against a fact layer that is called about ninety times more often.

    Attributes:
        output_file (Path | None): Annotation output path, or None to stream to stdout. ``STDOUT``
            (any case) and ``-`` mean stdout, as they do in Perl VEP.
        output_format (str): Annotation output format. Defaults to ``vcf``; Perl VEP defaults to ``vep``,
            which has no writer yet.
        force_overwrite (bool): Overwrite ``output_file`` if it already exists.
        no_progress (bool): Suppress interactive variant progress.
        no_stats (bool): Suppress the run summary entirely.
        stats_file (Path | None): Summary path, or None to derive ``<output_file>_summary.html``
            (no summary file when output goes to stdout).
        stats_format (str): Summary format.
    """

    section_name: ClassVar[str] = "run"

    output_file: Path | None = field(
        default=Path("variant_effect_output.txt"),
        converter=_as_output_file,
        metadata=option("annotation output file", short="o"),
    )
    output_format: str = field(
        default="vcf",
        metadata=option(
            "annotation output format (Perl VEP defaults to its own 'vep' format, "
            "which has no writer here yet)",
            constants={"vcf": "vcf", "json": "json"},
        ),
    )
    force_overwrite: bool = field(
        default=False,
        metadata=option("overwrite an existing output file", synonyms=("force",)),
    )
    no_progress: bool = field(default=False, metadata=option("do not show variant progress"))
    no_stats: bool = field(default=False, metadata=option("do not write a run summary"))
    stats_file: Path | None = field(
        default=None, metadata=option("run summary output file", synonyms=("sf",))
    )
    stats_format: str = field(
        default="html",
        metadata=option(
            "run summary format",
            choices=("html", "text", "json"),
            constants={"stats_html": "html", "stats_text": "text"},
        ),
    )
