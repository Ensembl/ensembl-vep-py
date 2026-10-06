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
"""Options for describing a synthetic corpus."""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

from attrs import define, field

from ensembl.vep.config.option_spec import option
from ensembl.vep.config.report import REPORT_FORMATS
from ensembl.vep.config.sections.base import ConfigSection

__all__ = ["CorpusDescribeConfig"]


@define(frozen=True, kw_only=True, auto_attribs=True)
class CorpusDescribeConfig(ConfigSection):
    """What to show about a corpus, and how.

    Attributes:
        select (str | None): List only the overlaps of the variant with this name, or None for
            every overlap.
        format (str): Table format.
        output_file (Path | None): Where to write, or None for stdout.
        annotate (bool): Also run the configured annotator over the corpus and count its consequence
            terms. Off by default, because it turns a listing into an annotation run.
    """

    section_name: ClassVar[str] = "corpus_describe"

    select: str | None = field(default=None, metadata=option("list only this variant's overlaps"))
    output_format: str | None = field(
        default=None,
        metadata=option("report format (default: markdown)", choices=REPORT_FORMATS),
    )
    output_file: Path | None = field(default=None, metadata=option("write output here"))
    annotate: bool = field(default=False, metadata=option("also annotate the corpus and count its terms"))
