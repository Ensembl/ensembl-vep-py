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
"""Options for measuring fact costs."""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

from attrs import define, field

from ensembl.vep.config.option_spec import option
from ensembl.vep.config.report import REPORT_FORMATS
from ensembl.vep.config.sections.base import ConfigSection

__all__ = ["CatalogueCostsConfig"]


@define(frozen=True, kw_only=True, auto_attribs=True)
class CatalogueCostsConfig(ConfigSection):
    """How the measured costs are reported and persisted.

    Attributes:
        format (str): Table format for the comparison.
        output_file (Path | None): Write the measured cost table here, or None to only report it.
    """

    section_name: ClassVar[str] = "catalogue_costs"

    output_format: str | None = field(
        default=None,
        metadata=option("report format (default: text)", choices=REPORT_FORMATS),
    )
    output_file: Path | None = field(default=None, metadata=option("write the measured cost table here"))
