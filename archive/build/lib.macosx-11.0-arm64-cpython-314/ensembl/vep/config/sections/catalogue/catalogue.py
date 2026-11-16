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
"""Which catalogue to act on, and whether to write it."""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

from attrs import define, field

from ensembl.vep.config.option_spec import option
from ensembl.vep.config.sections.base import ConfigSection

__all__ = ["CatalogueConfig"]


@define(frozen=True, kw_only=True, auto_attribs=True)
class CatalogueConfig(ConfigSection):
    """Which catalogue to act on. Every verb reads this one.

    Attributes:
        catalogue_file (Path | None): Catalogue JSON to read or write, or None for the packaged one.
    """

    section_name: ClassVar[str] = "catalogue"

    catalogue_file: Path | None = field(
        default=None, metadata=option("catalogue JSON to read or write", flag="catalogue")
    )
