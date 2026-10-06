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
"""Permission to change the catalogue on disk."""

from __future__ import annotations

from typing import ClassVar

from attrs import define, field

from ensembl.vep.config.option_spec import option
from ensembl.vep.config.sections.base import ConfigSection

__all__ = ["CatalogueWriteConfig"]


@define(frozen=True, kw_only=True, auto_attribs=True)
class CatalogueWriteConfig(ConfigSection):
    """Read only by the verbs that change the catalogue, so the listing verbs never offer it.

    Attributes:
        catalogue_overwrite (bool): Permit writing the packaged catalogue in place. Without this,
            or an explicit path, a write is refused: an installed ``vep`` must not silently rewrite
            its own wheel.
        dry_run (bool): Report the change and its validation without writing.
    """

    section_name: ClassVar[str] = "catalogue_write"

    catalogue_overwrite: bool = field(default=False, metadata=option("write the packaged catalogue in place"))
    dry_run: bool = field(default=False, metadata=option("report the change without writing it"))
