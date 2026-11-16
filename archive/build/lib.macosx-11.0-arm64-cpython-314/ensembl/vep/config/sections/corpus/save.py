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
"""Options for exporting a registered corpus to disk."""

from __future__ import annotations

from typing import ClassVar

from attrs import define, field

from ensembl.vep.config.option_spec import option
from ensembl.vep.config.sections.base import ConfigSection
from ensembl.vep.corpus.exporter import available_exporters

__all__ = ["CorpusSaveConfig"]


@define(frozen=True, kw_only=True, auto_attribs=True)
class CorpusSaveConfig(ConfigSection):
    """How to export a corpus. The destination is ``corpus.corpus_dir``.

    Attributes:
        exporter (str): Registered corpus exporter to use.
        force (bool): Overwrite the corpus files if the directory already holds them.
    """

    section_name: ClassVar[str] = "corpus_save"

    exporter: str = field(
        default="file", metadata=option("registered corpus exporter", choices=available_exporters)
    )
    force: bool = field(default=False, metadata=option("overwrite an existing saved corpus"))
