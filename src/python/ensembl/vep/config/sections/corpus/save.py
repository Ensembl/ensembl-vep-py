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
"""Options for writing a synthetic corpus to disk."""

from __future__ import annotations

from typing import ClassVar

from attrs import define, field

from ensembl.vep.config.option_spec import option
from ensembl.vep.config.sections.base import ConfigSection

__all__ = ["CorpusSaveConfig"]


@define(frozen=True, kw_only=True, auto_attribs=True)
class CorpusSaveConfig(ConfigSection):
    """Permission to replace a corpus already on disk. The destination is ``corpus.corpus_dir``.

    Attributes:
        force (bool): Overwrite the corpus files if the directory already holds them.
    """

    section_name: ClassVar[str] = "corpus_save"

    force: bool = field(default=False, metadata=option("overwrite an existing saved corpus"))
