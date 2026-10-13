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
"""The corpus exporter contract."""

from __future__ import annotations

from pathlib import Path
from typing import Protocol

from ensembl.vep.corpus.base import Corpus
from ensembl.vep.corpus.exporter.manifest import CorpusManifest

__all__ = ["CorpusExporter"]


class CorpusExporter(Protocol):
    """Write a corpus in one output format."""

    def export(self, corpus: Corpus, directory: Path, *, force: bool = False) -> CorpusManifest:
        """Write a corpus and return its manifest.

        Args:
            corpus (Corpus): Corpus to export.
            directory (Path): Destination directory.
            force (bool): Whether existing format files may be replaced.

        Returns:
            CorpusManifest: Description and identity of the exported representation.
        """
