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
"""The manifest written beside an exported corpus."""

from __future__ import annotations

import json

from attrs import asdict, frozen

__all__ = ["CorpusManifest"]


@frozen(kw_only=True, auto_attribs=True)
class CorpusManifest:
    """Describe an exported corpus.

    Attributes:
        kind (str): Which corpus kind was written.
        spec (dict[str, object]): Complete registered construction specification.
        fuzz (int | None): Compatibility fuzz count for synthetic corpora.
        step (int | None): Compatibility sweep step for synthetic corpora.
        seeds (dict[str, dict[str, int]]): The seeds each written corpus was built from, per kind.
        parameters (str): The spec's parameter hash.
        chromosomes (tuple[str, ...]): Sequence region names written.
        transcripts (int): Transcripts written.
        variants (int): Distinct variants written.
        overlaps (int): Overlaps the corpus pairs those variants into.
        identity (str): A hash of the corpus files.
        files (tuple[str, ...]): The file names written beside this manifest.
    """

    kind: str
    spec: dict[str, object]
    fuzz: int | None
    step: int | None
    seeds: dict[str, dict[str, int]]
    parameters: str
    chromosomes: tuple[str, ...]
    transcripts: int
    variants: int
    overlaps: int
    identity: str
    files: tuple[str, ...]

    def to_json(self) -> str:
        """Return the manifest as indented JSON."""
        return json.dumps(asdict(self), indent=2) + "\n"
