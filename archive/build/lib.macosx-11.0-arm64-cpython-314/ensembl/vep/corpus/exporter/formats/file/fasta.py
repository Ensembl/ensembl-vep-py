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
"""FASTA rendering for corpora."""

from __future__ import annotations

from collections.abc import Iterator, Sequence

from ensembl.vep.corpus.base import Corpus

__all__ = ["FASTA_NAME", "fasta_lines"]

FASTA_NAME = "corpus.fa"

_LINE_WIDTH = 60


def fasta_lines(corpora: Sequence[Corpus]) -> Iterator[str]:
    """Yield one FASTA record per corpus chromosome, wrapped to a fixed width."""
    for corpus in corpora:
        for sequence_region in corpus.sequence_regions:
            sequence = sequence_region.seq or ""
            yield f">{sequence_region.seq_region_name}"
            for offset in range(0, len(sequence), _LINE_WIDTH):
                yield sequence[offset : offset + _LINE_WIDTH]
