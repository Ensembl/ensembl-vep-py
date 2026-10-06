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
"""Synthetic corpus generation options."""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

from attrs import define, field

from ensembl.vep.config.option_spec import option
from ensembl.vep.config.sections.base import ConfigSection

__all__ = ["CorpusConfig"]

CORPUS_KINDS = ("default", "heavy", "both")


@define(frozen=True, kw_only=True, auto_attribs=True)
class CorpusConfig(ConfigSection):
    """Which synthetic corpus to build, and how.

    Consumed by ensembl.vep.corpus. Seeds are per stream rather than one global seed because the
    generator draws the chromosome sequence, the boundary variants and the fuzz variants from
    independent ``random.Random`` instances.

    Attributes:
        kind (str): ``"default"`` for the four-transcript boundary+fuzz corpus, ``"heavy"`` for the
            36-transcript exon-sharing locus, ``"both"`` for the two together. They exercise
            different things -- only the heavy locus shares exons between transcripts -- so a
            whole-catalogue check wants both.
        fuzz (int): Number of random variants to add to the default corpus.
        step (int): Position increment for the heavy-locus variant sweep.
        genome_seed (int): Seed for the default corpus's chromosome sequence and transcript layout.
        boundary_seed (int): Seed for the adversarial boundary variants.
        fuzz_seed (int): Seed for the random fuzz variants.
        heavy_genome_seed (int): Seed for the heavy locus's chromosome. Separate from
            ``genome_seed`` because the two corpora have always started from different seeds.
        sweep_seed (int): Seed for the heavy-locus variant sweep.
        corpus_dir (Path | None): Directory to load a saved corpus from, or write one to.
    """

    section_name: ClassVar[str] = "corpus"

    kind: str = field(default="both", metadata=option("synthetic corpus to use", choices=CORPUS_KINDS))
    fuzz: int = field(default=4000, metadata=option("fuzz variant count", metavar="N"))
    step: int = field(default=13, metadata=option("heavy-locus variant sweep step", metavar="N"))
    genome_seed: int = field(default=7, metadata=option("seed for the synthetic genome"))
    boundary_seed: int = field(default=11, metadata=option("seed for the boundary variants"))
    fuzz_seed: int = field(default=23, metadata=option("seed for the fuzz variants"))
    heavy_genome_seed: int = field(default=31, metadata=option("seed for the heavy-locus genome"))
    sweep_seed: int = field(default=41, metadata=option("seed for the heavy-locus variant sweep"))
    corpus_dir: Path | None = field(default=None, metadata=option("directory to load or save a corpus"))
