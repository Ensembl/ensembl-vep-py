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
"""Compatibility exports for corpus file saving."""

from __future__ import annotations

from pathlib import Path

from ensembl.vep.corpus.base import CorpusLike, as_corpus
from ensembl.vep.corpus.exporter.formats.file.base import (
    CORPUS_FILES,
    MANIFEST_NAME,
    distinct_variants,
)
from ensembl.vep.corpus.exporter.formats.file.fasta import FASTA_NAME, fasta_lines
from ensembl.vep.corpus.exporter.formats.file.gff3 import GFF3_NAME, gff3_lines
from ensembl.vep.corpus.exporter.formats.file.vcf import VCF_NAME, vcf_alleles, vcf_lines
from ensembl.vep.corpus.exporter.manifest import CorpusManifest
from ensembl.vep.corpus.exporter.registry import get_exporter

__all__ = [
    "CORPUS_FILES",
    "CorpusManifest",
    "FASTA_NAME",
    "GFF3_NAME",
    "MANIFEST_NAME",
    "VCF_NAME",
    "distinct_variants",
    "fasta_lines",
    "gff3_lines",
    "save_corpus",
    "vcf_alleles",
    "vcf_lines",
]


def save_corpus(corpus: CorpusLike, directory: Path, *, force: bool = False) -> CorpusManifest:
    """Write a corpus as FASTA, GFF3 and VCF, with a manifest describing it."""
    return get_exporter("file").export(as_corpus(corpus), directory, force=force)
