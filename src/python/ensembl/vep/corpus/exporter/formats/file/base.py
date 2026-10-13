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
"""A coordinated FASTA, GFF3 and VCF corpus exporter."""

from __future__ import annotations

import hashlib
from pathlib import Path
from tempfile import TemporaryDirectory

from ensembl.vep.corpus.base import Corpus
from ensembl.vep.corpus.exporter.formats.file.fasta import FASTA_NAME, fasta_lines
from ensembl.vep.corpus.exporter.formats.file.gff3 import GFF3_NAME, gff3_lines
from ensembl.vep.corpus.exporter.formats.file.vcf import VCF_NAME, distinct_variants, vcf_lines
from ensembl.vep.corpus.exporter.manifest import CorpusManifest
from ensembl.vep.corpus.exporter.registry import register_exporter

__all__ = [
    "CORPUS_FILES",
    "FileCorpusExporter",
    "MANIFEST_NAME",
    "distinct_variants",
]

MANIFEST_NAME = "manifest.json"
CORPUS_FILES = (FASTA_NAME, GFF3_NAME, VCF_NAME, MANIFEST_NAME)


def _optional_integer(value: object) -> int | None:
    """Return an optional integer manifest value."""
    return value if isinstance(value, int) and not isinstance(value, bool) else None


def _identity(contents: dict[str, str]) -> str:
    """Return a short digest over corpus files, in a fixed order."""
    digest = hashlib.sha256()
    for file_name in sorted(contents):
        digest.update(file_name.encode("utf-8"))
        digest.update(contents[file_name].encode("utf-8"))
    return digest.hexdigest()[:16]


@register_exporter("file")
class FileCorpusExporter:
    """Write a corpus as one coordinated file bundle."""

    def export(self, corpus: Corpus, directory: Path, *, force: bool = False) -> CorpusManifest:
        """Write a corpus as FASTA, GFF3 and VCF, with a manifest describing it.

        Args:
            corpus (Corpus): Corpus to render.
            directory (Path): Directory that receives the coordinated file bundle.
            force (bool): Whether existing corpus files may be replaced.

        Returns:
            CorpusManifest: Counts, parameters, and content identity for the written files.

        Raises:
            FileExistsError: If a corpus file exists and ``force`` was not given.
        """
        directory.mkdir(parents=True, exist_ok=True)
        present = [name for name in CORPUS_FILES if (directory / name).exists()]
        if present and not force:
            raise FileExistsError(
                f"{directory} already holds {', '.join(present)}; pass --force to overwrite"
            )

        corpora = (corpus,)
        contents = {
            FASTA_NAME: "\n".join(fasta_lines(corpora)) + "\n",
            GFF3_NAME: "\n".join(gff3_lines(corpora)) + "\n",
            VCF_NAME: "\n".join(vcf_lines(corpora)) + "\n",
        }
        spec = corpus.spec
        manifest = CorpusManifest(
            kind=spec.kind,
            spec=spec.as_dict(),
            fuzz=_optional_integer(spec.parameter("fuzz")),
            step=_optional_integer(spec.parameter("step")),
            seeds=spec.seeds_by_kind(),
            parameters=spec.parameter_hash(),
            chromosomes=tuple(region.seq_region_name for region in corpus.sequence_regions),
            transcripts=len(corpus.transcripts),
            variants=len(distinct_variants(corpus)),
            overlaps=len(corpus.named_overlaps),
            identity=_identity(contents),
            files=tuple(contents),
        )
        with TemporaryDirectory(prefix=".corpus-export-", dir=directory) as staging_name:
            staging = Path(staging_name)
            for file_name, text in contents.items():
                (staging / file_name).write_text(text, encoding="utf-8")
            (staging / MANIFEST_NAME).write_text(manifest.to_json(), encoding="utf-8")
            for file_name in contents:
                (staging / file_name).replace(directory / file_name)
            # Commit the manifest last, so readers never see a new identity before its data.
            (staging / MANIFEST_NAME).replace(directory / MANIFEST_NAME)
        return manifest
