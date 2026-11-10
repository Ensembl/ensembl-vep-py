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
"""The corpus on disk: written, indexed, annotated against, and removed again.

An end-to-end benchmark needs what a real run needs -- a VCF, a bgzip-compressed and tabix-indexed
GFF3, and a FASTA with a ``.fai`` beside it -- because an unindexed source takes a different code
path, and timing that path would answer a different question.

Everything this module writes is recorded as it is written, so ``clean()`` can remove exactly what
this invocation created and nothing else. A corpus already in the directory is reused and left
alone, which is what makes ``--corpus_dir`` safe to point at a corpus worth keeping.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path

from attrs import define, field

from ensembl.vep.corpus.base import Corpus
from ensembl.vep.corpus.exporter import get_exporter
from ensembl.vep.corpus.exporter.formats.file.base import CORPUS_FILES
from ensembl.vep.corpus.exporter.formats.file.fasta import FASTA_NAME
from ensembl.vep.corpus.exporter.formats.file.gff3 import GFF3_NAME
from ensembl.vep.corpus.exporter.formats.file.vcf import VCF_NAME

__all__ = ["INDEXING_TOOLS", "MissingTool", "Workspace"]

logger = logging.getLogger(__name__)

INDEXING_TOOLS = ("bgzip", "tabix")
"""The htslib executables the corpus has to be indexed with."""

GFF3_BGZ_NAME = "corpus.gff3.bgz"
"""The compressed, coordinate-sorted GFF3 the benchmark annotates against."""

_TABIX_INDEX_SUFFIX = ".tbi"
"""Suffix ``tabix -p gff`` writes beside the compressed GFF3."""

_FASTA_INDEX_SUFFIX = ".fai"
"""Suffix of the FASTA index the Rust reader needs for random access."""


class MissingTool(RuntimeError):
    """An external executable the benchmark cannot run without is not on PATH."""


def require_indexing_tools() -> None:
    """Check that every indexing executable is on PATH.

    Raises:
        MissingTool: If any of them is missing, naming all of them at once so one install fixes it.
    """
    missing = [name for name in INDEXING_TOOLS if shutil.which(name) is None]
    if missing:
        raise MissingTool(
            f"{', '.join(missing)} not found on PATH; the end-to-end benchmark annotates against an "
            "indexed corpus, so htslib is required (for example 'brew install htslib', or put an "
            "existing htslib bin directory on PATH)"
        )


def fasta_index_lines(fasta_text: str) -> list[str]:
    """Return ``.fai`` records for a FASTA whose records are wrapped to a fixed width.

    Written here rather than shelled out to ``samtools faidx``: the corpus writer controls the line
    width, so the offsets are arithmetic, and it keeps the tool requirement down to htslib.

    Args:
        fasta_text (str): The FASTA file's complete text.

    Returns:
        list[str]: One tab-separated ``name length offset linebases linewidth`` record per
        sequence, in file order.
    """
    records: list[str] = []
    name: str | None = None
    offset = 0
    length = 0
    linebases = 0
    linewidth = 0
    position = 0

    for line in fasta_text.splitlines(keepends=True):
        if line.startswith(">"):
            if name is not None:
                records.append(f"{name}\t{length}\t{offset}\t{linebases}\t{linewidth}")
            name = line[1:].strip().split(maxsplit=1)[0]
            position += len(line)
            offset = position
            length = 0
            linebases = 0
            linewidth = 0
            continue
        position += len(line)
        if name is None:
            continue
        sequence = line.strip()
        if not linebases:
            linebases = len(sequence)
            linewidth = len(line)
        length += len(sequence)

    if name is not None:
        records.append(f"{name}\t{length}\t{offset}\t{linebases}\t{linewidth}")
    return records


def sorted_gff3_text(gff3_text: str) -> str:
    """Return a GFF3 sorted the way tabix requires: headers first, then rows by region and start.

    The corpus writes its rows one transcript at a time, which tabix rejects. The sort is stable,
    so rows sharing a start keep the order the corpus emitted them in -- and with it the order the
    CSQ entries come out in.

    Args:
        gff3_text (str): The GFF3 file's complete text.

    Returns:
        str: The same rows, sorted, with a trailing newline.
    """
    headers: list[str] = []
    rows: list[tuple[str, int, str]] = []
    for line in gff3_text.splitlines():
        if not line:
            continue
        if line.startswith("#"):
            headers.append(line)
            continue
        columns = line.split("\t")
        rows.append((columns[0], int(columns[3]), line))

    rows.sort(key=lambda row: (row[0], row[1]))
    return "\n".join(headers + [row[2] for row in rows]) + "\n"


@define(kw_only=True, auto_attribs=True)
class Workspace:
    """A directory holding the corpus, its indexes and the benchmark's own output files.

    Attributes:
        directory (Path): Where everything lives.
        clean (bool): Remove what this invocation created when ``close`` is called.
        created (list[Path]): Files this invocation wrote, in creation order.
        created_directory (bool): Whether ``directory`` itself was created here.
    """

    directory: Path
    clean: bool = True
    created: list[Path] = field(factory=list, init=False)
    created_directory: bool = field(default=False, init=False)

    @property
    def fasta_path(self) -> Path:
        """Return the corpus FASTA path."""
        return self.directory / FASTA_NAME

    @property
    def gff3_path(self) -> Path:
        """Return the uncompressed, unsorted corpus GFF3 path."""
        return self.directory / GFF3_NAME

    @property
    def indexed_gff3_path(self) -> Path:
        """Return the compressed, sorted, tabix-indexed GFF3 path."""
        return self.directory / GFF3_BGZ_NAME

    @property
    def vcf_path(self) -> Path:
        """Return the corpus VCF input path."""
        return self.directory / VCF_NAME

    def open(self) -> None:
        """Create the directory, remembering whether it had to be created."""
        if not self.directory.exists():
            self.created_directory = True
        self.directory.mkdir(parents=True, exist_ok=True)

    def materialise(self, corpus: Corpus) -> None:
        """Write the corpus as files, or leave a corpus already in the directory alone.

        Args:
            corpus (Corpus): The corpus to write.
        """
        present = [name for name in CORPUS_FILES if (self.directory / name).exists()]
        if len(present) == len(CORPUS_FILES):
            logger.info("reusing the corpus already in %s", self.directory)
            return

        manifest = get_exporter("file").export(corpus, self.directory, force=True)
        for name in CORPUS_FILES:
            self.track(self.directory / name)
        logger.info(
            "wrote the corpus to %s: %d transcripts, %d variants",
            self.directory,
            manifest.transcripts,
            manifest.variants,
        )

    def variant_count(self) -> int:
        """Return how many variant rows the corpus VCF holds."""
        rows = 0
        with self.vcf_path.open(encoding="utf-8") as handle:
            for line in handle:
                if not line.startswith("#") and line.strip():
                    rows += 1
        return rows

    def index(self) -> None:
        """Write the FASTA index, then the sorted, compressed and tabix-indexed GFF3.

        Raises:
            MissingTool: If bgzip or tabix is not on PATH.
            subprocess.CalledProcessError: If either of them fails.
        """
        require_indexing_tools()
        self._write_fasta_index()
        self._write_indexed_gff3()

    def track(self, path: Path) -> None:
        """Record a path this invocation created, so ``close`` can remove it."""
        if path not in self.created:
            self.created.append(path)

    def close(self) -> None:
        """Remove what this invocation created, if cleaning was asked for."""
        if not self.clean:
            logger.info("keeping %d file(s) in %s", len(self.created), self.directory)
            return

        for path in reversed(self.created):
            path.unlink(missing_ok=True)
        if self.created_directory and self.directory.exists() and not any(self.directory.iterdir()):
            self.directory.rmdir()
        logger.info("removed %d file(s) from %s", len(self.created), self.directory)

    def _write_fasta_index(self) -> None:
        """Write ``corpus.fa.fai`` beside the corpus FASTA."""
        index_path = self.fasta_path.with_name(self.fasta_path.name + _FASTA_INDEX_SUFFIX)
        if index_path.exists():
            return
        records = fasta_index_lines(self.fasta_path.read_text(encoding="utf-8"))
        index_path.write_text("\n".join(records) + "\n", encoding="utf-8")
        self.track(index_path)

    def _write_indexed_gff3(self) -> None:
        """Compress the sorted GFF3 with bgzip and index it with tabix."""
        index_path = self.indexed_gff3_path.with_name(self.indexed_gff3_path.name + _TABIX_INDEX_SUFFIX)
        if self.indexed_gff3_path.exists() and index_path.exists():
            return

        sorted_text = sorted_gff3_text(self.gff3_path.read_text(encoding="utf-8"))
        compressed = subprocess.run(
            ["bgzip", "-c"],
            input=sorted_text.encode("utf-8"),
            capture_output=True,
            check=True,
        )
        self.indexed_gff3_path.write_bytes(compressed.stdout)
        self.track(self.indexed_gff3_path)

        subprocess.run(
            ["tabix", "-p", "gff", "-f", str(self.indexed_gff3_path)],
            capture_output=True,
            check=True,
        )
        self.track(index_path)
