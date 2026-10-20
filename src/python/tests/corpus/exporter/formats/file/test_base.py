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
"""Tests for coordinated file-bundle export."""

from __future__ import annotations

from pathlib import Path

import pytest

from ensembl.vep.corpus import Corpus
from ensembl.vep.corpus.exporter import get_exporter
from ensembl.vep.corpus.exporter.formats.file.base import CORPUS_FILES
from ensembl.vep.corpus.exporter.formats.file.vcf import VCF_NAME
from ensembl.vep.corpus.exporter.manifest import CorpusManifest


class TestFileBundle:
    """The coordinated FASTA, GFF3, VCF and manifest write."""

    def test_manifest_counts_and_files_match_the_corpus(
        self, source: Corpus, saved: tuple[Path, CorpusManifest]
    ) -> None:
        """The manifest reports the same aggregate values written by the exporter."""
        directory, manifest = saved
        assert manifest.kind == "synthetic.all"
        assert manifest.chromosomes == tuple(region.seq_region_name for region in source.sequence_regions)
        assert manifest.transcripts == len(source.transcripts)
        assert manifest.variants == len(source.named_variants)
        assert manifest.overlaps == len(source.named_overlaps)
        assert all((directory / name).is_file() for name in CORPUS_FILES)

    def test_same_corpus_exports_byte_identically(self, source: Corpus, tmp_path: Path) -> None:
        """A repeatable corpus has repeatable files and identity."""
        exporter = get_exporter("file")
        first = exporter.export(source, tmp_path / "one")
        second = exporter.export(source, tmp_path / "two")
        assert first.identity == second.identity
        for name in CORPUS_FILES:
            assert (tmp_path / "one" / name).read_bytes() == (tmp_path / "two" / name).read_bytes()

    def test_existing_bundle_requires_force(self, source: Corpus, tmp_path: Path) -> None:
        """A partial or complete existing bundle is never overwritten accidentally."""
        exporter = get_exporter("file")
        exporter.export(source, tmp_path)
        with pytest.raises(FileExistsError, match=VCF_NAME):
            exporter.export(source, tmp_path)
        assert exporter.export(source, tmp_path, force=True).identity
