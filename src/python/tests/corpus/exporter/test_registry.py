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
"""Tests for corpus exporter registration."""

from __future__ import annotations

from pathlib import Path

import pytest

from ensembl.vep.corpus.exporter import available_exporters, get_exporter, register_exporter
from ensembl.vep.corpus.exporter.formats.file.base import CORPUS_FILES, FileCorpusExporter
from ensembl.vep.corpus.exporter.formats.file.gff3 import escape_attribute_value
from ensembl.vep.corpus.registry import build_corpus
from ensembl.vep.corpus.spec import CorpusSpec


class TestFileExporter:
    """The built-in coordinated file exporter."""

    def test_is_registered(self) -> None:
        """The file exporter resolves from the registry."""
        assert "file" in available_exporters()
        assert isinstance(get_exporter("file"), FileCorpusExporter)

    def test_unknown_exporter_lists_the_available_names(self) -> None:
        """Unknown exporter names describe the registered alternatives."""
        with pytest.raises(KeyError, match="Unknown corpus exporter 'not.a.real.exporter'.*file"):
            get_exporter("not.a.real.exporter")

    def test_duplicate_exporter_registration_is_rejected(self) -> None:
        """Registration cannot silently replace the exporter consumers already resolved."""
        with pytest.raises(ValueError, match="already registered"):
            register_exporter("file")(FileCorpusExporter)

    def test_rejects_existing_corpus_files_without_force(self, tmp_path: Path) -> None:
        """One existing bundle file prevents a partial overwrite."""
        (tmp_path / CORPUS_FILES[0]).write_text("previous", encoding="utf-8")

        with pytest.raises(FileExistsError, match="pass --force to overwrite"):
            get_exporter("file").export(build_corpus(CorpusSpec(kind="default", fuzz=0)), tmp_path)

    def test_force_replaces_an_existing_file_bundle(self, tmp_path: Path) -> None:
        """Force writes every artifact after a previous partial bundle."""
        (tmp_path / CORPUS_FILES[0]).write_text("previous", encoding="utf-8")

        manifest = get_exporter("file").export(
            build_corpus(CorpusSpec(kind="default", fuzz=0)), tmp_path, force=True
        )

        assert manifest.files == CORPUS_FILES[:-1]
        assert all((tmp_path / file_name).is_file() for file_name in CORPUS_FILES)

    def test_gff3_attribute_values_escape_delimiters(self) -> None:
        """Loaded identifiers cannot inject another GFF3 attribute."""
        assert escape_attribute_value("gene;Name=other% value") == "gene%3BName%3Dother%25%20value"
