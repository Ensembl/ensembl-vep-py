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
"""Tests for the on-disk corpus a end-to-end benchmark annotates against."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from ensembl.vep.bench.workspace import (
    INDEXING_TOOLS,
    MissingTool,
    Workspace,
    fasta_index_lines,
    require_indexing_tools,
    sorted_gff3_text,
)
from ensembl.vep.corpus import Corpus, build_corpus
from ensembl.vep.corpus.exporter.formats.file.base import CORPUS_FILES
from ensembl.vep.corpus.spec import CorpusSpec

INDEXING_TOOLS_PRESENT = all(shutil.which(name) is not None for name in INDEXING_TOOLS)


@pytest.fixture(name="source")
def fixture_source() -> Corpus:
    """Return the smallest synthetic corpus, so writing it stays quick."""
    return build_corpus(CorpusSpec(kind="synthetic.default", fuzz=4))


class TestFastaIndexLines:
    """Deriving a .fai from a fixed-width FASTA."""

    def test_records_name_length_and_offset(self) -> None:
        """One record per sequence, carrying its length, first-base offset and line widths."""
        text = ">chr1\nACGTA\nCGTAC\nGT\n>chr2\nTTTT\n"

        records = fasta_index_lines(text)

        assert records == ["chr1\t12\t6\t5\t6", "chr2\t4\t27\t4\t5"]

    def test_ignores_sequence_before_any_header(self) -> None:
        """Text before the first header belongs to no sequence."""
        assert fasta_index_lines("ACGT\n>chr1\nAC\n") == ["chr1\t2\t11\t2\t3"]

    def test_empty_text_has_no_records(self) -> None:
        """A FASTA with nothing in it indexes to nothing."""
        assert not fasta_index_lines("")


class TestSortedGff3Text:
    """Reordering the corpus GFF3 into the order tabix requires."""

    def test_headers_come_first_and_rows_sort_by_region_and_start(self) -> None:
        """Every ## line is kept above rows ordered by sequence region, then start."""
        text = "\n".join(
            (
                "##gff-version 3",
                "chr2\tsrc\tgene\t50\t60\t.\t+\t.\tID=b",
                "chr1\tsrc\tgene\t90\t99\t.\t+\t.\tID=c",
                "chr1\tsrc\tgene\t10\t20\t.\t+\t.\tID=a",
            )
        )

        rows = sorted_gff3_text(text).splitlines()

        assert rows[0] == "##gff-version 3"
        assert [row.split("\t")[8] for row in rows[1:]] == ["ID=a", "ID=c", "ID=b"]

    def test_blank_lines_are_dropped(self) -> None:
        """A blank line is no row, and tabix would reject it."""
        text = "##gff-version 3\n\nchr1\tsrc\tgene\t10\t20\t.\t+\t.\tID=a\n"

        assert sorted_gff3_text(text).splitlines() == [
            "##gff-version 3",
            "chr1\tsrc\tgene\t10\t20\t.\t+\t.\tID=a",
        ]

    def test_rows_sharing_a_start_keep_their_order(self) -> None:
        """The sort is stable, so the corpus's own row order survives a tie."""
        text = "\n".join(
            (
                "chr1\tsrc\ttranscript\t10\t99\t.\t+\t.\tID=first",
                "chr1\tsrc\texon\t10\t20\t.\t+\t.\tID=second",
            )
        )

        rows = sorted_gff3_text(text).splitlines()

        assert [row.split("\t")[8] for row in rows] == ["ID=first", "ID=second"]


class TestRequireIndexingTools:
    """Refusing to run without htslib."""

    def test_missing_tools_are_named_together(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """One error names every missing executable, so one install fixes it."""
        monkeypatch.setattr(shutil, "which", lambda name: None)

        with pytest.raises(MissingTool) as raised:
            require_indexing_tools()

        assert "bgzip" in str(raised.value)
        assert "tabix" in str(raised.value)

    def test_present_tools_raise_nothing(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Nothing is raised once both executables resolve."""
        monkeypatch.setattr(shutil, "which", lambda name: f"/usr/bin/{name}")

        require_indexing_tools()


class TestWorkspaceMaterialise:
    """Writing the corpus, or reusing one already there."""

    def test_writes_the_corpus_and_tracks_every_file(self, tmp_path: Path, source: Corpus) -> None:
        """A corpus written here is recorded as this invocation's to remove."""
        workspace = Workspace(directory=tmp_path / "work")
        workspace.open()

        workspace.materialise(source)

        for name in CORPUS_FILES:
            assert (workspace.directory / name).exists()
            assert workspace.directory / name in workspace.created

    def test_reused_corpus_is_not_tracked(self, tmp_path: Path, source: Corpus) -> None:
        """A corpus that was already there belongs to whoever put it there."""
        first = Workspace(directory=tmp_path, clean=False)
        first.open()
        first.materialise(source)

        second = Workspace(directory=tmp_path)
        second.open()
        second.materialise(source)

        assert second.created == []

    def test_variant_count_matches_the_corpus_vcf(self, tmp_path: Path, source: Corpus) -> None:
        """The count skips headers and blank lines."""
        workspace = Workspace(directory=tmp_path)
        workspace.open()
        workspace.materialise(source)

        rows = [
            line
            for line in workspace.vcf_path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")
        ]

        assert workspace.variant_count() == len(rows)


class TestWorkspaceClose:
    """What survives the benchmark."""

    def test_clean_removes_only_what_was_created(self, tmp_path: Path, source: Corpus) -> None:
        """A file that was already in the directory is left alone."""
        existing = tmp_path / "keep.txt"
        existing.write_text("mine", encoding="utf-8")
        workspace = Workspace(directory=tmp_path)
        workspace.open()
        workspace.materialise(source)

        workspace.close()

        assert existing.exists()
        assert not (tmp_path / "corpus.vcf").exists()

    def test_clean_removes_a_directory_it_created(self, tmp_path: Path, source: Corpus) -> None:
        """The directory goes too, but only when this invocation made it."""
        directory = tmp_path / "made-here"
        workspace = Workspace(directory=directory)
        workspace.open()
        workspace.materialise(source)

        workspace.close()

        assert not directory.exists()

    def test_clean_keeps_a_directory_that_was_already_there(self, tmp_path: Path, source: Corpus) -> None:
        """A directory the caller chose outlives the benchmark, empty or not."""
        workspace = Workspace(directory=tmp_path)
        workspace.open()
        workspace.materialise(source)

        workspace.close()

        assert tmp_path.exists()

    def test_no_clean_keeps_everything(self, tmp_path: Path, source: Corpus) -> None:
        """Nothing is removed when cleaning was turned off."""
        workspace = Workspace(directory=tmp_path, clean=False)
        workspace.open()
        workspace.materialise(source)

        workspace.close()

        assert workspace.vcf_path.exists()

    def test_tracking_a_path_twice_records_it_once(self, tmp_path: Path) -> None:
        """Re-tracking an output between repeats does not duplicate the removal."""
        workspace = Workspace(directory=tmp_path)
        output = tmp_path / "out.vcf"

        workspace.track(output)
        workspace.track(output)

        assert workspace.created == [output]

    def test_close_tolerates_a_file_already_gone(self, tmp_path: Path) -> None:
        """A tracked file someone else removed does not break the cleanup."""
        workspace = Workspace(directory=tmp_path)
        workspace.track(tmp_path / "never-written.vcf")

        workspace.close()


@pytest.mark.integration
@pytest.mark.skipif(not INDEXING_TOOLS_PRESENT, reason="bgzip and tabix are required")
class TestWorkspaceIndex:
    """Indexing the corpus with htslib."""

    def test_writes_a_compressed_gff3_and_its_index(self, tmp_path: Path, source: Corpus) -> None:
        """Both the bgzip output and the tabix index are written and tracked."""
        workspace = Workspace(directory=tmp_path)
        workspace.open()
        workspace.materialise(source)

        workspace.index()

        index_path = workspace.indexed_gff3_path.with_suffix(workspace.indexed_gff3_path.suffix + ".tbi")
        assert workspace.indexed_gff3_path.exists()
        assert index_path.exists()
        assert workspace.indexed_gff3_path in workspace.created
        assert index_path in workspace.created

    def test_writes_a_fasta_index(self, tmp_path: Path, source: Corpus) -> None:
        """The FASTA gets the .fai the reader needs for random access."""
        workspace = Workspace(directory=tmp_path)
        workspace.open()
        workspace.materialise(source)

        workspace.index()

        assert workspace.fasta_path.with_suffix(".fa.fai").exists()

    def test_indexing_twice_writes_nothing_new(self, tmp_path: Path, source: Corpus) -> None:
        """A second call over the same workspace reuses what is already indexed."""
        workspace = Workspace(directory=tmp_path)
        workspace.open()
        workspace.materialise(source)
        workspace.index()
        created = list(workspace.created)

        workspace.index()

        assert workspace.created == created
