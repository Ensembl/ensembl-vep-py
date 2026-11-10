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
"""Tests for the corpus runners."""

from __future__ import annotations

from pathlib import Path

import pytest

import ensembl.vep.annotators  # noqa: F401  pylint: disable=unused-import

# ^ keeps the annotators-before-caches import order.
from ensembl.vep.cli import main
from ensembl.vep.corpus.exporter.formats.file.base import CORPUS_FILES
from ensembl.vep.runner.base import ExitCode

SMALL_CORPUS = ["--kind", "synthetic.default", "--fuzz", "5"]


class TestDescribe:
    """``vep corpus describe``."""

    def test_it_lists_the_shape_and_the_overlaps(self, capsys: pytest.CaptureFixture[str]) -> None:
        """Every table is present, and the listing carries a selectable name."""
        assert main(["corpus", "describe", *SMALL_CORPUS, "--output_format", "text"]) == ExitCode.OK
        out = capsys.readouterr().out
        for heading in ("Chromosomes", "Transcripts", "Variant classes", "Overlaps"):
            assert heading in out
        assert "ENST0001_cds3_del1" in out

    def test_html_is_a_page(self, capsys: pytest.CaptureFixture[str]) -> None:
        """The HTML format renders the tables into one document."""
        assert main(["corpus", "describe", *SMALL_CORPUS, "--output_format", "html"]) == ExitCode.OK
        out = capsys.readouterr().out
        assert "<table" in out
        assert "Corpus: synthetic.default" in out

    def test_select_narrows_to_one_variant(self, capsys: pytest.CaptureFixture[str]) -> None:
        """Only the chosen variant's overlaps are listed."""
        assert main(["corpus", "describe", *SMALL_CORPUS, "--select", "ENST0001_cds3_del1"]) == ExitCode.OK
        out = capsys.readouterr().out
        assert "ENST0001_cds3_del1" in out
        assert "ENST0001_intron0_donor+1" not in out

    def test_an_unknown_selection_is_a_usage_error(self) -> None:
        """A name that matches nothing is the caller's mistake."""
        assert main(["corpus", "describe", *SMALL_CORPUS, "--select", "no_such_variant"]) == ExitCode.USAGE

    def test_it_can_write_to_a_file(self, tmp_path: Path) -> None:
        """--output_file takes the rendering instead of stdout."""
        destination = tmp_path / "corpus.md"
        assert main(["corpus", "describe", *SMALL_CORPUS, "--output_file", str(destination)]) == ExitCode.OK
        assert "Overlaps" in destination.read_text(encoding="utf-8")

    def test_annotate_adds_the_term_counts(self, capsys: pytest.CaptureFixture[str]) -> None:
        """The corpus's terms, counted by the configured annotator."""
        assert (
            main(["corpus", "describe", *SMALL_CORPUS, "--output_format", "text", "--annotate"])
            == ExitCode.OK
        )
        out = capsys.readouterr().out
        assert "Most severe consequences" in out
        assert "Coding consequences" in out
        assert "inframe_deletion" in out

    def test_without_annotate_nothing_is_annotated(self, capsys: pytest.CaptureFixture[str]) -> None:
        """A listing stays a listing unless asked."""
        assert main(["corpus", "describe", *SMALL_CORPUS, "--output_format", "text"]) == ExitCode.OK
        assert "Most severe consequences" not in capsys.readouterr().out

    def test_registered_file_corpus_can_be_selected_from_the_command_line(self, tmp_path: Path) -> None:
        """Registry names and the bundle directory reach the generic file loader."""
        assert main(["corpus", "save", *SMALL_CORPUS, "--corpus_dir", str(tmp_path)]) == ExitCode.OK
        assert (
            main(
                [
                    "corpus",
                    "describe",
                    "--kind",
                    "loaded.files",
                    "--corpus_dir",
                    str(tmp_path),
                    "--select",
                    "ENST0001_cds3_del1",
                ]
            )
            == ExitCode.OK
        )

    def test_malformed_factory_parameter_is_a_usage_error(self) -> None:
        """A malformed repeatable parameter is reported without a traceback."""
        assert (
            main(
                [
                    "corpus",
                    "describe",
                    "--kind",
                    "loaded.files",
                    "--corpus_parameter",
                    "broken",
                ]
            )
            == ExitCode.USAGE
        )

    def test_missing_loaded_bundle_is_a_usage_error(self, tmp_path: Path) -> None:
        """Unreadable external corpus input is reported without a traceback."""
        missing = tmp_path / "missing"
        assert (
            main(
                [
                    "corpus",
                    "describe",
                    "--kind",
                    "loaded.files",
                    "--corpus_dir",
                    str(missing),
                ]
            )
            == ExitCode.USAGE
        )


class TestSave:
    """``vep corpus save``."""

    def test_it_writes_the_corpus_files(self, tmp_path: Path) -> None:
        """Every corpus file lands in the directory."""
        assert main(["corpus", "save", *SMALL_CORPUS, "--corpus_dir", str(tmp_path)]) == ExitCode.OK
        for file_name in CORPUS_FILES:
            assert (tmp_path / file_name).exists()

    def test_a_second_save_is_refused_without_force(self, tmp_path: Path) -> None:
        """Overwriting is explicit."""
        assert main(["corpus", "save", *SMALL_CORPUS, "--corpus_dir", str(tmp_path)]) == ExitCode.OK
        assert main(["corpus", "save", *SMALL_CORPUS, "--corpus_dir", str(tmp_path)]) == ExitCode.FAILURE
        assert (
            main(["corpus", "save", *SMALL_CORPUS, "--corpus_dir", str(tmp_path), "--force"]) == ExitCode.OK
        )

    def test_a_missing_directory_is_a_usage_error(self) -> None:
        """There is no default destination to guess."""
        assert main(["corpus", "save", *SMALL_CORPUS]) == ExitCode.USAGE

    def test_the_registered_file_exporter_can_be_selected(self, tmp_path: Path) -> None:
        """The command exposes the exporter registry without changing its default."""
        assert (
            main(
                [
                    "corpus",
                    "save",
                    *SMALL_CORPUS,
                    "--corpus_dir",
                    str(tmp_path),
                    "--exporter",
                    "file",
                ]
            )
            == ExitCode.OK
        )


class TestGolden:
    """``vep corpus golden``."""

    SPEC = ["--kind", "synthetic.default", "--fuzz", "0"]

    def test_record_then_check_is_ok(self, tmp_path: Path) -> None:
        """A check straight after recording finds nothing, reading the corpus from the file."""
        golden_file = tmp_path / "golden.json"
        assert main(
            ["corpus", "golden", *self.SPEC, "--golden_file", str(golden_file), "--update_golden"]
        ) == (ExitCode.OK)
        assert main(["corpus", "golden", "--golden_file", str(golden_file)]) == ExitCode.OK

    def test_a_tampered_file_fails(self, tmp_path: Path) -> None:
        """A moved term is a regression, and the command says so with its exit code."""
        golden_file = tmp_path / "golden.json"
        main(["corpus", "golden", *self.SPEC, "--golden_file", str(golden_file), "--update_golden"])
        text = golden_file.read_text(encoding="utf-8")
        golden_file.write_text(
            text.replace('["inframe_deletion"]', '["frameshift_variant"]', 1), encoding="utf-8"
        )
        assert main(["corpus", "golden", "--golden_file", str(golden_file)]) == ExitCode.FAILURE

    def test_another_strategy_checks_against_linears_terms(self, tmp_path: Path) -> None:
        """An explicit --annotator is checked against the file; a conforming strategy passes."""
        golden_file = tmp_path / "golden.json"
        main(["corpus", "golden", *self.SPEC, "--golden_file", str(golden_file), "--update_golden"])
        assert (
            main(["corpus", "golden", "--golden_file", str(golden_file), "--annotator", "consequence.staged"])
            == ExitCode.OK
        )

    def test_a_missing_file_is_a_usage_error(self, tmp_path: Path) -> None:
        """Nothing to check against, and the message says how to record one."""
        assert main(["corpus", "golden", "--golden_file", str(tmp_path / "absent.json")]) == ExitCode.USAGE

    def test_no_file_is_a_usage_error(self) -> None:
        """There is no default golden file to guess."""
        assert main(["corpus", "golden"]) == ExitCode.USAGE
