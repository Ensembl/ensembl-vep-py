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
"""Tests for the default annotation runner."""

from __future__ import annotations

import json
import stat
import sys
from io import StringIO
from pathlib import Path

import pytest

import ensembl.vep.annotators  # noqa: F401  pylint: disable=unused-import

# ^ registers the annotator names the --annotator option validates against.
from ensembl.vep.cli import main
from ensembl.vep.config import Config, RunConfig, SourceConfig
from ensembl.vep.config.parser import check_no_flag_collisions
from ensembl.vep.output.factory.write import JSONOutputWriter, register_output_writer
from ensembl.vep.output.factory.write import registry as writer_registry
from ensembl.vep.runner import ExitCode, available_runners, get_runner
from ensembl.vep.runner.run import default as run_module
from ensembl.vep.runner.run.default import Runner
from ensembl.vep.sources.file.fasta.sequence import FastaSequenceSource
from ensembl.vep.stats import RunSummaryCollector, RunSummaryDelta

GFF3 = """\
chr1\tensembl\tgene\t1000\t5000\t.\t+\t.\tID=gene:G1;Name=GENE1;biotype=protein_coding
chr1\tensembl\tmRNA\t1000\t5000\t.\t+\t.\tID=transcript:ENST001;Parent=gene:G1;biotype=protein_coding
chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tID=exon:E1;Parent=transcript:ENST001
chr1\tensembl\tCDS\t1100\t2000\t.\t+\t0\tID=cds:P1;Parent=transcript:ENST001;protein_id=P1
chr1\tensembl\texon\t3000\t5000\t.\t+\t.\tID=exon:E2;Parent=transcript:ENST001
chr1\tensembl\tCDS\t3000\t4900\t.\t+\t0\tID=cds:P1;Parent=transcript:ENST001;protein_id=P1
"""

VCF = """\
##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t1500\trs_coding\tA\tC\t.\tPASS\t.
chr1\t30000\trs_far\tA\tC\t.\tPASS\t.
chr1\t1700\trs_multi\tA\tC,G\t.\tPASS\t.
chr1\t1200\trs_indel\tC\tCAT\t.\tPASS\t.
"""
"""One variant in a transcript, one beyond every transcript, one multi-allelic, and one indel
carrying the anchor base a VCF adds."""


class _Terminal(StringIO):
    """String stream that represents an interactive terminal."""

    def isatty(self) -> bool:
        """Return true so progress is enabled."""
        return True


def _inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    """Write the VCF, GFF3 and FASTA a run needs, and return their paths."""
    variants = tmp_path / "variants.vcf"
    variants.write_text(VCF, encoding="utf-8")
    transcripts = tmp_path / "transcripts.gff3"
    transcripts.write_text(GFF3, encoding="utf-8")
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\n" + "ACGT" * 10000 + "\n", encoding="utf-8")
    fasta.with_suffix(".fa.fai").write_text("chr1\t40000\t6\t40000\t40001\n", encoding="utf-8")
    return variants, transcripts, fasta


def _run(tmp_path: Path, *extra: str) -> tuple[ExitCode, Path]:
    """Annotate the fixtures into a file, returning the exit code and that file."""
    variants, transcripts, fasta = _inputs(tmp_path)
    output = tmp_path / "annotated.vcf"
    code = main(
        [
            "run",
            "--input_file",
            str(variants),
            "--gff",
            str(transcripts),
            "--fasta",
            str(fasta),
            "--output_file",
            str(output),
            *extra,
        ]
    )
    return code, output


class TestRunUsage:
    """What the run path refuses, and why it is the caller's mistake rather than a failure."""

    def test_annotating_needs_an_input_file(self) -> None:
        """There is no sensible default input; a corpus is explored with another command."""
        assert main(["run"]) == ExitCode.USAGE

    def test_the_python_api_refuses_to_annotate_without_an_input_file(self) -> None:
        """A Python caller reaching ``annotate`` directly is told, rather than reading an empty stream."""
        with pytest.raises(ValueError, match="give --input_file"):
            Runner(Config()).annotate()

    def test_an_input_format_with_no_reader_is_rejected(self) -> None:
        """Refuse an input format nothing can read before any worker opens the file.

        The command line refuses an unregistered ``--format`` itself, so this is the Python
        caller's path into the same check.
        """
        config = Config(source=SourceConfig(input_file=Path("variants.bed"), format="bed"))

        assert Runner(config).run() == ExitCode.USAGE

    def test_an_output_format_with_no_writer_is_rejected(self) -> None:
        """Reject a configured format that has no writer implementation."""
        assert main(["run", "--input_file", "variants.vcf", "--output_format", "tab"]) == ExitCode.USAGE

    def test_an_existing_output_file_is_not_overwritten_silently(self, tmp_path: Path) -> None:
        """Annotation output is expensive to reproduce, so replacing one has to be asked for."""
        existing = tmp_path / "out.vcf"
        existing.write_text("keep me\n", encoding="utf-8")

        code = main(["run", "--input_file", "variants.vcf", "--output_file", str(existing)])

        assert code == ExitCode.USAGE
        assert existing.read_text(encoding="utf-8") == "keep me\n"

    @pytest.mark.parametrize("explicit_summary", [False, True])
    def test_an_existing_summary_is_not_overwritten_silently(
        self,
        explicit_summary: bool,
        tmp_path: Path,
    ) -> None:
        """Every file a run plans to replace requires explicit permission."""

        output = tmp_path / "new.vcf"
        summary = tmp_path / "report.txt" if explicit_summary else tmp_path / "new.vcf_summary.txt"
        summary.write_text("keep me\n", encoding="utf-8")
        arguments = [
            "run",
            "--input_file",
            "variants.vcf",
            "--output_file",
            str(output),
            "--stats_format",
            "text",
        ]
        if explicit_summary:
            arguments.extend(("--stats_file", str(summary)))

        code = main(arguments)

        assert code == ExitCode.USAGE
        assert not output.exists()
        assert summary.read_text(encoding="utf-8") == "keep me\n"

    def test_it_holds_the_configuration_it_was_given(self) -> None:
        """The runner is the composition root; nothing below it sees the config."""
        config = Config()
        assert Runner(config).config is config


class TestRunOutput:  # pylint: disable=too-many-public-methods
    """What a run writes."""

    @staticmethod
    def _data_lines(output: Path) -> list[str]:
        """The output's records, without its headers."""
        return [line for line in output.read_text(encoding="utf-8").splitlines() if not line.startswith("#")]

    def test_it_writes_a_csq_vcf(self, tmp_path: Path) -> None:
        """The output is a VCF carrying consequences in a CSQ INFO field it declares."""
        code, output = _run(tmp_path)

        assert code == ExitCode.OK
        text = output.read_text(encoding="utf-8")
        assert text.startswith("##fileformat=VCFv4")
        assert "##INFO=<ID=CSQ," in text
        assert "CSQ=" in text

    def test_a_synonyms_file_pairs_an_accession_input_with_chr_named_sources(self, tmp_path: Path) -> None:
        """A VCF spelling the region as a RefSeq accession annotates against chr-named sources."""
        _variants, transcripts, fasta = _inputs(tmp_path)
        variants = tmp_path / "accession.vcf"
        variants.write_text(
            "##fileformat=VCFv4.3\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "NC_000001.11\t1500\trs_coding\tA\tC\t.\tPASS\t.\n",
            encoding="utf-8",
        )
        synonyms = tmp_path / "synonyms.txt"
        synonyms.write_text("chr1\tNC_000001.11\n", encoding="utf-8")
        output = tmp_path / "annotated.vcf"

        code = main(
            [
                "run",
                "--input_file",
                str(variants),
                "--gff",
                str(transcripts),
                "--fasta",
                str(fasta),
                "--synonyms",
                str(synonyms),
                "--output_file",
                str(output),
                "--no_stats",
            ]
        )

        assert code == ExitCode.OK
        assert "ENST001" in self._data_lines(output)[0]

    def test_an_accession_input_without_synonyms_annotates_nothing(self, tmp_path: Path) -> None:
        """Without the map the accession is its own region, which no transcript sits on."""
        _variants, transcripts, fasta = _inputs(tmp_path)
        variants = tmp_path / "accession.vcf"
        variants.write_text(
            "##fileformat=VCFv4.3\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "NC_000001.11\t1500\trs_coding\tA\tC\t.\tPASS\t.\n",
            encoding="utf-8",
        )
        output = tmp_path / "annotated.vcf"

        code = main(
            [
                "run",
                "--input_file",
                str(variants),
                "--gff",
                str(transcripts),
                "--fasta",
                str(fasta),
                "--output_file",
                str(output),
                "--no_stats",
            ]
        )

        assert code == ExitCode.OK
        assert "ENST001" not in self._data_lines(output)[0]

    def test_it_writes_ensembl_vep_json_lines(self, tmp_path: Path) -> None:
        """Write one Ensembl VEP JSON object containing every consequence for each variant."""
        code, output = _run(tmp_path, "--output_format", "json")

        records = [json.loads(line) for line in output.read_text(encoding="utf-8").splitlines()]

        assert code == ExitCode.OK
        assert len(records) == 4
        assert records[0]["id"] == "rs_coding"
        assert records[0]["seq_region_name"] == "chr1"
        assert records[0]["most_severe_consequence"]
        assert records[0]["transcript_consequences"]
        assert records[1]["most_severe_consequence"] == "intergenic_variant"

    def test_cli_accepts_a_registered_output_writer(self, tmp_path: Path) -> None:
        """Allow an extension to add a CLI output format without changing run options."""

        @register_output_writer("test-json")
        def build_writer(_command_line: str) -> JSONOutputWriter:
            """Build the temporary test writer."""
            return JSONOutputWriter()

        try:
            code, output = _run(tmp_path, "--output_format", "test-json")
        finally:
            writer_registry._OUTPUT_WRITER_REGISTRY.pop("test-json", None)  # pylint: disable=protected-access

        assert code == ExitCode.OK
        assert json.loads(output.read_text(encoding="utf-8").splitlines()[0])["id"] == "rs_coding"

    def test_python_api_streams_format_neutral_variants(self, tmp_path: Path) -> None:
        """Expose detached annotations without creating the configured output file."""
        variants, transcripts, fasta = _inputs(tmp_path)
        unused_output = tmp_path / "unused.vcf"
        config = Config(
            source=SourceConfig(
                input_file=variants,
                gff=transcripts,
                fasta=fasta,
            ),
            run=RunConfig(
                output_file=unused_output,
                no_progress=True,
                no_stats=True,
            ),
        )

        annotations = Runner(config).annotate()
        with annotations:
            annotated_variants = list(annotations)

        assert [variant.id for variant in annotated_variants] == [
            "rs_coding",
            "rs_far",
            "rs_multi",
            "rs_indel",
        ]
        assert annotated_variants[0].annotations
        assert annotations.summary.variants == 4
        assert not unused_output.exists()

    def test_interactive_progress_counts_input_variants(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Report ordered completed records as they are annotated.

        The input is a plain VCF, so no index declares how many records it holds and the run
        does not read it twice to find out: the count is reported without a percentage.
        """
        terminal = _Terminal()
        monkeypatch.setattr(sys, "stderr", terminal)

        code, _output = _run(tmp_path, "--no_stats")

        rendered = terminal.getvalue()
        assert code == ExitCode.OK
        assert "Annotating variants: 4 |" in rendered
        assert "%" not in rendered

    @pytest.mark.parametrize("option", ["--no_progress", "--quiet"])
    def test_progress_can_be_suppressed(
        self,
        option: str,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Keep explicitly disabled and quiet runs free of progress output."""
        terminal = _Terminal()
        monkeypatch.setattr(sys, "stderr", terminal)

        code, _output = _run(tmp_path, "--no_stats", option)

        assert code == ExitCode.OK
        assert "Annotating variants:" not in terminal.getvalue()

    def test_default_output_is_written_in_the_working_directory(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """A run without ``--output_file`` follows Perl VEP's default filename."""
        variants, transcripts, fasta = _inputs(tmp_path)
        monkeypatch.chdir(tmp_path)

        code = main(
            [
                "run",
                "--input_file",
                str(variants),
                "--gff",
                str(transcripts),
                "--fasta",
                str(fasta),
                "--no_stats",
            ]
        )

        assert code == ExitCode.OK
        assert (tmp_path / "variant_effect_output.txt").exists()

    def test_new_output_uses_normal_creation_permissions(self, tmp_path: Path) -> None:
        """A staged result retains the mode a direct new output file would receive."""
        control = tmp_path / "ordinary.txt"
        control.write_text("control\n", encoding="utf-8")
        expected_mode = stat.S_IMODE(control.stat().st_mode)

        code, output = _run(tmp_path, "--no_stats")

        assert code == ExitCode.OK
        assert stat.S_IMODE(output.stat().st_mode) == expected_mode

    def test_forced_output_retains_existing_permissions(self, tmp_path: Path) -> None:
        """Replacing an existing result does not broaden or narrow its access mode."""
        variants, transcripts, fasta = _inputs(tmp_path)
        output = tmp_path / "annotated.vcf"
        output.write_text("previous output\n", encoding="utf-8")
        output.chmod(0o640)

        code = main(
            [
                "run",
                "--input_file",
                str(variants),
                "--gff",
                str(transcripts),
                "--fasta",
                str(fasta),
                "--output_file",
                str(output),
                "--force_overwrite",
                "--no_stats",
            ]
        )

        assert code == ExitCode.OK
        assert stat.S_IMODE(output.stat().st_mode) == 0o640

    @pytest.mark.parametrize("explicit_summary", [False, True])
    def test_force_overwrite_replaces_an_existing_summary(
        self,
        explicit_summary: bool,
        tmp_path: Path,
    ) -> None:
        """Explicit replacement permission applies to both transaction files."""

        variants, transcripts, fasta = _inputs(tmp_path)
        output = tmp_path / "annotated.vcf"
        summary = (
            tmp_path / "report.txt" if explicit_summary else output.with_name(f"{output.name}_summary.txt")
        )
        summary.write_text("previous summary\n", encoding="utf-8")
        arguments = [
            "run",
            "--input_file",
            str(variants),
            "--gff",
            str(transcripts),
            "--fasta",
            str(fasta),
            "--output_file",
            str(output),
            "--stats_format",
            "text",
            "--force_overwrite",
        ]
        if explicit_summary:
            arguments.extend(("--stats_file", str(summary)))

        code = main(arguments)

        assert code == ExitCode.OK
        assert output.exists()
        assert summary.read_text(encoding="utf-8") != "previous summary\n"

    @pytest.mark.parametrize("previous", [None, "previous output\n"])
    def test_failed_run_leaves_no_partial_output(
        self,
        previous: str | None,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """A failed forced run preserves an old result or leaves no new result behind."""
        variants = tmp_path / "variants.vcf"
        variants.write_text(VCF, encoding="utf-8")
        output = tmp_path / "annotated.vcf"
        if previous is not None:
            output.write_text(previous, encoding="utf-8")
        config = Config(
            source=SourceConfig(input_file=variants),
            run=RunConfig(output_file=output, force_overwrite=True, no_stats=True),
        )

        def fail_after_partial_write(_runner: Runner, destination: Path | None) -> int:
            assert destination is not None
            assert destination != output
            destination.write_text("partial output\n", encoding="utf-8")
            raise OSError("annotation failed")

        monkeypatch.setattr(Runner, "_run_buffered", fail_after_partial_write)

        assert Runner(config).run() == ExitCode.USAGE
        if previous is None:
            assert not output.exists()
        else:
            assert output.read_text(encoding="utf-8") == previous
        assert not list(tmp_path.glob(f".{output.name}.*.tmp"))

    def test_failed_summary_merge_preserves_existing_output(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """A late collection failure occurs before the completed file replaces its predecessor."""
        variants, transcripts, fasta = _inputs(tmp_path)
        output = tmp_path / "annotated.vcf"
        output.write_text("previous output\n", encoding="utf-8")

        def fail_merge(_collector: RunSummaryCollector, _delta: RunSummaryDelta) -> None:
            raise ValueError("summary merge failed")

        monkeypatch.setattr(RunSummaryCollector, "merge_ordered_delta", fail_merge)

        code = main(
            [
                "run",
                "--input_file",
                str(variants),
                "--gff",
                str(transcripts),
                "--fasta",
                str(fasta),
                "--output_file",
                str(output),
                "--force_overwrite",
                "--no_stats",
            ]
        )

        assert code == ExitCode.USAGE
        assert output.read_text(encoding="utf-8") == "previous output\n"
        assert not list(tmp_path.glob(f".{output.name}.*.tmp"))

    def test_failed_summary_render_preserves_existing_files(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Rendering must complete before either previous transaction file is displaced."""
        variants, transcripts, fasta = _inputs(tmp_path)
        output = tmp_path / "annotated.vcf"
        summary = output.with_name(f"{output.name}_summary.txt")
        output.write_text("previous output\n", encoding="utf-8")
        summary.write_text("previous summary\n", encoding="utf-8")

        def fail_render(_summary: object, _format_name: str) -> str:
            raise ValueError("summary rendering failed")

        monkeypatch.setattr(run_module, "render_reportable", fail_render)

        code = main(
            [
                "run",
                "--input_file",
                str(variants),
                "--gff",
                str(transcripts),
                "--fasta",
                str(fasta),
                "--output_file",
                str(output),
                "--stats_format",
                "text",
                "--force_overwrite",
            ]
        )

        assert code == ExitCode.USAGE
        assert output.read_text(encoding="utf-8") == "previous output\n"
        assert summary.read_text(encoding="utf-8") == "previous summary\n"
        assert not list(tmp_path.glob(".*.tmp"))

    def test_failed_summary_commit_restores_existing_files(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """A second replacement failure rolls the first destination back to the same generation."""
        variants, transcripts, fasta = _inputs(tmp_path)
        output = tmp_path / "annotated.vcf"
        summary = output.with_name(f"{output.name}_summary.txt")
        output.write_text("previous output\n", encoding="utf-8")
        summary.write_text("previous summary\n", encoding="utf-8")
        output.chmod(0o640)
        summary.chmod(0o600)
        replace = Path.replace

        def fail_summary_commit(path: Path, target: str | Path) -> Path:
            if path.suffix == ".tmp" and Path(target) == summary:
                raise OSError("summary replacement failed")
            return replace(path, target)

        monkeypatch.setattr(Path, "replace", fail_summary_commit)

        code = main(
            [
                "run",
                "--input_file",
                str(variants),
                "--gff",
                str(transcripts),
                "--fasta",
                str(fasta),
                "--output_file",
                str(output),
                "--stats_format",
                "text",
                "--force_overwrite",
            ]
        )

        assert code == ExitCode.USAGE
        assert output.read_text(encoding="utf-8") == "previous output\n"
        assert summary.read_text(encoding="utf-8") == "previous summary\n"
        assert stat.S_IMODE(output.stat().st_mode) == 0o640
        assert stat.S_IMODE(summary.stat().st_mode) == 0o600
        assert not list(tmp_path.glob(".*.tmp"))
        assert not list(tmp_path.glob(".*.backup"))

    @pytest.mark.parametrize("raced_destination", ["output", "summary"])
    def test_non_forced_commit_does_not_replace_a_file_created_during_the_run(
        self,
        raced_destination: str,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Atomically preserve a destination that appears after preflight validation."""

        output = tmp_path / "annotated.vcf"
        summary = output.with_name(f"{output.name}_summary.txt")
        raced_path = output if raced_destination == "output" else summary
        config = Config(
            source=SourceConfig(input_file=tmp_path / "variants.vcf"),
            run=RunConfig(output_file=output, stats_format="text"),
        )

        def create_racing_file(
            _runner: Runner,
            temporary_output: Path | None,
        ) -> tuple[int, RunSummaryCollector]:
            assert temporary_output is not None
            temporary_output.write_text("new output\n", encoding="utf-8")
            raced_path.write_text("racing writer\n", encoding="utf-8")
            return 0, RunSummaryCollector()

        monkeypatch.setattr(Runner, "_run_buffered", create_racing_file)

        code = Runner(config).run()

        assert code == ExitCode.USAGE
        assert raced_path.read_text(encoding="utf-8") == "racing writer\n"
        if raced_path == summary:
            assert not output.exists()
        assert not list(tmp_path.glob(".*.tmp"))
        assert not list(tmp_path.glob(".*.backup"))

    def test_summary_and_output_must_have_distinct_destinations(self, tmp_path: Path) -> None:
        """One path cannot represent both members of the output transaction."""
        destination = tmp_path / "result.txt"

        code = main(
            [
                "run",
                "--input_file",
                "variants.vcf",
                "--output_file",
                str(destination),
                "--stats_file",
                str(destination),
            ]
        )

        assert code == ExitCode.USAGE
        assert not destination.exists()

    def test_stdout_run_still_reports_its_summary(
        self,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Streaming the VCF does not discard the completed collector's log summary."""
        messages: list[str] = []
        config = Config(
            source=SourceConfig(input_file=Path("variants.vcf")),
            run=RunConfig(output_file=None),
        )

        def finish_without_output(
            _runner: Runner,
            destination: Path | None,
        ) -> tuple[int, RunSummaryCollector]:
            assert destination is None
            return 0, RunSummaryCollector()

        def capture_log(message: str, *arguments: object) -> None:
            messages.append(message % arguments)

        monkeypatch.setattr(Runner, "_run_buffered", finish_without_output)
        monkeypatch.setattr(run_module.logger, "info", capture_log)

        assert Runner(config).run() == ExitCode.OK
        assert any(message.startswith("run summary:") for message in messages)

    @pytest.mark.parametrize("spelling", ["STDOUT", "stdout", "Stdout", "-"])
    def test_stdout_spellings_stream_instead_of_writing_a_file(
        self,
        spelling: str,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """``--output_file STDOUT`` (any case) and ``-`` stream the VCF and create no file."""
        variants, transcripts, fasta = _inputs(tmp_path)
        working_directory = tmp_path / "cwd"
        working_directory.mkdir()
        monkeypatch.chdir(working_directory)

        code = main(
            [
                "run",
                "--input_file",
                str(variants),
                "--gff",
                str(transcripts),
                "--fasta",
                str(fasta),
                "--output_file",
                spelling,
            ]
        )

        assert code == ExitCode.OK
        captured = capsys.readouterr()
        # Logs go to stderr, so a piped stream carries the VCF and nothing else.
        assert "[INFO]" not in captured.out
        records = [line for line in captured.out.splitlines() if not line.startswith("#")]
        assert len(records) == 4
        assert not list(working_directory.iterdir())

    def test_stdout_with_an_explicit_stats_file_still_writes_the_summary(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Streaming the annotations does not stop a requested summary going to its file."""
        summary = tmp_path / "summary.txt"

        code, _output = _run(tmp_path, "--output_file", "STDOUT", "--stats_file", str(summary))

        assert code == ExitCode.OK
        assert summary.is_file()
        assert "#CHROM" in capsys.readouterr().out

    def test_every_input_variant_is_written(self, tmp_path: Path) -> None:
        """A variant that overlaps nothing still has an answer, so the output cannot be shorter
        than the input. Dropping it would lose exactly the intergenic variants."""
        _code, output = _run(tmp_path)

        assert len(self._data_lines(output)) == 4

    def test_annotation_windows_preserve_rendered_records(self, tmp_path: Path) -> None:
        """Smaller in-worker windows change memory lifetime, not VCF annotations."""
        default_path = tmp_path / "default"
        windowed_path = tmp_path / "windowed"
        default_path.mkdir()
        windowed_path.mkdir()
        _default_code, default_output = _run(default_path)
        _windowed_code, windowed_output = _run(
            windowed_path,
            "--annotation_window",
            "1",
        )

        assert self._data_lines(windowed_output) == self._data_lines(default_output)

    def test_a_variant_overlapping_nothing_is_intergenic(self, tmp_path: Path) -> None:
        """Reported, not omitted."""
        _code, output = _run(tmp_path)

        far = next(line for line in self._data_lines(output) if "rs_far" in line)
        assert "intergenic_variant" in far

    def test_the_input_record_is_echoed_rather_than_rebuilt(self, tmp_path: Path) -> None:
        """Rebuilding re-derives the padded REF from the reference sequence, which can disagree
        with what the input declared -- the output would contradict its own input."""
        _code, output = _run(tmp_path)

        indel = next(line for line in self._data_lines(output) if "rs_indel" in line)
        columns = indel.split("\t")
        assert columns[3] == "C"
        assert columns[4] == "CAT"
        assert columns[6] == "PASS"

    def test_format_and_sample_columns_are_preserved(self, tmp_path: Path) -> None:
        """Annotation changes INFO without normalising FORMAT or sample values."""
        variants, transcripts, fasta = _inputs(tmp_path)
        variants.write_text(
            "##fileformat=VCFv4.3\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A\tSAMPLE_B\n"
            "chr1\t1500\trs_genotyped\tA\tC\t10.00\tPASS\tSOURCE=raw\tGT:DP\t0|1:007\t1/1:042\n",
            encoding="utf-8",
        )
        output = tmp_path / "annotated.vcf"

        code = main(
            [
                "run",
                "--input_file",
                str(variants),
                "--gff",
                str(transcripts),
                "--fasta",
                str(fasta),
                "--output_file",
                str(output),
                "--no_stats",
            ]
        )

        assert code == ExitCode.OK
        text = output.read_text(encoding="utf-8")
        assert "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A\tSAMPLE_B" in text
        columns = self._data_lines(output)[0].split("\t")
        assert columns[:7] == ["chr1", "1500", "rs_genotyped", "A", "C", "10.00", "PASS"]
        assert columns[7].startswith("SOURCE=raw;CSQ=")
        assert columns[8:] == ["GT:DP", "0|1:007", "1/1:042"]

    def test_unindexed_fasta_still_produces_coding_consequences(self, tmp_path: Path) -> None:
        """Sequential FASTA fallback derives contig names and lengths from loaded records."""
        variants, transcripts, fasta = _inputs(tmp_path)
        Path(f"{fasta}.fai").unlink()
        output = tmp_path / "annotated.vcf"

        code = main(
            [
                "run",
                "--input_file",
                str(variants),
                "--gff",
                str(transcripts),
                "--fasta",
                str(fasta),
                "--output_file",
                str(output),
                "--no_stats",
            ]
        )

        coding = next(line for line in self._data_lines(output) if "rs_coding" in line)
        assert code == ExitCode.OK
        assert "CSQ=" in coding
        assert "intergenic_variant" not in coding

    def test_gff_and_fasta_region_spellings_are_resolved_independently(self, tmp_path: Path) -> None:
        """Match VCF/GFF name ``1`` on FASTA slice ``chr1`` and echo the VCF spelling."""
        variants, transcripts, fasta = _inputs(tmp_path)
        variants.write_text(
            "##fileformat=VCFv4.3\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "1\t1500\trs_coding\tA\tC\t.\tPASS\t.\n",
            encoding="utf-8",
        )
        transcripts.write_text(
            "##gff-version 3\n##sequence-region 1 1 40000\n" + GFF3.replace("chr1", "1"),
            encoding="utf-8",
        )
        output = tmp_path / "annotated.vcf"

        code = main(
            [
                "run",
                "--input_file",
                str(variants),
                "--gff",
                str(transcripts),
                "--fasta",
                str(fasta),
                "--output_file",
                str(output),
                "--stats_format",
                "json",
            ]
        )

        coding = self._data_lines(output)[0]
        summary_path = output.with_name(f"{output.name}_summary.json")
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        chromosomes = next(
            section for section in summary["sections"] if section["identifier"] == "chromosomes"
        )
        assert code == ExitCode.OK
        assert coding.split("\t")[0] == "1"
        assert "CSQ=" in coding
        assert "intergenic_variant" not in coding
        assert chromosomes["slices"] == [{"label": "1", "value": 1}]

    def test_a_multi_allelic_record_stays_one_line(self, tmp_path: Path) -> None:
        """One record in, one record out, with an annotation per alternate allele."""
        _code, output = _run(tmp_path)

        multi = next(line for line in self._data_lines(output) if "rs_multi" in line)
        consequences = multi.split("CSQ=")[1].split(",")
        assert len(consequences) == 2

    def test_a_variant_on_a_contig_the_fasta_lacks_is_written_through(self, tmp_path: Path) -> None:
        """It must be reported, not dropped, and not annotated.

        ClinVar names seven RefSeq scaffolds the Ensembl toplevel FASTA does not carry. Those 51
        variants used to abort a run of ~2.9 million. They now pass through carrying no CSQ, because
        nothing can be annotated against a region we have no sequence for -- and calling them
        intergenic would assert something we did not check.
        """
        variants = tmp_path / "variants.vcf"
        variants.write_text(
            "##fileformat=VCFv4.3\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t1500\trs_known\tA\tC\t.\tPASS\t.\n"
            "NT_113889.1\t54175\trs_absent\tG\tA\t.\tPASS\t.\n",
            encoding="utf-8",
        )
        transcripts = tmp_path / "transcripts.gff3"
        transcripts.write_text(GFF3, encoding="utf-8")
        fasta = tmp_path / "ref.fa"
        fasta.write_text(">chr1\n" + "ACGT" * 10000 + "\n", encoding="utf-8")
        fasta.with_suffix(".fa.fai").write_text("chr1\t40000\t6\t40000\t40001\n", encoding="utf-8")
        output = tmp_path / "annotated.vcf"

        code = main(
            [
                "run",
                "--input_file",
                str(variants),
                "--gff",
                str(transcripts),
                "--fasta",
                str(fasta),
                "--output_file",
                str(output),
                "--no_stats",
            ]
        )

        assert code == ExitCode.OK
        lines = self._data_lines(output)
        assert len(lines) == 2, "the absent-contig variant must not be dropped"

        absent = next(line for line in lines if "rs_absent" in line)
        assert "CSQ=" not in absent
        assert absent.split("\t")[:5] == ["NT_113889.1", "54175", "rs_absent", "G", "A"]

        known = next(line for line in lines if "rs_known" in line)
        assert "CSQ=" in known, "the ordinary variant must still be annotated"

    def test_absent_reference_region_never_enters_overlap_annotation(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Do not request placeholder-length sequence for a GXF transcript on an absent contig."""
        variants = tmp_path / "variants.vcf"
        variants.write_text(
            "##fileformat=VCFv4.3\n"
            '##INFO=<ID=CSQ,Number=.,Type=String,Description="Old incompatible schema">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "NT_113889.1\t54175\trs_absent\tG\tA\t.\tPASS\tSOURCE=raw;CSQ=T|OLD;END=55000\n",
            encoding="utf-8",
        )
        transcripts = tmp_path / "transcripts.gff3"
        transcripts.write_text(
            "NT_113889.1\tensembl\tgene\t54000\t55000\t.\t+\t.\tID=gene:absent\n"
            "NT_113889.1\tensembl\tmRNA\t54000\t55000\t.\t+\t.\t"
            "ID=transcript:absent;Parent=gene:absent;biotype=protein_coding\n"
            "NT_113889.1\tensembl\texon\t54000\t55000\t.\t+\t.\t"
            "ID=exon:absent;Parent=transcript:absent\n"
            "NT_113889.1\tensembl\tCDS\t54000\t55000\t.\t+\t0\t"
            "ID=cds:absent;Parent=transcript:absent;protein_id=protein:absent\n",
            encoding="utf-8",
        )
        fasta = tmp_path / "ref.fa"
        fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
        fasta.with_suffix(".fa.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
        output = tmp_path / "annotated.vcf"

        def reject_sequence_fetch(
            _source: object,
            sequence_region: str,
            _start: int,
            _end: int,
        ) -> str:
            raise AssertionError(f"unexpected sequence fetch for {sequence_region}")

        monkeypatch.setattr(FastaSequenceSource, "_fetch_sequence_chunk", reject_sequence_fetch)

        code = main(
            [
                "run",
                "--input_file",
                str(variants),
                "--gff",
                str(transcripts),
                "--fasta",
                str(fasta),
                "--output_file",
                str(output),
                "--no_stats",
            ]
        )

        line = self._data_lines(output)[0]
        text = output.read_text(encoding="utf-8")
        assert code == ExitCode.OK
        assert "CSQ=" not in line
        assert "Old incompatible schema" not in text
        assert "Consequence annotations from Ensembl VEP" in text
        assert line.split("\t")[7] == "SOURCE=raw;END=55000"
        assert line.split("\t")[:5] == ["NT_113889.1", "54175", "rs_absent", "G", "A"]

    def test_recurrent_sequence_regions_keep_input_order(self, tmp_path: Path) -> None:
        """A valid region order A, B, A completes without sorting or partial output."""
        variants, transcripts, fasta = _inputs(tmp_path)
        variants.write_text(
            "##fileformat=VCFv4.3\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t1500\tfirst\tA\tC\t.\tPASS\t.\n"
            "chr2\t20\tsecond\tA\tG\t.\tPASS\t.\n"
            "chr1\t1700\tthird\tA\tT\t.\tPASS\t.\n",
            encoding="utf-8",
        )
        output = tmp_path / "annotated.vcf"

        code = main(
            [
                "run",
                "--input_file",
                str(variants),
                "--gff",
                str(transcripts),
                "--fasta",
                str(fasta),
                "--output_file",
                str(output),
                "--stats_format",
                "text",
            ]
        )

        identifiers = [line.split("\t")[2] for line in self._data_lines(output)]
        assert code == ExitCode.OK
        assert identifiers == ["first", "second", "third"]
        assert output.with_name(f"{output.name}_summary.txt").exists()

    def test_the_summary_is_written_beside_the_output(self, tmp_path: Path) -> None:
        """Default-on, as Perl VEP has it, and named after the output it describes."""
        _code, output = _run(tmp_path, "--stats_format", "text")

        assert output.with_name(f"{output.name}_summary.txt").exists()

    def test_no_stats_suppresses_the_summary(self, tmp_path: Path) -> None:
        """Cheap as it is, a caller piping output may not want a second file appearing."""
        _code, output = _run(tmp_path, "--no_stats", "--stats_format", "text")

        assert not output.with_name(f"{output.name}_summary.txt").exists()


class TestOutputFileSpelling:
    """How every configuration layer reads the output destination."""

    @pytest.mark.parametrize("spelling", ["STDOUT", "stdout", "Stdout", "-", None])
    def test_stdout_spellings_become_none(self, spelling: str | None) -> None:
        """A direct caller gets the same stdout mapping as the command line."""
        assert RunConfig(output_file=spelling).output_file is None  # type: ignore[arg-type]

    def test_a_path_is_kept(self) -> None:
        """Anything else is an ordinary path, including a name that merely contains stdout."""
        assert RunConfig(output_file=Path("stdout.vcf")).output_file == Path("stdout.vcf")

    def test_a_directory_destination_is_reported_as_one(self, tmp_path: Path) -> None:
        """``--output_file .`` names a directory, which is a mistake, not a request to stream."""
        variants, transcripts, fasta = _inputs(tmp_path)

        code = main(
            [
                "run",
                "--input_file",
                str(variants),
                "--gff",
                str(transcripts),
                "--fasta",
                str(fasta),
                "--output_file",
                str(tmp_path),
                "--no_stats",
            ]
        )

        assert code == ExitCode.USAGE

    def test_an_empty_destination_is_not_stdout(self) -> None:
        """Only the two documented spellings stream; an empty value is a path the run will reject.

        An empty option is a mistake rather than a request to pipe, and the working directory is
        the path it names.
        """
        assert RunConfig(output_file="").output_file == Path(".")  # type: ignore[arg-type]

    def test_the_default_is_a_file_in_the_working_directory(self) -> None:
        """Omitting the option still writes variant_effect_output.txt."""
        assert RunConfig().output_file == Path("variant_effect_output.txt")

    def test_the_environment_layer_maps_stdout_to_none(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """An exported ENSEMBL_VEP_RUN_OUTPUT_FILE=STDOUT streams rather than writing a file."""
        variants, transcripts, fasta = _inputs(tmp_path)
        working_directory = tmp_path / "cwd"
        working_directory.mkdir()
        monkeypatch.chdir(working_directory)
        monkeypatch.setenv("ENSEMBL_VEP_RUN_OUTPUT_FILE", "STDOUT")

        code = main(["run", "--input_file", str(variants), "--gff", str(transcripts), "--fasta", str(fasta)])

        assert code == ExitCode.OK
        assert "CSQ=" in capsys.readouterr().out
        assert not list(working_directory.iterdir())


class TestRunnerContract:
    """Properties every registered runner must hold."""

    @pytest.mark.parametrize("name", list(available_runners()))
    def test_a_runner_declares_sections_that_exist(self, name: str) -> None:
        """A runner naming a missing section would fail only at invocation."""
        runner = get_runner(name)
        assert runner.sections
        Config.section_types_for(runner.sections)

    @pytest.mark.parametrize("name", list(available_runners()))
    def test_a_runner_reads_the_logging_section(self, name: str) -> None:
        """setup_logging is applied from it before the runner is constructed."""
        assert "logging" in get_runner(name).sections

    @pytest.mark.parametrize("name", list(available_runners()))
    def test_a_runner_describes_itself(self, name: str) -> None:
        """The description is the subcommand's help text."""
        assert get_runner(name).description

    @pytest.mark.parametrize("name", list(available_runners()))
    def test_a_runner_name_matches_its_registration(self, name: str) -> None:
        """A mismatch would make the usage listing disagree with the dispatch."""
        assert get_runner(name).name == name

    @pytest.mark.parametrize("name", list(available_runners()))
    def test_a_runner_reads_no_two_sections_that_share_a_flag(self, name: str) -> None:
        """The guard keeping the concern-based section split honest."""
        check_no_flag_collisions(Config.section_types_for(get_runner(name).sections))
