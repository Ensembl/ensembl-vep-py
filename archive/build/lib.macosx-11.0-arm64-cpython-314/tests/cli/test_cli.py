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
"""Tests for the ``vep`` command line."""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

import pytest

import ensembl.vep.annotators  # noqa: F401  pylint: disable=unused-import

# ^ registers the annotator names the --annotator option validates against.
from ensembl.vep.cli import DEFAULT_COMMAND, OTHER_SECTION, main, resolve_command
from ensembl.vep.config import Config, ConfigSource
from ensembl.vep.runner import ExitCode, available_runners, get_runner, register_runner


class TestCommandResolution:
    """The first argument selects the command; anything else falls back to the default."""

    @pytest.mark.parametrize(
        "argv,expected",
        [
            ([], (DEFAULT_COMMAND, [])),
            (["--input_file", "x.vcf"], (DEFAULT_COMMAND, ["--input_file", "x.vcf"])),
            (["-v"], (DEFAULT_COMMAND, ["-v"])),
            (["run"], ("run.default", [])),
            (["run", "--input_file", "x.vcf"], ("run.default", ["--input_file", "x.vcf"])),
            (["run", "default", "--input_file", "x.vcf"], ("run.default", ["--input_file", "x.vcf"])),
        ],
    )
    def test_resolution(self, argv: list[str], expected: tuple[str, list[str]]) -> None:
        """Only a leading flag, or nothing at all, means the default command."""
        assert resolve_command(argv) == expected

    @pytest.mark.parametrize("argv", [["nonsense"], ["frobnicate"], ["frobnicate", "-h"], ["Run"]])
    def test_an_unknown_leading_word_is_not_the_default_command(self, argv: list[str]) -> None:
        """It can only be a command attempt, so it must not quietly become ``vep run``."""
        assert resolve_command(argv) == ("", argv)

    def test_the_argument_list_is_not_mutated(self) -> None:
        """main() hands its own list in; resolution must not consume from it."""
        argv = ["run", "--input_file", "x.vcf"]
        resolve_command(argv)
        assert argv == ["run", "--input_file", "x.vcf"]


class TestUnknownCommand:
    """A leading word that is not a command is a usage error, not a silent fallback."""

    def test_it_exits_with_a_usage_code(self) -> None:
        """Distinguishable from a command that ran and failed."""
        assert main(["frobnicate", "-h"]) == ExitCode.USAGE

    def test_it_names_the_word_that_was_not_understood(self, capsys: pytest.CaptureFixture) -> None:
        """Reporting it against ``vep run`` would name a command the user never typed."""
        main(["nonsense"])
        assert "vep: 'nonsense' is not a vep command. See 'vep --help'." in capsys.readouterr().err

    def test_it_lists_the_available_commands(self, capsys: pytest.CaptureFixture) -> None:
        """The next thing anyone wants to know, as the same page ``vep -h`` prints."""
        main(["nonsense", "-h"])
        error = capsys.readouterr().err
        main(["-h"])
        assert capsys.readouterr().out in error

    def test_a_group_without_a_verb_lists_its_verbs(self, capsys: pytest.CaptureFixture) -> None:
        """A group is not itself runnable when it has several verbs."""
        for verb in ("alpha", "beta"):
            register_runner(f"sample.{verb}")(_sample_runner(f"sample.{verb}"))
        assert main(["sample"]) == ExitCode.USAGE
        error = capsys.readouterr().err
        assert "vep sample: needs a verb." in error
        assert "usage: vep sample <verb> [<options>]" in error
        assert "   alpha" in error
        assert "   beta" in error

    def test_a_group_verb_resolves(self) -> None:
        """``vep sample alpha`` reaches the runner registered under sample.alpha."""
        register_runner("sample.alpha")(_sample_runner("sample.alpha"))
        assert resolve_command(["sample", "alpha", "-x"]) == ("sample.alpha", ["-x"])

    def test_a_group_with_one_verb_is_runnable_by_its_domain(self) -> None:
        """``vep sample`` is unambiguous when the group holds a single runner."""
        register_runner("sample.only")(_sample_runner("sample.only"))
        assert resolve_command(["sample"]) == ("sample.only", [])


class TestEntryPoint:
    """main() dispatches and returns the runner's exit code."""

    def test_run_is_reached_without_a_subcommand(self) -> None:
        """``vep --input_file x.vcf`` still parses, as Perl VEP users expect.

        It reaches the run path and stops there for want of transcripts to pair against, which is
        a usage error: the dispatch is what this asserts, not the annotation.
        """
        assert main(["--input_file", "x.vcf"]) == ExitCode.USAGE

    def test_a_command_help_exits_zero(self) -> None:
        """Asking for a command's options is not an error."""
        with pytest.raises(SystemExit) as raised:
            main(["--input_file", "x.vcf", "--help"])
        assert raised.value.code == ExitCode.OK

    def test_an_options_page_names_the_command_as_typed(self, capsys: pytest.CaptureFixture) -> None:
        """``vep run default -h`` is titled with the words that reach it, not the registry name."""
        with pytest.raises(SystemExit):
            main(["run", "default", "-h"])
        assert capsys.readouterr().out.startswith("usage: vep run default ")


class TestHelp:
    """The command listings ``vep -h`` and ``vep <group> -h`` print, laid out as ``git -h`` does."""

    @pytest.mark.parametrize("flag", ["-h", "--help"])
    def test_top_level_help_lists_commands_not_run_options(
        self, flag: str, capsys: pytest.CaptureFixture
    ) -> None:
        """``vep -h`` is about vep, not about the default command it would otherwise reach."""
        assert main([flag]) == ExitCode.OK
        output = capsys.readouterr().out
        assert output.startswith("usage: vep [-h | --help] [--config <file>] <command> [<verb>] [<options>]")
        assert "--level" not in output

    def test_commands_are_listed_under_headings_and_aligned(self, capsys: pytest.CaptureFixture) -> None:
        """Each heading is flush left; its commands are indented three spaces into one column."""
        main(["-h"])
        lines = capsys.readouterr().out.splitlines()
        heading = lines.index("annotate variants")
        default_line, haplo_line = lines[heading + 1], lines[heading + 2]
        assert default_line.startswith("   run default ")
        assert haplo_line.startswith("   run haplo ")
        description_column = default_line.index("annotate variants against")
        assert haplo_line.index("annotate phased") == description_column

    def test_every_command_is_listed(self, capsys: pytest.CaptureFixture) -> None:
        """Nothing registered is left out of the listing."""
        main(["-h"])
        output = capsys.readouterr().out
        for name in available_runners():
            assert f"   {name.replace('.', ' ', 1)} " in output

    def test_a_group_help_lists_only_its_verbs(self, capsys: pytest.CaptureFixture) -> None:
        """``vep corpus -h`` is the same page narrowed to one group, instead of an error."""
        assert main(["corpus", "-h"]) == ExitCode.OK
        output = capsys.readouterr().out
        assert output.startswith("usage: vep corpus <verb> [<options>]")
        assert "   describe " in output
        assert "run default" not in output

    def test_a_group_with_a_default_verb_says_so(self, capsys: pytest.CaptureFixture) -> None:
        """``vep run -h`` lists the verbs, and names the one ``vep run`` reaches without one."""
        assert main(["run", "--help"]) == ExitCode.OK
        assert "Without a verb 'vep run' runs 'vep run default'." in capsys.readouterr().out

    def test_a_group_no_heading_covers_is_listed_under_other_commands(
        self, capsys: pytest.CaptureFixture
    ) -> None:
        """A newly registered group still appears, under a catch-all heading."""
        for verb in ("alpha", "beta"):
            register_runner(f"sample.{verb}")(_sample_runner(f"sample.{verb}"))
        main(["-h"])
        lines = capsys.readouterr().out.splitlines()
        assert lines[lines.index(OTHER_SECTION) + 1].startswith("   sample alpha ")

    def test_a_single_command_help_still_shows_its_options(self) -> None:
        """``explain`` is a command, not a group, so ``-h`` after it reaches argparse."""
        with pytest.raises(SystemExit) as raised:
            main(["explain", "-h"])
        assert raised.value.code == ExitCode.OK

    def test_an_unknown_option_is_a_usage_error(self) -> None:
        """argparse's own exit code, kept as the usage convention."""
        with pytest.raises(SystemExit) as raised:
            main(["--not-an-option"])
        assert raised.value.code == ExitCode.USAGE

    def test_an_invalid_section_value_is_a_usage_error(self, capsys: pytest.CaptureFixture) -> None:
        """Section validation is reported at the command boundary without a traceback."""
        assert main(["--processes", "0"]) == ExitCode.USAGE

        error = capsys.readouterr().err
        assert "processes must be at least one" in error
        assert "Traceback" not in error

    def test_a_config_file_is_read(self, tmp_path: Path) -> None:
        """--config is taken before the other layers are merged."""
        path = tmp_path / "vep.toml"
        path.write_text("[logging]\nlevel = 'error'\n", encoding="utf-8")
        assert main(["--config", str(path), "--input_file", "x.vcf"]) == ExitCode.USAGE

    def test_a_config_file_may_be_given_with_an_equals_sign(self, tmp_path: Path) -> None:
        """--config=FILE is the other spelling argparse users expect."""
        path = tmp_path / "vep.toml"
        path.write_text("[logging]\nlevel = 'error'\n", encoding="utf-8")
        assert main([f"--config={path}", "--input_file", "x.vcf"]) == ExitCode.USAGE


class TestConfigFlag:
    """``--config`` is taken off the whole command line, so it may stand before the command."""

    @pytest.mark.parametrize(
        "command,table,expected",
        [
            pytest.param(
                ["catalogue", "show"], "[catalogue_show]\nsubject = 'facts'\n", "overlap_cds", id="show"
            ),
            pytest.param(["corpus", "describe"], "[logging]\nlevel = 'error'\n", "# Corpus", id="describe"),
        ],
    )
    @pytest.mark.parametrize("spelling", ["separate", "equals"])
    def test_it_may_precede_the_command(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture,
        command: list[str],
        table: str,
        expected: str,
        spelling: str,
    ) -> None:
        """``vep --config f.toml <group> <verb>`` reaches the verb and reads the file."""
        path = tmp_path / "vep.toml"
        path.write_text(table, encoding="utf-8")
        flag = ["--config", str(path)] if spelling == "separate" else [f"--config={path}"]

        assert main([*flag, *command]) == ExitCode.OK
        assert expected in capsys.readouterr().out

    def test_it_may_follow_the_command(self, tmp_path: Path, capsys: pytest.CaptureFixture) -> None:
        """The position after the command keeps working."""
        path = tmp_path / "vep.toml"
        path.write_text("[catalogue_show]\nsubject = 'facts'\n", encoding="utf-8")
        assert main(["catalogue", "show", "--config", str(path)]) == ExitCode.OK
        assert "overlap_cds" in capsys.readouterr().out

    def test_it_is_not_mistaken_for_a_leading_flag_of_the_default_command(self, tmp_path: Path) -> None:
        """The old dispatch sent ``--config f corpus describe`` to ``run default`` and exited 2."""
        path = tmp_path / "vep.toml"
        path.write_text("[logging]\nlevel = 'error'\n", encoding="utf-8")
        assert main(["--config", str(path), "corpus", "describe"]) != ExitCode.USAGE

    @pytest.mark.parametrize("argv", [["--config"], ["corpus", "describe", "--config"], ["--config="]])
    def test_a_flag_without_a_file_is_a_usage_error(
        self, argv: list[str], capsys: pytest.CaptureFixture
    ) -> None:
        """Reported cleanly rather than crashing or being passed on to a command."""
        assert main(argv) == ExitCode.USAGE
        error = capsys.readouterr().err
        assert "--config: expected one argument" in error
        assert "Traceback" not in error

    @pytest.mark.parametrize("position", ["before", "after"])
    def test_a_missing_file_is_a_usage_error_naming_the_path(
        self, tmp_path: Path, capsys: pytest.CaptureFixture, position: str
    ) -> None:
        """A path that does not exist is the caller's mistake, not a traceback."""
        missing = tmp_path / "absent.toml"
        flag = ["--config", str(missing)]
        argv = [*flag, "corpus", "describe"] if position == "before" else ["corpus", "describe", *flag]

        assert main(argv) == ExitCode.USAGE
        error = capsys.readouterr().err
        assert f"cannot read config file '{missing}'" in error
        assert "Traceback" not in error

    def test_a_directory_is_a_usage_error(self, tmp_path: Path, capsys: pytest.CaptureFixture) -> None:
        """Any unreadable path is reported the same way as a missing one."""
        assert main(["--config", str(tmp_path), "corpus", "describe"]) == ExitCode.USAGE
        assert f"cannot read config file '{tmp_path}'" in capsys.readouterr().err

    def test_help_after_the_flag_lists_the_commands(
        self, tmp_path: Path, capsys: pytest.CaptureFixture
    ) -> None:
        """The flag is removed before the help check, so ``-h`` is still leading."""
        assert main(["--config", str(tmp_path / "vep.toml"), "-h"]) == ExitCode.OK
        assert capsys.readouterr().out.startswith("usage: vep [-h | --help]")


class TestPackageGuards:
    """Invariants that keep one option to one declaration."""

    def test_add_argument_appears_only_where_the_parser_is_derived(self) -> None:
        """One declaration per option means one place that builds arguments."""
        package = Path(__file__).resolve().parents[2] / "ensembl" / "vep"
        allowed = {package / "config" / "parser" / "cli.py", package / "cli.py"}
        offenders = [
            path
            for path in package.rglob("*.py")
            if path not in allowed and "add_argument(" in path.read_text(encoding="utf-8")
        ]
        assert not offenders, f"argparse used outside the derivation: {offenders}"

    def test_no_module_calls_basic_config(self) -> None:
        """It configures the root logger process-wide; setup_logging is the entry point."""
        package = Path(__file__).resolve().parents[2] / "ensembl" / "vep"
        offenders = [
            path for path in package.rglob("*.py") if "basicConfig" in path.read_text(encoding="utf-8")
        ]
        assert not offenders, f"logging.basicConfig used in {offenders}"

    def test_only_the_expected_sections_are_unread(self) -> None:
        """A section no command reads is dead configuration.

        The set shrinks as each runner lands; update it there rather than here.
        """
        read = {section for name in available_runners() for section in get_runner(name).sections}
        assert set(Config.section_types) - read == set()


class TestConfigPrecedenceThroughMain:
    """What reaches a runner's Config, and which layer it came from.

    Asserting on a command's output cannot tell a value read from the file apart from the same
    value defaulted, so these capture the Config the runner is handed.
    """

    @staticmethod
    def _capturing_runner(captured: list[Config]) -> type:
        """Return a runner that records the Config it is given and does nothing else."""

        class Capturing:
            """A runner registered only by these tests."""

            sections: ClassVar[tuple[str, ...]] = ("logging", "execution")
            description: ClassVar[str] = "records the configuration it is handed"

            def __init__(self, config: Config) -> None:
                captured.append(config)

            def run(self) -> ExitCode:
                """Do nothing successfully."""
                return ExitCode.OK

        Capturing.name = "sample.capturing"  # type: ignore[attr-defined]
        return Capturing

    def _config_from(self, argv: list[str], captured: list[Config]) -> Config:
        """Run the capturing runner and return the Config it received."""
        register_runner("sample.capturing")(self._capturing_runner(captured))
        assert main(argv) == ExitCode.OK
        assert len(captured) == 1
        return captured[0]

    def test_a_file_value_reaches_the_runner_and_is_attributed_to_the_file(self, tmp_path: Path) -> None:
        """The value must arrive, and ``source_of`` must say the file supplied it."""
        path = tmp_path / "vep.toml"
        path.write_text("[logging]\nlevel = 'error'\n", encoding="utf-8")
        captured: list[Config] = []

        config = self._config_from(["--config", str(path), "sample", "capturing"], captured)

        assert config.logging.log_level == "error"
        assert config.source_of("logging", "log_level") is ConfigSource.FILE

    def test_the_command_line_beats_the_file(self, tmp_path: Path) -> None:
        """Precedence is default < file < environment < command line."""
        path = tmp_path / "vep.toml"
        path.write_text("[logging]\nlevel = 'error'\n", encoding="utf-8")
        captured: list[Config] = []

        config = self._config_from(
            ["--config", str(path), "sample", "capturing", "--log_level", "debug"], captured
        )

        assert config.logging.log_level == "debug"
        assert config.source_of("logging", "log_level") is ConfigSource.COMMAND_LINE

    def test_the_environment_beats_the_file_and_loses_to_the_command_line(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """The middle layer, which no other test pins through main()."""
        path = tmp_path / "vep.toml"
        path.write_text("[execution]\nprocesses = 2\n", encoding="utf-8")
        monkeypatch.setenv("ENSEMBL_VEP_EXECUTION_PROCESSES", "3")
        captured: list[Config] = []

        config = self._config_from(["--config", str(path), "sample", "capturing"], captured)

        assert config.execution.processes == 3
        assert config.source_of("execution", "processes") is ConfigSource.ENVIRONMENT

        captured.clear()
        config = self._config_from(
            ["--config", str(path), "sample", "capturing", "--processes", "4"], captured
        )

        assert config.execution.processes == 4
        assert config.source_of("execution", "processes") is ConfigSource.COMMAND_LINE

    def test_a_default_is_attributed_to_no_layer(self) -> None:
        """An option nobody set reports DEFAULT, which is what makes the assertions above mean something."""
        captured: list[Config] = []

        config = self._config_from(["sample", "capturing"], captured)

        assert config.logging.log_level == "info"
        assert config.source_of("logging", "log_level") is ConfigSource.DEFAULT


def _sample_runner(name: str) -> type:
    """Return a runner class registered only by these tests."""

    class Sample:
        """A runner registered only by these tests."""

        sections: ClassVar[tuple[str, ...]] = ("logging",)
        description: ClassVar[str] = "a runner registered only by these tests"

        def __init__(self, config: Config) -> None:
            """Hold the configuration."""
            self.config = config

        def run(self) -> ExitCode:
            """Do nothing successfully."""
            return ExitCode.OK

    Sample.name = name  # type: ignore[attr-defined]
    return Sample
