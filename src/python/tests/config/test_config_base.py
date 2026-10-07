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
"""Tests for configuration loading."""

from __future__ import annotations

from importlib.util import find_spec
from pathlib import Path
from typing import Any

import pytest

from ensembl.vep.config import Config, ConfigSource, CorpusConfig, load
from ensembl.vep.config.option_spec import section_options

RUNNER_AVAILABLE = find_spec("ensembl.vep.runner") is not None


def get_runner(name: str) -> Any:
    """Return a runner when that later stack component is available."""
    runner = pytest.importorskip("ensembl.vep.runner")
    return runner.get_runner(name)


@pytest.fixture(name="toml_file")
def fixture_toml_file(tmp_path: Path) -> Path:
    """A TOML file setting two corpus options."""
    path = tmp_path / "vep.toml"
    path.write_text("[corpus]\nfuzz = 111\nstep = 222\n", encoding="utf-8")
    return path


class TestConfig:
    """The configuration root."""

    def test_build_section_reads_the_declared_options_once(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """The option table is computed once per section, not once per value."""
        calls: list[type] = []

        def counting(section_type: type) -> dict:
            calls.append(section_type)
            return section_options(section_type)

        monkeypatch.setattr("ensembl.vep.config.base.section_options", counting)
        built = Config.build_section(CorpusConfig, {"fuzz": 3, "step": 4, "nope": 5})
        assert isinstance(built, CorpusConfig)
        assert (built.fuzz, built.step) == (3, 4)
        assert len(calls) == 1

    def test_a_bad_environment_choice_fails_the_load(self) -> None:
        """An out-of-choice environment value is a ValueError out of load."""
        with pytest.raises(ValueError, match="logging.log_level: 'verbose'"):
            load(section_names=("logging",), environ={"ENSEMBL_VEP_LOGGING_LOG_LEVEL": "verbose"})

    def test_the_fields_are_the_section_registry(self) -> None:
        """The section registry matches the declared section attributes."""
        declared = {section.section_name: type(section) for section in Config().sections()}
        assert Config.section_types == declared

    def test_sections_are_reachable_by_name(self) -> None:
        """Section names are the lookup key."""
        assert Config().section("corpus").section_name == "corpus"

    def test_an_unknown_section_is_rejected(self) -> None:
        """Unknown section names fail fast."""
        with pytest.raises(KeyError):
            Config().section("nope")

    def test_as_mapping_covers_every_declared_option(self) -> None:
        """The mapping includes every declared option."""
        mapping = Config().as_mapping()
        for section in Config().sections():
            assert set(mapping[section.section_name]) == set(section_options(type(section)))

    def test_an_unset_origin_reads_as_default(self) -> None:
        """Unset values read as the default source."""
        assert Config().source_of("corpus", "fuzz") is ConfigSource.DEFAULT

    def test_the_command_is_recorded(self) -> None:
        """The command name is preserved."""
        assert load([], command="run", section_names=("corpus",), environ={}).command == "run"


@pytest.mark.skipif(not RUNNER_AVAILABLE, reason="runner component is not available")
class TestPerlSpellings:
    """The option names Perl VEP uses, as a run would be given them."""

    @staticmethod
    def _run(*argv: str) -> Config:
        """Load a run configuration from command-line arguments."""
        return load(list(argv), command="run.default", section_names=get_runner("run.default").sections)

    @pytest.mark.parametrize(
        ("argv", "section", "field_name", "expected"),
        [
            pytest.param(["--level", "debug"], "logging", "log_level", "debug", id="level"),
            pytest.param(
                ["--warning_file", "w.txt"], "logging", "log_file", Path("w.txt"), id="warning_file"
            ),
            pytest.param(["--input_format", "vcf"], "source", "format", "vcf", id="input_format"),
            pytest.param(["--flank", "250"], "source", "distance", (250,), id="flank"),
            pytest.param(["--force"], "run", "force_overwrite", True, id="force"),
            pytest.param(["--vcf"], "run", "output_format", "vcf", id="vcf"),
            pytest.param(["--json"], "run", "output_format", "json", id="json"),
            pytest.param(["--stats_html"], "run", "stats_format", "html", id="stats_html"),
            pytest.param(["--stats_text"], "run", "stats_format", "text", id="stats_text"),
            pytest.param(["--sf", "s.html"], "run", "stats_file", Path("s.html"), id="sf"),
            pytest.param(["-i", "in.vcf"], "source", "input_file", Path("in.vcf"), id="short_input"),
            pytest.param(["-o", "out.vcf"], "run", "output_file", Path("out.vcf"), id="short_output"),
        ],
    )
    def test_a_perl_spelling_sets_our_option(
        self, argv: list[str], section: str, field_name: str, expected: object
    ) -> None:
        """Each spelling Perl VEP documents reaches the option behind it."""
        config = self._run(*argv)

        assert getattr(config.section(section), field_name) == expected

    def test_the_canonical_spelling_still_works(self) -> None:
        """A synonym is an addition, so our own names keep working."""
        config = self._run("--log_level", "warning", "--distance", "10", "--force_overwrite")

        assert config.logging.log_level == "warning"
        assert config.source.distance == (10,)
        assert config.run.force_overwrite

    def test_a_legacy_config_file_is_read_with_perl_names(self, tmp_path: Path) -> None:
        """A config file from an existing VEP install needs no translation."""
        path = tmp_path / "vep.conf"
        path.write_text("force_overwrite\nvcf\nstats_text\ndistance 2500\n", encoding="utf-8")

        config = load(
            [],
            command="run.default",
            section_names=get_runner("run.default").sections,
            config_file=path,
        )

        assert config.run.force_overwrite
        assert config.run.output_format == "vcf"
        assert config.run.stats_format == "text"
        assert config.source.distance == (2500,)


@pytest.mark.skipif(not RUNNER_AVAILABLE, reason="runner component is not available")
class TestForCommand:
    """Building a configuration for one command from flat option names."""

    def test_each_option_reaches_the_section_that_declares_it(self) -> None:
        """A caller names options as the command line names them, not by section."""
        config = Config.for_command(
            "run.default",
            input_file="variants.vcf",
            fasta="reference.fa",
            output_file="annotated.vcf",
            no_stats=True,
        )

        assert config.command == "run.default"
        assert config.source.input_file == Path("variants.vcf")
        assert config.source.fasta == Path("reference.fa")
        assert config.run.output_file == Path("annotated.vcf")
        assert config.run.no_stats

    def test_values_are_coerced_as_the_command_line_coerces_them(self) -> None:
        """A path given as text becomes a Path, so an API caller need not import one."""
        config = Config.for_command("corpus.save", fuzz="7", corpus_dir="corpus", force="yes")

        assert config.corpus.fuzz == 7
        assert config.corpus.corpus_dir == Path("corpus")
        assert config.corpus_save.force is True

    def test_a_perl_synonym_and_shorthand_are_accepted(self) -> None:
        """The names a Perl VEP user already knows work here too."""
        config = Config.for_command("run.default", input_file="variants.vcf", flank=100, force=True, vcf=True)

        assert config.source.distance == (100,)
        assert config.source.upstream_distance == 100
        assert config.run.force_overwrite is True
        assert config.run.output_format == "vcf"

    def test_a_flag_spelling_is_accepted(self) -> None:
        """The flag a command advertises resolves to the field behind it."""
        config = Config.for_command("explain", list="true")

        assert config.explain.list_overlaps is True

    def test_sections_the_command_does_not_read_keep_their_defaults(self) -> None:
        """Only the command's own sections are built from the given options."""
        config = Config.for_command("run.default", input_file="variants.vcf")

        assert config.corpus == CorpusConfig()

    def test_an_option_the_command_does_not_read_is_rejected(self) -> None:
        """Silently dropping it would annotate with something other than what was asked for."""
        with pytest.raises(ValueError, match="'processes' is not an option of 'explain'"):
            Config.for_command("explain", processes=4)

    def test_an_unknown_option_is_rejected(self) -> None:
        """A typo fails at construction rather than at use."""
        with pytest.raises(ValueError, match="not an option of"):
            Config.for_command("run.default", fasta_file="reference.fa")

    def test_an_unknown_command_lists_the_registered_ones(self) -> None:
        """The registry's own message is the useful half of the error."""
        with pytest.raises(ValueError, match="Unknown command 'run.nope'"):
            Config.for_command("run.nope")


class TestSectionTypesFor:
    """Mapping section names to classes."""

    def test_names_resolve_in_the_order_given(self) -> None:
        """Names resolve in the order they were given."""
        resolved = Config.section_types_for(("corpus", "logging"))
        assert [section_type.section_name for section_type in resolved] == ["corpus", "logging"]

    def test_an_unknown_name_is_rejected(self) -> None:
        """Unknown names fail during resolution."""
        with pytest.raises(KeyError):
            Config.section_types_for(("nope",))


class TestForkSpelling:
    """``--fork`` is a spelling of ``processes``, so the usual layer precedence applies to it."""

    def test_the_command_line_fork_overrides_a_file_processes(self, tmp_path: Path) -> None:
        """A file's ``processes`` does not make ``--fork 4`` an error; the command line wins."""
        path = tmp_path / "vep.toml"
        path.write_text("[execution]\nprocesses = 2\n", encoding="utf-8")

        config = load(["--fork", "4"], section_names=("execution",), config_file=path, environ={})

        assert config.execution.processes == 4

    def test_the_command_line_fork_overrides_an_environment_processes(self) -> None:
        """The environment's ``processes`` is overridden by ``--fork``."""
        config = load(
            ["--fork", "4"],
            section_names=("execution",),
            environ={"ENSEMBL_VEP_EXECUTION_PROCESSES": "2"},
        )

        assert config.execution.processes == 4

    def test_the_command_line_processes_overrides_a_file_fork(self, tmp_path: Path) -> None:
        """The reverse also holds: a file's ``fork`` loses to ``--processes``."""
        path = tmp_path / "vep.toml"
        path.write_text("[execution]\nfork = 2\n", encoding="utf-8")

        config = load(["--processes", "3"], section_names=("execution",), config_file=path, environ={})

        assert config.execution.processes == 3

    def test_fork_in_a_legacy_config_file_sets_processes(self, tmp_path: Path) -> None:
        """A Perl config line ``fork 4`` sets the worker count."""
        path = tmp_path / "vep.conf"
        path.write_text("fork 4\n", encoding="utf-8")

        config = load([], section_names=("execution",), config_file=path, environ={})

        assert config.execution.processes == 4


class TestBooleanLayers:
    """A negative-sense switch set by a lower layer can be turned off from the command line."""

    def test_stats_turns_off_a_file_no_stats(self, tmp_path: Path) -> None:
        """``--stats`` overrides ``no_stats = true`` from a file."""
        path = tmp_path / "vep.toml"
        path.write_text("[run]\nno_stats = true\n", encoding="utf-8")

        kept = load([], section_names=("run",), config_file=path, environ={})
        turned_on = load(["--stats"], section_names=("run",), config_file=path, environ={})

        assert kept.run.no_stats
        assert not turned_on.run.no_stats

    def test_progress_turns_off_an_environment_no_progress(self) -> None:
        """``--progress`` overrides ``no_progress`` set in the environment."""
        config = load(["--progress"], section_names=("run",), environ={"ENSEMBL_VEP_RUN_NO_PROGRESS": "1"})

        assert not config.run.no_progress


class TestPrecedence:
    """Precedence order."""

    def test_defaults_win_when_nothing_else_sets_a_value(self) -> None:
        """Defaults apply when nothing else sets a value."""
        config = load([], section_names=("corpus",), environ={})
        assert (config.corpus.fuzz, config.source_of("corpus", "fuzz")) == (4000, ConfigSource.DEFAULT)

    def test_file_beats_defaults(self, toml_file: Path) -> None:
        """A file beats defaults."""
        config = load([], section_names=("corpus",), config_file=toml_file, environ={})
        assert (config.corpus.fuzz, config.source_of("corpus", "fuzz")) == (111, ConfigSource.FILE)

    def test_environment_beats_file(self, toml_file: Path) -> None:
        """Environment variables beat files."""
        config = load(
            [], section_names=("corpus",), config_file=toml_file, environ={"ENSEMBL_VEP_CORPUS_FUZZ": "555"}
        )
        assert (config.corpus.fuzz, config.source_of("corpus", "fuzz")) == (555, ConfigSource.ENVIRONMENT)

    def test_command_line_beats_environment(self, toml_file: Path) -> None:
        """Command-line values win."""
        config = load(
            ["--fuzz", "999"],
            section_names=("corpus",),
            config_file=toml_file,
            environ={"ENSEMBL_VEP_CORPUS_FUZZ": "555"},
        )
        assert (config.corpus.fuzz, config.source_of("corpus", "fuzz")) == (999, ConfigSource.COMMAND_LINE)

    def test_an_unset_flag_does_not_clobber_a_file_value(self, toml_file: Path) -> None:
        """Unset flags do not clobber file values."""
        config = load(["--step", "888"], section_names=("corpus",), config_file=toml_file, environ={})
        assert (config.corpus.step, config.corpus.fuzz) == (888, 111)

    def test_no_argv_skips_the_command_line_layer(self) -> None:
        """Passing no argv skips command-line parsing."""
        assert load(None, section_names=("corpus",), environ={}).corpus.fuzz == 4000

    def test_omitting_section_names_loads_every_section(self) -> None:
        """Omitting section names loads everything."""
        assert len(load(None, environ={}).sections()) == len(Config.section_types)
