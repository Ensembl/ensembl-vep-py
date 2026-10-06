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

from pathlib import Path

import pytest

from ensembl.vep.config import Config, ConfigSource, CorpusConfig, load
from ensembl.vep.config.option_spec import section_options


@pytest.fixture(name="toml_file")
def fixture_toml_file(tmp_path: Path) -> Path:
    """A TOML file setting two corpus options."""
    path = tmp_path / "vep.toml"
    path.write_text("[corpus]\nfuzz = 111\nstep = 222\n", encoding="utf-8")
    return path


class TestConfig:
    """The configuration root."""

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

    def test_evolve_replaces_one_section_and_leaves_the_rest(self) -> None:
        """Evolve replaces one section and leaves the rest alone."""
        config = Config().evolve(corpus=CorpusConfig(fuzz=1))
        assert config.corpus.fuzz == 1
        assert config.bench.repeats == Config().bench.repeats

    def test_an_unset_origin_reads_as_default(self) -> None:
        """Unset values read as the default source."""
        assert Config().source_of("corpus", "fuzz") is ConfigSource.DEFAULT

    def test_the_command_is_recorded(self) -> None:
        """The command name is preserved."""
        assert load([], command="run", section_names=("corpus",), environ={}).command == "run"


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

    def test_overrides_are_the_api_entry_point(self) -> None:
        """Overrides are the API entry point."""
        config = load(None, section_names=("corpus",), overrides={"corpus": {"fuzz": 12}}, environ={})
        assert config.corpus.fuzz == 12

    def test_no_argv_skips_the_command_line_layer(self) -> None:
        """Passing no argv skips command-line parsing."""
        assert load(None, section_names=("corpus",), environ={}).corpus.fuzz == 4000

    def test_omitting_section_names_loads_every_section(self) -> None:
        """Omitting section names loads everything."""
        assert len(load(None, environ={}).sections()) == len(Config.section_types)
