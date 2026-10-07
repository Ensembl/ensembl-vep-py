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
"""Tests for the derived argparse surface."""

from __future__ import annotations

import argparse
from enum import Enum
from pathlib import Path
from typing import ClassVar

import pytest
from attrs import define, field

from ensembl.vep.config import Config, load
from ensembl.vep.config.option_spec import option
from ensembl.vep.config.sections.base import ConfigSection
from ensembl.vep.config.parser.common import resolve_option, resolve_option_key
from ensembl.vep.config.parser.cli import (
    add_section_arguments,
    build_parser,
    check_no_flag_collisions,
    parse_args,
)
from ensembl.vep.config.sections import CorpusConfig


class Colour(Enum):
    """Enum option used only by these tests."""

    RED = "red"
    BLUE = "blue"


@define(frozen=True, kw_only=True, auto_attribs=True)
class SampleConfig(ConfigSection):
    """Section used only by these tests."""

    section_name: ClassVar[str] = "sample"

    count: int = field(default=7, metadata=option("how many", short="n"))
    label: str | None = field(default=None, metadata=option("what to call it"))
    enabled: bool = field(default=False, metadata=option("switch it on"))
    no_frills: bool = field(default=False, metadata=option("a negative-sense switch"))
    sizes: tuple[int, ...] = field(default=(1, 2), metadata=option("the sizes"))
    specs: tuple[str, ...] = field(default=(), metadata=option("whole specs", repeatable=True))
    path: Path | None = field(default=None, metadata=option("where"))
    ratio: float = field(default=0.5, metadata=option("a fraction"))
    mode: str = field(default="a", metadata=option("pick one", choices=("a", "b")))
    renamed: int = field(default=0, metadata=option("has its own flag", flag="other_name"))
    concealed: int = field(default=0, metadata=option("kept out of help", hidden=True))
    colour: Colour = field(default=Colour.RED, metadata=option("an enum-typed option"))
    aliased: str | None = field(default=None, metadata=option("has synonyms", synonyms=("alias", "second")))
    shaped: str = field(
        default="plain",
        metadata=option("set by name or by shorthand", constants={"fancy": "ornate", "bare": "plain"}),
    )


@define(frozen=True, kw_only=True, auto_attribs=True)
class OtherConfig(ConfigSection):
    """A second section declaring the same field name, as source and stats both declare format."""

    section_name: ClassVar[str] = "other"

    mode: str = field(default="z", metadata=option("the other mode", flag="other_mode"))


@pytest.fixture(name="parser")
def fixture_parser() -> argparse.ArgumentParser:
    """Parser derived from the sample section."""
    return build_parser("vep sample", "", (SampleConfig,))


class TestSparseParsing:
    """Sparse parsing."""

    def test_nothing_given_yields_nothing(self, parser: argparse.ArgumentParser) -> None:
        """An empty command line sets nothing."""
        assert not parse_args(parser, [], (SampleConfig,))

    def test_only_what_was_given_is_reported(self, parser: argparse.ArgumentParser) -> None:
        """Defaults belong to the section, not the parser."""
        assert parse_args(parser, ["--count", "9"], (SampleConfig,)) == {"sample": {"count": 9}}

    def test_several_options_at_once(self, parser: argparse.ArgumentParser) -> None:
        """Each value lands under its section."""
        parsed = parse_args(parser, ["--count", "1", "--label", "x"], (SampleConfig,))
        assert parsed == {"sample": {"count": 1, "label": "x"}}

    def test_a_section_with_nothing_set_is_absent(self) -> None:
        """An empty section is omitted."""
        types = Config.section_types_for(("corpus", "logging"))
        parser = build_parser("vep x", "", types)
        assert set(parse_args(parser, ["--fuzz", "2"], types)) == {"corpus"}

    def test_omitted_annotation_options_use_section_defaults(self) -> None:
        """Keep the benchmark-selected annotation defaults off the command line."""

        config = load(argv=[], section_names=("annotation",))

        assert config.annotation.annotator == "consequence.columnar-staged"
        assert not config.annotation.zone_map

    def test_zone_map_and_its_byte_budget_are_parsed_together(self) -> None:
        """--zone_map and --zone_map_bytes both reach the annotation section."""
        config = load(argv=["--zone_map", "--zone_map_bytes", "4096"], section_names=("annotation",))

        assert config.annotation.zone_map
        assert config.annotation.zone_map_bytes == 4096


class TestTypeDerivation:
    """Type-driven parsing."""

    def test_booleans_get_both_senses(self, parser: argparse.ArgumentParser) -> None:
        """Boolean flags support both senses."""
        assert parse_args(parser, ["--enabled"], (SampleConfig,)) == {"sample": {"enabled": True}}
        assert parse_args(parser, ["--no-enabled"], (SampleConfig,)) == {"sample": {"enabled": False}}

    def test_a_positive_boolean_takes_the_underscore_negative(self, parser: argparse.ArgumentParser) -> None:
        """``--no_enabled`` is the canonical negative, in Perl VEP's underscore style."""
        assert parse_args(parser, ["--no_enabled"], (SampleConfig,)) == {"sample": {"enabled": False}}

    def test_a_negative_sense_boolean_is_a_switch_with_a_counterpart(
        self, parser: argparse.ArgumentParser
    ) -> None:
        """A negative-sense flag keeps its Perl spelling and gains a positive that turns it off."""
        assert parse_args(parser, ["--no_frills"], (SampleConfig,)) == {"sample": {"no_frills": True}}
        assert parse_args(parser, ["--frills"], (SampleConfig,)) == {"sample": {"no_frills": False}}

    def test_the_last_boolean_spelling_wins(self, parser: argparse.ArgumentParser) -> None:
        """Giving both senses leaves the later one in force."""
        parsed = parse_args(parser, ["--no_frills", "--frills", "--enabled", "--no_enabled"], (SampleConfig,))

        assert parsed == {"sample": {"no_frills": False, "enabled": False}}

    def test_zone_map_accepts_every_negative_spelling(self) -> None:
        """Both ``--no_zone_map`` and the hyphenated ``--no-zone_map`` turn the zone map off."""
        types = Config.section_types_for(("annotation",))
        annotation_parser = build_parser("vep x", "", types)

        for spelling in ("--no_zone_map", "--no-zone_map"):
            parsed = parse_args(annotation_parser, ["--zone_map", spelling], types)
            assert parsed == {"annotation": {"zone_map": False}}

    def test_help_shows_the_underscore_negative_and_hides_the_counterparts(
        self, parser: argparse.ArgumentParser
    ) -> None:
        """Help lists ``--no_enabled`` but neither ``--no-enabled`` nor the ``--frills`` counterpart."""
        help_text = parser.format_help()

        assert "--enabled, --no_enabled" in help_text
        assert "--no-enabled" not in help_text
        assert "\n  --frills" not in help_text
        assert "(--frills turns it back on)" in help_text

    def test_a_repeatable_option_keeps_each_occurrence_whole(self, parser: argparse.ArgumentParser) -> None:
        """Each ``--specs`` is one item; commas inside it are not split."""
        parsed = parse_args(
            parser, ["--specs", "a.vcf,ClinVar,vcf", "--specs", "b.bed,Other"], (SampleConfig,)
        )

        assert parsed == {"sample": {"specs": ("a.vcf,ClinVar,vcf", "b.bed,Other")}}

    def test_an_unrepeatable_tuple_still_splits_on_commas(self, parser: argparse.ArgumentParser) -> None:
        """Only an option declared repeatable stops comma splitting."""
        assert parse_args(parser, ["--sizes", "1,2"], (SampleConfig,)) == {"sample": {"sizes": (1, 2)}}

    def test_tuples_are_comma_separated(self, parser: argparse.ArgumentParser) -> None:
        """Tuples are comma separated."""
        assert parse_args(parser, ["--sizes", "4,8"], (SampleConfig,)) == {"sample": {"sizes": (4, 8)}}

    def test_an_empty_tuple_is_accepted(self, parser: argparse.ArgumentParser) -> None:
        """An empty tuple is accepted."""
        assert parse_args(parser, ["--sizes", ""], (SampleConfig,)) == {"sample": {"sizes": ()}}

    def test_paths_are_converted(self, parser: argparse.ArgumentParser) -> None:
        """Path values are converted."""
        assert parse_args(parser, ["--path", "/tmp/x"], (SampleConfig,)) == {
            "sample": {"path": Path("/tmp/x")}
        }

    def test_floats_are_converted(self, parser: argparse.ArgumentParser) -> None:
        """Float values are converted."""
        assert parse_args(parser, ["--ratio", "0.25"], (SampleConfig,)) == {"sample": {"ratio": 0.25}}

    def test_a_short_flag_is_offered_when_declared(self, parser: argparse.ArgumentParser) -> None:
        """Short flags are offered when declared."""
        assert parse_args(parser, ["-n", "3"], (SampleConfig,)) == {"sample": {"count": 3}}

    def test_choices_are_enforced(self, parser: argparse.ArgumentParser) -> None:
        """Invalid choices fail."""
        with pytest.raises(SystemExit):
            parse_args(parser, ["--mode", "z"], (SampleConfig,))

    def test_a_bad_int_is_a_usage_error(self, parser: argparse.ArgumentParser) -> None:
        """Invalid ints fail at parse time."""
        with pytest.raises(SystemExit):
            parse_args(parser, ["--count", "abc"], (SampleConfig,))

    def test_an_enum_is_named_by_its_lowercased_member(self, parser: argparse.ArgumentParser) -> None:
        """Enum members are named by lowercased member name."""
        assert parse_args(parser, ["--colour", "blue"], (SampleConfig,)) == {
            "sample": {"colour": Colour.BLUE}
        }

    def test_an_unknown_enum_member_is_a_usage_error(self, parser: argparse.ArgumentParser) -> None:
        """Unknown enum members fail."""
        with pytest.raises(SystemExit):
            parse_args(parser, ["--colour", "green"], (SampleConfig,))


class TestSourceShapes:
    """The Perl-shaped ``--cache`` and ``--custom`` options."""

    def test_cache_is_a_switch_and_dir_cache_takes_the_directory(self) -> None:
        """``--cache`` takes no value; ``--dir_cache`` carries the path."""
        config = load(["--cache", "--dir_cache", "/data/vep"], section_names=("source",), environ={})

        assert config.source.cache is True
        assert config.source.dir_cache == Path("/data/vep")

    def test_custom_given_twice_keeps_each_spec_whole(self) -> None:
        """Each ``--custom`` is one spec, commas included."""
        config = load(
            ["--custom", "a.vcf.gz,ClinVar,vcf,exact", "--custom", "b.bed,Other,bed"],
            section_names=("source",),
            environ={},
        )

        assert config.source.custom == ("a.vcf.gz,ClinVar,vcf,exact", "b.bed,Other,bed")


class TestHelp:
    """Help output."""

    def test_a_hidden_option_is_absent_from_help(self, parser: argparse.ArgumentParser) -> None:
        """Hidden options stay out of help."""
        assert "--concealed" not in parser.format_help()

    def test_a_hidden_option_still_parses(self, parser: argparse.ArgumentParser) -> None:
        """Hidden options still parse."""
        assert parse_args(parser, ["--concealed", "4"], (SampleConfig,)) == {"sample": {"concealed": 4}}

    def test_help_text_comes_from_the_declaration(self, parser: argparse.ArgumentParser) -> None:
        """Help text comes from the declaration."""
        assert "how many" in parser.format_help()

    def test_enum_choices_use_the_accepted_command_line_names(self, parser: argparse.ArgumentParser) -> None:
        """Show enum member names rather than Python enum representations."""

        help_text = parser.format_help()

        assert "--colour {red,blue}" in help_text
        assert "Colour.RED" not in help_text

    def test_options_are_grouped_by_section(self) -> None:
        """Sections are grouped in help output."""
        types = Config.section_types_for(("corpus", "logging"))
        assert "corpus options:" in build_parser("vep x", "", types).format_help()


class TestSynonyms:
    """Further long spellings of one option, for the names Perl VEP uses."""

    def test_each_synonym_sets_the_option(self, parser: argparse.ArgumentParser) -> None:
        """Every spelling reaches the same field."""
        for spelling in ("--aliased", "--alias", "--second"):
            layer = parse_args(parser, [spelling, "given"], (SampleConfig,))
            assert layer["sample"]["aliased"] == "given"

    def test_a_synonym_is_shown_with_the_option(self, parser: argparse.ArgumentParser) -> None:
        """A caller reading --help sees the spelling they may already know."""
        assert "--alias" in parser.format_help()

    def test_a_synonym_resolves_from_a_flat_key(self) -> None:
        """A config file may use the synonym as its key."""
        assert resolve_option_key("alias", (SampleConfig,)) == ("sample", "aliased")


class TestConstants:
    """Flags that set an option to a fixed value, as Perl VEP's --vcf does."""

    def test_a_constant_flag_sets_its_value(self, parser: argparse.ArgumentParser) -> None:
        """The flag carries the value, so none is given with it."""
        layer = parse_args(parser, ["--fancy"], (SampleConfig,))
        assert layer["sample"]["shaped"] == "ornate"

    def test_the_option_itself_still_takes_a_value(self, parser: argparse.ArgumentParser) -> None:
        """The shorthand is an addition, not a replacement."""
        layer = parse_args(parser, ["--shaped", "ornate"], (SampleConfig,))
        assert layer["sample"]["shaped"] == "ornate"

    def test_a_constant_flag_is_described_in_help(self, parser: argparse.ArgumentParser) -> None:
        """The help says which option it stands for."""
        assert "shorthand for --shaped ornate" in parser.format_help()

    def test_a_constant_resolves_from_a_flat_key_with_its_value(self) -> None:
        """A legacy config file writes "vcf 1"; the flag's own value is what it means."""
        assert resolve_option("fancy", "1", (SampleConfig,)) == ("sample", "shaped", "ornate")


class TestSectionNamespacing:
    """Two sections of one command may declare the same field name."""

    def test_a_value_reaches_only_the_section_it_was_declared_in(self) -> None:
        """``--other_mode`` is the other section's; ``--mode`` is the sample's."""
        parser = build_parser("vep two", "", (SampleConfig, OtherConfig))

        layer = parse_args(parser, ["--other_mode", "z"], (SampleConfig, OtherConfig))

        assert layer["other"]["mode"] == "z"
        assert "mode" not in layer.get("sample", {})

    def test_the_value_name_in_help_follows_the_flag(self) -> None:
        """The namespaced destination must not leak into the usage line."""
        parser = build_parser("vep two", "", (SampleConfig, OtherConfig))

        assert "--other_mode OTHER_MODE" in parser.format_help()
        assert "SAMPLE__" not in parser.format_help()


class TestFlagCollisions:
    """Flag collisions."""

    def test_a_collision_names_both_sections(self) -> None:
        """The error names both sections."""
        with pytest.raises(ValueError, match="output_file"):
            check_no_flag_collisions(Config.section_types_for(("run", "stats")))

    def test_disjoint_sections_are_accepted(self) -> None:
        """Disjoint sections are accepted."""
        check_no_flag_collisions(Config.section_types_for(("corpus", "logging")))

    def test_building_a_parser_checks_first(self) -> None:
        """Collisions are checked before argparse runs."""
        with pytest.raises(ValueError, match="Sections declare colliding flags"):
            build_parser("vep x", "", Config.section_types_for(("run", "stats")))

    def test_every_registered_command_has_a_usable_section_set(self) -> None:
        """Every runner uses a valid section set."""
        runner = pytest.importorskip("ensembl.vep.runner")

        for name in runner.available_runners():
            check_no_flag_collisions(Config.section_types_for(runner.get_runner(name).sections))


class TestAddSectionArguments:
    """Adding one section to a parser."""

    def test_arguments_land_in_a_named_group(self) -> None:
        """Arguments appear in a named group."""
        parser = argparse.ArgumentParser()
        add_section_arguments(parser, CorpusConfig)
        assert "corpus options:" in parser.format_help()
