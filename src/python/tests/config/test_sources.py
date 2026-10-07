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
"""Tests for configuration sources and coercion."""

from __future__ import annotations

import logging
from enum import Enum
from pathlib import Path
from typing import ClassVar

import pytest
from attrs import define, field

from ensembl.vep.config import Config
from ensembl.vep.config.option_spec import option
from ensembl.vep.config.sections.base import ConfigSection
from ensembl.vep.config.sources import (
    ConfigLayer,
    ConfigSource,
    coerce_layer,
    coerce_value,
    merge_layers,
)


class Colour(Enum):
    """Enum used only by these tests."""

    RED = "red"


@define(frozen=True, kw_only=True, auto_attribs=True)
class ToySection(ConfigSection):
    """Section with a tuple option restricted to a few choices."""

    section_name: ClassVar[str] = "toy"

    picks: tuple[str, ...] = field(default=(), metadata=option("picks", choices=("a", "b")))


class TestConfigSource:
    """The four layers, lowest precedence first."""

    def test_every_layer_has_a_readable_label(self) -> None:
        """The sources appear in a stable order."""
        assert [source.value for source in ConfigSource] == [
            "default",
            "file",
            "environment",
            "command line",
        ]


class TestCoercion:
    """Value coercion."""

    @pytest.mark.parametrize(
        "raw,annotation,expected",
        [
            ("12", int, 12),
            (12, int, 12),
            ("0.5", float, 0.5),
            ("yes", bool, True),
            ("on", bool, True),
            ("1", bool, True),
            ("off", bool, False),
            ("no", bool, False),
            (True, bool, True),
            ("a,b", tuple[str, ...], ("a", "b")),
            ("1, 2", tuple[int, ...], (1, 2)),
            ([1, 2], tuple[int, ...], (1, 2)),
            ("", tuple[int, ...], ()),
            ("/tmp/x", Path, Path("/tmp/x")),
            ("", str | None, None),
            (None, str | None, None),
            ("x", str | None, "x"),
        ],
    )
    def test_coercion(self, raw: object, annotation: object, expected: object) -> None:
        """Each declared type has a matching reader."""
        assert coerce_value(raw, annotation) == expected

    def test_a_bad_boolean_is_rejected(self) -> None:
        """Invalid booleans fail."""
        with pytest.raises(ValueError):
            coerce_value("maybe", bool)

    @pytest.mark.parametrize("raw", ["red", "RED", Colour.RED])
    def test_an_enum_reads_from_its_member_name(self, raw: object) -> None:
        """Enum names are matched case-insensitively."""
        assert coerce_value(raw, Colour) is Colour.RED

    def test_an_unknown_enum_name_is_a_value_error_listing_the_members(self) -> None:
        """A bad enum name is a ValueError, not a KeyError."""
        with pytest.raises(ValueError, match="'green' is not one of red"):
            coerce_value("green", Colour)

    def test_an_unrecognised_type_passes_through(self) -> None:
        """Unknown types pass through."""
        marker = object()
        assert coerce_value(marker, object) is marker


class TestCoerceLayer:
    """Coercing a whole layer."""

    def test_unknown_sections_and_fields_are_dropped(self) -> None:
        """Unknown sections and fields are dropped."""
        types = Config.section_types_for(("corpus",))
        layer: ConfigLayer = {"corpus": {"fuzz": "5", "nope": 1}, "absent": {"x": 1}}
        assert coerce_layer(layer, types) == {"corpus": {"fuzz": 5}}

    def test_values_are_coerced_to_the_declared_type(self) -> None:
        """Values are coerced to their declared type."""
        types = Config.section_types_for(("corpus",))
        assert coerce_layer({"corpus": {"fuzz": "7"}}, types) == {"corpus": {"fuzz": 7}}


class TestLayerChoices:
    """Declared choices are enforced on every layer."""

    def test_a_value_outside_the_choices_names_section_and_field(self) -> None:
        """The message names the option and lists the choices."""
        types = Config.section_types_for(("logging",))
        with pytest.raises(ValueError, match="logging.log_level: 'loud' is not one of debug, info,"):
            coerce_layer({"logging": {"log_level": "loud"}}, types)

    def test_a_value_inside_the_choices_passes(self) -> None:
        """Permitted values are unchanged."""
        types = Config.section_types_for(("logging",))
        assert coerce_layer({"logging": {"log_level": "warning"}}, types) == {
            "logging": {"log_level": "warning"}
        }

    def test_a_none_for_an_optional_choice_field_passes(self) -> None:
        """An unset optional field is not checked."""
        types = Config.section_types_for(("source",))
        assert coerce_layer({"source": {"format": ""}}, types) == {"source": {"format": None}}

    def test_every_member_of_a_tuple_is_checked(self) -> None:
        """A tuple value is checked item by item."""
        with pytest.raises(ValueError, match="toy.picks: 'c' is not one of a, b"):
            coerce_layer({"toy": {"picks": "a,c"}}, (ToySection,))

    def test_a_conversion_failure_names_section_and_field(self) -> None:
        """A bad number names the option."""
        types = Config.section_types_for(("corpus",))
        with pytest.raises(ValueError, match="corpus.fuzz: invalid literal"):
            coerce_layer({"corpus": {"fuzz": "x"}}, types)

    def test_an_unknown_field_in_a_known_section_is_warned_about(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        """The dropped key is named in a warning."""
        types = Config.section_types_for(("corpus",))
        with caplog.at_level(logging.WARNING, logger="ensembl.vep"):
            coerce_layer({"corpus": {"nope": 1}}, types)
        assert "unknown option 'nope' in section 'corpus'" in caplog.text


class TestMerge:
    """Layer merging."""

    def test_later_layers_win(self) -> None:
        """Later layers win."""
        assert merge_layers({"a": {"x": 1, "y": 2}}, {"a": {"y": 3}}) == {"a": {"x": 1, "y": 3}}

    def test_sections_are_combined_rather_than_replaced(self) -> None:
        """Sections are merged, not replaced."""
        assert merge_layers({"a": {"x": 1}}, {"b": {"y": 2}}) == {"a": {"x": 1}, "b": {"y": 2}}

    def test_merging_nothing_yields_nothing(self) -> None:
        """Merging nothing yields nothing."""
        assert not merge_layers()
