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
"""Tests for option declarations."""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

import pytest
from attrs import define, field

from ensembl.vep.config.sections.base import ConfigSection
from ensembl.vep.config.option_spec import (
    OptionSpec,
    option,
    section_options,
    unwrap_optional_type,
)


@define(frozen=True, kw_only=True, auto_attribs=True)
class SampleConfig(ConfigSection):
    """Section used only by these tests."""

    section_name: ClassVar[str] = "sample"

    count: int = field(default=7, metadata=option("how many"))
    label: str | None = field(default=None, metadata=option("what to call it"))
    enabled: bool = field(default=False, metadata=option("switch it on"))
    sizes: tuple[int, ...] = field(default=(1, 2), metadata=option("the sizes"))
    path: Path | None = field(default=None, metadata=option("where"))
    undeclared: int = 3


class TestOptionSpec:
    """Individual option specs."""

    def test_choices_may_be_a_literal_sequence(self) -> None:
        """Literal choices are returned unchanged."""
        assert OptionSpec(help="h", choices=("a", "b")).resolved_choices() == ("a", "b")

    def test_choices_may_be_a_callable(self) -> None:
        """Callable choices are resolved when needed."""
        assert OptionSpec(help="h", choices=lambda: ["x"]).resolved_choices() == ("x",)

    def test_no_choices_means_unrestricted(self) -> None:
        """Missing choices stay unrestricted."""
        assert OptionSpec(help="h").resolved_choices() is None

    def test_option_returns_metadata_under_one_key(self) -> None:
        """The spec is stored on the field metadata."""
        metadata = option("how many", short="n", metavar="N")
        spec = next(iter(metadata.values()))
        assert (spec.help, spec.short, spec.metavar) == ("how many", "n", "N")


class TestOptionSpecs:
    """Reading declarations back from a section."""

    def test_only_fields_carrying_option_metadata_are_options(self) -> None:
        """Only fields with option metadata count."""
        assert set(section_options(SampleConfig)) == {"count", "label", "enabled", "sizes", "path"}

    def test_declaration_order_is_preserved(self) -> None:
        """Declaration order is preserved."""
        assert list(section_options(SampleConfig)) == ["count", "label", "enabled", "sizes", "path"]

    def test_field_types_are_resolved_rather_than_left_as_strings(self) -> None:
        """Field types are resolved, not left as strings."""
        specs = section_options(SampleConfig)
        assert specs["count"][0].type is int
        assert specs["sizes"][0].type == tuple[int, ...]

    def test_a_section_with_no_options_is_empty_rather_than_an_error(self) -> None:
        """The base section declares no options."""
        assert not section_options(ConfigSection)


class TestUnwrapOptional:
    """Optional type handling."""

    @pytest.mark.parametrize(
        "annotation,expected",
        [(int, (int, False)), (int | None, (int, True)), (str | None, (str, True))],
    )
    def test_unwrapping(self, annotation: object, expected: tuple[object, bool]) -> None:
        """Bare and optional annotations are separated."""
        assert unwrap_optional_type(annotation) == expected

    def test_a_multi_member_union_is_left_alone(self) -> None:
        """Multi-member unions are left alone."""
        annotation = int | str | None
        assert unwrap_optional_type(annotation) == (annotation, False)
