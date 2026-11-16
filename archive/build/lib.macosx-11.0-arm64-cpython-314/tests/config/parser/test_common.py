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
"""Tests for parser helpers shared across readers."""

from __future__ import annotations

from typing import ClassVar

import pytest
from attrs import define, field

from ensembl.vep.config import Config
from ensembl.vep.config.option_spec import option, section_options
from ensembl.vep.config.sections.base import ConfigSection
from ensembl.vep.config.parser.common import flag_for, resolve_option_key


@define(frozen=True, kw_only=True, auto_attribs=True)
class SampleConfig(ConfigSection):
    """Section used only by these tests."""

    section_name: ClassVar[str] = "sample"

    count: int = field(default=7, metadata=option("how many"))
    renamed: int = field(default=0, metadata=option("has its own flag", flag="other_name"))


class TestFlagFor:
    """Flag names come from the declaration."""

    def test_the_field_name_is_the_default_spelling(self) -> None:
        """Underscores are kept."""
        _field, spec = section_options(SampleConfig)["count"]
        assert flag_for("count", spec) == "count"

    def test_an_explicit_flag_wins(self) -> None:
        """An explicit flag name is respected."""
        _field, spec = section_options(SampleConfig)["renamed"]
        assert flag_for("renamed", spec) == "other_name"


class TestResolveOptionKey:
    """Resolve flat keys from legacy config files."""

    @pytest.mark.parametrize("key", ["fuzz", "--fuzz", "  fuzz  "])
    def test_a_declared_key_resolves_to_its_section(self, key: str) -> None:
        """Leading dashes and surrounding space are ignored."""
        types = Config.section_types_for(("corpus",))
        assert resolve_option_key(key, types) == ("corpus", "fuzz")

    def test_a_dashed_spelling_resolves(self) -> None:
        """Dashes and underscores both resolve."""
        types = Config.section_types_for(("run",))
        assert resolve_option_key("output-file", types) == ("run", "output_file")

    def test_an_explicit_flag_name_resolves_to_its_field(self) -> None:
        """An explicit flag resolves to its field."""
        types = Config.section_types_for(("explain",))
        assert resolve_option_key("list", types) == ("explain", "list_overlaps")

    def test_an_unknown_key_is_unresolved(self) -> None:
        """Unknown keys stay unresolved."""
        types = Config.section_types_for(("corpus",))
        assert resolve_option_key("sift", types) is None
