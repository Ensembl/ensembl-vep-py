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
"""Tests for reading configuration from the environment."""

from __future__ import annotations

import pytest

from ensembl.vep.config import Config
from ensembl.vep.config.parser.env import default_env_name, read_environment


class TestEnvironment:
    """Environment loading."""

    def test_names_are_derived_from_section_and_field(self) -> None:
        """Variable names are derived."""
        assert default_env_name("corpus", "fuzz") == "ENSEMBL_VEP_CORPUS_FUZZ"

    def test_only_variables_that_are_set_are_reported(self) -> None:
        """Unset variables are ignored."""
        types = Config.section_types_for(("corpus",))
        assert not read_environment(types, {})

    def test_a_declared_env_name_overrides_the_derived_one(self) -> None:
        """Legacy names still resolve."""
        types = Config.section_types_for(("annotation",))
        layer = read_environment(types, {"ENSEMBL_VEP_FACT_COST_OVERRIDE": "1"})
        assert layer == {"annotation": {"fact_cost_override": True}}

    def test_values_are_coerced_from_text(self) -> None:
        """Environment values are coerced."""
        types = Config.section_types_for(("corpus",))
        assert read_environment(types, {"ENSEMBL_VEP_CORPUS_FUZZ": "9"}) == {"corpus": {"fuzz": 9}}

    def test_a_value_outside_its_choices_is_an_error_naming_the_option(self) -> None:
        """Choices are enforced for environment values."""
        types = Config.section_types_for(("logging",))
        with pytest.raises(ValueError, match="logging.log_level: 'verbose' is not one of debug, info"):
            read_environment(types, {"ENSEMBL_VEP_LOGGING_LOG_LEVEL": "verbose"})

    def test_a_bad_enum_name_is_an_error_naming_the_option(self) -> None:
        """An unknown enum name is a ValueError listing the accepted names, not a KeyError."""
        types = Config.section_types_for(("execution",))
        with pytest.raises(ValueError, match="execution.start_method: 'bogus' is not one of auto, "):
            read_environment(types, {"ENSEMBL_VEP_EXECUTION_START_METHOD": "bogus"})

    def test_a_bad_number_is_an_error_naming_the_option(self) -> None:
        """A conversion failure names the option."""
        types = Config.section_types_for(("corpus",))
        with pytest.raises(ValueError, match="corpus.fuzz:"):
            read_environment(types, {"ENSEMBL_VEP_CORPUS_FUZZ": "many"})
