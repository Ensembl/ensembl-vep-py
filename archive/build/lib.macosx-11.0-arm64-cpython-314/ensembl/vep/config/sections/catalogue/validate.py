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
"""Validation options for the catalogue rule set."""

from __future__ import annotations

from typing import ClassVar

from attrs import define, field

from ensembl.vep.config.option_spec import option
from ensembl.vep.config.sections.base import ConfigSection

__all__ = ["CatalogueValidateConfig"]

SOLVERS = ("none", "z3")


@define(frozen=True, kw_only=True, auto_attribs=True)
class CatalogueValidateConfig(ConfigSection):
    """Validation tiers to run.

    Attributes:
        solver (str): Whole-ruleset SAT backend. ``z3`` needs the ``validation`` extra.
        verify_requires (bool): Check every declared ``requires`` claim against a real corpus.
    """

    section_name: ClassVar[str] = "catalogue_validate"

    solver: str = field(default="none", metadata=option("whole-ruleset SAT backend", choices=SOLVERS))
    verify_requires: bool = field(
        default=False, metadata=option("check every declared requires against a corpus")
    )
