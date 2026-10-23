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
"""Fixtures shared by the consequence rule tests."""

from __future__ import annotations

import pytest

from ensembl.vep.annotators.consequence.rules.rulesets import FactRuleSet, build_fact_ruleset


@pytest.fixture(scope="session")
def transcript_ruleset() -> FactRuleSet:
    """The compiled transcript rule set the real catalogue produces."""
    return build_fact_ruleset()
