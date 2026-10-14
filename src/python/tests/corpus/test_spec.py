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
"""Tests for corpus construction specifications."""

from __future__ import annotations

import pytest

from ensembl.vep.corpus.spec import CorpusSpec


class TestCorpusSpec:
    """Registered corpus names and factory-owned parameters."""

    def test_accepts_generic_factory_parameters(self) -> None:
        """A corpus factory owns the meaning and validation of its parameters."""
        specification = CorpusSpec(kind="custom", parameters={"window": 250, "strict": True})
        assert specification.parameters == {"strict": True, "window": 250}

    def test_default_uses_the_canonical_combined_name(self) -> None:
        """The default specification expands to both canonical synthetic corpora."""
        specification = CorpusSpec()
        assert specification.kind == "synthetic.all"
        assert specification.kinds == ("synthetic.default", "synthetic.heavy")

    def test_parameters_are_returned_as_a_copy(self) -> None:
        """Callers cannot mutate a specification through the public mapping."""
        specification = CorpusSpec(kind="custom", parameters={"seed": 7})
        parameters = specification.parameters
        parameters["seed"] = 8
        assert specification.parameters == {"seed": 7}

    def test_rejects_non_scalar_parameters(self) -> None:
        """Manifest-safe scalar values form the registry boundary."""
        with pytest.raises(TypeError, match="must be a scalar"):
            CorpusSpec(kind="custom", parameters={"items": [1, 2]})  # type: ignore[dict-item]

    def test_json_representation_round_trips(self) -> None:
        """Specifications retain identity through manifests."""
        specification = CorpusSpec(kind="custom", parameters={"seed": 7})
        assert CorpusSpec.from_dict(specification.as_dict()) == specification

    def test_parameter_hash_is_stable(self) -> None:
        """Equivalent mappings have the same construction identity."""
        first = CorpusSpec(kind="custom", parameters={"seed": 7, "window": 250})
        second = CorpusSpec(kind="custom", parameters={"window": 250, "seed": 7})
        assert first.parameter_hash() == second.parameter_hash()
