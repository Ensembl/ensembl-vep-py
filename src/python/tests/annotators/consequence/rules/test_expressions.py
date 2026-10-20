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

# pylint: disable=redefined-outer-name
"""Tests for consequence rule expressions."""

from __future__ import annotations

import pytest

from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.annotators.consequence.rules.expressions import Atom, Not, all_of, any_of, expr_from_json, no
from ensembl.vep.common import Strand
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.models.variation import VariationAllele


def _create_snv_fact_context(transcript: Transcript, *, position: int) -> FactContext:
    """Create a fact context for a single-nucleotide transcript variation allele."""
    variation_allele = VariationAllele(
        start=position,
        end=position,
        strand=Strand.FORWARD,
        slice=transcript.slice,
        reference_allele="A",
        alternate_allele="T",
    )
    transcript_variation_allele = TranscriptVariationAllele(
        feature=transcript, overlapping_feature=variation_allele
    )
    return FactContext(overlap=transcript_variation_allele)


class TestAtom:
    """Test Atom expressions."""

    def test_repr_is_the_bare_name(self) -> None:
        """Test that repr returns the atom name."""
        assert repr(Atom("is_snv")) == "is_snv"

    def test_atoms_returns_the_atom_name(self) -> None:
        """Test that atoms returns the atom name."""
        assert Atom("is_snv").atoms() == {"is_snv"}

    def test_evaluate_reads_the_fact(self, synthetic_transcripts_by_stable_id: dict[str, Transcript]) -> None:
        """Test that evaluate reads the named fact from the context."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        context = _create_snv_fact_context(transcript, position=transcript.start)
        assert Atom("is_snv").evaluate(context) is True

    def test_evaluate_partial_returns_known_values(self) -> None:
        """Test that partial evaluation returns known values and preserves unknown values."""
        assert Atom("is_snv").evaluate_partial({"is_snv": True}) is True
        assert Atom("is_snv").evaluate_partial({"is_snv": False}) is False
        assert Atom("is_snv").evaluate_partial({}) is None


class TestNot:
    """Test Not expressions."""

    def test_repr_prefixes_the_child_expression(self) -> None:
        """Test that repr prefixes the child expression with an exclamation mark."""
        assert repr(Not(Atom("is_snv"))) == "!is_snv"

    def test_atoms_returns_the_child_atoms(self) -> None:
        """Test that atoms returns the child expression's atoms."""
        assert Not(Atom("is_snv")).atoms() == {"is_snv"}

    def test_evaluate_partial_settles_only_when_child_settles(self) -> None:
        """Test that partial evaluation settles only when the child settles."""
        expr = Not(Atom("is_snv"))
        assert expr.evaluate_partial({"is_snv": True}) is False
        assert expr.evaluate_partial({"is_snv": False}) is True
        assert expr.evaluate_partial({}) is None


class TestAndOr:
    """Test And and Or expressions."""

    @pytest.mark.parametrize(
        ("known", "expected"),
        [
            ({"a": True, "b": True}, True),
            ({"a": True, "b": False}, False),
            ({"a": False}, False),  # short-circuits: b unknown but a alone settles it
            ({"a": True}, None),
        ],
    )
    def test_and_evaluate_partial_kleene(self, known: dict[str, bool], expected: bool | None) -> None:
        """Test partial evaluation of And expressions with three-valued logic."""
        expr = all_of("a", "b")
        assert expr.evaluate_partial(known) is expected

    @pytest.mark.parametrize(
        ("known", "expected"),
        [
            ({"a": False, "b": False}, False),
            ({"a": True, "b": False}, True),
            ({"a": True}, True),  # short-circuits: a alone settles it
            ({"a": False}, None),
        ],
    )
    def test_or_evaluate_partial_kleene(self, known: dict[str, bool], expected: bool | None) -> None:
        """Test partial evaluation of Or expressions with three-valued logic."""
        expr = any_of("a", "b")
        assert expr.evaluate_partial(known) is expected

    def test_and_atoms_union_children(self) -> None:
        """Test that And expressions return the union of their child atoms."""
        assert all_of("a", "b", no("c")).atoms() == {"a", "b", "c"}

    def test_operator_overloads_build_the_same_shape_as_helpers(self) -> None:
        """Test that operator overloads and helper functions build equivalent expressions."""
        assert repr(Atom("a") & Atom("b")) == repr(all_of("a", "b"))
        assert repr(Atom("a") | Atom("b")) == repr(any_of("a", "b"))
        assert repr(~Atom("a")) == repr(no("a"))

    def test_nested_expression_evaluates_against_a_real_context(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """Test that a nested expression evaluates against a fact context."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        context = _create_snv_fact_context(transcript, position=transcript.start)
        expr = all_of("is_snv", no("is_insertion"), any_of("is_snv", "is_deletion"))
        assert expr.evaluate(context) is True


class TestExprFromJson:
    """Test constructing expressions from JSON-compatible values."""

    def test_bare_string_is_an_atom(self) -> None:
        """Test that a string creates an Atom expression."""
        assert repr(expr_from_json("is_snv")) == "is_snv"

    def test_not_mapping_creates_a_not_expression(self) -> None:
        """Test that a not mapping creates a Not expression."""
        assert repr(expr_from_json({"not": "is_snv"})) == "!is_snv"

    def test_all_mapping_creates_an_and_expression(self) -> None:
        """Test that an all mapping creates an And expression."""
        assert repr(expr_from_json({"all": ["a", "b"]})) == repr(all_of("a", "b"))

    def test_any_mapping_creates_an_or_expression(self) -> None:
        """Test that an any mapping creates an Or expression."""
        assert repr(expr_from_json({"any": ["a", "b"]})) == repr(any_of("a", "b"))

    def test_nested_mapping_creates_a_nested_expression(self) -> None:
        """Test that nested mappings create a nested expression."""
        payload = {"all": ["overlap_cds", {"not": "partial_last_codon"}, {"any": ["is_snv", "is_deletion"]}]}
        expr = expr_from_json(payload)
        assert expr.atoms() == {"overlap_cds", "partial_last_codon", "is_snv", "is_deletion"}

    @pytest.mark.parametrize("bad_payload", [42, None, [], {"all": ["a"], "any": ["b"]}, {"unknown": ["a"]}])
    def test_rejects_malformed_payloads(self, bad_payload: object) -> None:
        """Test that malformed payloads raise ValueError."""
        with pytest.raises(ValueError):
            expr_from_json(bad_payload)
