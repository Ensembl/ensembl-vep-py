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
"""Tests for fact base utilities."""

from __future__ import annotations

import importlib

import pytest

from ensembl.vep.annotators.consequence import facts as facts_package
from ensembl.vep.annotators.consequence.facts.base import (
    REGISTRY,
    Fact,
    apply_fact_cost_overrides,
    apply_packaged_fact_cost_overrides,
    atom_names,
    closure_cost,
    get_fact,
    implied_false,
    load_fact_cost_overrides_from_package,
    register_fact,
    total_cost,
    verify_requires,
)
from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.trace import FactOrigin, TraceRecorder
from ensembl.vep.trace import reconstruct
from ensembl.vep.common import CostTier, Strand
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.models.variation import VariationAllele


def _create_transcript_variation_allele(transcript: Transcript) -> TranscriptVariationAllele:
    """Create a single-nucleotide transcript variation allele."""
    variation_allele = VariationAllele(
        start=transcript.start,
        end=transcript.start,
        strand=Strand.FORWARD,
        slice=transcript.slice,
        reference_allele="A",
        alternate_allele="T",
    )
    return TranscriptVariationAllele(feature=transcript, overlapping_feature=variation_allele)


def _create_non_overlapping_fact_context(transcript: Transcript) -> FactContext:
    """A context whose variant misses the transcript entirely, so overlap_feature is False."""
    far_outside_transcript_position = transcript.start - 100_000
    variation_allele = VariationAllele(
        start=far_outside_transcript_position,
        end=far_outside_transcript_position,
        strand=Strand.FORWARD,
        slice=transcript.slice,
        reference_allele="A",
        alternate_allele="T",
    )
    overlap = TranscriptVariationAllele(feature=transcript, overlapping_feature=variation_allele)
    return FactContext(overlap=overlap)


class TestFactDecorator:
    """The fact decorator registers a Fact and enforces its requires."""

    def test_the_decorator_enforces_requires_rather_than_trusting_the_body(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
    ) -> None:
        """A body that ignores its own ``requires`` is still False when the requirement fails.

        The body below unconditionally returns True; the decorator's guard is the only thing
        standing between that and a wrong answer, and it never runs the body at all.

        Removed from REGISTRY again in a finally block so this test cannot leak state into others.
        """
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]

        @register_fact(CostTier.FREE, requires=("overlap_feature",))
        def _always_true_despite_requires(_fact_context: FactContext) -> bool:
            """A fact that ignores overlap_feature entirely and always returns True."""
            return True

        try:
            fact_context = _create_non_overlapping_fact_context(transcript)
            assert get_fact(fact_context, "overlap_feature") is False
            assert get_fact(fact_context, "_always_true_despite_requires") is False
        finally:
            del REGISTRY["_always_true_despite_requires"]

    def test_the_decorator_folds_requires_into_deps(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
    ) -> None:
        """An atom that must be *true* first must also be *computed* first, so deps absorbs it.

        Declaring the same name in both places was only ever a chance to declare them inconsistently,
        so the decorator does it and callers pass ``requires`` alone.
        """
        _ = synthetic_transcripts_by_stable_id

        @register_fact(CostTier.FREE, deps=("is_snv",), requires=("overlap_feature",))
        def _folds_requires_into_deps(_fact_context: FactContext) -> bool:
            """A fact declaring a dep and a separate requirement."""
            return True

        try:
            registered = REGISTRY["_folds_requires_into_deps"]
            assert set(registered.deps) == {"is_snv", "overlap_feature"}
            assert registered.requires == ("overlap_feature",)
        finally:
            del REGISTRY["_folds_requires_into_deps"]

    def test_fact_decorator_records_docstring_as_documentation(self) -> None:
        """Fact decorator records docstring as documentation."""
        assert isinstance(REGISTRY["is_snv"], Fact)
        assert REGISTRY["is_snv"].documentation


class TestGetFact:
    """get_fact returns the computed value of a fact, or False if its requires are not met."""

    def test_get_fact_returns_computed_value(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """get_fact returns the computed value of a fact."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        context = FactContext(overlap=_create_transcript_variation_allele(transcript))
        assert get_fact(context, "is_snv") is True

    def test_get_fact_returns_false_if_requires_not_met(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """get_fact returns False if requires are not met."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        context = FactContext(overlap=_create_transcript_variation_allele(transcript))
        assert get_fact(context, "overlap_cds") is False

    def test_get_fact_memoises(self, synthetic_transcripts_by_stable_id: dict[str, Transcript]) -> None:
        """get_fact computes a fact once and serves the memo thereafter."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        context = FactContext(
            overlap=_create_transcript_variation_allele(transcript), recorder=TraceRecorder()
        )
        assert get_fact(context, "is_short") is True
        assert get_fact(context, "is_short") is True
        computed = [
            event
            for event in reconstruct.memo_fact_events(context)
            if event.name == "is_short" and event.origin is FactOrigin.COMPUTED
        ]
        assert len(computed) == 1, "memoised, not recomputed"


class TestGetFactTracing:
    """Tracing records reads, not recomputation."""

    def test_no_recorder_by_default_and_untouched_by_get_fact(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """Tracing is fully opt-in."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        context = FactContext(overlap=_create_transcript_variation_allele(transcript))
        assert context.recorder is None
        get_fact(context, "is_snv")
        assert context.recorder is None

    def test_counts_every_read_including_memo_hits(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """Read counts include memo hits."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        context = FactContext(
            overlap=_create_transcript_variation_allele(transcript), recorder=TraceRecorder()
        )
        get_fact(context, "is_snv")
        get_fact(context, "is_snv")
        get_fact(context, "is_snv")
        assert context.recorder is not None
        assert context.recorder.read_counts["is_snv"] == 3

    def test_a_fact_read_many_times_is_one_event_carrying_its_read_count(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """One event per fact carries the read count."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        context = FactContext(
            overlap=_create_transcript_variation_allele(transcript), recorder=TraceRecorder()
        )
        get_fact(context, "is_snv")
        get_fact(context, "is_snv")

        matching_events = [event for event in reconstruct.memo_fact_events(context) if event.name == "is_snv"]
        assert len(matching_events) == 1
        assert matching_events[0].origin is FactOrigin.COMPUTED
        assert matching_events[0].reads == 2

    def test_a_computed_fact_carries_its_declared_cost(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """Computed facts carry their declared cost."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        context = FactContext(
            overlap=_create_transcript_variation_allele(transcript), recorder=TraceRecorder()
        )
        get_fact(context, "is_short")

        event = reconstruct.memo_fact_events(context)[0]
        assert event.name == "is_short"
        assert event.value is True
        assert event.cost == REGISTRY["is_short"].cost


class TestVerifyRequires:
    """verify_requires is the backstop for a Fact built bypassing the decorator."""

    def test_verify_requires_for_a_well_behaved_fact(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
    ) -> None:
        """verify_requires returns known problematic terms for a well-behaved Fact."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        problems = verify_requires([_create_non_overlapping_fact_context(transcript)])
        problematic_terms = [problem.split()[0] for problem in problems]

        # The known problematic terms are the facts that ignore their own requires and return True.
        # But that can be okay if we do not care about the fact body without the requires guard.
        assert not problematic_terms

    def test_verify_requires_catches_a_fact_that_ignores_its_requires(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
    ) -> None:
        """A body that ignores its own ``requires`` will be flagged by verify_requires."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]

        @register_fact(CostTier.FREE, requires=("overlap_feature",))
        def _always_true_despite_requires(_fact_context: FactContext) -> bool:
            """A fact that ignores overlap_feature entirely and always returns True."""
            return True

        try:
            fact_context = _create_non_overlapping_fact_context(transcript)
            assert get_fact(fact_context, "overlap_feature") is False
            assert get_fact(fact_context, "_always_true_despite_requires") is False
            problems = verify_requires([_create_non_overlapping_fact_context(transcript)])
            assert any("_always_true_despite_requires" in problem for problem in problems)
        finally:
            del REGISTRY["_always_true_despite_requires"]

    def test_verify_requires_catches_a_fact_registered_without_the_guard(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
    ) -> None:
        """verify_requires remains the backstop for a Fact built bypassing the decorator.

        Enforcement makes the check close to tautological for anything registered the normal way, so
        this constructs the Fact directly -- which is also what a bug *in the guard wrapper itself*
        would look like from verify_requires' point of view. That is the one failure it can still
        see, and the reason it is kept.
        """
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]

        REGISTRY["_unguarded_broken_requires"] = Fact(
            name="_unguarded_broken_requires",
            cost=CostTier.FREE,
            compute=lambda _fact_context: True,
            requires=("overlap_feature",),
        )
        try:
            problems = verify_requires([_create_non_overlapping_fact_context(transcript)])
            assert any("_unguarded_broken_requires" in problem for problem in problems)
        finally:
            del REGISTRY["_unguarded_broken_requires"]


class TestAtomNames:
    """Atom names are the keys of REGISTRY, minus any internal facts."""

    def test_atom_names_are_sorted(self) -> None:
        """Atom names are sorted."""
        names = atom_names()
        assert names == sorted(names)

    def test_atom_names_are_unique(self) -> None:
        """Atom names are unique."""
        names = atom_names()
        assert len(names) == len(set(names))

    def test_atom_names_excludes_internal_facts(self) -> None:
        """Atom names excludes internal facts."""
        names = atom_names()
        assert all(not name.startswith("_") for name in names)
        assert names == sorted(names)


class TestImpliedFalse:
    """Implied false facts are added to a fact context when a fact is False."""

    def test_implied_false_propagates_through_a_requires_chain(self) -> None:
        """Implied false propagates through a requires chain."""
        extended = implied_false({"overlap_feature": False})
        assert extended["overlap_exon"] is False
        assert extended["overlap_cds"] is False

    def test_implied_false_leaves_unrelated_facts_alone(self) -> None:
        """Implied false leaves unrelated facts alone."""
        extended = implied_false({"overlap_feature": True})
        assert "overlap_exon" not in extended


class TestClosureCost:
    """Closure cost is the sum of a fact's own cost and all its transitive dependencies' costs."""

    def test_closure_cost_is_at_least_own_cost(self) -> None:
        """Closure cost is at least own cost."""
        for name in atom_names():
            assert closure_cost(name) >= REGISTRY[name].cost

    def test_closure_cost_includes_transitive_dependencies(self) -> None:
        """Closure cost includes transitive dependencies."""
        assert closure_cost("overlap_cds") > REGISTRY["overlap_cds"].cost


class TestTotalCost:
    """Total cost is the sum of closure costs for a set of facts."""

    def test_total_cost_is_at_least_sum_of_own_costs(self) -> None:
        """Total cost is at least sum of own costs."""
        names = ["overlap_feature", "overlap_exon", "overlap_cds"]
        total = total_cost(names)
        own_sum = sum(REGISTRY[name].cost for name in names)
        assert total >= own_sum

    def test_total_cost_counts_shared_dependencies_once(self) -> None:
        """Total cost counts shared dependencies once."""
        combined = total_cost(["peptide_changed", "alt_peptide_has_ter"])
        assert combined < closure_cost("peptide_changed") + closure_cost("alt_peptide_has_ter")


class TestFactCostOverrides:
    """The optional, off-by-default measured-cost override (vep catalogue costs)."""

    def test_apply_fact_cost_overrides_replaces_declared_cost(self) -> None:
        """Apply fact cost overrides replaces declared cost."""
        original_cost = REGISTRY["is_snv"].cost
        try:
            apply_fact_cost_overrides({"is_snv": original_cost + 7})
            assert REGISTRY["is_snv"].cost == original_cost + 7
        finally:
            apply_fact_cost_overrides({"is_snv": original_cost})
        assert REGISTRY["is_snv"].cost == original_cost

    def test_apply_fact_cost_overrides_ignores_unknown_names(self) -> None:
        """An override for a fact that does not exist is silently ignored, not an error."""
        apply_fact_cost_overrides({"this_fact_does_not_exist": 99})
        assert "this_fact_does_not_exist" not in REGISTRY

    def test_apply_fact_cost_overrides_invalidates_closure_cost_cache(self) -> None:
        """A stale closure_cost must not survive an override of one of its dependencies."""
        original_cost = REGISTRY["overlap_feature"].cost
        stale_closure_cost = closure_cost("overlap_exon")  # populate the cache before overriding
        try:
            apply_fact_cost_overrides({"overlap_feature": original_cost + 100})
            assert closure_cost("overlap_exon") == stale_closure_cost + 100
        finally:
            apply_fact_cost_overrides({"overlap_feature": original_cost})

    def test_apply_packaged_fact_cost_overrides_applies_the_packaged_table(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """The packaged table replaces the declared cost of the facts it names."""
        original_cost = REGISTRY["is_snv"].cost
        monkeypatch.setattr(
            "ensembl.vep.annotators.consequence.facts.base.load_fact_cost_overrides_from_package",
            lambda: {"is_snv": original_cost + 123},
        )
        try:
            apply_packaged_fact_cost_overrides()
            assert REGISTRY["is_snv"].cost == original_cost + 123
        finally:
            apply_fact_cost_overrides({"is_snv": original_cost})

    def test_apply_packaged_fact_cost_overrides_ignores_the_environment(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Importing the facts package never reads ENSEMBL_VEP_FACT_COST_OVERRIDE; the config does."""
        monkeypatch.setenv("ENSEMBL_VEP_FACT_COST_OVERRIDE", "1")
        monkeypatch.setattr(
            "ensembl.vep.annotators.consequence.facts.base.load_fact_cost_overrides_from_package",
            lambda: {"is_snv": REGISTRY["is_snv"].cost + 123},
        )
        original_cost = REGISTRY["is_snv"].cost

        importlib.reload(facts_package)

        assert REGISTRY["is_snv"].cost == original_cost

    def test_load_fact_cost_overrides_from_package_returns_empty_when_resource_missing(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Load fact cost overrides from package returns empty when resource missing."""
        monkeypatch.setattr(
            "ensembl.vep.annotators.consequence.facts.base.FACT_COST_OVERRIDE_RESOURCE",
            "data/this_file_does_not_exist.json",
        )
        assert not load_fact_cost_overrides_from_package()
