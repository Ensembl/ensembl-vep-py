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
"""Tests for rebuilding a derivation from a context the run has finished with."""

from __future__ import annotations

from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.annotators.consequence.rules.expressions import And, Atom, Expr, Not, Or
from ensembl.vep.annotators.consequence.rules.rulesets import FactRule
from ensembl.vep.trace import reconstruct
from ensembl.vep.trace.events import FactOrigin, TraceRecorder


def _fact_rule(term: str, rank: int, when: Expr, subsumes: tuple[str, ...] = ()) -> FactRule:
    """Return a minimal rule, without loading the consequence catalogue."""
    return FactRule(
        term=term, accession=f"SO:{rank:07d}", rank=rank, impact=None, when=when, subsumes=subsumes
    )


def _context(recorder: TraceRecorder | None = None) -> FactContext:
    """Return a context with no real overlap behind it: reconstruction reads only the memo."""
    return FactContext(overlap=None, recorder=recorder)  # type: ignore[arg-type]


class TestMemoFactEvents:
    """Deriving each memoised fact's origin, cost and read count."""

    def test_a_memoised_fact_no_cache_claimed_reads_as_computed(self) -> None:
        """Computed is the default: it is what is left once the caches are accounted for."""
        context = _context(TraceRecorder(read_counts={"overlap_exon": 4}))
        context.memo["overlap_exon"] = True

        event = reconstruct.memo_fact_events(context)[0]
        assert event.origin is FactOrigin.COMPUTED
        assert event.value is True
        assert event.reads == 4
        assert event.cost > 0

    def test_a_served_fact_takes_the_caches_origin_and_charges_nothing(self) -> None:
        """A prefilled fact is uncharged, which is the whole point of the zone map."""
        recorder = TraceRecorder()
        recorder.note_served(["overlap_exon"], FactOrigin.PREFILLED)
        context = _context(recorder)
        context.memo["overlap_exon"] = True

        event = reconstruct.memo_fact_events(context)[0]
        assert event.origin is FactOrigin.PREFILLED
        assert event.cost == 0

    def test_events_follow_memo_insertion_order(self) -> None:
        """The memo's order is the order the derivation materialised its facts in."""
        context = _context(TraceRecorder())
        context.memo["overlap_feature"] = True
        context.memo["overlap_exon"] = False

        names = [event.name for event in reconstruct.memo_fact_events(context)]
        assert names == ["overlap_feature", "overlap_exon"]

    def test_a_context_without_a_recorder_still_reconstructs(self) -> None:
        """No recorder means no read counts and no cache attribution, not an error."""
        context = _context()
        context.memo["overlap_exon"] = True

        event = reconstruct.memo_fact_events(context)[0]
        assert event.origin is FactOrigin.COMPUTED
        assert event.reads == 0


class TestForcedFalseEvents:
    """Recovering the facts requires-propagation proved False without computing them."""

    def test_false_and_absent_from_the_memo_is_forced(self) -> None:
        """The staged assignment plus the memo is enough: nothing recorded it."""
        known = {"overlap_exon": False, "overlap_cds": False}
        events = reconstruct.forced_false_events(known, memo={"overlap_exon": False})

        assert [event.name for event in events] == ["overlap_cds"]
        assert events[0].origin is FactOrigin.FORCED_FALSE
        assert events[0].value is False

    def test_a_forced_fact_carries_the_required_atoms_that_forced_it(self) -> None:
        """The cause comes from the fact's own declared requires, checked against the assignment."""
        known = {"overlap_exon": False, "overlap_cds": False}
        events = reconstruct.forced_false_events(known, memo={"overlap_exon": False})
        assert "overlap_exon" in events[0].forced_by

    def test_a_forced_fact_carries_the_cost_it_avoided(self) -> None:
        """Its declared cost is reported as cost avoided rather than as cost charged."""
        events = reconstruct.forced_false_events({"overlap_cds": False}, memo={})
        assert events[0].cost > 0

    def test_a_true_or_computed_fact_is_not_forced(self) -> None:
        """Only a False name absent from the memo was settled by propagation alone."""
        known = {"overlap_exon": True, "overlap_cds": False}
        events = reconstruct.forced_false_events(known, memo={"overlap_cds": False})
        assert not events


class TestExplainPartialSettlement:
    """The recursive walk over Atom/Not/And/Or that names the settling atom."""

    def test_true_atom_reports_its_value(self) -> None:
        """A settled atom reports itself."""
        assert reconstruct.explain_partial_settlement(Atom("overlap_exon"), {"overlap_exon": True}) == [
            ("overlap_exon", True)
        ]

    def test_unknown_atom_returns_no_reason(self) -> None:
        """An atom the assignment says nothing about settles nothing."""
        assert reconstruct.explain_partial_settlement(Atom("overlap_exon"), {}) == []

    def test_not_recurses_into_child(self) -> None:
        """Negation reports the atom underneath it."""
        assert reconstruct.explain_partial_settlement(Not(Atom("overlap_exon")), {"overlap_exon": False}) == [
            ("overlap_exon", False)
        ]

    def test_and_reports_the_false_conjunct(self) -> None:
        """An And is settled False by its first False conjunct."""
        expr = And(Atom("overlap_cds"), Atom("partial_last_codon"))
        assert reconstruct.explain_partial_settlement(expr, {"overlap_cds": False}) == [
            ("overlap_cds", False)
        ]

    def test_or_reports_the_true_disjunct(self) -> None:
        """An Or is settled True by its first True disjunct."""
        expr = Or(Atom("is_upstream"), Atom("is_downstream"))
        known = {"is_upstream": False, "is_downstream": True}
        assert reconstruct.explain_partial_settlement(expr, known) == [("is_downstream", True)]

    def test_and_all_true_has_no_single_settling_atom(self) -> None:
        """An And settled True by every conjunct has no single atom to blame."""
        assert not reconstruct.explain_partial_settlement(And(Atom("a"), Atom("b")), {"a": True, "b": True})


class TestPlainRuleEvents:
    """Every rule gets its evaluated formula as the reason."""

    def test_records_fired_and_non_firing_rules(self) -> None:
        """Both verdicts are recorded, and neither settled early."""
        rules = (
            _fact_rule("missense_variant", 13, Atom("overlap_cds")),
            _fact_rule("synonymous_variant", 22, Atom("peptide_changed")),
        )
        events = {event.term: event for event in reconstruct.plain_rule_events(rules, ["missense_variant"])}
        assert events["missense_variant"].fired is True
        assert events["missense_variant"].settled_early is False
        assert events["missense_variant"].reason == "overlap_cds"
        assert events["synonymous_variant"].fired is False


class TestGuardedRuleEvents:
    """A cheap-guard rejection is told apart from a formula that was really evaluated."""

    def test_candidate_rule_gets_its_evaluated_formula(self) -> None:
        """A rule that passed the guard was evaluated, so its formula is the reason."""
        rules = (_fact_rule("missense_variant", 13, Atom("overlap_cds")),)
        event = reconstruct.guarded_rule_events(
            rules, {"missense_variant"}, ["missense_variant"], guard_flags={}
        )[0]
        assert event.fired is True
        assert event.settled_early is False
        assert event.reason == "overlap_cds"

    def test_non_candidate_rule_names_the_guard_atoms_that_failed(self) -> None:
        """A rule the guard rejected never had its formula evaluated, and does not claim to have."""
        rule = _fact_rule("stop_gained", 4, And(Atom("overlap_cds"), Atom("alt_peptide_has_ter")))
        event = reconstruct.guarded_rule_events((rule,), set(), [], guard_flags={"overlap_cds": False})[0]
        assert event.fired is False
        assert event.settled_early is True
        assert "guard failed" in event.reason


class TestPartialRuleEvents:
    """Settled-early versus fell-through-to-the-expensive-tier reasons."""

    def test_settled_early_names_the_settling_atom(self) -> None:
        """A rule the assignment settled reports which atom settled it."""
        rule = _fact_rule("stop_gained", 4, Atom("overlap_cds"))
        event = reconstruct.partial_rule_events((rule,), {"overlap_cds": False}, [])[0]
        assert event.fired is False
        assert event.settled_early is True
        assert event.reason == "settled early: overlap_cds=False"

    def test_undecided_falls_through(self) -> None:
        """A rule the assignment could not settle says so instead of quoting a formula."""
        rule = _fact_rule("missense_variant", 13, Atom("alt_peptide_has_ter"))
        event = reconstruct.partial_rule_events((rule,), {}, ["missense_variant"])[0]
        assert event.fired is True
        assert event.settled_early is False
        assert "expensive tier" in event.reason


class TestSubsumptionEvents:
    """Replaying which firing term dropped which, without touching prune()."""

    def test_only_co_firing_subsumption_is_recorded(self) -> None:
        """A term is dropped only when the term subsuming it also fired."""
        by_term = {
            "missense_variant": _fact_rule(
                "missense_variant", 13, Atom("overlap_cds"), subsumes=("coding_sequence_variant",)
            ),
            "coding_sequence_variant": _fact_rule("coding_sequence_variant", 23, Atom("overlap_cds")),
        }
        events = reconstruct.subsumption_events(by_term, ["missense_variant", "coding_sequence_variant"])
        assert len(events) == 1
        assert events[0].dropped_term == "coding_sequence_variant"
        assert events[0].subsumed_by == "missense_variant"

    def test_no_event_when_the_subsumed_term_did_not_fire(self) -> None:
        """Nothing is dropped when the subsumed term never fired in the first place."""
        by_term = {
            "missense_variant": _fact_rule(
                "missense_variant", 13, Atom("overlap_cds"), subsumes=("coding_sequence_variant",)
            )
        }
        assert not reconstruct.subsumption_events(by_term, ["missense_variant"])
