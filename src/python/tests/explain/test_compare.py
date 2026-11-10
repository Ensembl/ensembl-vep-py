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
"""Tests for comparing two derivations, or one against somebody else's terms."""

from __future__ import annotations

from ensembl.vep.trace.annotation import AnnotationTrace
from ensembl.vep.trace.events import FactEvent, FactOrigin, PruneEvent, RuleEvent
from ensembl.vep.explain.compare import diff_terms, diff_traces


def _rule_event(term: str, rank: int, fired: bool, reason: str = "formula") -> RuleEvent:
    """Return a rule event with the fields a diff reads."""
    return RuleEvent(term=term, rank=rank, fired=fired, settled_early=False, reason=reason)


def _fact_event(name: str, value: object) -> FactEvent:
    """Return a computed fact event carrying a value."""
    return FactEvent(name=name, value=value, cost=1, origin=FactOrigin.COMPUTED)


class TestDiffTraces:
    """Where two annotators' derivations of one overlap part company."""

    def test_agreeing_traces_report_agreement(self) -> None:
        """Two annotators that emit the same terms and the same verdicts agree everywhere."""
        events = (_rule_event("missense_variant", 13, True),)
        trace_a = AnnotationTrace(
            rule_events=events, emitted_terms=("missense_variant",), annotator_name="consequence.linear"
        )
        trace_b = AnnotationTrace(
            rule_events=events, emitted_terms=("missense_variant",), annotator_name="consequence.staged"
        )
        report = diff_traces(trace_a, trace_b)
        assert "emitted terms AGREE" in report
        assert "no rule-verdict divergence" in report

    def test_divergent_terms_are_reported_with_both_sides(self) -> None:
        """When the outputs differ, both term lists appear so the difference is readable."""
        trace_a = AnnotationTrace(emitted_terms=("missense_variant",), annotator_name="a")
        trace_b = AnnotationTrace(emitted_terms=("synonymous_variant",), annotator_name="b")
        report = diff_traces(trace_a, trace_b)
        assert "emitted terms DIVERGE" in report
        assert "missense_variant" in report
        assert "synonymous_variant" in report

    def test_the_first_divergent_rule_is_reported_by_rank(self) -> None:
        """Rank order is severity order, so the first disagreement is the one that matters."""
        trace_a = AnnotationTrace(
            rule_events=(_rule_event("stop_gained", 4, True), _rule_event("missense_variant", 13, True)),
            annotator_name="a",
        )
        trace_b = AnnotationTrace(
            rule_events=(_rule_event("stop_gained", 4, False), _rule_event("missense_variant", 13, False)),
            annotator_name="b",
        )
        report = diff_traces(trace_a, trace_b)
        assert "FIRST RULE-VERDICT DIVERGENCE: [4] stop_gained" in report

    def test_a_fact_both_sides_valued_differently_is_reported(self) -> None:
        """A disagreement upstream of the rules is the one worth finding first."""
        trace_a = AnnotationTrace(fact_events=(_fact_event("overlap_cds", True),), annotator_name="a")
        trace_b = AnnotationTrace(fact_events=(_fact_event("overlap_cds", False),), annotator_name="b")
        assert "FACT VALUE DIVERGENCE" in diff_traces(trace_a, trace_b)

    def test_a_fact_only_one_side_reached_is_not_a_divergence(self) -> None:
        """Staging exists to avoid facts, so touching fewer of them is correct, not a difference."""
        trace_a = AnnotationTrace(
            fact_events=(_fact_event("overlap_cds", True), _fact_event("peptide_changed", True)),
            annotator_name="a",
        )
        trace_b = AnnotationTrace(fact_events=(_fact_event("overlap_cds", True),), annotator_name="b")
        assert "FACT VALUE DIVERGENCE" not in diff_traces(trace_a, trace_b)

    def test_rules_only_one_side_recorded_are_listed(self) -> None:
        """A rule one annotator never reached is reported rather than silently dropped."""
        trace_a = AnnotationTrace(rule_events=(_rule_event("stop_gained", 4, False),), annotator_name="a")
        trace_b = AnnotationTrace(annotator_name="b")
        assert "rules traced only by a" in diff_traces(trace_a, trace_b)


class TestDiffTerms:
    """Explaining a derivation against a term list something else assigned."""

    def test_agreement_is_reported_plainly(self) -> None:
        """Nothing to explain when the two term lists match."""
        trace = AnnotationTrace(emitted_terms=("missense_variant",))
        assert "terms AGREE" in diff_terms(trace, ["missense_variant"], "VEP CSQ")

    def test_a_term_we_did_not_emit_reports_our_rules_verdict(self) -> None:
        """The rule ledger turns a bare disagreement into a specific FALSE verdict."""
        trace = AnnotationTrace(
            rule_events=(_rule_event("stop_gained", 4, False, reason="overlap_cds=False"),),
            emitted_terms=(),
        )
        report = diff_terms(trace, ["stop_gained"], "VEP CSQ")
        assert "our rule was FALSE -- overlap_cds=False" in report

    def test_a_term_we_pruned_reports_the_subsumption_edge(self) -> None:
        """A term we fired and then dropped is a different finding from one we never fired."""
        trace = AnnotationTrace(
            rule_events=(_rule_event("coding_sequence_variant", 23, True),),
            prune_events=(
                PruneEvent(dropped_term="coding_sequence_variant", subsumed_by="missense_variant"),
            ),
            emitted_terms=("missense_variant",),
        )
        report = diff_terms(trace, ["coding_sequence_variant", "missense_variant"], "VEP CSQ")
        assert "subsumed by missense_variant" in report

    def test_a_term_the_catalogue_does_not_model_is_unexplainable(self) -> None:
        """The catalogue is a documented subset, so this reports rather than blames a rule."""
        trace = AnnotationTrace(emitted_terms=())
        report = diff_terms(trace, ["mature_miRNA_variant"], "VEP CSQ")
        assert "not modelled by this catalogue" in report

    def test_a_term_only_we_emit_reports_the_rule_that_fired(self) -> None:
        """The other direction is answered from the same ledger."""
        trace = AnnotationTrace(
            rule_events=(_rule_event("intron_variant", 28, True, reason="overlap_intron"),),
            emitted_terms=("intron_variant",),
        )
        report = diff_terms(trace, [], "VEP CSQ")
        assert "only ours: intron_variant -- fired on overlap_intron" in report
