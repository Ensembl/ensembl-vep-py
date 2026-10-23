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
"""Tests for the run-scoped rollup of per-overlap derivations."""

from __future__ import annotations

from ensembl.vep.trace.annotation import AnnotationTrace
from ensembl.vep.trace.events import FactEvent, FactOrigin
from ensembl.vep.trace.run import RunTrace


def _trace(*events: FactEvent, terms: tuple[str, ...] = ()) -> AnnotationTrace:
    """Return a trace carrying only the fact events and terms a test cares about."""
    return AnnotationTrace(fact_events=events, emitted_terms=terms)


def _computed(name: str, cost: int, reads: int = 1) -> FactEvent:
    """Return a computed fact event."""
    return FactEvent(name=name, value=True, cost=cost, origin=FactOrigin.COMPUTED, reads=reads)


class TestRunTraceFold:
    """Accumulating many derivations into one view of the run."""

    def test_counts_overlaps_and_terms(self) -> None:
        """Every trace folded is one overlap, and every emitted term is counted once per overlap."""
        run_trace = RunTrace.fold(
            [
                _trace(_computed("overlap_cds", 20), terms=("missense_variant",)),
                _trace(_computed("overlap_cds", 20), terms=("missense_variant", "intron_variant")),
            ]
        )
        assert run_trace.overlaps == 2
        assert run_trace.term_counts == (("missense_variant", 2), ("intron_variant", 1))

    def test_sums_cost_and_reads_per_fact(self) -> None:
        """A fact reached by two overlaps charges twice and carries both overlaps' reads."""
        run_trace = RunTrace.fold(
            [_trace(_computed("overlap_cds", 20, reads=3)), _trace(_computed("overlap_cds", 20, reads=2))]
        )
        usage = run_trace.fact_usage[0]
        assert usage.name == "overlap_cds"
        assert usage.computed == 2
        assert usage.reads == 5
        assert usage.charged_cost == 40
        assert run_trace.charged_cost == 40

    def test_separates_the_uncharged_origins(self) -> None:
        """A cache hit and a forced-False fact are counted apart from a computed one."""
        run_trace = RunTrace.fold(
            [
                _trace(
                    FactEvent(name="overlap_exon", value=True, cost=0, origin=FactOrigin.PREFILLED),
                    FactEvent(name="_codon_alleles", value=(), cost=0, origin=FactOrigin.CACHED),
                    FactEvent(name="overlap_cds", value=False, cost=20, origin=FactOrigin.FORCED_FALSE),
                )
            ]
        )
        usage = {entry.name: entry for entry in run_trace.fact_usage}
        assert usage["overlap_exon"].prefilled == 1
        assert usage["_codon_alleles"].cached == 1
        assert usage["overlap_cds"].forced_false == 1
        assert run_trace.charged_cost == 0
        assert run_trace.avoided_cost == 20

    def test_orders_facts_by_what_they_cost(self) -> None:
        """The expensive facts come first: that is what a reader of this table is looking for."""
        run_trace = RunTrace.fold([_trace(_computed("cheap", 1), _computed("expensive", 60))])
        assert [usage.name for usage in run_trace.fact_usage] == ["expensive", "cheap"]

    def test_tier_reach_names_the_declared_tiers(self) -> None:
        """Tier reach says how far up the cost tiers the run actually went."""
        run_trace = RunTrace.fold([_trace(_computed("cheap", 0), _computed("expensive", 60))])
        assert dict(run_trace.tier_reach) == {"FREE": 1, "SEQUENCE": 1}

    def test_folding_nothing_gives_an_empty_view(self) -> None:
        """Every field is an empty collection, so a renderer needs no branch."""
        run_trace = RunTrace.fold([])
        assert run_trace.overlaps == 0
        assert run_trace.fact_usage == ()
        assert run_trace.term_counts == ()
        assert not run_trace.fact_rows()

    def test_touched_counts_every_origin(self) -> None:
        """A fact needed by an overlap is touched however it became known."""
        run_trace = RunTrace.fold(
            [
                _trace(_computed("overlap_cds", 20)),
                _trace(FactEvent(name="overlap_cds", value=True, cost=0, origin=FactOrigin.PREFILLED)),
            ]
        )
        assert run_trace.fact_usage[0].touched == 2
