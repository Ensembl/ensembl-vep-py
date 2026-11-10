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
"""Tests for fact-level profiling."""

from __future__ import annotations

import pytest

import ensembl.vep.annotators  # noqa: F401  pylint: disable=unused-import

# ^ keeps the annotators-before-caches import order.
from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.annotators.consequence.facts.base import Fact, REGISTRY, get_fact
from ensembl.vep.annotators.factory import build_annotator
from ensembl.vep.corpus.synthetic.default import default_corpus
from ensembl.vep.stats.profile import count_fact_computations, trace_run


class TestRegistryHoldsNoValues:
    """The invariant profiling leans on: a fact's value never lives in the registry."""

    def test_fact_has_no_value_field(self) -> None:
        """Only descriptors: name, cost, compute, deps, requires, documentation."""
        assert "value" not in {attribute.name for attribute in Fact.__attrs_attrs__}

    def test_a_value_lives_in_the_context_memo_only(self) -> None:
        """Two contexts over the same overlap compute independently."""
        _name, overlap = default_corpus(fuzz_count=0).named_overlaps[0]
        first, second = FactContext(overlap=overlap), FactContext(overlap=overlap)
        get_fact(first, "overlap_feature")
        assert "overlap_feature" in first.memo
        assert "overlap_feature" not in second.memo


class TestCountFactComputations:
    """The shim counts, and leaves the registry exactly as it found it."""

    def test_counts_computations_and_cost(self) -> None:
        """A fact computed once is counted once, at its declared cost."""
        _name, overlap = default_corpus(fuzz_count=0).named_overlaps[0]
        with count_fact_computations() as counters:
            get_fact(FactContext(overlap=overlap), "overlap_feature")
        assert counters["computed"] >= 1
        assert counters["cost"] >= REGISTRY["overlap_feature"].cost

    def test_every_compute_is_restored_by_identity(self) -> None:
        """Nothing instrumented survives the block."""
        before = {name: fact.compute for name, fact in REGISTRY.items()}
        with count_fact_computations():
            assert REGISTRY["overlap_feature"].compute is not before["overlap_feature"]
        assert {name: fact.compute for name, fact in REGISTRY.items()} == before

    def test_restored_even_when_the_body_raises(self) -> None:
        """The finally is the contract."""
        before = {name: fact.compute for name, fact in REGISTRY.items()}
        with pytest.raises(RuntimeError):
            with count_fact_computations():
                raise RuntimeError("partway")
        assert {name: fact.compute for name, fact in REGISTRY.items()} == before


class TestTraceRun:
    """Folding a traced run."""

    def test_a_python_annotator_yields_a_fold_over_every_overlap(self) -> None:
        """Every overlap is traced, and the fold names what it profiled."""
        corpus = default_corpus(fuzz_count=0)
        overlaps = [overlap for _name, overlap in corpus.named_overlaps[:20]]
        run_trace = trace_run(
            build_annotator(annotator="consequence.linear"),
            overlaps,
            annotator_name="consequence.linear",
            cost_table="declared tiers",
        )
        assert run_trace.overlaps == 20
        assert run_trace.charged_cost > 0
        assert "consequence.linear" in run_trace.title
        assert "declared tiers" in run_trace.subtitle
        assert run_trace.fact_usage
