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
"""Tests for the strategy x cache matrix."""

from __future__ import annotations

import pytest

import ensembl.vep.annotators  # noqa: F401  pylint: disable=unused-import

# ^ registers the annotators and the facts they read.
from ensembl.vep.bench.harness import TimingPolicy
from ensembl.vep.bench.suites import Configuration, build_configurations, parse_cache_set, run_matrix
from ensembl.vep.corpus.synthetic.default import default_corpus


class TestCacheSets:
    """Parsing a cache set."""

    def test_none_and_joined_sets(self) -> None:
        """``none`` is the empty set; names join with ``+``, spaces allowed."""
        assert parse_cache_set("none") == frozenset()
        assert parse_cache_set("zonemap+codon") == frozenset({"zonemap", "codon"})
        assert parse_cache_set("zonemap + mapper") == frozenset({"zonemap", "mapper"})

    def test_an_unknown_cache_is_named(self) -> None:
        """The error says what was not understood."""
        with pytest.raises(ValueError, match="bogus"):
            parse_cache_set("zonemap+bogus")


class TestConfigurations:
    """Crossing strategies with cache sets."""

    def test_a_label_names_the_caches_in_a_fixed_order(self) -> None:
        """The same set always reads the same way."""
        configuration = Configuration(strategy="consequence.staged", caches=frozenset({"codon", "zonemap"}))
        assert configuration.label == "consequence.staged +zonemap+codon"
        assert Configuration(strategy="consequence.linear").label == "consequence.linear"

    def test_strategies_are_outermost_and_repeats_dropped(self) -> None:
        """The first configuration is the first strategy with the first cache set."""
        configurations = build_configurations(
            ["consequence.linear", "consequence.staged"], ["none", "zonemap", "none"], 256
        )
        assert [c.label for c in configurations] == [
            "consequence.linear",
            "consequence.linear +zonemap",
            "consequence.staged",
            "consequence.staged +zonemap",
        ]

    def test_each_build_gets_fresh_caches(self) -> None:
        """No configuration is ever timed against a cache another one warmed."""
        configuration = Configuration(strategy="consequence.staged", caches=frozenset({"zonemap"}))
        first, second = configuration.build(), configuration.build()
        first_cache = first.annotator_context.zone_map_cache  # type: ignore[attr-defined]
        second_cache = second.annotator_context.zone_map_cache  # type: ignore[attr-defined]
        assert first_cache is not None
        assert first_cache is not second_cache


class TestRunMatrix:
    """Measuring a small matrix."""

    def test_every_configuration_gets_a_row(self) -> None:
        """Rows carry the corpus, positive timings and cost, and the caches' counters."""
        overlaps = [overlap for _name, overlap in default_corpus(fuzz_count=0).named_overlaps[:30]]
        configurations = [
            Configuration(strategy="consequence.linear"),
            Configuration(strategy="consequence.staged", caches=frozenset({"zonemap"})),
        ]
        rows = run_matrix("default", overlaps, configurations, TimingPolicy(repeats=1, warmup=0))
        assert [row.configuration for row in rows] == ["consequence.linear", "consequence.staged +zonemap"]
        for row in rows:
            assert row.corpus == "default"
            assert row.contexts == 30
            assert row.ms_per_1k > 0
            assert row.batch_ms_per_1k > 0
            assert row.cost_per_context > 0
        assert "zone map maps built" in rows[1].extras
        assert set(rows[0].to_json()) >= {"corpus", "configuration", "ms_per_1k", "cost_per_context"}
