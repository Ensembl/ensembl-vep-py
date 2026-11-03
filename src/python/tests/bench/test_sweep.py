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
"""Tests for sweeping a setting."""

from __future__ import annotations

import pytest

import ensembl.vep.annotators  # noqa: F401  pylint: disable=unused-import

# ^ registers the annotators and the facts they read.
from ensembl.vep.bench.harness import TimingPolicy
from ensembl.vep.bench.sweep import Sweep, sweep_strategy
from ensembl.vep.corpus.synthetic.default import default_corpus
from ensembl.vep.models.overlap import TranscriptVariationAllele


@pytest.fixture(name="overlaps", scope="module")
def fixture_overlaps() -> list[TranscriptVariationAllele]:
    """Forty boundary overlaps."""
    return [overlap for _name, overlap in default_corpus(fuzz_count=0).named_overlaps[:40]]


class TestSweep:
    """Sweeping one setting."""

    def test_an_unsweepable_setting_is_rejected(self) -> None:
        """Only numeric annotation settings can be swept."""
        with pytest.raises(ValueError, match="readthrough"):
            Sweep(knob="readthrough", values=(1,))

    def test_batch_sizes_are_checked_then_timed(self, overlaps: list[TranscriptVariationAllele]) -> None:
        """Every value agrees with the reference, so every value is timed."""
        points, failures = sweep_strategy(
            "default",
            "consequence.columnar-staged",
            overlaps,
            Sweep(knob="batch_size", values=(16, 64)),
            TimingPolicy(repeats=1, warmup=0),
        )
        assert not failures
        assert [point.configuration for point in points] == [
            "consequence.columnar-staged batch_size=16",
            "consequence.columnar-staged batch_size=64",
        ]
        assert all(point.ms_per_1k > 0 for point in points)

    def test_a_cache_setting_switches_its_cache_on(self, overlaps: list[TranscriptVariationAllele]) -> None:
        """Sweeping a cache's size is meaningless with the cache off, so it is switched on."""
        points, failures = sweep_strategy(
            "default",
            "consequence.staged",
            overlaps,
            Sweep(knob="codon_cache_entries", values=(8,)),
            TimingPolicy(repeats=1, warmup=0),
        )
        assert not failures
        assert points[0].to_json()["value"] == 8

    def test_zone_map_bytes_switches_the_zone_map_on(self, overlaps: list[TranscriptVariationAllele]) -> None:
        """Sweeping the zone map's budget is meaningless with the map off, so it is switched on."""
        points, failures = sweep_strategy(
            "default",
            "consequence.columnar-staged",
            overlaps,
            Sweep(knob="zone_map_bytes", values=(4096,)),
            TimingPolicy(repeats=1, warmup=0),
        )
        assert not failures
        assert points[0].to_json()["value"] == 4096
