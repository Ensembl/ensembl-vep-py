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
"""Run-level aggregation of annotation traces.

``AnnotationTrace`` describes one overlap. ``RunTrace`` aggregates fact use, charged and avoided
costs, reached cost tiers, and emitted terms across a run. It remains separate from the lightweight
term-level ``RunSummary`` because fact-level tracing is optional and more expensive.
"""

from __future__ import annotations

from collections.abc import Iterable

from attrs import define, field, frozen

from ensembl.vep.common import CostTier
from ensembl.vep.render.chart import Slice
from ensembl.vep.render.report import ChartKind, ReportSection
from ensembl.vep.trace.annotation import AnnotationTrace
from ensembl.vep.trace.events import FactOrigin

__all__ = ["FACT_HEADERS", "FactUsage", "RunTrace"]

FACT_HEADERS = ("Fact", "Computed", "Prefilled", "Cached", "Forced False", "Reads", "Cost", "Avoided")
"""Column headers matching ``RunTrace.fact_rows``."""


@frozen(kw_only=True, auto_attribs=True)
class FactUsage:
    """Aggregate one fact's use across a run.

    Attributes:
        name (str): The fact's name.
        computed (int): Overlaps that computed the fact.
        prefilled (int): Overlaps supplied by a zone map.
        cached (int): Overlaps supplied by the codon cache.
        forced_false (int): Overlaps that proved the fact false without computation.
        reads (int): Total ``get_fact`` requests.
        charged_cost (int): Total declared cost charged.
        avoided_cost (int): Total declared cost avoided.
    """

    name: str
    computed: int = field(default=0)
    prefilled: int = field(default=0)
    cached: int = field(default=0)
    forced_false: int = field(default=0)
    reads: int = field(default=0)
    charged_cost: int = field(default=0)
    avoided_cost: int = field(default=0)

    @property
    def touched(self) -> int:
        """Return the number of overlaps that used this fact."""
        return self.computed + self.prefilled + self.cached + self.forced_false


@define(kw_only=True, auto_attribs=True)
class _FactUsageBuilder:
    """Accumulate one fact's usage during a fold, then freeze it to a ``FactUsage``."""

    name: str
    computed: int = field(default=0)
    prefilled: int = field(default=0)
    cached: int = field(default=0)
    forced_false: int = field(default=0)
    reads: int = field(default=0)
    charged_cost: int = field(default=0)
    avoided_cost: int = field(default=0)
    declared_cost: int = field(default=0)

    def freeze(self) -> FactUsage:
        """Return the accumulated totals as a ``FactUsage``."""
        return FactUsage(
            name=self.name,
            computed=self.computed,
            prefilled=self.prefilled,
            cached=self.cached,
            forced_false=self.forced_false,
            reads=self.reads,
            charged_cost=self.charged_cost,
            avoided_cost=self.avoided_cost,
        )


@frozen(kw_only=True, auto_attribs=True)
class RunTrace:
    """Record fact-level totals for an annotation run.

    Attributes:
        overlaps (int): Number of traced overlaps.
        fact_usage (tuple[FactUsage, ...]): Per-fact totals ordered by charged cost.
        term_counts (tuple[tuple[str, int], ...]): Emitted-term counts ordered by frequency.
        tier_reach (tuple[tuple[str, int], ...]): Computed-fact counts by cost tier.
        charged_cost (int): Total declared cost charged.
        avoided_cost (int): Total declared cost avoided.
        annotator_name (str): Annotator represented by the trace.
        cost_table (str): Cost table used for declared fact costs.
    """

    overlaps: int = field(default=0)
    fact_usage: tuple[FactUsage, ...] = field(default=())
    term_counts: tuple[tuple[str, int], ...] = field(default=())
    tier_reach: tuple[tuple[str, int], ...] = field(default=())
    charged_cost: int = field(default=0)
    avoided_cost: int = field(default=0)
    annotator_name: str = field(default="")
    cost_table: str = field(default="")

    @property
    def title(self) -> str:
        """Return the report title."""
        return f"Fact profile: {self.annotator_name}" if self.annotator_name else "Fact profile"

    @property
    def subtitle(self) -> str:
        """Return the report totals and cost-table description."""
        table = f"; costs from {self.cost_table}" if self.cost_table else ""
        return (
            f"{self.overlaps} overlaps, cost {self.charged_cost} charged, "
            f"{self.avoided_cost} avoided{table}"
        )

    def sections(self) -> tuple[ReportSection, ...]:
        """Return the profile as report sections."""
        return (
            ReportSection(
                identifier="facts", title="Facts", headers=FACT_HEADERS, rows=tuple(self.fact_rows())
            ),
            ReportSection(
                identifier="tiers",
                title="Facts computed, by cost tier",
                slices=tuple(Slice(label=label, value=count) for label, count in self.tier_reach),
                chart=ChartKind.BARS,
            ),
            ReportSection(
                identifier="terms",
                title="Emitted terms",
                slices=tuple(Slice(label=term, value=count) for term, count in self.term_counts),
                chart=ChartKind.BARS,
            ),
        )

    @classmethod
    def fold(
        cls, traces: Iterable[AnnotationTrace], *, annotator_name: str = "", cost_table: str = ""
    ) -> RunTrace:
        """Build run-level totals from per-overlap traces.

        Args:
            traces (Iterable[AnnotationTrace]): Derivations to aggregate.
            annotator_name (str): Annotator that produced the traces.
            cost_table (str): Cost table used for declared fact costs.

        Returns:
            RunTrace: Aggregated run trace.
        """
        builders: dict[str, _FactUsageBuilder] = {}
        term_totals: dict[str, int] = {}
        overlaps = 0

        for trace in traces:
            overlaps += 1
            for term in trace.emitted_terms:
                term_totals[term] = term_totals.get(term, 0) + 1
            for event in trace.fact_events:
                builder = builders.get(event.name)
                if builder is None:
                    builder = _FactUsageBuilder(name=event.name)
                    builders[event.name] = builder
                builder.reads += event.reads
                if event.origin is FactOrigin.COMPUTED:
                    builder.computed += 1
                    builder.charged_cost += event.cost
                    builder.declared_cost = max(builder.declared_cost, event.cost)
                elif event.origin is FactOrigin.PREFILLED:
                    builder.prefilled += 1
                elif event.origin is FactOrigin.CACHED:
                    builder.cached += 1
                else:
                    builder.forced_false += 1
                    builder.avoided_cost += event.cost
                    builder.declared_cost = max(builder.declared_cost, event.cost)

        usage = [builder.freeze() for builder in builders.values()]
        usage.sort(key=lambda entry: (-entry.charged_cost, entry.name))

        tier_totals: dict[str, int] = {}
        for builder in builders.values():
            if builder.computed:
                label = CostTier.label(builder.declared_cost)
                tier_totals[label] = tier_totals.get(label, 0) + builder.computed

        terms = sorted(term_totals.items(), key=lambda entry: (-entry[1], entry[0]))
        tiers = sorted(tier_totals.items(), key=lambda entry: (-entry[1], entry[0]))
        return cls(
            overlaps=overlaps,
            fact_usage=tuple(usage),
            term_counts=tuple(terms),
            tier_reach=tuple(tiers),
            charged_cost=sum(entry.charged_cost for entry in usage),
            avoided_cost=sum(entry.avoided_cost for entry in usage),
            annotator_name=annotator_name,
            cost_table=cost_table,
        )

    def fact_rows(self) -> list[tuple[object, ...]]:
        """Return the per-fact totals as table rows, in ``FACT_HEADERS`` order."""
        rows: list[tuple[object, ...]] = []
        for usage in self.fact_usage:
            rows.append(
                (
                    usage.name,
                    usage.computed,
                    usage.prefilled,
                    usage.cached,
                    usage.forced_false,
                    usage.reads,
                    usage.charged_cost,
                    usage.avoided_cost,
                )
            )
        return rows
