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
"""Run-scoped rollup of per-overlap derivations.

``AnnotationTrace`` is strictly one overlap. This folds a stream of them into what a whole run
looks like: how often each fact was needed, what it charged, how much the caches and ``requires``
propagation saved, how far up the cost tiers the run actually reached, and which terms came out.

Kept a separate type from Phase 6's ``RunSummary`` on purpose. A summary is term-level and costs
about a millisecond a run; this is fact-level and only exists when tracing is on, and merging the
two would hide that difference behind one signature.
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
    """What one fact cost a run.

    Attributes:
        name (str): The fact's name.
        computed (int): Overlaps that computed it and were charged for it.
        prefilled (int): Overlaps a zone map served it to, uncharged.
        cached (int): Overlaps the codon cache served it to, uncharged.
        forced_false (int): Overlaps that proved it False without computing it.
        reads (int): Total ``get_fact`` calls for it across the run.
        charged_cost (int): Declared cost actually charged, summed.
        avoided_cost (int): Declared cost the caches and propagation kept it from charging.
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
        """How many overlaps needed this fact at all, by any origin."""
        return self.computed + self.prefilled + self.cached + self.forced_false


@define(kw_only=True, auto_attribs=True)
class _FactTally:
    """Mutable accumulator behind one FactUsage, so the fold is a single pass."""

    computed: int = field(default=0)
    prefilled: int = field(default=0)
    cached: int = field(default=0)
    forced_false: int = field(default=0)
    reads: int = field(default=0)
    charged_cost: int = field(default=0)
    avoided_cost: int = field(default=0)
    declared_cost: int = field(default=0)


@frozen(kw_only=True, auto_attribs=True)
class RunTrace:
    """Fact-level totals for a whole annotation run.

    Attributes:
        overlaps (int): How many overlaps were traced.
        fact_usage (tuple[FactUsage, ...]): Per-fact totals, most expensive first.
        term_counts (tuple[tuple[str, int], ...]): Emitted term -> overlaps that emitted it, most
            frequent first.
        tier_reach (tuple[tuple[str, int], ...]): Cost tier label -> facts computed at that tier.
        charged_cost (int): Declared cost the run charged in total.
        avoided_cost (int): Declared cost the caches and propagation avoided.
        annotator_name (str): The annotator that was traced, for the report.
        cost_table (str): Which cost table the declared costs came from, for the report -- a
            measured override rewrites them process-wide, so a profile must say which it used.
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
        """The document title."""
        return f"Fact profile: {self.annotator_name}" if self.annotator_name else "Fact profile"

    @property
    def subtitle(self) -> str:
        """The headline figures, and which cost table they are in."""
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
        """Accumulate a stream of per-overlap traces into one run-scoped view.

        Args:
            traces (Iterable[AnnotationTrace]): The derivations to fold, in any order.
            annotator_name (str): The annotator that produced them, for the report.
            cost_table (str): Which cost table their declared costs came from, for the report.

        Returns:
            RunTrace: The totals.
        """
        tallies: dict[str, _FactTally] = {}
        term_totals: dict[str, int] = {}
        overlaps = 0

        for trace in traces:
            overlaps += 1
            for term in trace.emitted_terms:
                term_totals[term] = term_totals.get(term, 0) + 1
            for event in trace.fact_events:
                tally = tallies.get(event.name)
                if tally is None:
                    tally = _FactTally()
                    tallies[event.name] = tally
                tally.reads += event.reads
                if event.origin is FactOrigin.COMPUTED:
                    tally.computed += 1
                    tally.charged_cost += event.cost
                    tally.declared_cost = max(tally.declared_cost, event.cost)
                elif event.origin is FactOrigin.PREFILLED:
                    tally.prefilled += 1
                elif event.origin is FactOrigin.CACHED:
                    tally.cached += 1
                else:
                    tally.forced_false += 1
                    tally.avoided_cost += event.cost
                    tally.declared_cost = max(tally.declared_cost, event.cost)

        usage = []
        for name, tally in tallies.items():
            usage.append(
                FactUsage(
                    name=name,
                    computed=tally.computed,
                    prefilled=tally.prefilled,
                    cached=tally.cached,
                    forced_false=tally.forced_false,
                    reads=tally.reads,
                    charged_cost=tally.charged_cost,
                    avoided_cost=tally.avoided_cost,
                )
            )
        usage.sort(key=lambda entry: (-entry.charged_cost, entry.name))

        tier_totals: dict[str, int] = {}
        for name, tally in tallies.items():
            if tally.computed:
                label = CostTier.label(tally.declared_cost)
                tier_totals[label] = tier_totals.get(label, 0) + tally.computed

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
