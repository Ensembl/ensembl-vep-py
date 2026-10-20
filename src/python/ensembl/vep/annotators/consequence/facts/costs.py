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
"""Measure fact evaluation costs.

Declared tiers provide relative ordering rather than wall-clock timings. Profiling measures each
fact's self-time and normalises it to the ``BOUNDS = 1`` scale used by declared costs.

Profiling temporarily replaces registered compute functions with timing wrappers and therefore
belongs only in explicit measurement runs.
"""

from __future__ import annotations

import json
import statistics
import time
from collections.abc import Callable
from importlib.resources import files
from pathlib import Path

from attrs import evolve

from ensembl.vep.common import CostTier
from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.annotators.consequence.facts.base import (
    FACT_COST_OVERRIDE_RESOURCE,
    Fact,
    REGISTRY,
    get_fact,
)

__all__ = [
    "declared_vs_measured",
    "fact_costs_path",
    "normalised_cost_table",
    "profile_self_nanoseconds",
    "write_cost_table",
]


def profile_self_nanoseconds(fact_contexts: list[FactContext]) -> dict[str, float]:
    """Profile every registered fact's mean self-time, in nanoseconds, over ``fact_contexts``.

    Temporarily wraps every Fact in REGISTRY with a timing shim (restored in a
    ``finally`` block, so a failure partway through still leaves the registry as it found it), then
    computes every fact for every context.

    Args:
        fact_contexts (list[FactContext]): Contexts to compute every fact over. Each is reset
            first.

    Returns:
        dict[str, float]: Fact name -> mean self-time in nanoseconds, over every context in which
            that fact was computed (every fact is computed for every context here, so the divisor
            is ``len(fact_contexts)`` for all of them).
    """
    self_nanoseconds_total: dict[str, float] = {}
    computation_counts: dict[str, int] = {}
    nested_time_stack: list[float] = []
    original_compute_functions = {name: registered_fact.compute for name, registered_fact in REGISTRY.items()}

    def make_timed_compute(
        fact_name: str, original_compute: Callable[[FactContext], object]
    ) -> Callable[[FactContext], object]:
        def timed_compute(fact_context: FactContext) -> object:
            nested_time_stack.append(0.0)
            started_at = time.perf_counter_ns()
            value = original_compute(fact_context)
            total_elapsed_nanoseconds = time.perf_counter_ns() - started_at
            nested_elapsed_nanoseconds = nested_time_stack.pop()
            self_elapsed_nanoseconds = total_elapsed_nanoseconds - nested_elapsed_nanoseconds
            if nested_time_stack:
                nested_time_stack[-1] += total_elapsed_nanoseconds
            self_nanoseconds_total[fact_name] = (
                self_nanoseconds_total.get(fact_name, 0.0) + self_elapsed_nanoseconds
            )
            computation_counts[fact_name] = computation_counts.get(fact_name, 0) + 1
            return value

        return timed_compute

    try:
        for fact_name, registered_fact in list(REGISTRY.items()):
            REGISTRY[fact_name] = _with_compute(
                registered_fact, make_timed_compute(fact_name, registered_fact.compute)
            )

        fact_names = list(REGISTRY)
        for fact_context in fact_contexts:
            fact_context.reset()
            for fact_name in fact_names:
                get_fact(fact_context, fact_name)
    finally:
        for fact_name, original_compute in original_compute_functions.items():
            REGISTRY[fact_name] = _with_compute(REGISTRY[fact_name], original_compute)

    return {
        fact_name: total_nanoseconds / computation_counts[fact_name]
        for fact_name, total_nanoseconds in self_nanoseconds_total.items()
    }


def _with_compute(registered_fact: Fact, compute: Callable[[FactContext], object]) -> Fact:
    """Return a copy of ``registered_fact`` with its ``compute`` replaced.

    A tiny local helper rather than importing ``attrs.evolve`` at call sites throughout
    profile_self_nanoseconds, since Fact is frozen.

    Args:
        registered_fact (Fact): The fact to copy.
        compute (Callable[[FactContext], object]): The replacement ``compute`` callable.

    Returns:
        Fact: The copy, with every other field unchanged.
    """
    return evolve(registered_fact, compute=compute)


def normalised_cost_table(measured_self_nanoseconds: dict[str, float]) -> dict[str, int]:
    """Normalise every fact's measured self-time to BOUNDS = 1 and write it to ``output_path``.

    The inverse of ``apply_fact_cost_overrides``: BOUNDS = 1 is the same normalisation the declared
    tier constants use, so a measured table is interchangeable with the declared ones.

    Args:
        measured_self_nanoseconds (dict[str, float]): Fact name -> mean self-time in nanoseconds.

    Returns:
        dict[str, int]: Fact name -> cost in BOUNDS units.

    Raises:
        RuntimeError: If no BOUNDS-tier fact was measured, leaving nothing to normalise against.
    """
    bounds_facts = [name for name in measured_self_nanoseconds if REGISTRY[name].cost == CostTier.BOUNDS]
    if not bounds_facts:
        raise RuntimeError("no BOUNDS-tier fact was measured; cannot normalise -- widen the corpus")
    bounds_mean_nanoseconds = statistics.mean(measured_self_nanoseconds[name] for name in bounds_facts)

    cost_table = {}
    for fact_name, self_nanoseconds in measured_self_nanoseconds.items():
        cost_table[fact_name] = max(0, round(self_nanoseconds / bounds_mean_nanoseconds))
    return cost_table


def fact_costs_path() -> Path:
    """Return the filesystem path of the packaged measured-cost table.

    Resolved the same way ``load_fact_cost_overrides_from_package`` reads it, so ``--write`` and
    the loader can no longer disagree about where the table lives.
    """
    resource = files("ensembl.vep.annotators.consequence").joinpath(FACT_COST_OVERRIDE_RESOURCE)
    return Path(str(resource))


def write_cost_table(cost_table: dict[str, int], output_path: Path | None = None) -> Path:
    """Write a measured cost table as JSON, and return where it was written.

    Args:
        cost_table (dict[str, int]): Fact name -> cost in BOUNDS units.
        output_path (Path | None): Where to write, or None for the packaged table.

    Returns:
        Path: The path written to.
    """
    destination = output_path or fact_costs_path()
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(cost_table, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return destination


def declared_vs_measured(measured_self_nanoseconds: dict[str, float]) -> list[tuple[str, str, float, int]]:
    """Return one row per fact comparing its declared tier against what was measured.

    Returns:
        list[tuple[str, str, float, int]]: ``(fact name, declared tier label, measured ns, measured
            cost in BOUNDS units)``, ordered by declared cost then name.
    """
    cost_table = normalised_cost_table(measured_self_nanoseconds)
    rows = []
    for fact_name, self_nanoseconds in measured_self_nanoseconds.items():
        declared = REGISTRY[fact_name].cost
        rows.append((fact_name, CostTier.label(declared), self_nanoseconds, cost_table[fact_name]))
    rows.sort(key=lambda row: (REGISTRY[row[0]].cost, row[0]))
    return rows
