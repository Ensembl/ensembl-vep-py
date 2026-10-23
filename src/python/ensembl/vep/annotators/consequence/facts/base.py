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
"""Lazy, memoised, cost-tagged facts for consequence annotation.

Facts form a dependency graph and are computed at most once per ``FactContext``. Declared costs
allow strategies to order work and measure avoided computation.

Registered functions have two forms:

* Boolean facts have public names and may appear in catalogue formulas.
* Internal value facts have private names and share non-Boolean intermediate values between facts.

Internal value facts must not declare ``requires`` because a failed guard returns ``False``.
Keep single-use calculations as plain helpers; registration only benefits values shared by
multiple facts.
"""

from __future__ import annotations

import json
from collections.abc import Callable, Iterable, Mapping
from functools import wraps
from importlib.resources import files
from typing import TYPE_CHECKING, Any

from attrs import Factory, evolve, field, frozen

if TYPE_CHECKING:
    from ensembl.vep.annotators.consequence.context.fact import FactContext

__all__ = [
    "FACT_COST_OVERRIDE_RESOURCE",
    "Fact",
    "REGISTRY",
    "apply_fact_cost_overrides",
    "apply_packaged_fact_cost_overrides",
    "atom_names",
    "build_requires_guard",
    "closure_cost",
    "get_fact",
    "implied_false",
    "load_fact_cost_overrides_from_package",
    "register_fact",
    "total_cost",
    "verify_requires",
]


@frozen(kw_only=True, auto_attribs=True)
class Fact:
    """Represent a named, cost-tagged property derived from a context.

    ``deps`` names values the implementation reads. ``requires`` names Boolean preconditions;
    these are also dependencies and allow evaluation to return ``False`` before running the fact.

    Attributes:
        name (str): Registry key and rule atom name.
        cost (int): Marginal evaluation cost, excluding dependencies.
        compute (Callable[[FactContext[Any]], object]): Guarded computation called by
            ``get_fact``.
        unguarded_compute (Callable[[FactContext[Any]], object]): Original computation used to
            validate ``requires`` declarations.
        deps (tuple[str, ...]): Facts read by the computation.
        requires (tuple[str, ...]): Boolean preconditions that gate the computation.
        documentation (str): Description copied from the registered function.
    """

    name: str
    cost: int
    compute: Callable[[FactContext[Any]], object]
    unguarded_compute: Callable[[FactContext[Any]], object] = field(
        default=Factory(lambda self: self.compute, takes_self=True)
    )
    deps: tuple[str, ...] = field(default=())
    requires: tuple[str, ...] = field(default=())
    documentation: str = field(default="")


REGISTRY: dict[str, Fact] = {}
"""Every registered fact, keyed by name. Populated by importing the ``facts`` subpackage via the decorator."""


def register_fact(
    cost: int, deps: tuple[str, ...] = (), requires: tuple[str, ...] = ()
) -> Callable[[Callable[[FactContext[Any]], object]], Callable[[FactContext[Any]], object]]:
    """Register a function as a named Fact in REGISTRY.

    ``requires`` names are folded into ``deps`` automatically because an atom that has to
    be *true* first necessarily has to be *computed* first, so declaring it twice was only ever
    an opportunity to declare it inconsistently.

    Args:
        cost (int): The fact's marginal cost.
        deps (tuple[str, ...]): Names of facts that must be computed before this one.
        requires (tuple[str, ...]): Atoms this fact is gated on.

    Returns:
        Callable[[Callable[[FactContext[Any]], object]], Callable[[FactContext[Any]], object]]: A
            decorator that registers the wrapped function and returns the original function
            unchanged, so a module-level reference to it still calls the unguarded body.
    """

    def register(compute: Callable[[FactContext[Any]], object]) -> Callable[[FactContext[Any]], object]:
        REGISTRY[compute.__name__] = Fact(
            name=compute.__name__,
            cost=cost,
            compute=build_requires_guard(compute, requires),
            unguarded_compute=compute,
            deps=deps + tuple(name for name in requires if name not in deps),
            documentation=compute.__doc__ or "",
            requires=requires,
        )
        return compute

    return register


def build_requires_guard(
    compute: Callable[[FactContext[Any]], object], requires: tuple[str, ...]
) -> Callable[[FactContext[Any]], object]:
    """Wrap ``Fact.compute`` so it short-circuits to False when any required atom is False.

    Args:
        compute (Callable[[FactContext[Any]], object]): The fact's raw computation.
        requires (tuple[str, ...]): Atoms that must all be true for ``compute`` to run.

    Returns:
        Callable[[FactContext[Any]], object]: ``compute`` itself when ``requires`` is empty,
            else a guarded wrapper around it.
    """
    if not requires:
        return compute

    # TBD: wrapping value-fact causes problem here (e.g. _nearby_introns, _cds_footprint)
    @wraps(compute)
    def guarded(fact_context: FactContext[Any]) -> object:
        for required_name in requires:
            if not get_fact(fact_context, required_name):
                return False
        return compute(fact_context)

    return guarded


def get_fact(fact_context: FactContext, name: str) -> object:
    """Fetch a fact's value for this context, computing it at most once.

    The hottest function in the system, so it does the least it can when tracing is off: one
    ``is not None`` check, shared by both paths, whose body is a single dict bump. Everything else
    a derivation trace shows -- origins, costs, forced-False facts -- is reconstructed afterwards
    from the memo table (``ensembl.vep.trace.reconstruct``), because it can be. The read count is
    the one thing that cannot: a replay sees only memo hits, so it never observes the original
    read pattern. Making the counter unconditional to drop the branch measures 1.30x here.

    Args:
        fact_context (FactContext): The per-overlap memoisation context to read from and populate.
        name (str): The name of a registered fact.

    Returns:
        object: The fact's value, from the memo table if already computed.
    """
    recorder = fact_context.recorder
    if recorder is not None:
        read_counts = recorder.read_counts
        read_counts[name] = read_counts.get(name, 0) + 1

    memo = fact_context.memo
    if name in memo:
        return memo[name]

    registered_fact = REGISTRY[name]
    value = registered_fact.compute(fact_context)
    memo[name] = value
    return value


def atom_names() -> list[str]:
    """Return the Boolean facts a rule may reference.

    Returns:
        list[str]: Sorted names of every registered fact that does not start with ``_``.
    """
    return sorted(name for name in REGISTRY if not name.startswith("_"))


def implied_false(known: Mapping[str, bool]) -> dict[str, bool]:
    """Extend known facts that are implied False by the current assignment.

    The result also explains itself later: a name this added is False here and absent from the
    memo, which is how ``trace.reconstruct.forced_false_events`` recovers the forced-False set
    and its cause without this function recording anything, and without the explanation
    duplicating the propagation.

    Args:
        known (Mapping[str, bool]): A partial assignment of fact names to boolean values.

    Returns:
        dict[str, bool]: ``known`` extended with every fact that ``requires`` proves False given
            the input, iterated to a fixed point.
    """
    extended = dict(known)
    changed = True
    while changed:
        changed = False
        for name, registered_fact in REGISTRY.items():
            if name in extended or not registered_fact.requires:
                continue
            if any(extended.get(required) is False for required in registered_fact.requires):
                extended[name] = False
                changed = True
    return extended


def verify_requires(fact_contexts: Iterable[FactContext[Any]]) -> list[str]:
    """Check every ``requires`` claim against real contexts.

    A wrong declaration would silently produce wrong consequences. This is a sanity check to make sure
    that never happens.

    Checks a fact by running a fact's compute function without its requires guard, and then checking
    when the result is True all of its required atoms are True. if one required atom is False, the claim
    is false and the fact's compute function is mis-declared.

    Running the **unguarded** body is the whole point: a guarded call can never disprove its own
    guard, since the guard returns False whenever a required atom is False. Checking through
    ``get_fact`` would therefore pass for every declaration, correct or not.

    A body that *raises* when its required atoms are false is reported too. That is the same class
    of defect seen from the other side: the body is relying on the guard having already
    short-circuited, rather than being independently consistent with what it declares.

    Args:
        fact_contexts (Iterable[FactContext]): Contexts to check the claims against. Each is
            reset before use.

    Returns:
        list[str]: A description of every ``requires`` claim found to be false, or an empty list
            if every claim held.
    """
    problems: list[str] = []
    for fact_context in fact_contexts:
        fact_context.reset()
        for name, registered_fact in REGISTRY.items():
            if not registered_fact.requires:
                continue
            required_all_true = all(get_fact(fact_context, required) for required in registered_fact.requires)
            if required_all_true:
                continue
            try:
                value = registered_fact.unguarded_compute(fact_context)
            except Exception as exc:  # pylint: disable=broad-except
                # A body that cannot even run unguarded meaning the fact fails to compute when any require is
                # False is itself a finding. It may mean the requirement declaration is not just about truth,
                # and the body silently depends on the guard having already short-circuited.
                problems.append(
                    f"{name} raised {type(exc).__name__}: {exc} for overlap "
                    f"{fact_context.overlap!r} when requires={registered_fact.requires} was not met"
                )
                continue
            if value:
                problems.append(
                    f"{name} is true for overlap {fact_context.overlap!r} despite "
                    f"requires={registered_fact.requires}"
                )
    return problems


_closure_cost_cache: dict[str, int] = {}


def closure_cost(name: str) -> int:
    """Return the cost of a fact and its dependencies.

    Args:
        name (str): The name of a registered fact.

    Returns:
        int: The total cost of computing ``name`` and every fact it (transitively) depends on,
            counting each shared dependency once.
    """
    if name in _closure_cost_cache:
        return _closure_cost_cache[name]
    seen: set[str] = set()

    def walk(fact_name: str) -> int:
        if fact_name in seen:
            return 0
        seen.add(fact_name)
        registered_fact = REGISTRY[fact_name]
        return registered_fact.cost + sum(walk(dep) for dep in registered_fact.deps)

    _closure_cost_cache[name] = walk(name)
    return _closure_cost_cache[name]


def total_cost(names: Iterable[str]) -> int:
    """Return the cost of facts and their shared dependencies.

    Args:
        names (Iterable[str]): Names of registered facts.

    Returns:
        int: The total cost, with every shared dependency (transitively) counted exactly once.
    """
    seen: set[str] = set()
    total = 0
    stack = list(names)
    while stack:
        name = stack.pop()
        if name in seen:
            continue
        seen.add(name)
        registered_fact = REGISTRY[name]
        total += registered_fact.cost
        stack.extend(registered_fact.deps)
    return total


# --------------------------------------------------------------------------
# Optional measured-cost override:
#
# ``vep catalogue costs`` profiles every fact's real self-time over a corpus and can persist the
# result to ``annotators/consequence/data/fact_costs.json`` -- see
# ensembl.vep.annotators.consequence.facts.costs.
#
# This section only ever overrides an already-registered fact's ``cost`` field, and only when the
# ``fact_cost_override`` configuration option asks for it (see annotators/factory.py). Could be
# useful if in certain platform the measured costs are significantly different from the declared
# ones.
# --------------------------------------------------------------------------

FACT_COST_OVERRIDE_RESOURCE = "data/fact_costs.json"


def apply_fact_cost_overrides(overrides: Mapping[str, int]) -> None:
    """Override REGISTRY entries' ``cost`` in place from a measured fact-cost table.

    Declared cost tiers are untouched for any fact not present in ``overrides``. Also closure_cost's
    memoisation cache is cleared so subsequent calls reflect the new costs rather than a stale value
    computed before the override was applied.

    Args:
        overrides (Mapping[str, int]): Fact name -> cost unit, typically loaded from
            ``annotators/consequence/data/fact_costs.json``.
    """
    for fact_name, cost in overrides.items():
        registered_fact = REGISTRY.get(fact_name)
        if registered_fact is not None:
            REGISTRY[fact_name] = evolve(registered_fact, cost=int(cost))
    _closure_cost_cache.clear()


def load_fact_cost_overrides_from_package() -> dict[str, int]:
    """Read the checked-in measured fact-cost table, if any, without applying it.

    Returns:
        dict[str, int]: The ``{fact_name: cost_units}`` table from ``annotators/data/
            fact_costs.json``, or an empty dict if that resource does not exist.
    """
    resource = files("ensembl.vep.annotators.consequence").joinpath(FACT_COST_OVERRIDE_RESOURCE)
    if not resource.is_file():
        return {}
    return dict(json.loads(resource.read_text(encoding="utf-8")))


def apply_packaged_fact_cost_overrides() -> None:
    """Apply the packaged measured-cost table, if the package ships one.

    Call after every fact submodule has registered, since REGISTRY is empty before that. Costs
    change evaluation order only, never the terms emitted.
    """
    overrides = load_fact_cost_overrides_from_package()
    if overrides:
        apply_fact_cost_overrides(overrides)
