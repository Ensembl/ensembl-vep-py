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
"""Runner registry keyed by ``<domain>.<verb>`` names."""

from __future__ import annotations

from collections.abc import Callable
from typing import TypeVar

from ensembl.vep.runner.base import Runner

__all__ = [
    "available_runners",
    "get_runner",
    "register_runner",
    "runner_default",
    "runner_group",
    "runner_groups",
]

RunnerFactory = type[Runner]
"""Runner classes."""

RunnerT = TypeVar("RunnerT", bound=RunnerFactory)
"""A particular runner class, so the decorator hands back what it was given."""

_RUNNER_REGISTRY: dict[str, RunnerFactory] = {}
_RUNNER_DEFAULTS: dict[str, str] = {}


def register_runner(name: str, *, default_for_group: bool = False) -> Callable[[RunnerT], RunnerT]:
    """Register a runner class under a name."""

    def decorator(factory: RunnerT) -> RunnerT:
        _RUNNER_REGISTRY[name] = factory
        if default_for_group:
            _RUNNER_DEFAULTS[runner_group(name)] = name
        return factory

    return decorator


def get_runner(name: str) -> RunnerFactory:
    """Return the registered runner class."""
    try:
        return _RUNNER_REGISTRY[name]
    except KeyError as error:
        available = ", ".join(sorted(_RUNNER_REGISTRY)) or "<none>"
        raise KeyError(f"Unknown command '{name}'. Available commands: {available}") from error


def runner_default(group: str) -> str | None:
    """Return a group's default runner, if it declares one."""
    return _RUNNER_DEFAULTS.get(group)


def available_runners() -> tuple[str, ...]:
    """Return the sorted set of registered runner names.

    Returns:
        tuple[str, ...]: Every registered runner name, sorted.
    """
    return tuple(sorted(_RUNNER_REGISTRY))


def runner_group(name: str) -> str:
    """Return the domain part of a runner name.

    Args:
        name (str): A runner name, e.g. "bench.matrix".

    Returns:
        str: The part before the first dot, e.g. "bench". A name with no dot is its own group.
    """
    return name.split(".", 1)[0]


def runner_groups() -> dict[str, tuple[str, ...]]:
    """Return the registered runners grouped by domain.

    Returns:
        dict[str, tuple[str, ...]]: Domain -> its runner names, sorted. A domain with exactly one
            runner named after it is a subcommand in its own right; one with several is a group
            whose verbs are the subcommands.
    """
    groups: dict[str, list[str]] = {}
    for name in available_runners():
        groups.setdefault(runner_group(name), []).append(name)
    return {group: tuple(sorted(names)) for group, names in sorted(groups.items())}
