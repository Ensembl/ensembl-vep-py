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
"""Annotator registry keyed by ``<domain>.<strategy>`` names."""

from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING, Any, TypeVar, cast

from ensembl.vep.annotators.base import Annotator

__all__ = ["available_annotators", "get_annotator", "register_annotator"]

# What a registered annotator produces.
if TYPE_CHECKING:
    AnnotatorResultT = TypeVar("AnnotatorResultT", default=Any)
else:
    AnnotatorResultT = TypeVar("AnnotatorResultT")

AnnotatorFactory = Callable[..., Annotator[AnnotatorResultT]]
"""Factory for annotators, parameterized by result type."""

_ANNOTATOR_REGISTRY: dict[str, AnnotatorFactory[Any]] = {}


def register_annotator(
    name: str, factory: AnnotatorFactory[AnnotatorResultT] | None = None
) -> (
    AnnotatorFactory[AnnotatorResultT]
    | Callable[[AnnotatorFactory[AnnotatorResultT]], AnnotatorFactory[AnnotatorResultT]]
):
    """Register an annotator factory under a name.

    Args:
        name (str): Registry name, e.g. ``"consequence.linear"``.
        factory (AnnotatorFactory | None): Factory to register, or ``None`` when used as a
            decorator.

    Returns:
        AnnotatorFactory | Callable[[AnnotatorFactory], AnnotatorFactory]: The registered factory
            or a decorator that registers one.
    """

    def decorator(
        registered_factory: AnnotatorFactory[AnnotatorResultT],
    ) -> AnnotatorFactory[AnnotatorResultT]:
        _ANNOTATOR_REGISTRY[name] = cast(AnnotatorFactory[Any], registered_factory)
        return registered_factory

    if factory is not None:
        return decorator(factory)
    return decorator


def annotator_factory(name: str) -> AnnotatorFactory[AnnotatorResultT]:
    """Return the registered factory without calling it.

    Args:
        name (str): Registry name, e.g. ``"consequence.linear"``.

    Returns:
        AnnotatorFactory: The registered factory, usually an annotator class.

    Raises:
        KeyError: If no annotator is registered under ``name``.
    """
    try:
        return cast(AnnotatorFactory[AnnotatorResultT], _ANNOTATOR_REGISTRY[name])
    except KeyError as error:
        available = ", ".join(sorted(_ANNOTATOR_REGISTRY)) or "<none>"
        raise KeyError(f"Unknown annotator '{name}'. Available annotators: {available}") from error


def get_annotator(name: str, /, *args: object, **kwargs: object) -> Annotator[AnnotatorResultT]:
    """Return a registered annotator instance.

    Args:
        name (str): Registry name, e.g. ``"consequence.linear"``.
        *args (object): Positional arguments forwarded to the registered factory.
        **kwargs (object): Keyword arguments forwarded to the registered factory.

    Returns:
        Annotator[AnnotatorResultT]: A new annotator instance.

    Raises:
        KeyError: If no annotator is registered under ``name``.
    """
    return annotator_factory(name)(*args, **kwargs)


def available_annotators() -> tuple[str, ...]:
    """Return the sorted set of registered annotator names."""
    return tuple(sorted(_ANNOTATOR_REGISTRY))
