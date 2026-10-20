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
"""Corpus factories and their registry."""

from __future__ import annotations

from collections.abc import Callable

from ensembl.vep.corpus.base import Corpus
from ensembl.vep.corpus.spec import CorpusSpec

__all__ = [
    "CorpusFactory",
    "CorpusOptionsError",
    "available_corpora",
    "build_corpus",
    "build_corpus_from_options",
    "corpus_factory",
    "corpus_spec_from_options",
    "register_corpus",
]


class CorpusOptionsError(ValueError):
    """A configured corpus name, parameter, or external input cannot be used."""


CorpusFactory = Callable[[CorpusSpec], Corpus]
"""Build a corpus from its specification."""

_CORPUS_FACTORIES: dict[str, CorpusFactory] = {}
_CORPUS_ALIASES: dict[str, str] = {}


def register_corpus(name: str, *, aliases: tuple[str, ...] = ()) -> Callable[[CorpusFactory], CorpusFactory]:
    """Register a corpus factory and any compatibility aliases.

    Args:
        name (str): Canonical registry name.
        aliases (tuple[str, ...]): Additional names resolving to the same factory.

    Returns:
        Callable[[CorpusFactory], CorpusFactory]: Factory decorator.

    Raises:
        ValueError: If a name is blank, repeated, or already registered.
    """
    requested_names = (name, *aliases)
    if any(not requested_name.strip() for requested_name in requested_names):
        raise ValueError("Corpus registry names must not be empty")

    def decorator(factory: CorpusFactory) -> CorpusFactory:
        if len(set(requested_names)) != len(requested_names):
            raise ValueError(f"Corpus registration repeats a name in {requested_names!r}")
        conflicts = [
            requested_name
            for requested_name in requested_names
            if requested_name in _CORPUS_FACTORIES or requested_name in _CORPUS_ALIASES
        ]
        if conflicts:
            raise ValueError(f"Corpus name {conflicts[0]!r} is already registered")

        _CORPUS_FACTORIES[name] = factory
        for alias in aliases:
            _CORPUS_ALIASES[alias] = name
        return factory

    return decorator


def corpus_factory(name: str) -> CorpusFactory:
    """Return the factory registered under ``name`` or one of its aliases.

    Args:
        name (str): Canonical name or alias.

    Returns:
        CorpusFactory: Registered corpus constructor.

    Raises:
        KeyError: If ``name`` is not registered.
    """
    canonical_name = _CORPUS_ALIASES.get(name, name)
    try:
        return _CORPUS_FACTORIES[canonical_name]
    except KeyError as error:
        available = ", ".join(available_corpora()) or "<none>"
        raise KeyError(f"Unknown corpus {name!r}. Available corpora: {available}") from error


def build_corpus(specification: CorpusSpec) -> Corpus:
    """Build the corpus described by ``specification``.

    Args:
        specification (CorpusSpec): Registered name and factory parameters.

    Returns:
        Corpus: Constructed corpus implementation.

    Raises:
        KeyError: If the requested corpus kind is not registered.
    """
    return corpus_factory(specification.kind)(specification)


def corpus_spec_from_options(options: object) -> CorpusSpec:
    """Build a specification from command options.

    Args:
        options (object): Configuration-like object carrying corpus attributes.

    Returns:
        CorpusSpec: Normalised construction specification.

    Raises:
        CorpusOptionsError: If an option cannot form a valid specification.
    """
    try:
        return CorpusSpec.from_options(options)
    except (TypeError, ValueError) as error:
        raise CorpusOptionsError(str(error)) from error


def build_corpus_from_options(options: object) -> Corpus:
    """Build a configured corpus with one user-facing error type.

    Args:
        options (object): Configuration-like object carrying corpus attributes.

    Returns:
        Corpus: Constructed corpus implementation.

    Raises:
        CorpusOptionsError: If configuration, lookup, or external input cannot be used.
    """
    try:
        return build_corpus(corpus_spec_from_options(options))
    except (KeyError, TypeError, ValueError, OSError) as error:
        if isinstance(error, CorpusOptionsError):
            raise
        message = error.args[0] if isinstance(error, KeyError) and error.args else str(error)
        raise CorpusOptionsError(str(message)) from error


def available_corpora() -> tuple[str, ...]:
    """Return registered corpus names in sorted order."""
    return tuple(sorted(_CORPUS_FACTORIES))
