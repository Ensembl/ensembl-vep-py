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
"""Registry for corpus exporters."""

from __future__ import annotations

from collections.abc import Callable
from typing import TypeVar

from ensembl.vep.corpus.exporter.base import CorpusExporter

__all__ = ["available_exporters", "get_exporter", "register_exporter"]

ExporterFactory = Callable[[], CorpusExporter]
"""Construct a corpus exporter."""

ExporterFactoryT = TypeVar("ExporterFactoryT", bound=ExporterFactory)
"""A particular exporter factory for the registration decorator."""

_EXPORTER_REGISTRY: dict[str, ExporterFactory] = {}


def register_exporter(name: str) -> Callable[[ExporterFactoryT], ExporterFactoryT]:
    """Register an exporter factory under a name."""

    def decorator(factory: ExporterFactoryT) -> ExporterFactoryT:
        if name in _EXPORTER_REGISTRY:
            raise ValueError(f"Corpus exporter {name!r} is already registered")
        _EXPORTER_REGISTRY[name] = factory
        return factory

    return decorator


def get_exporter(name: str) -> CorpusExporter:
    """Return a new registered exporter.

    Raises:
        KeyError: If no exporter is registered under ``name``.
    """
    try:
        return _EXPORTER_REGISTRY[name]()
    except KeyError as error:
        available = ", ".join(sorted(_EXPORTER_REGISTRY)) or "<none>"
        raise KeyError(f"Unknown corpus exporter '{name}'. Available exporters: {available}") from error


def available_exporters() -> tuple[str, ...]:
    """Return the sorted set of registered exporter names."""
    return tuple(sorted(_EXPORTER_REGISTRY))
