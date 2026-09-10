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
"""Buffer registry keyed by ``<domain>.<strategy>`` names."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any, TypeVar, cast

from ensembl.vep.buffers.buffer import Buffer
from ensembl.vep.models import Feature
from ensembl.vep.sources import FeatureSource

__all__ = ["available_buffers", "get_buffer", "register_buffer"]

BufferFeatureT = TypeVar("BufferFeatureT", bound=Feature)
"""What feature type a registered buffer stores."""

BufferSourceT = TypeVar("BufferSourceT", bound=FeatureSource[Any])
"""What source type a registered buffer reads from."""

BufferFactory = Callable[..., Buffer[BufferFeatureT, BufferSourceT]]
"""Factory for buffers, parameterized by feature and source type."""

_BUFFER_REGISTRY: dict[str, BufferFactory[Any, Any]] = {}


def register_buffer(
    name: str, factory: BufferFactory[BufferFeatureT, BufferSourceT] | None = None
) -> (
    BufferFactory[BufferFeatureT, BufferSourceT]
    | Callable[[BufferFactory[BufferFeatureT, BufferSourceT]], BufferFactory[BufferFeatureT, BufferSourceT]]
):
    """Register a buffer factory under a name.

    Args:
        name: Registry name, e.g. ``"input.fixed-size"``.
        factory: Factory to register, or ``None`` when used as a decorator.

    Returns:
        The registered factory or a decorator that registers one.
    """

    def decorator(
        registered_factory: BufferFactory[BufferFeatureT, BufferSourceT],
    ) -> BufferFactory[BufferFeatureT, BufferSourceT]:
        _BUFFER_REGISTRY[name] = cast(BufferFactory[Any, Any], registered_factory)
        return registered_factory

    if factory is not None:
        return decorator(factory)
    return decorator


def buffer_factory(name: str) -> BufferFactory[BufferFeatureT, BufferSourceT]:
    """Return the registered factory without calling it.

    Args:
        name: Registry name, e.g. ``"input.fixed-size"``.

    Returns:
        The registered factory, usually a buffer class.

    Raises:
        KeyError: If no buffer is registered under ``name``.
    """
    try:
        return cast(BufferFactory[BufferFeatureT, BufferSourceT], _BUFFER_REGISTRY[name])
    except KeyError as error:
        available = ", ".join(sorted(_BUFFER_REGISTRY)) or "<none>"
        raise KeyError(f"Unknown buffer '{name}'. Available buffers: {available}") from error


def get_buffer(name: str, /, *args: object, **kwargs: object) -> Buffer[BufferFeatureT, BufferSourceT]:
    """Return a registered buffer instance.

    Args:
        name: Registry name, e.g. ``"input.fixed-size"``.
        *args: Positional arguments forwarded to the registered factory.
        **kwargs: Keyword arguments forwarded to the registered factory.

    Returns:
        A new buffer instance.

    Raises:
        KeyError: If no buffer is registered under ``name``.
    """
    return buffer_factory(name)(*args, **kwargs)


def available_buffers() -> tuple[str, ...]:
    """Return the sorted set of registered buffer names.

    Returns:
        Every registered buffer name, sorted.
    """
    return tuple(sorted(_BUFFER_REGISTRY))
