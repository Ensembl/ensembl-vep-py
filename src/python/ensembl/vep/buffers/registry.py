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
"""Buffer registry keyed by ``<object_type>.<data_provider>.<size>`` names."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any, TypeVar, cast, overload

from ensembl.vep.buffers.buffer import Buffer

__all__ = ["available_buffers", "get_buffer", "register_buffer"]

BufferItemT = TypeVar("BufferItemT")
"""What item type a registered buffer stores."""

BufferFactory = Callable[..., Buffer[BufferItemT]]
"""Factory for buffers, parameterized by item type."""

_BUFFER_REGISTRY: dict[str, BufferFactory[Any]] = {}


@overload
def register_buffer(
    name: str, factory: None = None
) -> Callable[[BufferFactory[BufferItemT]], BufferFactory[BufferItemT]]: ...


@overload
def register_buffer(name: str, factory: BufferFactory[BufferItemT]) -> BufferFactory[BufferItemT]: ...


def register_buffer(
    name: str, factory: BufferFactory[BufferItemT] | None = None
) -> BufferFactory[BufferItemT] | Callable[[BufferFactory[BufferItemT]], BufferFactory[BufferItemT]]:
    """Register a buffer factory under a name.

    Args:
        name: Registry name, e.g. ``"variation.source.fixed"``.
        factory: Factory to register, or ``None`` when used as a decorator.

    Returns:
        The registered factory or a decorator that registers one.
    """

    def decorator(
        registered_factory: BufferFactory[BufferItemT],
    ) -> BufferFactory[BufferItemT]:
        _BUFFER_REGISTRY[name] = cast(BufferFactory[Any], registered_factory)
        return registered_factory

    if factory is not None:
        return decorator(factory)
    return decorator


def buffer_factory(name: str) -> BufferFactory[BufferItemT]:
    """Return the registered factory without calling it.

    Args:
        name: Registry name, e.g. ``"variation.source.fixed"``.

    Returns:
        The registered factory, usually a buffer class.

    Raises:
        KeyError: If no buffer is registered under ``name``.
    """
    try:
        return cast(BufferFactory[BufferItemT], _BUFFER_REGISTRY[name])
    except KeyError as error:
        available = ", ".join(sorted(_BUFFER_REGISTRY)) or "<none>"
        raise KeyError(f"Unknown buffer '{name}'. Available buffers: {available}") from error


def get_buffer(name: str, /, *args: object, **kwargs: object) -> Buffer[BufferItemT]:
    """Return a registered buffer instance.

    Args:
        name: Registry name, e.g. ``"variation.source.fixed"``.
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
