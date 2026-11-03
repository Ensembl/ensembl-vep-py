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
"""Benchmark cached-property invalidation strategies.

Each pair compares a cached-property implementation with and without an attrs ``on_setattr``
hook that clears the cached value when ``x`` changes.

Result:
    WithSetAttr           : init                    212.0 ns/call  443.3 ns/call
    WithSetAttr           : trivial get              10.4 ns/call   10.4 ns/call
    WithSetAttr           : invalidation            148.8 ns/call  335.4 ns/call
    WithProperty          : init                     81.7 ns/call  157.3 ns/call
    WithProperty          : trivial get              34.2 ns/call   34.0 ns/call
    WithProperty          : invalidation              8.4 ns/call  601.2 ns/call
    WithCachedProperty    : init                     85.8 ns/call  161.3 ns/call
    WithCachedProperty    : trivial get               8.2 ns/call    8.1 ns/call
    WithCachedProperty    : invalidation              8.5 ns/call  426.8 ns/call
    WithManualCache       : init                     87.7 ns/call  161.8 ns/call
    WithManualCache       : trivial get               8.2 ns/call    8.6 ns/call
    WithManualCache       : invalidation              8.8 ns/call  427.6 ns/call
    WithCachedPropertyMangling: init                     84.7 ns/call  159.1 ns/call
    WithCachedPropertyMangling: trivial get              48.9 ns/call   49.5 ns/call
    WithCachedPropertyMangling: invalidation              8.4 ns/call  526.8 ns/call

The first number is measured without invalidation and the second with invalidation. Invalidation
increases assignment and construction cost, so it is best reserved for mutable core attributes
whose derived values must remain coherent.
"""

from __future__ import annotations

import functools
import statistics
import timeit
from collections.abc import Callable
from typing import Any

from attrs import define, field


class MyCachedProperty(functools.cached_property):
    """Cached property that refuses direct assignment."""

    def __set__(self, obj: object, value: object) -> None:
        """Refuse assignment to the derived value."""
        raise AttributeError(
            f"'{self.attrname}' is a cached property — assign to its dependencies, not to it directly"
        )


def invalidates(*cached_names: str) -> Callable[[Any, Any, Any], Any]:
    """Return an attrs hook that clears named cached properties."""

    def hook(instance: Any, _attribute: Any, new_value: Any) -> Any:
        try:
            for name in cached_names:
                delattr(instance, name)
        except AttributeError as error:
            raise NotImplementedError("Feature is not implemented yet.") from error

        return new_value  # attrs requires the hook to return the value

    return hook


@define(kw_only=True, auto_attribs=True)
class WithSetAttr:
    """Protect a cached value without invalidating it.

    Attributes:
        x (int): Input used by the derived value.
    """

    x: int = field(default=1000)

    def __setattr__(self, name: str, value: object) -> None:
        """Refuse assignment to cached values."""
        if name in ("value", "tvalue"):
            raise AttributeError(
                f"'{name}' is a cached property — assign to its dependencies, not to it directly"
            )
        super().__setattr__(name, value)

    @functools.cached_property
    def tvalue(self) -> int:
        """Return the trivial derived value."""
        return self.x * 2  # trivial compute


@define(kw_only=True, auto_attribs=True)
class WithSetAttrInvalidates:
    """Protect and invalidate a cached value through ``__setattr__``.

    Attributes:
        x (int): Input used by the derived value.
    """

    x: int = field(default=1000)

    def __setattr__(self, name: str, value: object) -> None:
        """Refuse assignment to cached values and invalidate them on change."""
        if name in ("value", "tvalue"):
            raise AttributeError(
                f"'{name}' is a cached property — assign to its dependencies, not to it directly"
            )
        if name == "x":
            try:
                delattr(self, "tvalue")
            except AttributeError as error:
                raise NotImplementedError("Feature is not implemented yet.") from error
        super().__setattr__(name, value)

    @functools.cached_property
    def tvalue(self) -> int:
        """Return the trivial derived value."""
        return self.x * 2  # trivial compute


@define(kw_only=True, auto_attribs=True)
class WithProperty:
    """Expose uncached read-only properties.

    Attributes:
        x (int): Input used by the derived values.
    """

    x: int = field(default=1000)

    @property
    def value(self) -> int:
        """Return the expensive derived value."""
        return sum(range(self.x))

    @property
    def tvalue(self) -> int:
        """Return the trivial derived value."""
        return self.x * 2  # trivial compute


@define(kw_only=True, auto_attribs=True)
class WithPropertyInvalidates:
    """Expose read-only properties with an unused invalidation hook.

    Attributes:
        x (int): Input used by the derived values.
    """

    x: int = field(default=1000, on_setattr=invalidates("tvalue"))

    @property
    def value(self) -> int:
        """Return the expensive derived value."""
        return sum(range(self.x))

    @property
    def tvalue(self) -> int:
        """Return the trivial derived value."""
        return self.x * 2  # trivial compute


@define(kw_only=True, auto_attribs=True)
class WithCachedProperty:
    """Expose cached properties without invalidation.

    Attributes:
        x (int): Input used by the derived values.
    """

    x: int = field(default=1000)

    @functools.cached_property
    def value(self) -> int:
        """Return the expensive derived value."""
        return sum(range(self.x))

    @functools.cached_property
    def tvalue(self) -> int:
        """Return the trivial derived value."""
        return self.x * 2  # trivial compute


@define(kw_only=True, auto_attribs=True)
class WithCachedPropertyInvalidates:
    """Invalidate a standard cached property when ``x`` changes.

    Attributes:
        x (int): Input used by the derived values.
    """

    x: int = field(default=1000, on_setattr=invalidates("tvalue"))

    @functools.cached_property
    def value(self) -> int:
        """Return the expensive derived value."""
        return sum(range(self.x))

    @functools.cached_property
    def tvalue(self) -> int:
        """Return the trivial derived value."""
        return self.x * 2  # trivial compute


@define(kw_only=True, auto_attribs=True)
class WithManualCache:
    """Expose protected cached values without invalidation.

    Attributes:
        x (int): Input used by the derived values.
    """

    x: int = field(default=1000)

    @MyCachedProperty
    def value(self) -> int:
        """Return the expensive derived value."""
        return sum(range(self.x))

    @MyCachedProperty
    def tvalue(self) -> int:
        """Return the trivial derived value."""
        return self.x * 2  # trivial compute


@define(kw_only=True, auto_attribs=True)
class WithManualCacheInvalidates:
    """Invalidate a protected cached property when ``x`` changes.

    Attributes:
        x (int): Input used by the derived values.
    """

    x: int = field(default=1000, on_setattr=invalidates("tvalue"))

    @MyCachedProperty
    def value(self) -> int:
        """Return the expensive derived value."""
        return sum(range(self.x))

    @MyCachedProperty
    def tvalue(self) -> int:
        """Return the trivial derived value."""
        return self.x * 2  # trivial compute


@define(kw_only=True, auto_attribs=True)
class WithCachedPropertyMangling:
    """Hide standard cached properties behind read-only properties.

    Attributes:
        x (int): Input used by the derived values.
    """

    x: int = field(default=1000)

    @property
    def value(self) -> int:
        """Return the expensive derived value."""
        return self.__value

    @functools.cached_property
    def __value(self) -> int:
        return sum(range(self.x))

    @property
    def tvalue(self) -> int:
        """Return the trivial derived value."""
        return self.__tvalue

    @functools.cached_property
    def __tvalue(self) -> int:
        return self.x * 2  # trivial compute


@define(kw_only=True, auto_attribs=True)
class WithCachedPropertyManglingInvalidates:
    """Invalidate a hidden cached property when ``x`` changes.

    Attributes:
        x (int): Input used by the derived value.
    """

    x: int = field(default=1000, on_setattr=invalidates("__tvalue"))

    @property
    def tvalue(self) -> int:
        """Return the trivial derived value."""
        return self.__tvalue

    @functools.cached_property
    def __tvalue(self) -> int:
        return self.x * 2  # trivial compute


def bench(cls1: type[Any], cls2: type[Any], n: int = 1_000_000) -> None:
    """Compare construction, reads, and assignment for two classes."""
    obj1 = cls1()  # warm up
    obj2 = cls2()  # warm up
    times1 = timeit.repeat("cls1()", globals={"cls1": cls1}, number=n, repeat=5)
    times2 = timeit.repeat("cls2()", globals={"cls2": cls2}, number=n, repeat=5)
    mean_ns1 = statistics.mean(times1) / n * 1e9
    mean_ns2 = statistics.mean(times2) / n * 1e9
    print(f"{cls1.__name__:22}:", f"{'init':22}", f"{mean_ns1:6.1f} ns/call", f"{mean_ns2:6.1f} ns/call")

    _ = obj1.tvalue  # warm up
    _ = obj2.tvalue  # warm up
    times1 = timeit.repeat("obj1.tvalue", globals={"obj1": obj1}, number=n, repeat=5)
    times2 = timeit.repeat("obj2.tvalue", globals={"obj2": obj2}, number=n, repeat=5)
    mean_ns1 = statistics.mean(times1) / n * 1e9
    mean_ns2 = statistics.mean(times2) / n * 1e9
    print(
        f"{cls1.__name__:22}:", f"{'trivial get':22}", f"{mean_ns1:6.1f} ns/call", f"{mean_ns2:6.1f} ns/call"
    )

    obj1.x = 2000  # warm up invalidation
    obj2.x = 2000  # warm up invalidation
    times1 = timeit.repeat("obj1.x = 2000", globals={"obj1": obj1}, number=n, repeat=5)
    times2 = timeit.repeat("obj2.x = 2000", globals={"obj2": obj2}, number=n, repeat=5)
    mean_ns1 = statistics.mean(times1) / n * 1e9
    mean_ns2 = statistics.mean(times2) / n * 1e9
    print(
        f"{cls1.__name__:22}:", f"{'invalidation':22}", f"{mean_ns1:6.1f} ns/call", f"{mean_ns2:6.1f} ns/call"
    )


def main() -> None:
    """Run every invalidation benchmark pair."""
    bench(WithSetAttr, WithSetAttrInvalidates)
    bench(WithProperty, WithPropertyInvalidates)
    bench(WithCachedProperty, WithCachedPropertyInvalidates)
    bench(WithManualCache, WithManualCacheInvalidates)
    bench(WithCachedPropertyMangling, WithCachedPropertyManglingInvalidates)


if __name__ == "__main__":
    main()
