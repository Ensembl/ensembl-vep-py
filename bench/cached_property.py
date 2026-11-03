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
"""Benchmark cached-property implementations.

The benchmark compares a read-only property, the standard cached-property descriptor, a data
descriptor that refuses assignment, explicit ``__setattr__`` protection, and a name-mangled
cached property.

Result:
    WithSetAttr           : init                    216.7 ns/call
    WithSetAttr           : get                      21.2 ns/call
    WithSetAttr           : trivial get              21.2 ns/call
    WithProperty          : init                     91.1 ns/call
    WithProperty          : get                    9390.3 ns/call
    WithProperty          : trivial get              35.1 ns/call
    WithCachedProperty    : init                     97.6 ns/call
    WithCachedProperty    : get                      21.2 ns/call
    WithCachedProperty    : trivial get              21.5 ns/call
    WithManualCache       : init                     94.8 ns/call
    WithManualCache       : get                     108.4 ns/call
    WithManualCache       : trivial get             106.5 ns/call
    WithCachedPropertyMangling: init                     94.6 ns/call
    WithCachedPropertyMangling: get                      37.0 ns/call
    WithCachedPropertyMangling: trivial get              39.6 ns/call

The uncached property makes the expensive read much slower. Name mangling adds a small property
call to a standard cached-property read, while ``__setattr__`` protection increases construction
cost. The custom data descriptor is slower because its ``__get__`` method remains on every access.
"""

from __future__ import annotations

import functools
import statistics
import timeit
from typing import Any


class MyCachedProperty(functools.cached_property):
    """Cached property that refuses direct assignment."""

    def __set__(self, obj: object, value: object) -> None:
        """Refuse assignment to the derived value."""
        raise AttributeError(
            f"'{self.attrname}' is a cached property — assign to its dependencies, not to it directly"
        )


class WithSetAttr:
    """Protect cached values through ``__setattr__``.

    Attributes:
        x (int): Input used by the derived values.
    """

    def __init__(self) -> None:
        self.x = 1000

    def __setattr__(self, name: str, value: object) -> None:
        """Refuse assignment to either cached value."""
        if name in ("value", "tvalue"):
            raise AttributeError(
                f"'{name}' is a cached property — assign to its dependencies, not to it directly"
            )
        super().__setattr__(name, value)

    @functools.cached_property
    def value(self) -> int:
        """Return the expensive derived value."""
        return sum(range(self.x))

    @functools.cached_property
    def tvalue(self) -> int:
        """Return the trivial derived value."""
        return self.x * 2  # trivial compute


class WithProperty:
    """Expose uncached read-only properties.

    Attributes:
        x (int): Input used by the derived values.
    """

    def __init__(self) -> None:
        self.x = 1000

    @property
    def value(self) -> int:
        """Return the expensive derived value."""
        return sum(range(self.x))

    @property
    def tvalue(self) -> int:
        """Return the trivial derived value."""
        return self.x * 2  # trivial compute


class WithCachedProperty:
    """Expose assignable standard cached properties.

    Attributes:
        x (int): Input used by the derived values.
    """

    def __init__(self) -> None:
        self.x = 1000

    @functools.cached_property
    def value(self) -> int:
        """Return the expensive derived value."""
        return sum(range(self.x))

    @functools.cached_property
    def tvalue(self) -> int:
        """Return the trivial derived value."""
        return self.x * 2  # trivial compute


class WithManualCache:
    """Expose cached values through a data descriptor.

    Attributes:
        x (int): Input used by the derived values.
    """

    def __init__(self) -> None:
        self.x = 1000

    @MyCachedProperty
    def value(self) -> int:
        """Return the expensive derived value."""
        return sum(range(self.x))

    @MyCachedProperty
    def tvalue(self) -> int:
        """Return the trivial derived value."""
        return self.x * 2  # trivial compute


class WithCachedPropertyMangling:
    """Hide standard cached properties behind read-only properties.

    Attributes:
        x (int): Input used by the derived values.
    """

    def __init__(self) -> None:
        self.x = 1000

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


def bench(cls: type[Any], n: int = 10_000_000) -> None:
    """Benchmark construction and derived-value reads for one class."""
    obj = cls()  # warm up
    times = timeit.repeat("cls()", globals={"cls": cls}, number=n, repeat=5)
    mean_ns = statistics.mean(times) / n * 1e9
    print(f"{cls.__name__:22}:", f"{'init':22}", f"{mean_ns:6.1f} ns/call")

    _ = obj.value  # warm up
    times = timeit.repeat("obj.value", globals={"obj": obj}, number=n, repeat=5)
    mean_ns = statistics.mean(times) / n * 1e9
    print(f"{cls.__name__:22}:", f"{'get':22}", f"{mean_ns:6.1f} ns/call")

    _ = obj.tvalue  # warm up
    times = timeit.repeat("obj.tvalue", globals={"obj": obj}, number=n, repeat=5)
    mean_ns = statistics.mean(times) / n * 1e9
    print(f"{cls.__name__:22}:", f"{'trivial get':22}", f"{mean_ns:6.1f} ns/call")


def main() -> None:
    """Run every cached-property benchmark."""
    bench(WithSetAttr)
    bench(WithProperty)
    bench(WithCachedProperty)
    bench(WithManualCache)
    bench(WithCachedPropertyMangling)


if __name__ == "__main__":
    main()
