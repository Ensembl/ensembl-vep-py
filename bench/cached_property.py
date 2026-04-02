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
"""Benchmarks for class descriptors.

This is used to evaluate different strategies for implementing cached properties. The benchmarks include:
- `WithSetAttr`: A class that uses `__setattr__` to prevent direct assignment to the `functools.cached_property` attributes.
- `WithProperty`: A class that uses a regular `@property` without a setter to prevent direct assignment but no caching.
- `WithCachedProperty`: A class that uses `functools.cached_property` to cache but no prohibition on direct assignment.
- `WithManualCache`: A class that uses a custom desciptor from `functools.cached_property` but uses __set__ to prevent direct assignment.
- `WithCachedPropertyMangling`: A class that uses `functools.cached_property` with name mangling to make the cached value less accessible.

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

As expected, WithProperty is very slow on get due to no caching. WithCachedProperty and WithCachedPropertyMangling
are much faster due to caching, there is a small overhead in WithCachedPropertyMangling due to calling the property
method to access the cached value, but it is still very fast. WithSetAttr is also fast on get, but object
initialization is slow because of __setattr__ being called on object initialization. WithManualCache is slower than
WithCachedProperty due to the overhead of the __set__ method in the custom descriptor, it make the attribute a
data descriptor instead of a non-data descriptor, which means the descriptor's __get__ method will be called on every
access, and not directly from the instance's __dict__ after the first access.

WithCachedPropertyMangling seems to be the best option in terms of performance and safety, it is almost as
fast as WithCachedProperty but with the added benefit of making the cached value less accessible and preventing
direct assignment.
"""

import timeit
import functools
import statistics


class MyCachedProperty(functools.cached_property):
    def __set__(self, obj, value):
        raise AttributeError(
            f"'{self.attrname}' is a cached property — assign to its dependencies, not to it directly"
        )


class WithSetAttr:
    def __init__(self):
        self.x = 1000

    def __setattr__(self, name, value):
        if name in ("value", "tvalue"):
            raise AttributeError(
                f"'{name}' is a cached property — assign to its dependencies, not to it directly"
            )
        super().__setattr__(name, value)

    @functools.cached_property
    def value(self):
        return sum(range(self.x))

    @functools.cached_property
    def tvalue(self):
        return self.x * 2  # trivial compute


class WithProperty:
    def __init__(self):
        self.x = 1000

    @property
    def value(self):
        return sum(range(self.x))

    @property
    def tvalue(self):
        return self.x * 2  # trivial compute


class WithCachedProperty:
    def __init__(self):
        self.x = 1000

    @functools.cached_property
    def value(self):
        return sum(range(self.x))

    @functools.cached_property
    def tvalue(self):
        return self.x * 2  # trivial compute


class WithManualCache:
    def __init__(self):
        self.x = 1000

    @MyCachedProperty
    def value(self):
        return sum(range(self.x))

    @MyCachedProperty
    def tvalue(self):
        return self.x * 2  # trivial compute


class WithCachedPropertyMangling:
    def __init__(self):
        self.x = 1000

    @property
    def value(self):
        return self.__value

    @functools.cached_property
    def __value(self):
        return sum(range(self.x))

    @property
    def tvalue(self):
        return self.__tvalue

    @functools.cached_property
    def __tvalue(self):
        return self.x * 2  # trivial compute


def bench(cls, n=10_000_000):
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


bench(WithSetAttr)
bench(WithProperty)
bench(WithCachedProperty)
bench(WithManualCache)
bench(WithCachedPropertyMangling)
