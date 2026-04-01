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

This is used to evaluate different strategies for invalidating cached properties when their dependencies change. 
The benchmarks include:
- {calss}: A class that does not invalidate the cached property when its dependencies change.
- {class}Invalidates: A class that uses `invalidates` on_setattr hook to clear the cached property when 
    its dependencies change.

The classes here are different class used in `cached_property.py`.

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

The first result is for the classes without invalidation, and the second result is for the classes with invalidation. 
In all cases, the invalidation adds a significant overhead to the get operation, 50-70x increase. Furthermore, it also
increased object initialization time by 2x.

Therefore, it makes invalidation undesirable. Invalidation should be used only with core attributes and updating core 
attributes should be discouraged.

The object initialization is also a problem. WithCachedPropertyMangling is still the better option than WithSetAttr as
both is increasing 2x.
"""
import timeit
import functools
import statistics
from attr import define, field

class MyCachedProperty(functools.cached_property):
    def __set__(self, obj, value):
        raise AttributeError(
            f"'{self.attrname}' is a cached property — assign to its dependencies, not to it directly"
        )    

def invalidates(*cached_names):
    """
    Returns an attrs on_setattr hook that clears the named
    cached_property entries whenever this field is set.
    """
    def hook(instance, attrib, new_value):
        try:
            for name in cached_names:
                delattr(instance, name)
        except AttributeError:
            pass

        # cannot be used with __slots__
        # d = instance.__dict__
        # for name in cached_names:
        #     d.pop(name, None)  # Fast, no exception if not present

        return new_value           # attrs requires the hook to return the value
    return hook

@define
class WithSetAttr:
    x : int = field(default=1000)

    def __setattr__(self, name, value):
        if name in ("value", "tvalue"):
            raise AttributeError(
                f"'{name}' is a cached property — assign to its dependencies, not to it directly"
            )
        super().__setattr__(name, value)
    
    @functools.cached_property
    def tvalue(self):
        return self.x * 2       # trivial compute
    

@define
class WithSetAttrInvalidates:
    x : int = field(default=1000)

    def __setattr__(self, name, value):
        if name in ("value", "tvalue"):
            raise AttributeError(
                f"'{name}' is a cached property — assign to its dependencies, not to it directly"
            )
        if name == "x":
            try:
                delattr(self, "tvalue")
            except AttributeError:
                pass
        super().__setattr__(name, value)
    
    @functools.cached_property
    def tvalue(self):
        return self.x * 2       # trivial compute

@define
class WithProperty:
    x : int = field(default=1000)

    @property
    def value(self):
        return sum(range(self.x))
    
    @property
    def tvalue(self):
        return self.x * 2       # trivial compute
    
@define
class WithPropertyInvalidates():
    x : int = field(default=1000, on_setattr=invalidates("tvalue"))

    @property
    def value(self):
        return sum(range(self.x))
    
    @property
    def tvalue(self):
        return self.x * 2       # trivial compute

@define
class WithCachedProperty:
    x : int = field(default=1000)

    @functools.cached_property
    def value(self):
        return sum(range(self.x))
    
    @functools.cached_property
    def tvalue(self):
        return self.x * 2       # trivial compute

@define
class WithCachedPropertyInvalidates:
    x : int = field(default=1000, on_setattr=invalidates("tvalue"))

    @functools.cached_property
    def value(self):
        return sum(range(self.x))
    
    @functools.cached_property
    def tvalue(self):
        return self.x * 2       # trivial compute

@define
class WithManualCache:
    x : int = field(default=1000)

    @MyCachedProperty
    def value(self):
        return sum(range(self.x))
    
    @MyCachedProperty
    def tvalue(self):
        return self.x * 2       # trivial compute
    
@define
class WithManualCacheInvalidates:
    x : int = field(default=1000, on_setattr=invalidates("tvalue"))

    @MyCachedProperty
    def value(self):
        return sum(range(self.x))
    
    @MyCachedProperty
    def tvalue(self):
        return self.x * 2       # trivial compute

@define
class WithCachedPropertyMangling:
    x : int = field(default=1000)

    @property
    def value(self):
        self.__value

    @functools.cached_property
    def __value(self):
        return sum(range(self.x))
    
    @property
    def tvalue(self):
        self.__tvalue

    @functools.cached_property
    def __tvalue(self):
        return self.x * 2       # trivial compute
    
@define
class WithCachedPropertyManglingInvalidates:
    x : int = field(default=1000, on_setattr=invalidates("__tvalue"))
    
    @property
    def tvalue(self):
        self.__tvalue

    @functools.cached_property
    def __tvalue(self):
        return self.x * 2       # trivial compute

def bench(cls1, cls2, n=1_000_000):
    obj1 = cls1()                        # warm up
    obj2 = cls2()                        # warm up
    times1 = timeit.repeat("cls1()", globals={"cls1": cls1}, number=n, repeat=5)
    times2 = timeit.repeat("cls2()", globals={"cls2": cls2}, number=n, repeat=5)
    mean_ns1 = statistics.mean(times1) / n * 1e9
    mean_ns2 = statistics.mean(times2) / n * 1e9
    print(f"{cls1.__name__:22}:", f"{'init':22}", f"{mean_ns1:6.1f} ns/call", f"{mean_ns2:6.1f} ns/call")

    _ = obj1.tvalue                  # warm up
    _ = obj2.tvalue                  # warm up
    times1 = timeit.repeat("obj1.tvalue", globals={"obj1": obj1}, number=n, repeat=5)
    times2 = timeit.repeat("obj2.tvalue", globals={"obj2": obj2}, number=n, repeat=5)
    mean_ns1 = statistics.mean(times1) / n * 1e9
    mean_ns2 = statistics.mean(times2) / n * 1e9
    print(f"{cls1.__name__:22}:", f"{'trivial get':22}", f"{mean_ns1:6.1f} ns/call", f"{mean_ns2:6.1f} ns/call")

    obj1.x = 2000                    # warm up invalidation
    obj2.x = 2000                    # warm up invalidation
    times1 = timeit.repeat("obj1.x = 2000", globals={"obj1": obj1}, number=n, repeat=5)
    times2 = timeit.repeat("obj2.x = 2000", globals={"obj2": obj2}, number=n, repeat=5)
    mean_ns1 = statistics.mean(times1) / n * 1e9
    mean_ns2 = statistics.mean(times2) / n * 1e9
    print(f"{cls1.__name__:22}:", f"{'invalidation':22}", f"{mean_ns1:6.1f} ns/call", f"{mean_ns2:6.1f} ns/call")

bench(WithSetAttr, WithSetAttrInvalidates)
bench(WithProperty, WithPropertyInvalidates)
bench(WithCachedProperty, WithCachedPropertyInvalidates)
bench(WithManualCache, WithManualCacheInvalidates)
bench(WithCachedPropertyMangling, WithCachedPropertyManglingInvalidates)