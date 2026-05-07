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

# pylint: disable=redefined-outer-name
"""Unit tests for common mixins."""

from functools import cached_property
import pytest
from attrs import define, field


from ensembl.vep.common import invalidates
from ensembl.vep.common import UnsafeSetterMixin, CachedAttributeGuardMixin

@define
class DummyUnsafeSetter(UnsafeSetterMixin):
    """Dummy class for UnsafeSetterMixin tests."""
    x: int = field(on_setattr=invalidates("__y"))
    __y: int = field(init=False, default=99)

    @property
    def y(self) -> int:
        """Return the value of y."""
        return self.__y


@define
class DummyCachedAttributeGuard(CachedAttributeGuardMixin):
    """Dummy class for CachedAttributeGuardMixin tests."""
    a: int = field()
    _DERIVED = {"b"}
    _INVALIDATES = {"a": ("b",)}

    @cached_property
    def b(self) -> int:
        """Return b as a*2."""
        return self.a * 2


class TestUnsafeSetterMixin:
    """Tests for UnsafeSetterMixin."""

    def test_unsafe_setter_invalidation(self) -> None:
        """Setting x should invalidate y (delete y)."""
        d = DummyUnsafeSetter(x=1)
        d.x = 42
        assert d.x == 42
        assert not hasattr(d, "y")

    def test_unsafe_setter_no_invalidation(self) -> None:
        """In unsafe context, y should NOT be invalidated."""
        d = DummyUnsafeSetter(x=1)
        with d.unsafe():
            d.x = 100
        assert d.x == 100
        assert d.y == 99

class TestCachedAttributeGuardMixin:
    """Tests for CachedAttributeGuardMixin."""

    def test_invalidate_and_guard(self) -> None:
        """Test derived property invalidation and guard."""
        d = DummyCachedAttributeGuard(a=1)
        assert d.a == 1
        assert d.b == 2  # b is derived from a

        # Test that setting a derived attribute raises an error
        with pytest.raises(AttributeError):
            d.b = 10

        # Test invalidation
        d.a = 5
        assert d.b == 10    # b gets invalidated and recalculated
