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

"""Unit tests for common mixins."""

from functools import cached_property
from typing import ClassVar

import pytest
from attrs import define, field

from ensembl.vep.common import CachedAttributeGuardMixin, UnsafeSetterMixin, invalidates


@define(kw_only=True, auto_attribs=True)
class _UnsafeSetterProbe(UnsafeSetterMixin):
    """Probe unsafe assignment and dependent attribute invalidation."""

    source_value: int = field(on_setattr=invalidates("_seeded_value"))
    _seeded_value: int = field(init=False, default=99)

    SEED_INVALIDATES: ClassVar[dict[str, tuple[str, ...]]] = {"seeded_value": ("_dependent_value",)}

    @property
    def seeded_value(self) -> int:
        """Return the seedable value."""
        return self._seeded_value

    @property
    def dependent_value(self) -> int:
        """Return the value derived from seeded_value."""
        return self._dependent_value

    @cached_property
    def _dependent_value(self) -> int:
        """Calculate the value derived from seeded_value."""
        return self._seeded_value + 1000


@define(kw_only=True, auto_attribs=True)
class _CachedAttributeGuardProbe(CachedAttributeGuardMixin):
    """Probe cached attribute invalidation and assignment guards."""

    source_value: int
    _DERIVED = {"derived_value"}
    _INVALIDATES = {"source_value": ("derived_value",)}

    @cached_property
    def derived_value(self) -> int:
        """Calculate a value from source_value."""
        return self.source_value * 2


class TestUnsafeSetterMixin:
    """Tests for the unsafe setter mixin."""

    def test_unsafe_setter_invalidation(self) -> None:
        """Test normal assignment invalidates the dependent value."""
        probe = _UnsafeSetterProbe(source_value=1)

        probe.source_value = 42

        assert probe.source_value == 42
        assert not hasattr(probe, "seeded_value")

    def test_unsafe_setter_no_invalidation(self) -> None:
        """Test unsafe assignment does not invalidate the dependent value."""
        probe = _UnsafeSetterProbe(source_value=1)

        with probe.unsafe():
            probe.source_value = 100

        assert probe.source_value == 100
        assert probe.seeded_value == 99

    def test_seed_attrs(self) -> None:
        """Test seed_attrs supplies a derived attribute value."""
        probe = _UnsafeSetterProbe(source_value=1)

        probe.seed_attrs(seeded_value=123)
        assert probe.seeded_value == 123

        with probe.unsafe():
            probe.source_value = 50
        assert probe.seeded_value == 123

        probe.source_value = 50
        assert not hasattr(probe, "seeded_value")

    def test_seed_attrs_with_dependencies(self) -> None:
        """Test seed_attrs preserves already resolved dependent attributes."""
        probe = _UnsafeSetterProbe(source_value=1)
        assert probe.seeded_value == 99
        assert probe.dependent_value == 1099

        probe.seed_attrs(seeded_value=123)

        assert probe.seeded_value == 123
        assert probe.dependent_value == 1099

    def test_seed_and_invalidate_attrs(self) -> None:
        """Test seed_and_invalidate_attrs recalculates dependent attributes."""
        probe = _UnsafeSetterProbe(source_value=1)
        assert probe.seeded_value == 99
        assert probe.dependent_value == 1099

        probe.seed_and_invalidate_attrs(seeded_value=123)

        assert probe.seeded_value == 123
        assert probe.dependent_value == 1123

    def test_seed_and_invalidate_attrs_invalidation(self) -> None:
        """Test seed_and_invalidate_attrs seed dependent attribute."""
        probe = _UnsafeSetterProbe(source_value=1)
        assert probe.seeded_value == 99

        probe.seed_and_invalidate_attrs(dependent_value=10)

        assert probe.seeded_value == 99
        assert probe.dependent_value == 10

    def test_seed_and_invalidate_attrs_avoid_invalidation(self) -> None:
        """Test seed_and_invalidate_attrs doesn't invalidate dependent attribute if already seeded."""
        probe = _UnsafeSetterProbe(source_value=1)
        assert probe.seeded_value == 99

        probe.seed_and_invalidate_attrs(dependent_value=10, seeded_value=123)

        assert probe.seeded_value == 123
        assert probe.dependent_value == 10


class TestCachedAttributeGuardMixin:
    """Tests for the cached attribute guard mixin."""

    def test_invalidate_and_guard(self) -> None:
        """Test derived attributes reject assignment and recalculate after invalidation."""
        probe = _CachedAttributeGuardProbe(source_value=1)
        assert probe.source_value == 1
        assert probe.derived_value == 2

        with pytest.raises(AttributeError):
            probe.derived_value = 10

        probe.source_value = 5

        assert probe.derived_value == 10
