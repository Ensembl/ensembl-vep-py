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
"""Common mixins. Mixins are for reusable behavior without data."""

from contextlib import contextmanager


class UnsafeSetterMixin:
    """
    Can be used to temporarily disable cache invalidation to allow faster but unsafe attribute update.

    Caveat:
    This can be very expansive in terms of storage as it stores the type(self) for every instance.
    It also create a temporary subclass under the unsafe context although that gets cleaned up
    after the context exits.
    """

    __slots__ = ("_safe_cls",)

    def __attrs_post_init__(self) -> None:
        """Set up the safe class reference for later use in the unsafe context manager."""
        object.__setattr__(self, "_safe_cls", type(self))

    @contextmanager
    def unsafe(self):  # type: ignore
        """Context manager that temporarily swaps the class of this instance to a version with an
        unsafe __setattr__ that bypasses attrs' on_setattr hooks."""
        safe_cls = self._safe_cls

        # create a throwaway subclass with a bare __setattr__
        unsafe_cls = type(
            f"{safe_cls.__name__}_unsafe",
            (safe_cls,),
            {
                "__slots__": (),
                "__setattr__": object.__setattr__,
            },
        )

        # point this instance at the unsafe subclass
        object.__setattr__(self, "__class__", unsafe_cls)
        try:
            yield self
        finally:
            # restore the original class
            object.__setattr__(self, "__class__", safe_cls)


class CachedAttributeGuardMixin:
    """Mixin that guards against setting attributes that are internal properties
    and/or invalidate cached properties when their dependencies change.

    Caveat:
    These can be very expansive as it will be called every time any attribute is updated.

    Currently not in use as a less preferable method. Instead, name mangling is used for guarding update to
    internal properties and attrs' on_setattr hook for invalidation.
    """

    __slots__ = ()

    def __setattr__(self, name, value):  # type: ignore
        """Override __setattr__ to guard against setting attributes that are internal properties
        and/or invalidate cached properties when their dependencies change.

        Args:
            name (str): The name of the attribute being set.
            value: The value to set the attribute to.

        Raises:
            AttributeError: If attempting to set a cached property directly.
        """
        if name not in type(self)._DERIVED and name not in type(self)._INVALIDATES:
            return super().__setattr__(name, value)

        if name in type(self)._DERIVED and name in self.__slots__:
            raise AttributeError(f"'{name}' is a cached property and cannot be reassigned")

        for invalidation_candidate in type(self)._INVALIDATES.get(name, ()):
            try:
                delattr(self, invalidation_candidate)
            except AttributeError as exc:
                raise NotImplementedError("Feature is not implemented yet.") from exc

        return super().__setattr__(name, value)
