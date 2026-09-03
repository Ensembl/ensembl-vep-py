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

from ensembl.vep.common.methods import invalidate


class UnsafeSetterMixin:
    """
    Can be used to temporarily disable cache invalidation to allow faster but unsafe attribute update.

    Caveat:
    This can be very expansive in terms of storage as it stores the type(self) for every instance.
    It also create a temporary subclass under the unsafe context although that gets cleaned up
    after the context exits.

    It also allows setting of cached properties directly using seed_attrs and seed_and_invalidate_attrs
    which can lead to inconsistent state if not used carefully.
    """

    __slots__ = ()

    @contextmanager
    def unsafe(self):  # type: ignore
        """Context manager that temporarily swaps the class of this instance to a version with an
        unsafe __setattr__ that bypasses attrs' on_setattr hooks."""
        safe_cls = type(self)

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

    def seed_attrs(self, **seeds: object) -> None:
        """(Pre-)populate derived and cached attributes, by name, without triggering their normal
        derivation.

        Writes directly into each named protected slot (``object.__setattr__``, bypassing ``on_setattr``
        hooks and any setter), exactly as a normal ``cached_property`` computation would populate it on
        first access.

        Seeded attribute can be cleared by the ``invalidates()`` hooks that clear a computed value.
        Therefore, Seeding a property derived from such a field must happen after that field is set, or
        the seed is invalidated before anything ever reads it.

        A derived attribute's computation may depend on other derived attributes. Seeding a derived attribute
        does not invalidate derived attributes that depend on it. Use seed_and_invalidate_attrs() if you want
        to seed and invalidate dependent derived attributes.

        Args:
            object: The object to seed.
            **seeds (object): Cached-property name (without the leading underscore) with value to
                seed, e.g. ``transcript.seed_attrs(cds=my_cds)``.
        """
        for name, value in seeds.items():
            object.__setattr__(self, f"_{name}", value)

    def seed_and_invalidate_attrs(self, **seeds: object) -> None:
        """(Pre-)populate derived and cached attributes, by name, without triggering their
        normal derivation and invalidate dependent derived attributes. Don't invalidate dependent
        derived attributes if they have also been seeded in the call to ``seed_and_invalidate_attrs``.

        This depends on the class having a SEED_INVALIDATES mapping that defines which derived
        attributes should be invalidated when a given attribute is seeded.

        Args:
            object: The object to seed.
            **seeds (object): Cached-property name (without the leading underscore) with value to seed,
                e.g. ``transcript.seed_and_invalidate_attrs(translation=my_translation, cds=my_cds)``.
        """
        seed_invalidations = getattr(type(self), "SEED_INVALIDATES", {})
        seeded_slots: set[str] = set()
        invalidation_candidates: set[str] = set()

        for name, value in seeds.items():
            object.__setattr__(self, f"_{name}", value)
            seeded_slots.add(f"_{name}")
            invalidation_candidates.update(seed_invalidations.get(name, ()))

        for invalidation_candidate in invalidation_candidates - seeded_slots:
            invalidate(self, invalidation_candidate)


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
            except AttributeError:
                pass  # If the attribute doesn't exist, we can ignore it

        return super().__setattr__(name, value)
