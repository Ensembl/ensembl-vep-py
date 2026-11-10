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
"""Tests for the buffer registry."""

from __future__ import annotations

from collections.abc import Iterator

import pytest
from attrs import define

from ensembl.vep.buffer import registry as registry_module
from ensembl.vep.buffer.base import Buffer
from ensembl.vep.buffer.input_buffer import FixedSizeInputBuffer, InputBuffer
from ensembl.vep.buffer.registry import available_buffers, get_buffer, register_buffer
from ensembl.vep.models.variation import Variation


@define(kw_only=True, auto_attribs=True)
class _EagerBuffer(InputBuffer[Variation]):
    """Minimal buffer used to test registration and lookup."""

    def next(self, items: Iterator[Variation]) -> None:
        """Take every remaining object."""

        self._fill(list(items), exhausted=True)


class TestAvailableBuffers:
    """What the registry knows about."""

    def test_includes_the_built_ins(self) -> None:
        """Register the built-in buffers on import."""
        assert "variation.fixed" in available_buffers()


class TestLookup:
    """Looking a buffer up by name."""

    def test_unknown_name_raises_key_error_with_the_available_list(self) -> None:
        """Name the registered buffers when asked for one that is not there."""
        with pytest.raises(KeyError, match="Unknown buffer 'not.a.real.buffer'"):
            get_buffer("not.a.real.buffer")

    def test_builtin_fixed_size_input_buffer_can_be_resolved(self) -> None:
        """Build the built-in fixed-size input buffer, forwarding its options."""
        buffer: Buffer[Variation] = get_buffer("variation.fixed", buffer_size=10)

        assert isinstance(buffer, FixedSizeInputBuffer)
        assert buffer.buffer_size == 10


class TestRegistration:
    """Registering a buffer of one's own."""

    def test_register_buffer_as_a_plain_call(self) -> None:
        """Register a factory directly."""
        register_buffer("test.eager", _EagerBuffer)
        try:
            assert "test.eager" in available_buffers()
            assert isinstance(get_buffer("test.eager"), _EagerBuffer)
        finally:
            registry_module._BUFFER_REGISTRY.pop("test.eager", None)  # pylint: disable=protected-access

    def test_register_buffer_as_a_decorator(self) -> None:
        """Register a factory as a decorator."""

        @register_buffer("test.decorated")  # type: ignore[arg-type]
        @define(kw_only=True, auto_attribs=True)
        class _DecoratedBuffer(_EagerBuffer):
            """Minimal decorated buffer used to test registration."""

        try:
            assert "test.decorated" in available_buffers()
            assert isinstance(get_buffer("test.decorated"), _DecoratedBuffer)
        finally:
            registry_module._BUFFER_REGISTRY.pop("test.decorated", None)  # pylint: disable=protected-access
