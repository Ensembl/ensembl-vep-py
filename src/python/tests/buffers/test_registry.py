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
"""Unit tests for ensembl.vep.buffers.registry."""

from __future__ import annotations

from attrs import define
import pytest

from ensembl.vep.buffers import registry as registry_module
from ensembl.vep.buffers.buffer import Buffer, InputBuffer
from ensembl.vep.buffers.fixed_size_input_buffer import FixedSizeInputBuffer
from ensembl.vep.buffers.registry import available_buffers, get_buffer, register_buffer
from ensembl.vep.common import SourceType
from ensembl.vep.models.variation.variation import Variation
from tests.buffers.conftest import DummyVariationSource


@define(kw_only=True, auto_attribs=True)
class _ConcreteBuffer(InputBuffer[Variation]):
    """Minimal buffer used to test registration and lookup."""

    def next(self) -> None:
        """Fill the buffer from the source."""
        self._buffer = list(self.source.iter_all_features())
        self._is_empty = True


class TestAvailableBuffers:
    """The registry exposes the built-in buffers."""

    def test_includes_the_built_ins(self) -> None:
        """Built-ins are registered on import."""
        assert "variation.source.fixed" in available_buffers()


class TestLookup:
    """Lookup behavior for registered and unregistered buffers."""

    def test_unknown_name_raises_key_error_with_the_available_list(self) -> None:
        """Unknown names raise KeyError with the registry summary in the message."""
        with pytest.raises(KeyError, match="Unknown buffer 'not.a.real.buffer'"):
            get_buffer("not.a.real.buffer")

    def test_builtin_fixed_size_input_buffer_can_be_resolved(self, variation_list: list[Variation]) -> None:
        """The registry can build the built-in fixed-size input buffer."""
        source = DummyVariationSource(source_type=SourceType.FILE, features=variation_list)

        buffer: Buffer[Variation] = get_buffer("variation.source.fixed", source=source)

        assert isinstance(buffer, FixedSizeInputBuffer)
        assert buffer.source is source


class TestRegistration:
    """Registration works as a plain call or as a decorator."""

    def test_register_buffer_as_a_plain_call(self, variation_list: list[Variation]) -> None:
        """A factory can be registered directly."""
        register_buffer("test.concrete", _ConcreteBuffer)
        try:
            source = DummyVariationSource(source_type=SourceType.FILE, features=variation_list)

            assert "test.concrete" in available_buffers()
            assert isinstance(get_buffer("test.concrete", source=source), _ConcreteBuffer)
        finally:
            registry_module._BUFFER_REGISTRY.pop("test.concrete", None)  # pylint: disable=protected-access

    def test_register_buffer_as_a_decorator(self, variation_list: list[Variation]) -> None:
        """A factory can also be registered as a decorator."""

        @register_buffer("test.decorated")  # type: ignore[arg-type]
        @define(kw_only=True, auto_attribs=True)
        class _DecoratedBuffer(_ConcreteBuffer):
            """Minimal decorated buffer used to test registration."""

        try:
            source = DummyVariationSource(source_type=SourceType.FILE, features=variation_list)

            assert "test.decorated" in available_buffers()
            assert isinstance(get_buffer("test.decorated", source=source), _DecoratedBuffer)
        finally:
            registry_module._BUFFER_REGISTRY.pop("test.decorated", None)  # pylint: disable=protected-access
