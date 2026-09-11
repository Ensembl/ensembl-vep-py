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
"""Buffer interface."""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections import deque
from functools import cached_property
from typing import Generic, Iterator, TypeVar
from attrs import define, field

# from ensembl.vep.config.base import Config
from ensembl.vep.models.variation.variation import Variation
from ensembl.vep.sources import FeatureSource

T = TypeVar("T")
V = TypeVar("V", bound=Variation)


@define(kw_only=True, auto_attribs=True)
class Config:
    """Stand in object for Config class"""

    buffer_size: int = field(default=5, init=False)


@define(kw_only=True, auto_attribs=True)
class Buffer(ABC, Generic[T]):
    """Buffer interface.

    Attributes:
        config: A Config object for configuring the buffer
        buffer: A list of objects in the buffer
    """

    config: Config = field(factory=Config)
    _buffer: list[T] = field(factory=list, init=False)
    _pre_buffer: deque[T] = field(factory=deque, init=False)
    _is_empty: bool = field(default=False, init=False)

    @property
    def buffer(self) -> list[T]:
        """Public accessor for _buffer.

        Returns:
            A list of Feature objects
        """
        return self._buffer

    @property
    def is_empty(self) -> bool:
        """Public accessor for _is_empty.

        Returns:
            True if the buffer is empty and the source is exhausted, False otherwise.
        """
        return self._is_empty

    @abstractmethod
    def next(self) -> None:
        """Reset and re-fill the buffer from the source.

        Returns:
            None
        """

        raise NotImplementedError


@define(kw_only=True, auto_attribs=True)
class InputBuffer(Buffer[V], Generic[V]):
    """InputBuffer interface.

    Attributes:
        config: A Config object for configuring the buffer
        buffer: A queue of ``Variation`` objects obtained from the source
        source: The source of the ``Variation`` objects.
    """

    source: FeatureSource[V]
    _buffer: list[V] = field(factory=list, init=False)
    _pre_buffer: deque[V] = field(factory=deque, init=False)

    @property
    def variation_iterator(self) -> Iterator[V]:
        """Public accessor for _variation_iterator.

        Returns:
            An iterator over the source variations.
        """
        return self._variation_iterator

    @cached_property
    def _variation_iterator(self) -> Iterator[V]:
        """Returns the iterator over the source variations.

        Returns:
            An iterator over the source variations.
        """
        return self.source.iter_all_features()

    @abstractmethod
    def next(self) -> None:
        """Reset and re-fill the buffer from the source.

        Returns:
            None
        """

        raise NotImplementedError
