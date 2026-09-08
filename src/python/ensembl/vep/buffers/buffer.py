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
from ensembl.vep.models import Feature
from ensembl.vep.sources import FeatureSource

F = TypeVar("F", bound=Feature)
S = TypeVar("S", bound=FeatureSource)


@define(kw_only=True, auto_attribs=True)
class Config:
    """Stand in object for Config class"""

    buffer_size: int = field(default=5, init=False)


@define(kw_only=True, auto_attribs=True)
class Buffer(ABC, Generic[F, S]):
    """Buffer interface.

    Attributes:
        source: The source of the ``Features``
        config: A Config object for configuring the buffer
        buffer: A queue of ``Feature`` objects obtained from the source
    """

    source: S
    config: Config = field(factory=Config)
    _buffer: list[F] = field(factory=list, init=False)
    _pre_buffer: deque[F] = field(factory=deque, init=False)
    _is_empty: bool = field(default=False, init=False)

    @property
    def buffer(self) -> list[F]:
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

    @property
    def feature_iterator(self) -> Iterator[F]:
        """Public accessor for _feature_iterator.

        Returns:
            An iterator over the source features.
        """
        return self._feature_iterator

    @cached_property
    def _feature_iterator(self) -> Iterator[F]:
        """Returns the iterator over the source features.

        Returns:
            An iterator over the source features.
        """
        return self.source.iter_all_features()

    @abstractmethod
    def next(self) -> None:
        """Reset and re-fill the buffer from the source.

        Returns:
            None
        """

        raise NotImplementedError
