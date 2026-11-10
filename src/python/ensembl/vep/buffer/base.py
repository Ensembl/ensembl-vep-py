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
"""The buffer contract.

A buffer holds one batch of input objects. It is handed the objects to draw from, rather than
holding a source of its own, so a buffer stays picklable and reading stays with the input source.
"""

from __future__ import annotations

from collections.abc import Iterator
from typing import Protocol, TypeVar, runtime_checkable

__all__ = ["Buffer"]

ItemT = TypeVar("ItemT")
"""What item type a buffer holds."""


@runtime_checkable
class Buffer(Protocol[ItemT]):
    """Refillable batch of input items."""

    @property
    def objects(self) -> tuple[ItemT, ...]:
        """Return the items in the current batch, in input order."""
        raise NotImplementedError

    @property
    def is_empty(self) -> bool:
        """Return whether the input the batch is drawn from is spent."""
        raise NotImplementedError

    def next(self, items: Iterator[ItemT]) -> None:
        """Reset the batch and refill it from ``items``."""
        raise NotImplementedError

    def clear(self) -> None:
        """Drop whatever was indexed for the current batch."""
        raise NotImplementedError
