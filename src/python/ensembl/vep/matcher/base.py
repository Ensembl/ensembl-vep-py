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
"""Matcher contracts and policies."""

from __future__ import annotations

from enum import Enum
from typing import Protocol

__all__ = ["AlleleLike", "InsertionMatchPolicy", "IntervalLike"]


class IntervalLike(Protocol):
    """Inclusive coordinates used for interval matching."""

    @property
    def start(self) -> int:
        """Return the inclusive start coordinate."""
        raise NotImplementedError

    @property
    def end(self) -> int:
        """Return the inclusive end coordinate."""
        raise NotImplementedError


class AlleleLike(IntervalLike, Protocol):
    """Coordinates and allele strings used for colocation."""

    @property
    def reference_allele(self) -> str | None:
        """Return the reference allele string."""
        raise NotImplementedError

    @property
    def alternate_allele(self) -> str | None:
        """Return the alternate allele string."""
        raise NotImplementedError


class InsertionMatchPolicy(str, Enum):
    """Policy for associating an interbase insertion with a feature."""

    CONTAINED = "contained"
    TOUCHING = "touching"
    THREE_PRIME = "three_prime"
