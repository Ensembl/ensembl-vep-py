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
"""Test class definitions for common methods' unit testing."""

from functools import cached_property
from attrs import define, field

from ensembl.vep.common import invalidates


@define(kw_only=True, auto_attribs=True, slots=False)
class Foo:
    """Class for testing cached property invalidation."""

    start: int = field(on_setattr=invalidates("_length"))
    end: int = field(on_setattr=invalidates("_length"))

    @cached_property
    def _length(self) -> int:
        """
        Returns the length.

        Returns:
            int: The length.
        """
        return self.end - self.start + 1

    @cached_property
    def _range(self) -> str:
        """
        Returns the range notation.

        Returns:
            str: The range notation.
        """
        return f"{self.start}..{self.end}"

    @property
    def length(self) -> int:
        """
        Public accessor for the _invalidated_length.

        Returns:
            int: The length.
        """
        return self._length

    @property
    def range(self) -> str:
        """
        Public accessor to get the range notation

        Returns:
            int: The length.
        """
        return self._range


class Bar(Foo):
    """Class for testing cached property invalidation on inheritance."""
