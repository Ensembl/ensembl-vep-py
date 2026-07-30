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
"""CircularSlice class."""

from attrs import define

from ensembl.vep.models.sequence.slice import Slice


@define(kw_only=True, auto_attribs=True)

__all__ = ["CircularSlice"]


class CircularSlice(Slice):
    """
    Represents a circular genomic region.
    """

    def is_circular(self) -> bool:
        """Indicates whether the slice is circular.

        Returns:
            bool: True, indicating the slice is circular.
        """
        return True
