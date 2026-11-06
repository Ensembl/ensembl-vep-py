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
"""The contract a writer of annotated output implements."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Protocol

from ensembl.vep.output.factory.base import AnnotatedVariant

__all__ = ["OutputWriter"]


class OutputWriter(Protocol):
    """Serializer for one annotated-variant output format."""

    def headers(self, input_headers: Sequence[str] = ()) -> list[str]:
        """Return the lines written before the annotated variants."""

    def line(self, variant: AnnotatedVariant) -> str:
        """Return one serialized annotated variant."""
