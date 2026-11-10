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
"""The contract for a format adapter over annotated variants."""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from typing import Protocol, TypeVar

from ensembl.vep.output.factory.base import AnnotatedVariant

__all__ = ["AnnotationStreamFactory"]

RecordT_co = TypeVar("RecordT_co", covariant=True)


class AnnotationStreamFactory(Protocol[RecordT_co]):
    """Record factory for one annotation output representation."""

    def records(self, variants: Iterable[AnnotatedVariant]) -> Iterator[RecordT_co]:
        """Yield one output record for each annotated variant."""
