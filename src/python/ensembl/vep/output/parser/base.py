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
"""The reader contract."""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Iterator
from pathlib import Path

from ensembl.vep.models.variation import Variation

__all__ = ["OutputParser"]


class OutputParser(ABC):
    """Parser for annotated files containing variation consequences.

    A reader needs no transcript or sequence source: the consequences are seeded straight onto
    each allele, leaving the feature overlaps empty, and everything a summary reads still answers.
    """

    @abstractmethod
    def parse(self, path: str | Path) -> Iterator[Variation]:
        """Yield one variation per record, with its alleles' consequences attached."""
        raise NotImplementedError
