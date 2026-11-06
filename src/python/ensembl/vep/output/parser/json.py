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
"""JSON output. A declared seam; the writer it mirrors does not exist yet."""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

from ensembl.vep.models.variation import Variation
from ensembl.vep.output.parser.base import OutputParser

__all__ = ["JSONOutputParser"]


class JSONOutputParser(OutputParser):
    """Parser for VEP JSON output."""

    def parse(self, path: str | Path) -> Iterator[Variation]:
        """Not implemented yet.

        Raises:
            NotImplementedError: Always; the JSON format has no reader yet.
        """
        raise NotImplementedError(f"reading JSON output ({path}) is not implemented yet")
