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
"""Write annotated variants as JSON lines."""

from __future__ import annotations

import json
from collections.abc import Iterable, Iterator, Sequence

from ensembl.vep.output.factory.base import AnnotatedVariant
from ensembl.vep.output.factory.stream.json import JSONStreamFactory

__all__ = ["JSONOutputWriter"]


class JSONOutputWriter:
    """Write one Ensembl VEP-schema JSON object per line."""

    def __init__(self, stream_factory: JSONStreamFactory | None = None) -> None:
        """Initialize the writer.

        Args:
            stream_factory (JSONStreamFactory | None): The record mapper, or None for the default.
        """
        self.stream_factory = stream_factory or JSONStreamFactory()

    def headers(self, input_headers: Sequence[str] = ()) -> list[str]:
        """Return no headers for JSON lines."""
        del input_headers
        return []

    def line(self, variant: AnnotatedVariant) -> str:
        """Return one compact JSON object for an annotated variant."""
        record = self.stream_factory.record(variant)
        return json.dumps(record, separators=(",", ":"), ensure_ascii=False)

    def lines(self, variants: Iterable[AnnotatedVariant]) -> Iterator[str]:
        """Yield one compact JSON object per annotated variant."""
        for variant in variants:
            yield self.line(variant)
