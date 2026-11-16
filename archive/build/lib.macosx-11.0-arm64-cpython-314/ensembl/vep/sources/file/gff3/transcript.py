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
"""GFF3-backed transcript feature source."""

from __future__ import annotations

from typing import ClassVar

from attrs import define, field

from ensembl.vep.parsers import GFF3Parser, GFF3Record
from ensembl.vep.sources.file.gxf.transcript import GXFTranscriptSource


@define(kw_only=True, auto_attribs=True)
class GFF3TranscriptSource(GXFTranscriptSource[GFF3Parser]):
    """Transcript source backed by a GFF3 file.

    Attributes:
        adaptor: The parser adaptor used to stream GFF3 records.
    """

    transcript_feature_types: ClassVar[frozenset[str]] = frozenset(
        {"mrna", "transcript", "lnc_rna", "pseudogenic_transcript", "processed_transcript"}
    )
    adaptor: GFF3Parser = field(init=False, factory=GFF3Parser)

    def _query_typed_records(
        self,
        region: str,
        *,
        feature_types: frozenset[str] | None = None,
    ) -> list[GFF3Record]:
        """Seek the persistent GFF3 query parser and read its records."""

        query_parser = self._open_query_parser()
        if not query_parser.seek(region):
            return []
        return query_parser.records(limit=None, feature_types=feature_types)

    def get_record_id(self, record: GFF3Record) -> str | None:
        """Return the GFF3 record identifier.

        Args:
            record: The parser record to inspect.

        Returns:
            str | None: The record ID, or the `Name` fallback when the record
            has no `ID` attribute.
        """

        _id = self._first_attr(record, "ID") or self._first_attr(record, "Name")
        if _id is None:
            _id = f"{record.feature_type.lower()}:{record.seqid}:{record.start}-{record.end}"
        return _id

    def get_parent_id(self, record: GFF3Record) -> str | None:
        """Return the GFF3 parent identifier, if present.

        Args:
            record: The parser record to inspect.

        Returns:
            str | None: The parent identifier, or None when the record has no
            parent.
        """

        return self._first_attr(record, "Parent")
