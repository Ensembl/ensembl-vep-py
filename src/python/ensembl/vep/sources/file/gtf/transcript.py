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
"""GTF-backed transcript feature source."""

from __future__ import annotations

from typing import ClassVar

from attrs import define, field

from ensembl.vep.parsers import GTFParser, GTFRecord
from ensembl.vep.sources.file.gxf.transcript import GXFTranscriptSource


@define(kw_only=True, auto_attribs=True)
class GTFTranscriptSource(GXFTranscriptSource[GTFParser]):
    """Transcript source backed by a GTF file.

    Attributes:
        adaptor: The parser adaptor used to stream GTF records.
        _PARENTS: Mapping of feature types to their parent feature types.
    """

    transcript_feature_types: ClassVar[frozenset[str]] = frozenset({"transcript"})
    coding_feature_types: ClassVar[frozenset[str]] = frozenset({"cds", "stop_codon"})
    adaptor: GTFParser = field(init=False, factory=GTFParser)
    _PARENTS: ClassVar[dict[str, str]] = {
        "exon": "transcript",
        "cds": "transcript",
        "stop_codon": "transcript",
        "transcript": "gene",
    }

    def _query_typed_records(
        self,
        region: str,
        *,
        feature_types: frozenset[str] | None = None,
    ) -> list[GTFRecord]:
        """Seek the persistent GTF query parser and read its records."""

        query_parser = self._open_query_parser()
        if not query_parser.seek(region):
            return []
        return query_parser.records(limit=None, feature_types=feature_types)

    def get_record_id(self, record: GTFRecord) -> str:
        """Return the GTF record identifier for a record.

        GTF stores identifiers in attributes named after the feature type, such
        as `gene_id`, `transcript_id`, or `exon_id`.

        Args:
            record: The record to inspect.

        Returns:
            str: The feature-specific identifier.
        """

        _id = self._first_attr(record, f"{record.feature_type.lower()}_id")
        if _id is None:
            _id = f"{record.feature_type.lower()}:{record.seqid}:{record.start}-{record.end}"
        return _id

    def get_parent_ids(self, record: GTFRecord) -> tuple[str, ...]:
        """Return the GTF parent identifier for a record, if one exists.

        Args:
            record: The record to inspect.

        Returns:
            tuple[str, ...]: A tuple of parent identifiers, or an empty tuple when the record has no
            parent feature.
        """

        parent_feature_type = self._PARENTS.get(record.feature_type.lower())
        if parent_feature_type is None:
            return ()
        return record.attributes.get(f"{parent_feature_type}_id", ())
