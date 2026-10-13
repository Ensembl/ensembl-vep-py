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
"""Tests for GFF3 corpus export."""

from __future__ import annotations

from pathlib import Path

import pytest

from ensembl.vep.common import Strand
from ensembl.vep.corpus import Corpus
from ensembl.vep.corpus.exporter.formats.file.gff3 import GFF3_NAME
from ensembl.vep.corpus.exporter.manifest import CorpusManifest
from ensembl.vep.parsers.gff3 import GFF3Parser

READER_FEATURE_TYPES = {"gene", "transcript", "exon", "cds", "five_prime_utr", "three_prime_utr"}


@pytest.fixture(name="records", scope="module")
def fixture_records(saved: tuple[Path, CorpusManifest]) -> list:
    """Parse the exported GFF3 once through the production parser."""
    directory, _manifest = saved
    return GFF3Parser.parse_file(directory / GFF3_NAME, limit=None)


class TestGff3:
    """Transcript annotations reconstructed from GFF3."""

    def test_uses_only_reader_supported_feature_types(self, records: list) -> None:
        """Every exported row is recognised by the transcript reader."""
        assert {record.feature_type.lower() for record in records} <= READER_FEATURE_TYPES

    def test_transcript_rows_preserve_coordinates_and_strand(self, source: Corpus, records: list) -> None:
        """Transcript rows retain their identifiers, spans and orientation."""
        rows = {
            record.attributes["ID"][0]: record for record in records if record.feature_type == "transcript"
        }
        for transcript in source.transcripts:
            row = rows[f"transcript:{transcript.stable_id}"]
            assert (row.start, row.end) == (transcript.start, transcript.end)
            assert row.strand_symbol == ("+" if transcript.strand == Strand.FORWARD else "-")

    def test_exon_rows_preserve_spliced_exons(self, source: Corpus, records: list) -> None:
        """Exon spans remain grouped under their parent transcript."""
        spans_by_parent: dict[str, set[tuple[int, int]]] = {}
        for record in records:
            if record.feature_type == "exon":
                spans_by_parent.setdefault(record.attributes["Parent"][0], set()).add(
                    (record.start, record.end)
                )
        for transcript in source.transcripts:
            expected = {(exon.start, exon.end) for exon in transcript.spliced_exons}
            assert spans_by_parent[f"transcript:{transcript.stable_id}"] == expected

    def test_cds_rows_preserve_spans_and_initial_phase(self, source: Corpus, records: list) -> None:
        """Coding blocks are exact and begin in phase zero in reading order."""
        rows_by_parent: dict[str, list] = {}
        for record in records:
            if record.feature_type == "CDS":
                rows_by_parent.setdefault(record.attributes["Parent"][0], []).append(record)
        for transcript in source.transcripts:
            rows = rows_by_parent.get(f"transcript:{transcript.stable_id}", [])
            spans = transcript.cds_genomic_spans or []
            assert {(row.start, row.end) for row in rows} == set(spans)
            if spans:
                assert next(row for row in rows if row.start == spans[0][0]).phase_value == 0

    def test_gene_parentage_and_biotype_are_not_invented(self, records: list) -> None:
        """Parentless transcripts stay parentless, while genes retain their biotype."""
        transcript_rows = {
            record.attributes["ID"][0]: record for record in records if record.feature_type == "transcript"
        }
        genes = {record.attributes["ID"][0]: record for record in records if record.feature_type == "gene"}
        assert "Parent" not in transcript_rows["transcript:ENST0001"].attributes
        assert transcript_rows["transcript:ENST0002"].attributes["Parent"] == ("gene:ENSG0002",)
        assert genes["gene:ENSG0004"].attributes["biotype"] == ("nonsense_mediated_decay",)
