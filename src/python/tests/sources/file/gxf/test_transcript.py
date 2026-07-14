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
"""Unit tests for the shared GXF transcript source."""

from attrs import define, field

from ensembl.vep.parsers import GFF3Parser
from ensembl.vep.sources.file.gxf.transcript import GXFTranscriptSource, _GXFTranscriptGroup

_GFF3_ONE_TRANSCRIPT = """\
chr1\tensembl\tgene\t1000\t5000\t.\t+\t.\tID=gene:ENSG001;Name=BRCA1
chr1\tensembl\tmRNA\t1000\t5000\t.\t+\t.\tID=transcript:ENST001;Parent=gene:ENSG001;biotype=protein_coding
chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tID=exon:ENSE001;Parent=transcript:ENST001
chr1\tensembl\tCDS\t1100\t2000\t.\t+\t0\tParent=transcript:ENST001
chr1\tensembl\texon\t3000\t5000\t.\t+\t.\tID=exon:ENSE002;Parent=transcript:ENST001
chr1\tensembl\tCDS\t3000\t4900\t.\t+\t0\tParent=transcript:ENST001
"""

_GFF3_TWO_TRANSCRIPTS = _GFF3_ONE_TRANSCRIPT + """\
chr1\tensembl\tmRNA\t1000\t4000\t.\t+\t.\tID=transcript:ENST002;Parent=gene:ENSG001;biotype=protein_coding
chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tID=exon:ENSE003;Parent=transcript:ENST002
chr1\tensembl\tCDS\t1100\t2000\t.\t+\t0\tParent=transcript:ENST002
chr1\tensembl\texon\t3000\t4000\t.\t+\t.\tID=exon:ENSE004;Parent=transcript:ENST002
chr1\tensembl\tCDS\t3000\t3900\t.\t+\t0\tParent=transcript:ENST002
"""

_GFF3_NO_GENE_OR_CDS = """\
chr1\tensembl\tmRNA\t1000\t5000\t.\t+\t.\tID=transcript:ENST001;Parent=gene:ENSG001
chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tID=exon:ENSE001;Parent=transcript:ENST001
"""


class DummySequenceSource:
    """Minimal sequence source fixture for slice-length assertions."""

    def __init__(self, lengths: dict[str, int], sequences: dict[str, str] | None = None) -> None:
        self.lengths = lengths
        self.sequences = sequences or {}
        self.closed = False

    def get_sequence_by_start_end_strand(self, seq_region_name, start, end, strand):  # pragma: no cover - helper
        return self.sequences.get(seq_region_name)

    def get_seq_region_length(self, seq_region_name: str) -> int | None:
        return self.lengths.get(seq_region_name)

    def close(self) -> None:
        self.closed = True


@define(kw_only=True, auto_attribs=True)
class DummyGXFTranscriptSource(GXFTranscriptSource[object]):
    """Concrete source used to exercise the shared GXF assembly helpers."""

    transcript_feature_types = frozenset({"mrna", "transcript", "lnc_rna", "pseudogenic_transcript", "processed_transcript"})
    adaptor: object = field(init=False, default=None)

    def get_record_id(self, record):
        return record.attributes.get("ID", (None,))[0] or record.attributes.get("transcript_id", (None,))[0]

    def get_parent_id(self, record):
        return record.attributes.get("Parent", (None,))[0] or record.attributes.get("gene_id", (None,))[0]

    def iter_all_features(self):  # pragma: no cover - fixture helper
        return iter(())

    def iter_features_by_region(self, region=None, *, seq_region_name=None, start=None, end=None):  # pragma: no cover
        return iter(())


class TestGroupRecords:
    """Tests for the shared record-grouping helper."""

    def test_group_records_uses_shared_grouping(self) -> None:
        """Test grouping delegates to the shared per-transcript grouping logic."""

        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        source = DummyGXFTranscriptSource()

        groups = source._group_records(records)  # pylint: disable=protected-access

        assert len(groups) == 1
        assert isinstance(groups[0], _GXFTranscriptGroup)

    def test_group_records_with_multiple_transcripts(self) -> None:
        """Test multiple transcripts are grouped into separate groups."""

        records = GFF3Parser.parse_text(_GFF3_TWO_TRANSCRIPTS)
        source = DummyGXFTranscriptSource()

        groups = source._group_records(records)  # pylint: disable=protected-access

        assert len(groups) == 2
        assert isinstance(groups[0], _GXFTranscriptGroup)
        assert isinstance(groups[1], _GXFTranscriptGroup)

    def test_group_records_without_gene_or_cds(self) -> None:
        """Test a transcript with no gene or CDS records is still grouped."""

        records = GFF3Parser.parse_text(_GFF3_NO_GENE_OR_CDS)
        source = DummyGXFTranscriptSource()

        groups = source._group_records(records)  # pylint: disable=protected-access

        assert len(groups) == 1
        assert isinstance(groups[0], _GXFTranscriptGroup)


class TestResolveGeneRecord:
    """Tests for resolving a transcript group's gene record."""

    def test_returns_group_gene_record_when_present(self) -> None:
        """Test the group's own gene record is preferred when present."""

        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        source = DummyGXFTranscriptSource()
        group = source._group_records(records)[0]  # pylint: disable=protected-access

        gene_record = source._resolve_gene_record(group, {})  # pylint: disable=protected-access

        assert gene_record is group.gene_record

    def test_falls_back_to_previously_seen_gene_records(self) -> None:
        """Test a gene seen earlier in the stream is resolved via the parent id."""

        records = GFF3Parser.parse_text(_GFF3_NO_GENE_OR_CDS)
        source = DummyGXFTranscriptSource()
        group = source._group_records(records)[0]  # pylint: disable=protected-access
        assert group.gene_record is None

        gene_record = source._resolve_gene_record(  # pylint: disable=protected-access
            group, {"gene:ENSG001": records[0]}
        )

        assert gene_record is records[0]

    def test_returns_none_when_gene_record_cannot_be_found(self) -> None:
        """Test None is returned when no gene record is known for the parent id."""

        records = GFF3Parser.parse_text(_GFF3_NO_GENE_OR_CDS)
        source = DummyGXFTranscriptSource()
        group = source._group_records(records)[0]  # pylint: disable=protected-access

        gene_record = source._resolve_gene_record(group, {})  # pylint: disable=protected-access

        assert gene_record is None


class TestResolveGene:
    """Tests for attaching a Gene to a Transcript."""

    def test_no_op_when_gene_record_is_none(self) -> None:
        """Test nothing is attached when there is no gene record."""

        records = GFF3Parser.parse_text(_GFF3_NO_GENE_OR_CDS)
        source = DummyGXFTranscriptSource()
        group = source._group_records(records)[0]  # pylint: disable=protected-access
        transcript = source._build_transcript_from_record_group(group, {})  # pylint: disable=protected-access

        source._resolve_gene(None, transcript, {}, {})  # pylint: disable=protected-access

        assert transcript.gene is None

    def test_reuses_cached_gene_for_repeated_gene_id(self) -> None:
        """Test the same Gene instance is reused across transcripts sharing a gene id."""

        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        source = DummyGXFTranscriptSource()
        group = source._group_records(records)[0]  # pylint: disable=protected-access
        slice_cache: dict[str, object] = {}
        gene_cache: dict[str, object] = {}

        transcript_1 = source._build_transcript_from_record_group(group, slice_cache)  # pylint: disable=protected-access
        transcript_2 = source._build_transcript_from_record_group(group, slice_cache)  # pylint: disable=protected-access
        source._resolve_gene(group.gene_record, transcript_1, gene_cache, slice_cache)  # pylint: disable=protected-access
        source._resolve_gene(group.gene_record, transcript_2, gene_cache, slice_cache)  # pylint: disable=protected-access

        assert transcript_1.gene is transcript_2.gene
        assert transcript_1.gene.transcripts == [transcript_1, transcript_2]


class TestBuildTranscriptFromRecordGroup:
    """Tests for assembling a Transcript from one grouped record block."""

    def test_build_transcript_from_group_builds_full_graph(self) -> None:
        """Test a full transcript graph is assembled from grouped records."""

        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        source = DummyGXFTranscriptSource(sequence_source=DummySequenceSource({"chr1": 9000}))
        group = _GXFTranscriptGroup(
            transcript_record=records[1],
            gene_record=records[0],
            exon_records=[records[2], records[4]],
            cds_records=[records[3], records[5]],
        )

        transcript = source._build_transcript_from_record_group(group, {})  # pylint: disable=protected-access

        assert transcript.stable_id == "transcript:ENST001"
        assert transcript.slice.end == 9000
        assert len(transcript.spliced_exons) == 2
        assert len(transcript.introns) == 1
        assert transcript.translation is not None
        assert transcript.translation.start == 1100
        assert transcript.translation.end == 4900

    def test_build_transcript_without_cds_records_has_no_translation(self) -> None:
        """Test a transcript with no CDS records has no translation."""

        records = GFF3Parser.parse_text(_GFF3_NO_GENE_OR_CDS)
        source = DummyGXFTranscriptSource()
        group = source._group_records(records)[0]  # pylint: disable=protected-access

        transcript = source._build_transcript_from_record_group(group, {})  # pylint: disable=protected-access

        assert transcript.translation is None

    def test_exons_and_transcript_share_slice_instance(self) -> None:
        """Test the transcript and its exons share the same Slice instance."""

        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        source = DummyGXFTranscriptSource()
        group = source._group_records(records)[0]  # pylint: disable=protected-access

        transcript = source._build_transcript_from_record_group(group, {})  # pylint: disable=protected-access

        assert transcript.slice is transcript.spliced_exons[0].slice
        assert transcript.slice is transcript.spliced_exons[1].slice


class TestBuildTranslation:
    """Tests for building a Translation from CDS records."""

    def test_returns_none_for_empty_cds_records(self) -> None:
        """Test no translation is built when there are no CDS records."""

        source = DummyGXFTranscriptSource()

        assert source._build_translation([]) is None  # pylint: disable=protected-access

    def test_translation_spans_min_start_to_max_end(self) -> None:
        """Test the translation spans the full range of the CDS records."""

        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        source = DummyGXFTranscriptSource()
        cds_records = [records[3], records[5]]

        translation = source._build_translation(cds_records)  # pylint: disable=protected-access

        assert translation is not None
        assert translation.start == 1100
        assert translation.end == 4900


class TestGetOrBuildSlice:
    """Tests for the shared, cached Slice builder."""

    def test_get_or_build_slice_uses_sequence_source_length(self) -> None:
        """Test the built slice uses the sequence source's region length."""

        source = DummyGXFTranscriptSource(
            sequence_source=DummySequenceSource({"chr1": 42}, {"chr1": "ACGT"})
        )

        slice_cache: dict[str, object] = {}
        seq_slice = source._get_or_build_slice("chr1", slice_cache)  # pylint: disable=protected-access

        assert seq_slice.end == 42
        assert seq_slice.seq_region_length == 42
        assert seq_slice.seq == "ACGT"

    def test_get_or_build_slice_without_sequence_source_uses_placeholder_end(self) -> None:
        """Test the built slice falls back to a placeholder end without a sequence source."""

        source = DummyGXFTranscriptSource()

        seq_slice = source._get_or_build_slice("chr1", {})  # pylint: disable=protected-access

        assert seq_slice.start == 1
        assert seq_slice.seq_region_length is None
        assert seq_slice.end == 2**31 - 1

    def test_get_or_build_slice_reuses_cached_instance(self) -> None:
        """Test the same Slice instance is returned for the same region."""

        source = DummyGXFTranscriptSource()
        slice_cache: dict[str, object] = {}

        first = source._get_or_build_slice("chr1", slice_cache)  # pylint: disable=protected-access
        second = source._get_or_build_slice("chr1", slice_cache)  # pylint: disable=protected-access

        assert first is second

    def test_get_or_build_slice_different_regions_get_different_instances(self) -> None:
        """Test different regions produce distinct Slice instances."""

        source = DummyGXFTranscriptSource()
        slice_cache: dict[str, object] = {}

        chr1_slice = source._get_or_build_slice("chr1", slice_cache)  # pylint: disable=protected-access
        chr2_slice = source._get_or_build_slice("chr2", slice_cache)  # pylint: disable=protected-access

        assert chr1_slice is not chr2_slice
        assert chr1_slice.seq_region_name == "chr1"
        assert chr2_slice.seq_region_name == "chr2"


class TestFirstAttr:
    """Tests for the shared attribute lookup helper."""

    def test_returns_first_value_for_known_key(self) -> None:
        """Test the first attribute value is returned for a known key."""

        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)

        assert DummyGXFTranscriptSource._first_attr(records[0], "Name") == "BRCA1"  # pylint: disable=protected-access

    def test_returns_none_for_missing_key(self) -> None:
        """Test None is returned when the attribute key is absent."""

        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)

        assert DummyGXFTranscriptSource._first_attr(records[0], "missing") is None  # pylint: disable=protected-access
