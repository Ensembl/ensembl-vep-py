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

from pathlib import Path
import subprocess
from collections.abc import Iterator, Mapping
from typing import Any, Protocol

from attrs import define, field
import pytest

from ensembl.vep.caches.feature_pool import TranscriptPool
from ensembl.vep.common import CoordSystemLevel, Strand, is_resolved
from ensembl.vep.models.genic import Exon, Gene, Translation
from ensembl.vep.models.genic.transcript import _CODING_DERIVATIONS, Transcript
from ensembl.vep.models.sequence.coord_system import CoordSystem
from ensembl.vep.parsers import GFF3Parser, GFF3Record
from ensembl.vep.sources.sequence_source import SequenceSource
from ensembl.vep.sources.file.gff3.transcript import GFF3TranscriptSource
from ensembl.vep.models.sequence import Slice
from ensembl.vep.sources.file.gxf.transcript import GXFTranscriptSource, _GXFRecord, _GXFTranscriptGroup
from ensembl.vep.sources.synonym import SequenceRegionResolver


class _HtslibToolchain(Protocol):
    """Structural stand-in for the session ``htslib_toolchain`` fixture."""

    bgzip: str
    tabix: str


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

_GFF3_DUPLICATE_TRANSCRIPT_ROW_LOSES_PARENT = """\
chr1\tensembl\tgene\t1000\t5000\t.\t+\t.\tID=gene:ENSG001;Name=BRCA1
chr1\tensembl\tmRNA\t1000\t5000\t.\t+\t.\tID=transcript:ENST001;Parent=gene:ENSG001
chr1\tensembl\tmRNA\t1000\t5000\t.\t+\t.\tID=transcript:ENST001
"""


def _filter_records_by_feature_type(
    records: list[GFF3Record],
    feature_types: frozenset[str] | None,
) -> list[GFF3Record]:
    """Mirror the parser-side feature type filter used by indexed stage-one queries."""

    if feature_types is None:
        return records
    return [record for record in records if record.feature_type.lower() in feature_types]


@define(auto_attribs=True)
class DummySequenceSource(SequenceSource):
    """Minimal sequence source fixture for slice-length assertions."""

    lengths: dict[str, int]
    sequences: dict[str, str] = field(factory=dict)
    closed: bool = field(default=False, init=False)

    def get_sequence_by_start_end_strand(  # pylint: disable=unused-argument
        self, seq_region_name: str, start: int, end: int, *, strand: Strand = Strand.FORWARD
    ) -> str | None:
        """Return the whole stored sequence, ignoring the requested span."""
        return self.sequences.get(seq_region_name)

    def close(self) -> None:
        """Record that the source was closed."""
        self.closed = True

    def _fetch_sequence_chunk(self, seq_region_name: str, start: int, end: int) -> str:
        """Return Ns; spans are served by ``get_sequence_by_start_end_strand``."""
        return "N" * (end - start + 1)

    def _get_seq_region_length(self, seq_region_name: str) -> int | None:
        """Return the configured length for a sequence region."""
        return self.lengths.get(seq_region_name)


@define(kw_only=True, auto_attribs=True)
class DummyGXFTranscriptSource(GXFTranscriptSource[GFF3Parser]):
    """Concrete source used to exercise the shared GXF assembly helpers."""

    transcript_feature_types = frozenset(
        {"mrna", "transcript", "lnc_rna", "pseudogenic_transcript", "processed_transcript"}
    )
    adaptor: GFF3Parser = field(init=False, factory=GFF3Parser)

    def get_record_id(self, record: _GXFRecord) -> str | None:
        """Return the GFF3 ``ID`` or GTF ``transcript_id``."""
        return record.attributes.get("ID", (None,))[0] or record.attributes.get("transcript_id", (None,))[0]

    def get_parent_ids(self, record: _GXFRecord) -> tuple[str, ...]:
        """Return the GFF3 ``Parent`` or GTF ``gene_id``."""
        gff_parent_ids = record.attributes.get("Parent")
        if gff_parent_ids is not None:
            return gff_parent_ids
        return record.attributes.get("gene_id", ())

    def _query_typed_records(  # pragma: no cover
        self, region: str, *, feature_types: frozenset[str] | None = None
    ) -> list[Any]:
        return []

    def iter_all_features(self) -> Iterator[Transcript]:  # pragma: no cover - fixture helper
        """Yield nothing."""
        return iter(())


class TestGroupRecords:
    """Tests for the shared record-grouping helper."""

    def test_group_records_uses_shared_grouping(self) -> None:
        """Test grouping delegates to the shared per-transcript grouping logic."""

        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        source = DummyGXFTranscriptSource(file_path="")

        groups = source._group_records(records)  # pylint: disable=protected-access

        assert len(groups) == 1
        assert isinstance(groups[0], _GXFTranscriptGroup)

    def test_duplicate_transcript_row_keeps_an_earlier_gene_id(self) -> None:
        """A later duplicate row with no Parent does not clear the gene id an earlier row gave."""

        records = GFF3Parser.parse_text(_GFF3_DUPLICATE_TRANSCRIPT_ROW_LOSES_PARENT)
        source = DummyGXFTranscriptSource(file_path="")

        (group,) = source._group_records(records)  # pylint: disable=protected-access

        assert group.gene_id == "gene:ENSG001"
        assert group.gene_record is records[0]
        # The record itself still comes from the last duplicate row seen.
        assert group.transcript_record is records[2]


class TestCompleteRegionFeatures:
    """Complete transcript assembly from indexed regions."""

    def test_distant_exons_and_cds_are_loaded_from_the_full_transcript_span(
        self,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """A narrow hit reloads the whole transcript before building its model."""
        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        source = DummyGXFTranscriptSource(
            file_path="",
            sequence_source=DummySequenceSource({"chr1": 9000}),
        )
        calls: list[tuple[int, int]] = []

        def records_for_region(
            _source: DummyGXFTranscriptSource,
            _region: str,
            start: int,
            end: int,
            *,
            feature_types: frozenset[str] | None = None,
        ) -> list[GFF3Record]:
            calls.append((start, end))
            if (start, end) == (1500, 1600):
                return _filter_records_by_feature_type(records[:4], feature_types)
            return _filter_records_by_feature_type(records, feature_types)

        monkeypatch.setattr(DummyGXFTranscriptSource, "_records_for_region", records_for_region)

        transcripts = source.get_features_by_region(seq_region_name="chr1", start=1500, end=1600)

        assert len(transcripts) == 1
        assert len(transcripts[0].spliced_exons) == 2
        assert transcripts[0].translation is not None
        assert calls == [(1500, 1600), (1000, 5000)]

    def test_overlapping_candidate_spans_are_merged_for_stage_two_query(
        self,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Overlapping candidate transcript spans are reloaded with one combined full-span query."""
        records = GFF3Parser.parse_text(_GFF3_TWO_TRANSCRIPTS)
        source = DummyGXFTranscriptSource(
            file_path="",
            sequence_source=DummySequenceSource({"chr1": 9000}),
        )
        stage_two_calls: list[tuple[int | None, int | None]] = []

        def records_for_region(
            _source: DummyGXFTranscriptSource,
            _region: str,
            start: int | None,
            end: int | None,
            *,
            feature_types: frozenset[str] | None = None,
        ) -> list[GFF3Record]:
            if feature_types is None:
                stage_two_calls.append((start, end))
            return _filter_records_by_feature_type(records, feature_types)

        monkeypatch.setattr(DummyGXFTranscriptSource, "_records_for_region", records_for_region)

        transcripts = source.get_features_by_region(seq_region_name="chr1", start=1500, end=1600)

        assert [transcript.stable_id for transcript in transcripts] == [
            "transcript:ENST001",
            "transcript:ENST002",
        ]
        assert stage_two_calls == [(1000, 5000)]

    def test_stage_one_query_is_filtered_to_transcript_feature_types(
        self,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Stage 1 passes ``transcript_feature_types``; stage 2's span query passes no filter."""
        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        source = DummyGXFTranscriptSource(
            file_path="",
            sequence_source=DummySequenceSource({"chr1": 9000}),
        )
        calls: list[tuple[tuple[int, int], frozenset[str] | None]] = []

        def records_for_region(
            _source: DummyGXFTranscriptSource,
            _region: str,
            start: int,
            end: int,
            *,
            feature_types: frozenset[str] | None = None,
        ) -> list[GFF3Record]:
            calls.append(((start, end), feature_types))
            return _filter_records_by_feature_type(records, feature_types)

        monkeypatch.setattr(DummyGXFTranscriptSource, "_records_for_region", records_for_region)

        source.get_features_by_region(seq_region_name="chr1", start=1500, end=1600)

        assert calls == [
            ((1500, 1600), DummyGXFTranscriptSource.transcript_feature_types),
            ((1000, 5000), None),
        ]

    def test_resolves_reference_name_to_an_advertised_source_name(
        self,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Use the GXF name for indexed lookup when slices use a FASTA name."""
        source = DummyGXFTranscriptSource(
            file_path="",
            source_region_resolver=SequenceRegionResolver(source_names=frozenset({"1"})),
        )
        calls: list[str] = []

        def query_records(
            _source: DummyGXFTranscriptSource,
            sequence_region: str,
            _start: int,
            _end: int,
            *,
            feature_types: frozenset[str] | None = None,
        ) -> list[GFF3Record]:
            del feature_types
            calls.append(sequence_region)
            return []

        monkeypatch.setattr(DummyGXFTranscriptSource, "_query_records", query_records)

        source._records_for_region("chr1", 100, 200)  # pylint: disable=protected-access

        assert calls == ["1"]

    def test_retries_chr_spelling_when_the_source_has_no_region_metadata(
        self,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Give GTF and headerless GFF sources a raw-name fallback."""
        source = DummyGXFTranscriptSource(
            file_path="",
            source_region_resolver=SequenceRegionResolver(source_names=frozenset()),
        )
        calls: list[str] = []

        def query_records(
            _source: DummyGXFTranscriptSource,
            sequence_region: str,
            _start: int,
            _end: int,
            *,
            feature_types: frozenset[str] | None = None,
        ) -> list[GFF3Record]:
            del feature_types
            calls.append(sequence_region)
            return [object()] if sequence_region == "1" else []

        monkeypatch.setattr(DummyGXFTranscriptSource, "_query_records", query_records)

        records = source._records_for_region("chr1", 100, 200)  # pylint: disable=protected-access

        assert calls == ["chr1", "1"]
        assert len(records) == 1

    def test_retries_chr_spelling_when_region_metadata_is_partial(
        self,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Treat an undeclared query name like a headerless source region."""

        source = DummyGXFTranscriptSource(
            file_path="",
            source_region_resolver=SequenceRegionResolver(source_names=frozenset({"1"})),
        )
        calls: list[str] = []

        def query_records(
            _source: DummyGXFTranscriptSource,
            sequence_region: str,
            _start: int,
            _end: int,
            *,
            feature_types: frozenset[str] | None = None,
        ) -> list[GFF3Record]:
            del feature_types
            calls.append(sequence_region)
            return [object()] if sequence_region == "2" else []

        monkeypatch.setattr(DummyGXFTranscriptSource, "_query_records", query_records)

        records = source._records_for_region("chr2", 100, 200)  # pylint: disable=protected-access

        assert calls == ["chr2", "2"]
        assert len(records) == 1

    def test_group_records_with_multiple_transcripts(self) -> None:
        """Test multiple transcripts are grouped into separate groups."""

        records = GFF3Parser.parse_text(_GFF3_TWO_TRANSCRIPTS)
        source = DummyGXFTranscriptSource(file_path="")

        groups = source._group_records(records)  # pylint: disable=protected-access

        assert len(groups) == 2
        assert isinstance(groups[0], _GXFTranscriptGroup)
        assert isinstance(groups[1], _GXFTranscriptGroup)

    def test_group_records_without_gene_or_cds(self) -> None:
        """Test a transcript with no gene or CDS records is still grouped."""

        records = GFF3Parser.parse_text(_GFF3_NO_GENE_OR_CDS)
        source = DummyGXFTranscriptSource(file_path="")

        groups = source._group_records(records)  # pylint: disable=protected-access

        assert len(groups) == 1
        assert isinstance(groups[0], _GXFTranscriptGroup)


class TestResolveGeneRecord:
    """Tests for resolving a transcript group's gene record."""

    def test_returns_group_gene_record_when_present(self) -> None:
        """Test the group's own gene record is preferred when present."""

        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        source = DummyGXFTranscriptSource(file_path="")
        group = source._group_records(records)[0]  # pylint: disable=protected-access

        gene_record = source._resolve_gene_record(group, {})  # pylint: disable=protected-access

        assert gene_record is group.gene_record

    def test_falls_back_to_previously_seen_gene_records(self) -> None:
        """Test a gene seen earlier in the stream is resolved via the parent id."""

        records = GFF3Parser.parse_text(_GFF3_NO_GENE_OR_CDS)
        source = DummyGXFTranscriptSource(file_path="")
        group = source._group_records(records)[0]  # pylint: disable=protected-access
        assert group.gene_record is None

        gene_record = source._resolve_gene_record(  # pylint: disable=protected-access
            group, {"gene:ENSG001": records[0]}
        )

        assert gene_record is records[0]

    def test_returns_none_when_gene_record_cannot_be_found(self) -> None:
        """Test None is returned when no gene record is known for the parent id."""

        records = GFF3Parser.parse_text(_GFF3_NO_GENE_OR_CDS)
        source = DummyGXFTranscriptSource(file_path="")
        group = source._group_records(records)[0]  # pylint: disable=protected-access

        gene_record = source._resolve_gene_record(group, {})  # pylint: disable=protected-access

        assert gene_record is None


class TestResolveGene:
    """Tests for attaching a Gene to a Transcript."""

    def test_no_op_when_gene_record_is_none(self) -> None:
        """Test nothing is attached when there is no gene record."""

        records = GFF3Parser.parse_text(_GFF3_NO_GENE_OR_CDS)
        source = DummyGXFTranscriptSource(file_path="", sequence_source=DummySequenceSource({"chr1": 9000}))
        group = source._group_records(records)[0]  # pylint: disable=protected-access
        transcript = source._build_transcript_from_record_group(group)  # pylint: disable=protected-access

        source._resolve_gene(None, None, transcript, {})  # pylint: disable=protected-access

        assert transcript.gene is None

    def test_reuses_cached_gene_for_repeated_gene_id(self) -> None:
        """Test the same Gene instance is reused across transcripts sharing a gene id."""

        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        source = DummyGXFTranscriptSource(file_path="", sequence_source=DummySequenceSource({"chr1": 9000}))
        group = source._group_records(records)[0]  # pylint: disable=protected-access
        # slice_cache: dict[str, Slice] = {}
        gene_cache: dict[str, Gene] = {}

        # pylint: disable-next=protected-access
        transcript_1 = source._build_transcript_from_record_group(group)
        # pylint: disable-next=protected-access
        transcript_2 = source._build_transcript_from_record_group(group)
        # pylint: disable-next=protected-access
        source._resolve_gene(group.gene_record, group.gene_id, transcript_1, gene_cache)
        # pylint: disable-next=protected-access
        source._resolve_gene(group.gene_record, group.gene_id, transcript_2, gene_cache)

        assert transcript_1.gene is transcript_2.gene
        assert transcript_1.gene is not None
        assert transcript_1.gene.transcripts == [transcript_1, transcript_2]


class TestBuildTranscriptFromRecordGroup:
    """Tests for assembling a Transcript from one grouped record block."""

    def test_build_transcript_from_group_builds_full_graph(self) -> None:
        """Test a full transcript graph is assembled from grouped records."""

        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        source = DummyGXFTranscriptSource(file_path="", sequence_source=DummySequenceSource({"chr1": 9000}))
        group = _GXFTranscriptGroup(
            transcript_id="transcript:ENST001",
            gene_id="gene:ENSG001",
            transcript_record=records[1],
            gene_record=records[0],
            exon_records=[records[2], records[4]],
            cds_records=[records[3], records[5]],
        )

        transcript = source._build_transcript_from_record_group(group)  # pylint: disable=protected-access

        assert transcript.stable_id == "transcript:ENST001"
        assert transcript.slice is not None
        assert transcript.slice.end == 9000
        assert len(transcript.spliced_exons) == 2
        assert len(transcript.introns) == 1
        assert transcript.translation is not None
        # Offsets from the transcript's start, not genomic coordinates: CDS 1100-4900 on a
        # transcript starting at 1000.
        assert transcript.translation.start == 101
        assert transcript.translation.end == 3901

    def test_build_transcript_without_cds_records_has_no_translation(self) -> None:
        """Test a transcript with no CDS records has no translation."""

        records = GFF3Parser.parse_text(_GFF3_NO_GENE_OR_CDS)
        source = DummyGXFTranscriptSource(file_path="", sequence_source=DummySequenceSource({"chr1": 9000}))
        group = source._group_records(records)[0]  # pylint: disable=protected-access

        transcript = source._build_transcript_from_record_group(group)  # pylint: disable=protected-access

        assert transcript.translation is None

    def test_exons_and_transcript_share_slice_instance(self) -> None:
        """Test the transcript and its exons share the same Slice instance."""

        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        source = DummyGXFTranscriptSource(file_path="", sequence_source=DummySequenceSource({"chr1": 9000}))
        group = source._group_records(records)[0]  # pylint: disable=protected-access

        transcript = source._build_transcript_from_record_group(group)  # pylint: disable=protected-access

        assert transcript.slice is transcript.spliced_exons[0].slice
        assert transcript.slice is transcript.spliced_exons[1].slice


class TestBuildTranslation:
    """Tests for building a Translation from CDS records."""

    @staticmethod
    def _transcript() -> Transcript:
        """The transcript the offsets are relative to, built the way the source builds it."""

        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        source = DummyGXFTranscriptSource(file_path="", sequence_source=DummySequenceSource({"chr1": 9000}))
        group = source._group_records(records)[0]  # pylint: disable=protected-access
        return source._build_transcript_from_record_group(group)  # pylint: disable=protected-access

    def test_returns_none_for_empty_cds_records(self) -> None:
        """Test no translation is built when there are no CDS records."""

        source = DummyGXFTranscriptSource(file_path="", sequence_source=DummySequenceSource({"chr1": 9000}))

        # pylint: disable-next=protected-access
        assert source._build_translation([], self._transcript()) is None

    def test_translation_is_relative_to_the_transcript_start(self) -> None:
        """The CDS spans 1100-4900 on a transcript starting at 1000, so the translation runs from
        offset 101 to 3901.

        Translation.start and end are offsets from the transcript's start, and Transcript re-maps
        them with `self.start + start - 1`. Asserting the raw genomic span here is what let a real
        bug through: the coding boundaries landed outside the transcript, `cds` resolved to None,
        and every coding consequence vanished while positional ones survived.
        """

        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        source = DummyGXFTranscriptSource(file_path="", sequence_source=DummySequenceSource({"chr1": 9000}))
        cds_records = [records[3], records[5]]

        # pylint: disable-next=protected-access
        translation = source._build_translation(cds_records, self._transcript())

        assert translation is not None
        assert translation.start == 101
        assert translation.end == 3901


class TestFirstAttr:
    """Tests for the shared attribute lookup helper."""

    def test_returns_first_value_for_known_key(self) -> None:
        """Test the first attribute value is returned for a known key."""

        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        first_attr = DummyGXFTranscriptSource._first_attr  # pylint: disable=protected-access

        assert first_attr(records[0], "Name") == "BRCA1"

    def test_returns_none_for_missing_key(self) -> None:
        """Test None is returned when the attribute key is absent."""

        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        first_attr = DummyGXFTranscriptSource._first_attr  # pylint: disable=protected-access

        assert first_attr(records[0], "missing") is None

    def test_returns_first_value_for_multi_valued_key(self) -> None:
        """Test only the first value is returned for a multi-valued key."""

        records = GFF3Parser.parse_text("chr1\tensembl\tgene\t1\t10\t.\t+\t.\tAlias=a,b,c\n")
        first_attr = DummyGXFTranscriptSource._first_attr  # pylint: disable=protected-access

        assert first_attr(records[0], "Alias") == "a"

    def test_works_against_a_plain_python_implementation(self) -> None:
        """Test ``_first_attr`` delegates correctly to any ``_GXFRecord``, not only the Rust type."""

        @define(kw_only=True, auto_attribs=True)
        class PlainRecord:
            """A minimal ``_GXFRecord`` implementing ``first_attribute`` in plain Python."""

            feature_type: str = "gene"
            seqid: str = "chr1"
            start: int = 1
            end: int = 10
            strand_symbol: str = "+"
            phase_value: int | None = None
            attributes: dict[str, tuple[str, ...]] = field(factory=dict)

            def first_attribute(self, key: str) -> str | None:
                """Return the first attribute value for a key, or None when absent."""
                values = self.attributes.get(key)
                return values[0] if values else None

        record = PlainRecord(attributes={"Name": ("BRCA1",)})
        first_attr = DummyGXFTranscriptSource._first_attr  # pylint: disable=protected-access

        assert first_attr(record, "Name") == "BRCA1"
        assert first_attr(record, "missing") is None


class TestBuilderFastPath:
    """The builder's fast path (object.__setattr__, seed_attrs) matches ordinary attribute
    assignment (attrs on_setattr hooks, seed_and_invalidate_attrs) field-for-field."""

    def test_fast_built_transcript_matches_hooked_assignment(self) -> None:
        """Compare the builder's output against a transcript assembled the slow, hooked way."""

        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        source = DummyGXFTranscriptSource(file_path="", sequence_source=DummySequenceSource({"chr1": 9000}))
        group = source._group_records(records)[0]  # pylint: disable=protected-access

        fast = source._build_transcript_from_record_group(group)  # pylint: disable=protected-access

        reference_slice = Slice(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME, default=True),
            seq_region_name="chr1",
            start=1,
            end=9000,
            strand=Strand.FORWARD,
            sequence_source=DummySequenceSource({"chr1": 9000}),
        )
        reference = Transcript(
            stable_id="transcript:ENST001", start=1000, end=5000, strand=Strand.FORWARD, slice=reference_slice
        )
        # Ordinary attribute assignment: goes through spliced_exons' on_setattr hook.
        reference.biotype = "protein_coding"
        reference.spliced_exons = [
            Exon(
                stable_id="exon:ENSE001", start=1000, end=2000, strand=Strand.FORWARD, slice=reference_slice
            ),
            Exon(
                stable_id="exon:ENSE002", start=3000, end=5000, strand=Strand.FORWARD, slice=reference_slice
            ),
        ]
        reference_translation = Translation(start=101, end=3901, stable_id="")
        reference_translation.transcript = reference
        reference.seed_and_invalidate_attrs(translation=reference_translation)

        assert fast.stable_id == reference.stable_id
        assert fast.biotype == reference.biotype
        assert [(exon.start, exon.end) for exon in fast.spliced_exons] == [
            (exon.start, exon.end) for exon in reference.spliced_exons
        ]
        assert [(intron.start, intron.end) for intron in fast.introns] == [
            (intron.start, intron.end) for intron in reference.introns
        ]
        assert (fast.coding_start, fast.coding_end) == (reference.coding_start, reference.coding_end)
        assert fast.cds_genomic_spans == reference.cds_genomic_spans
        assert fast.cds_mapper is not None and reference.cds_mapper is not None
        assert fast.cds_mapper.length == reference.cds_mapper.length
        assert fast.cdna_mapper is not None and reference.cdna_mapper is not None
        assert fast.cdna_mapper.length == reference.cdna_mapper.length
        assert fast.translation is not None and reference.translation is not None
        assert (fast.translation.start, fast.translation.end) == (
            reference.translation.start,
            reference.translation.end,
        )

    def test_fast_build_resolves_no_coding_derivations(self) -> None:
        """Nothing in the coding derivation closure is computed as a side effect of building."""

        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        source = DummyGXFTranscriptSource(file_path="", sequence_source=DummySequenceSource({"chr1": 9000}))
        group = source._group_records(records)[0]  # pylint: disable=protected-access

        transcript = source._build_transcript_from_record_group(group)  # pylint: disable=protected-access

        for name in _CODING_DERIVATIONS:
            assert not is_resolved(transcript, name), f"{name} was resolved as a side effect of building"
        assert not is_resolved(transcript, "_introns")
        assert not is_resolved(transcript, "_spliced_seq")


class TestCandidateFiltering:
    """Non-candidate transcripts pulled in by the merged-span query are dropped before any model
    is built for them -- not built, then filtered from the result afterwards."""

    def test_non_candidate_groups_are_not_built(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Only the transcript found by the narrow stage-1 query is built."""

        records = GFF3Parser.parse_text(_GFF3_TWO_TRANSCRIPTS)
        source = DummyGXFTranscriptSource(file_path="", sequence_source=DummySequenceSource({"chr1": 9000}))
        build_calls: list[str] = []
        # pylint: disable-next=protected-access
        original_build = GXFTranscriptSource._build_transcript_from_record_group

        def counting_build(self, group):  # type: ignore[no-untyped-def]
            build_calls.append(group.transcript_id)
            return original_build(self, group)

        monkeypatch.setattr(DummyGXFTranscriptSource, "_build_transcript_from_record_group", counting_build)

        def records_for_region(
            _source: DummyGXFTranscriptSource,
            _region: str,
            start: int,
            end: int,
            *,
            feature_types: frozenset[str] | None = None,
        ) -> list[GFF3Record]:
            del feature_types
            if (start, end) == (1500, 1600):
                # Stage 1: the narrow query only ever touches the ENST001 transcript row.
                return [records[1]]
            # Stage 2: the merged-span query returns the whole file, including ENST002 and its
            # shared gene -- neither of which was ever a candidate.
            return records

        monkeypatch.setattr(DummyGXFTranscriptSource, "_records_for_region", records_for_region)

        transcripts = source.get_features_by_region(seq_region_name="chr1", start=1500, end=1600)

        assert [transcript.stable_id for transcript in transcripts] == ["transcript:ENST001"]
        assert build_calls == ["transcript:ENST001"]

    def test_gene_transcripts_lists_only_kept_transcripts(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """A shared gene's ``transcripts`` list reflects only the candidate transcripts that were
        actually built, not every transcript the merged-span query happened to return.

        This is a deliberate behaviour change: filtering by candidate id now happens while
        grouping, before any model (including the shared Gene) is built, rather than on the
        finished transcript list afterwards. Nothing in the production code reads
        ``Gene.transcripts`` back (only tests do), so this is safe.
        """

        records = GFF3Parser.parse_text(_GFF3_TWO_TRANSCRIPTS)
        source = DummyGXFTranscriptSource(file_path="", sequence_source=DummySequenceSource({"chr1": 9000}))

        def records_for_region(
            _source: DummyGXFTranscriptSource,
            _region: str,
            start: int,
            end: int,
            *,
            feature_types: frozenset[str] | None = None,
        ) -> list[GFF3Record]:
            del feature_types
            if (start, end) == (1500, 1600):
                return [records[1]]
            return records

        monkeypatch.setattr(DummyGXFTranscriptSource, "_records_for_region", records_for_region)

        (transcript,) = source.get_features_by_region(seq_region_name="chr1", start=1500, end=1600)

        assert transcript.gene is not None
        assert transcript.gene.transcripts == [transcript]


class TestWholeRegionQuery:
    """A region without coordinates queries the whole sequence region."""

    def test_region_without_bounds_queries_the_whole_region(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Test stage 1 passes no bounds through, so the query covers the contig."""

        queried: list[tuple[str, int | None, int | None]] = []

        def records_for_region(
            _self: DummyGXFTranscriptSource,
            seq_region_name: str,
            start: int | None = None,
            end: int | None = None,
            *,
            feature_types: frozenset[str] | None = None,
        ) -> list[Any]:
            del feature_types
            queried.append((seq_region_name, start, end))
            return []

        monkeypatch.setattr(DummyGXFTranscriptSource, "_records_for_region", records_for_region)
        source = DummyGXFTranscriptSource(file_path="")

        assert source.get_features_by_region("chr1") == []
        assert queried == [("chr1", None, None)]


class TestStreamingIteration:
    """Streaming transcript assembly from sequential parser records."""

    def test_parentless_transcript_streams_as_its_own_gene_block(self, tmp_path: Path) -> None:
        """A transcript row without Parent still forms a streamable synthetic block."""
        path = tmp_path / "parentless.gff3"
        path.write_text(
            "##gff-version 3\n"
            "chr1\tensembl\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript:ORPHAN\n"
            "chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tID=exon:ORPHAN;Parent=transcript:ORPHAN\n"
        )
        source = GFF3TranscriptSource(
            file_path=str(path),
            sequence_source=DummySequenceSource({"chr1": 9000}),
        )

        transcripts = list(source.iter_all_features())

        assert [transcript.stable_id for transcript in transcripts] == ["transcript:ORPHAN"]
        assert transcripts[0].gene is None
        assert [(exon.start, exon.end) for exon in transcripts[0].spliced_exons] == [(1000, 2000)]


class TestTranscriptPoolReuse:
    """Reusing pooled transcripts across ``get_features_by_region`` calls."""

    def test_pool_hit_skips_stage_two_query(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """A transcript already in the pool is reused without a merged-span query."""
        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        pool = TranscriptPool()
        source = DummyGXFTranscriptSource(
            file_path="",
            sequence_source=DummySequenceSource({"chr1": 9000}),
            transcript_pool=pool,
        )

        def records_for_region(
            _source: DummyGXFTranscriptSource,
            _region: str,
            _start: int,
            _end: int,
            *,
            feature_types: frozenset[str] | None = None,
        ) -> list[GFF3Record]:
            return _filter_records_by_feature_type(records, feature_types)

        monkeypatch.setattr(DummyGXFTranscriptSource, "_records_for_region", records_for_region)

        (transcript,) = source.get_features_by_region(seq_region_name="chr1", start=1500, end=1600)
        assert len(pool) == 1

        def stage_two_forbidden(
            _source: DummyGXFTranscriptSource,
            _region: str,
            _start: int,
            _end: int,
            *,
            feature_types: frozenset[str] | None = None,
        ) -> list[GFF3Record]:
            if feature_types is None:
                raise AssertionError("stage 2 must not run on a pool hit")
            return _filter_records_by_feature_type(records, feature_types)

        monkeypatch.setattr(DummyGXFTranscriptSource, "_records_for_region", stage_two_forbidden)

        (cached,) = source.get_features_by_region(seq_region_name="chr1", start=1500, end=1600)

        assert cached is transcript
        assert pool.hit_count == 1

    def test_pool_miss_builds_and_stores(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """A cold pool is populated by the first call that builds a transcript."""
        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        pool = TranscriptPool()
        source = DummyGXFTranscriptSource(
            file_path="",
            sequence_source=DummySequenceSource({"chr1": 9000}),
            transcript_pool=pool,
        )
        monkeypatch.setattr(
            DummyGXFTranscriptSource,
            "_records_for_region",
            lambda _source, _region, start, end, *, feature_types=None: _filter_records_by_feature_type(
                records, feature_types
            ),
        )

        (transcript,) = source.get_features_by_region(seq_region_name="chr1", start=1500, end=1600)

        assert pool.miss_count == 1
        assert pool.get(("chr1", transcript.stable_id, transcript.start, transcript.end)) is transcript

    def test_pooled_conflicting_span_falls_back_to_unpooled_rebuild(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A pooled ID resurfacing at a different span poisons the pool and rebuilds all candidates."""
        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        (conflicting_record,) = GFF3Parser.parse_text(
            "chr1\tensembl\tmRNA\t2000\t6000\t.\t+\t.\tID=transcript:ENST001;Parent=gene:ENSG001\n"
        )
        pool = TranscriptPool()
        source = DummyGXFTranscriptSource(
            file_path="",
            sequence_source=DummySequenceSource({"chr1": 9000}),
            transcript_pool=pool,
        )
        # pylint: disable-next=protected-access
        (cached,) = source._build_transcripts(records, gene_records={}, gene_cache={})
        pool.put(("chr1", cached.stable_id, cached.start, cached.end), cached)
        candidates = {
            "transcript:ENST001": (cached.start, cached.end),
            "transcript:ENST002": (2000, 6000),
        }
        fallback_calls: list[tuple[str, dict[str, tuple[int, int]]]] = []

        def records_for_region(
            _source: DummyGXFTranscriptSource,
            _region: str,
            start: int | None = None,
            end: int | None = None,
            *,
            feature_types: frozenset[str] | None = None,
        ) -> list[GFF3Record]:
            del start, end, feature_types
            return [conflicting_record]

        def fallback_build_candidates(
            _source: DummyGXFTranscriptSource,
            seq_region_name: str,
            rebuild_candidates: Mapping[str, tuple[int, int]],
        ) -> list[Transcript]:
            fallback_calls.append((seq_region_name, dict(rebuild_candidates)))
            return [cached]

        monkeypatch.setattr(DummyGXFTranscriptSource, "_records_for_region", records_for_region)
        monkeypatch.setattr(DummyGXFTranscriptSource, "_build_candidates", fallback_build_candidates)

        # pylint: disable-next=protected-access
        rebuilt = source._build_candidates_with_pool("chr1", candidates, pool)

        assert rebuilt == [cached]
        assert source.transcript_pool_poisoned
        assert len(pool) == 0
        assert fallback_calls == [("chr1", candidates)]

    def test_disabled_pool_behaves_like_no_pool(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """``max_bytes <= 0`` disables reuse even when a pool instance is configured."""
        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        pool = TranscriptPool(max_bytes=0)
        source = DummyGXFTranscriptSource(
            file_path="",
            sequence_source=DummySequenceSource({"chr1": 9000}),
            transcript_pool=pool,
        )
        stage_two_calls: list[tuple[int, int]] = []

        def records_for_region(
            _source: DummyGXFTranscriptSource,
            _region: str,
            start: int,
            end: int,
            *,
            feature_types: frozenset[str] | None = None,
        ) -> list[GFF3Record]:
            if feature_types is None:
                stage_two_calls.append((start, end))
            return _filter_records_by_feature_type(records, feature_types)

        monkeypatch.setattr(DummyGXFTranscriptSource, "_records_for_region", records_for_region)

        source.get_features_by_region(seq_region_name="chr1", start=1500, end=1600)
        source.get_features_by_region(seq_region_name="chr1", start=1500, end=1600)

        assert stage_two_calls == [(1000, 5000), (1000, 5000)]
        assert len(pool) == 0

    def test_duplicate_candidate_id_falls_back_to_non_pooled_path(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A transcript id repeated at a different span in stage 1 disables pooling for that call.

        The second row must have a *different* span from the first -- the same row turning up
        twice (a collateral re-fetch, not a duplicate) must never disable the pool; see
        ``TestWideTranscriptIsNotMistakenForADuplicate``.
        """
        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        (other_span_record,) = GFF3Parser.parse_text(
            "chr1\tensembl\tmRNA\t2000\t6000\t.\t+\t.\tID=transcript:ENST001;Parent=gene:ENSG001\n"
        )
        duplicated = [records[0], records[1], other_span_record, records[2], records[3]]
        pool = TranscriptPool()
        source = DummyGXFTranscriptSource(
            file_path="",
            sequence_source=DummySequenceSource({"chr1": 9000}),
            transcript_pool=pool,
        )

        def records_for_region(
            _source: DummyGXFTranscriptSource,
            _region: str,
            start: int,
            end: int,
            *,
            feature_types: frozenset[str] | None = None,
        ) -> list[GFF3Record]:
            if (start, end) == (1500, 1600):
                return _filter_records_by_feature_type(duplicated, feature_types)
            return _filter_records_by_feature_type(records, feature_types)

        monkeypatch.setattr(DummyGXFTranscriptSource, "_records_for_region", records_for_region)

        transcripts = source.get_features_by_region(seq_region_name="chr1", start=1500, end=1600)

        assert [transcript.stable_id for transcript in transcripts] == ["transcript:ENST001"]
        assert len(pool) == 0

    def test_result_order_matches_stage_one_regardless_of_build_order(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Combining a pool hit and a fresh build still returns stage-1 (file) order.

        This is the hard CSQ-order invariant: the result must not depend on whether a candidate
        came from the pool or was just built, nor on the order records were grouped in.
        """
        records = GFF3Parser.parse_text(_GFF3_TWO_TRANSCRIPTS)
        pool = TranscriptPool()
        source = DummyGXFTranscriptSource(
            file_path="",
            sequence_source=DummySequenceSource({"chr1": 9000}),
        )
        groups = source._group_records(records)  # pylint: disable=protected-access
        enst002_group = next(group for group in groups if group.transcript_id == "transcript:ENST002")
        enst002 = source._build_transcript_from_record_group(  # pylint: disable=protected-access
            enst002_group
        )
        pool.put(("chr1", "transcript:ENST002", enst002.start, enst002.end), enst002)
        source.transcript_pool = pool

        def records_for_region(
            _source: DummyGXFTranscriptSource,
            _region: str,
            start: int,
            end: int,
            *,
            feature_types: frozenset[str] | None = None,
        ) -> list[GFF3Record]:
            # A real tabix query only returns records overlapping the queried span; ENST002 is
            # pooled and so must never legitimately turn up when only ENST001's span is queried.
            overlapping = [record for record in records if record.start <= end and record.end >= start]
            return _filter_records_by_feature_type(overlapping, feature_types)

        monkeypatch.setattr(DummyGXFTranscriptSource, "_records_for_region", records_for_region)

        transcripts = source.get_features_by_region(seq_region_name="chr1", start=1000, end=5000)

        assert [transcript.stable_id for transcript in transcripts] == [
            "transcript:ENST001",
            "transcript:ENST002",
        ]
        assert transcripts[1] is enst002


class TestDuplicateTranscriptIdDisablesPool:
    """A repeated transcript id makes pool-on diverge from pool-off unless pooling is disabled.

    ``dup.gff3.gz`` (id ``X`` at 4000-4600 and again, unrelated, at 5000-8000; id ``Y`` at
    4500-4700) reproduces the failure: querying 6000-6100 first pools the second ``X`` alone (1
    exon), then querying 4650-5100 hits that pooled entry by (region, id, start, end) -- the same
    key both times -- while a fresh build for that same call also picks up the first ``X``'s exon
    collaterally through ``Y``'s span query, giving 2 exons and a different order.
    """

    _GFF3_DUPLICATE_ID = """\
chr1\tt\tgene\t4000\t8000\t.\t+\t.\tID=gene:G
chr1\tt\tmRNA\t4000\t4600\t.\t+\t.\tID=transcript:X;Parent=gene:G
chr1\tt\texon\t4000\t4600\t.\t+\t.\tParent=transcript:X
chr1\tt\tmRNA\t4500\t4700\t.\t+\t.\tID=transcript:Y;Parent=gene:G
chr1\tt\texon\t4500\t4700\t.\t+\t.\tParent=transcript:Y
chr1\tt\tmRNA\t5000\t8000\t.\t+\t.\tID=transcript:X;Parent=gene:G
chr1\tt\texon\t5000\t8000\t.\t+\t.\tParent=transcript:X
"""

    @pytest.fixture(name="indexed_duplicate_id_path")
    def fixture_indexed_duplicate_id_path(
        self, tmp_path_factory: pytest.TempPathFactory, htslib_toolchain: _HtslibToolchain
    ) -> Path:
        """The duplicate-id fixture above, bgzipped and tabix-indexed with a real HTSlib toolchain."""
        directory = tmp_path_factory.mktemp("duplicate_transcript_id")
        gff3_path = directory / "dup.gff3"
        gff3_path.write_text(self._GFF3_DUPLICATE_ID)

        subprocess.run([htslib_toolchain.bgzip, "-f", str(gff3_path)], check=True)
        gz_path = gff3_path.with_name(gff3_path.name + ".gz")
        subprocess.run([htslib_toolchain.tabix, "-f", "-p", "gff", str(gz_path)], check=True)
        return gz_path

    @staticmethod
    def _query_both_windows(source: GFF3TranscriptSource) -> list[tuple[str, int, int, int]]:
        """Query the two windows that reproduce the divergence, and describe the second result."""
        source.get_features_by_region(seq_region_name="chr1", start=6000, end=6100)
        second = source.get_features_by_region(seq_region_name="chr1", start=4650, end=5100)
        return [
            (transcript.stable_id, transcript.start, transcript.end, len(transcript.spliced_exons))
            for transcript in second
        ]

    @pytest.mark.integration
    def test_pool_on_matches_pool_off_after_the_first_duplicate_is_seen(
        self, indexed_duplicate_id_path: Path
    ) -> None:
        """Pool-on and pool-off must describe the second window identically."""
        baseline_source = GFF3TranscriptSource(file_path=str(indexed_duplicate_id_path))
        baseline = self._query_both_windows(baseline_source)
        baseline_source.close()

        pooled_source = GFF3TranscriptSource(
            file_path=str(indexed_duplicate_id_path), transcript_pool=TranscriptPool()
        )
        pooled = self._query_both_windows(pooled_source)
        pooled_source.close()

        assert pooled == baseline

    @pytest.mark.integration
    def test_pool_is_disabled_and_cleared_once_a_duplicate_is_seen(
        self, indexed_duplicate_id_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """The pool is cleared and marked poisoned, with exactly one warning logged."""
        pool = TranscriptPool()
        source = GFF3TranscriptSource(file_path=str(indexed_duplicate_id_path), transcript_pool=pool)

        with caplog.at_level("WARNING"):
            self._query_both_windows(source)
            self._query_both_windows(source)

        source.close()

        assert len(pool) == 0
        assert source._pool_poisoned  # pylint: disable=protected-access
        duplicate_warnings = [
            record for record in caplog.records if "repeated transcript id" in record.message
        ]
        assert len(duplicate_warnings) == 1


class TestWideTranscriptIsNotMistakenForADuplicate:
    """A single transcript row must never disable the pool.

    Reproduces the false positive found on real GRCh38 data: ``WIDE`` (1000-9000) is pooled by an
    early query that touches only its middle. A later query asks for two *other* candidates
    (``A`` at 1000-1100, ``B`` at 8900-9000) far enough apart that their merged spans stay
    separate -- each span is queried on its own, and each one also overlaps ``WIDE``, so its one
    row is fetched twice, once per query, with the exact same (start, end) both times. That must
    read as the same row seen again, not two different rows sharing an id.
    """

    _GFF3_WIDE_TRANSCRIPT = """\
chr1\tt\tgene\t1000\t9000\t.\t+\t.\tID=gene:G
chr1\tt\tmRNA\t1000\t9000\t.\t+\t.\tID=transcript:WIDE;Parent=gene:G
chr1\tt\texon\t1000\t9000\t.\t+\t.\tParent=transcript:WIDE
chr1\tt\tmRNA\t1000\t1100\t.\t+\t.\tID=transcript:A;Parent=gene:G
chr1\tt\texon\t1000\t1100\t.\t+\t.\tParent=transcript:A
chr1\tt\tmRNA\t8900\t9000\t.\t+\t.\tID=transcript:B;Parent=gene:G
chr1\tt\texon\t8900\t9000\t.\t+\t.\tParent=transcript:B
"""

    @pytest.fixture(name="indexed_wide_transcript_path")
    def fixture_indexed_wide_transcript_path(
        self, tmp_path_factory: pytest.TempPathFactory, htslib_toolchain: _HtslibToolchain
    ) -> Path:
        """The wide-transcript fixture above, bgzipped and tabix-indexed."""
        directory = tmp_path_factory.mktemp("wide_transcript_not_a_duplicate")
        gff3_path = directory / "wide.gff3"
        gff3_path.write_text(self._GFF3_WIDE_TRANSCRIPT)

        subprocess.run([htslib_toolchain.bgzip, "-f", str(gff3_path)], check=True)
        gz_path = gff3_path.with_name(gff3_path.name + ".gz")
        subprocess.run([htslib_toolchain.tabix, "-f", "-p", "gff", str(gz_path)], check=True)
        return gz_path

    @pytest.mark.integration
    def test_pool_stays_enabled_when_a_pooled_wide_transcript_resurfaces_twice(
        self, indexed_wide_transcript_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """Querying A and B after WIDE is already pooled must not poison the pool."""
        pool = TranscriptPool()
        source = GFF3TranscriptSource(file_path=str(indexed_wide_transcript_path), transcript_pool=pool)

        # Pools WIDE alone: this window touches only its middle, not A or B.
        first = source.get_features_by_region(seq_region_name="chr1", start=4000, end=5000)
        assert [transcript.stable_id for transcript in first] == ["transcript:WIDE"]

        # WIDE is now a pool hit; A and B are missing and queried separately (their spans do not
        # merge), each collaterally re-fetching WIDE's own row at its own, unchanged span.
        with caplog.at_level("WARNING"):
            second = source.get_features_by_region(seq_region_name="chr1", start=1000, end=9000)

        source.close()

        assert [transcript.stable_id for transcript in second] == [
            "transcript:WIDE",
            "transcript:A",
            "transcript:B",
        ]
        assert second[0] is first[0]  # the pooled WIDE instance was reused, not rebuilt
        assert not source._pool_poisoned  # pylint: disable=protected-access
        assert pool.hit_count == 1
        assert not [record for record in caplog.records if "repeated transcript id" in record.message]
