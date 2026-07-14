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
"""Unit tests for the GFF3 transcript source."""

import os
from pathlib import Path

from ensembl.vep.common import ExonPhase, Strand
from ensembl.vep.parsers import GFF3Parser
from ensembl.vep.models.genic import Transcript
from ensembl.vep.sources.file.gff3.transcript import GFF3TranscriptSource


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


_GFF3_ONE_TRANSCRIPT = """\
chr1\tensembl\tgene\t1000\t5000\t.\t+\t.\tID=gene:ENSG001;Name=BRCA1;biotype=protein_coding
chr1\tensembl\tmRNA\t1000\t5000\t.\t+\t.\tID=transcript:ENST001;Parent=gene:ENSG001
chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tID=exon:ENSE001;Parent=transcript:ENST001
chr1\tensembl\tCDS\t1100\t2000\t.\t+\t0\tID=cds:ENSP001;Parent=transcript:ENST001;protein_id=ENSP001
chr1\tensembl\texon\t3000\t5000\t.\t+\t.\tID=exon:ENSE002;Parent=transcript:ENST001
chr1\tensembl\tCDS\t3000\t4900\t.\t+\t0\tID=cds:ENSP001;Parent=transcript:ENST001;protein_id=ENSP001
"""

_GFF3_TWO_TRANSCRIPTS = """\
chr1\tensembl\tmRNA\t1000\t5000\t.\t+\t.\tID=transcript:ENST001;Parent=gene:ENSG001
chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tID=exon:ENSE001;Parent=transcript:ENST001
chr1\tensembl\tmRNA\t2000\t8000\t.\t-\t.\tID=transcript:ENST002;Parent=gene:ENSG002
chr1\tensembl\texon\t2000\t3000\t.\t-\t.\tID=exon:ENSE002;Parent=transcript:ENST002
chr1\tensembl\texon\t6000\t8000\t.\t-\t.\tID=exon:ENSE003;Parent=transcript:ENST002
"""

_GFF3_REVERSE_STRAND = """\
chr2\tensembl\tmRNA\t500\t1100\t.\t-\t.\tID=transcript:ENST003;Parent=gene:ENSG003
chr2\tensembl\texon\t500\t600\t.\t-\t.\tID=exon:ENSE010;Parent=transcript:ENST003
chr2\tensembl\texon\t800\t900\t.\t-\t.\tID=exon:ENSE011;Parent=transcript:ENST003
chr2\tensembl\texon\t1000\t1100\t.\t-\t.\tID=exon:ENSE012;Parent=transcript:ENST003
"""

_GFF3_EXON_WITH_PHASE = """\
chr1\tensembl\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript:ENST001;Parent=gene:ENSG001
chr1\tensembl\texon\t1000\t2000\t.\t+\t1\tID=exon:ENSE001;Parent=transcript:ENST001
"""

_GFF3_NO_ID_OR_NAME = """\
chr1\tensembl\tmRNA\t1000\t2000\t.\t+\t.\tParent=gene:ENSG001
"""

HEADER_LINES = ["##gff-version 3"]
CHR1_TRANSCRIPT_ROW = "chr1\tensembl\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript:ENST001;Parent=gene:ENSG001"
CHR1_EXON_ROW = "chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tID=exon:ENSE001;Parent=transcript:ENST001"
CHR2_TRANSCRIPT_ROW = "chr2\tensembl\tmRNA\t100\t200\t.\t+\t.\tID=transcript:ENST002;Parent=gene:ENSG002"


def _write_gff3(path: Path, data: str = _GFF3_ONE_TRANSCRIPT) -> Path:
    """Write a GFF3 fixture file."""

    path.write_text(data)
    return path


def _write_query_gff3(path: Path, data: str) -> Path:
    """Write a query-mode GFF3 fixture with a placeholder tabix index."""

    path.write_text(data)
    path.with_suffix(path.suffix + ".tbi").write_text("")
    return path


def _write_query_mock_tabix(tmp_path: Path, monkeypatch) -> Path:
    """Write a mock tabix executable on PATH that returns deterministic GFF3 records.

    GFF3TranscriptSource's adaptor is constructed without access to the source's own
    tabix_path field, so the mock is put on PATH under the literal name "tabix" rather
    than passed via the tabix_path constructor argument.
    """

    tabix = tmp_path / "tabix"
    tabix.write_text(
        f"""#!/bin/sh
if [ "$1" = "--version" ]; then
    exit 0
fi

case "$2" in
    chr1)
        printf '%s\\n' '{CHR1_TRANSCRIPT_ROW}'
        printf '%s\\n' '{CHR1_EXON_ROW}'
        ;;
    chr1:1000-2000)
        printf '%s\\n' '{CHR1_TRANSCRIPT_ROW}'
        printf '%s\\n' '{CHR1_EXON_ROW}'
        ;;
    chr2:100-200)
        printf '%s\\n' '{CHR2_TRANSCRIPT_ROW}'
        ;;
esac
"""
    )
    tabix.chmod(0o755)
    monkeypatch.setenv("PATH", f"{tmp_path}{os.pathsep}{os.environ['PATH']}")
    return tabix


class TestGetRecordId:
    """Tests for the GFF3-specific record identifier resolution."""

    def test_get_record_id_uses_id_attribute(self) -> None:
        """Test the ID attribute is used verbatim as the record identifier."""

        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        source = GFF3TranscriptSource()

        assert source.get_record_id(records[1]) == "transcript:ENST001"

    def test_get_record_id_falls_back_to_name(self) -> None:
        """Test Name is used when ID is absent."""

        record = GFF3Parser.parse_text("chr1\tensembl\tgene\t1\t10\t.\t+\t.\tName=BRCA1\n")[0]
        source = GFF3TranscriptSource()

        assert source.get_record_id(record) == "BRCA1"

    def test_get_record_id_synthesizes_id_when_missing(self) -> None:
        """Test a synthetic identifier is built when ID and Name are both absent."""

        record = GFF3Parser.parse_text(_GFF3_NO_ID_OR_NAME)[0]
        source = GFF3TranscriptSource()

        assert source.get_record_id(record) == "mrna:chr1:1000-2000"


class TestGetParentId:
    """Tests for the GFF3-specific parent identifier resolution."""

    def test_get_parent_id_returns_parent_attribute(self) -> None:
        """Test the Parent attribute is returned as the parent identifier."""

        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        source = GFF3TranscriptSource()

        assert source.get_parent_id(records[1]) == "gene:ENSG001"

    def test_get_parent_id_returns_none_when_absent(self) -> None:
        """Test None is returned when the record has no Parent attribute."""

        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        source = GFF3TranscriptSource()

        assert source.get_parent_id(records[0]) is None


class TestGetAllFeatures:
    """Integration-style tests for GFF3-backed transcript building."""

    def test_get_all_features_builds_transcript(self, tmp_path: Path) -> None:
        """Test that a transcript object is built from a GFF3 file."""

        path = _write_gff3(tmp_path / "sample.gff3")
        source = GFF3TranscriptSource(file_path=str(path))

        transcripts = source.get_all_features()

        assert len(transcripts) == 1
        transcript = transcripts[0]
        assert isinstance(transcript, Transcript)
        assert transcript.stable_id == "transcript:ENST001"
        assert transcript.gene is not None
        assert transcript.gene.stable_id == "gene:ENSG001"
        assert transcript.gene.name == "BRCA1"
        assert transcript.gene.biotype == "protein_coding"
        assert len(transcript.spliced_exons) == 2
        assert len(transcript.introns) == 1
        assert transcript.translation is not None
        assert transcript.translation.start == 1100
        assert transcript.translation.end == 4900
        assert transcript.translation.stable_id == "ENSP001"

    def test_get_all_features_returns_one_transcript_per_group(self, tmp_path: Path) -> None:
        """Test that multiple transcripts are built from a GFF3 file with multiple transcript groups."""

        path = _write_gff3(tmp_path / "sample.gff3", _GFF3_TWO_TRANSCRIPTS)
        source = GFF3TranscriptSource(file_path=str(path))

        transcripts = source.get_all_features()

        assert len(transcripts) == 2
        assert {t.stable_id for t in transcripts} == {"transcript:ENST001", "transcript:ENST002"}

    def test_transcript_without_cds_has_no_translation(self, tmp_path: Path) -> None:
        """Test that a transcript without a CDS feature has no translation object."""

        path = _write_gff3(tmp_path / "sample.gff3", _GFF3_TWO_TRANSCRIPTS)
        source = GFF3TranscriptSource(file_path=str(path))

        transcripts = source.get_all_features()

        assert all(t.translation is None for t in transcripts)

    def test_exon_phase_is_set_from_record(self, tmp_path: Path) -> None:
        """Test that the exon phase is correctly set from the GFF3 record."""

        path = _write_gff3(tmp_path / "sample.gff3", _GFF3_EXON_WITH_PHASE)
        source = GFF3TranscriptSource(file_path=str(path))

        transcript = source.get_all_features()[0]

        assert transcript.spliced_exons[0].phase == ExonPhase.ONE

    def test_reverse_strand_transcript_orders_exons_and_introns(self, tmp_path: Path) -> None:
        """Test that a reverse strand transcript orders exons and introns correctly."""

        path = _write_gff3(tmp_path / "sample.gff3", _GFF3_REVERSE_STRAND)
        source = GFF3TranscriptSource(file_path=str(path))

        transcript = source.get_all_features()[0]

        assert transcript.strand == Strand.REVERSE
        assert len(transcript.spliced_exons) == 3
        assert transcript.spliced_exons[0].start == 1000
        assert transcript.spliced_exons[0].end == 1100
        assert transcript.spliced_exons[1].start == 800
        assert transcript.spliced_exons[1].end == 900
        assert transcript.spliced_exons[2].start == 500
        assert transcript.spliced_exons[2].end == 600
        assert len(transcript.introns) == 2
        assert transcript.introns[0].start == 901
        assert transcript.introns[0].end == 999
        assert transcript.introns[1].start == 601
        assert transcript.introns[1].end == 799

    def test_slice_end_uses_sequence_source_length(self, tmp_path: Path) -> None:
        """Test that the slice end and seq_region_length are set from the attached sequence source."""

        path = _write_gff3(tmp_path / "sample.gff3")
        source = GFF3TranscriptSource(
            file_path=str(path), sequence_source=DummySequenceSource({"chr1": 42}, {"chr1": "ACGT"})
        )

        transcript = source.get_all_features()[0]

        assert transcript.slice.end == 42
        assert transcript.slice.seq_region_length == 42
        assert transcript.slice.seq == "ACGT"

    def test_transcript_and_exon_share_slice_instance(self, tmp_path: Path) -> None:
        """Test that the transcript and its exons share the same Slice instance."""

        path = _write_gff3(tmp_path / "sample.gff3")
        source = GFF3TranscriptSource(file_path=str(path))

        transcript = source.get_all_features()[0]

        assert transcript.slice is transcript.spliced_exons[0].slice
        assert transcript.slice is transcript.gene.slice


class TestIterAllFeatures:
    """Integration-style tests for GFF3-backed transcript building via iteration."""

    def test_iter_all_features_yields_transcript(self, tmp_path: Path) -> None:
        """Test that iter_all_features yields a transcript object from a GFF3 file."""

        path = _write_gff3(tmp_path / "sample.gff3")
        source = GFF3TranscriptSource(file_path=str(path))

        transcripts = list(source.iter_all_features())

        assert len(transcripts) == 1
        transcript = transcripts[0]
        assert transcript.stable_id == "transcript:ENST001"
        assert transcript.gene is not None
        assert transcript.gene.stable_id == "gene:ENSG001"
        assert transcript.gene.name == "BRCA1"
        assert transcript.gene.biotype == "protein_coding"
        assert len(transcript.spliced_exons) == 2
        assert len(transcript.introns) == 1
        assert transcript.translation is not None
        assert transcript.translation.start == 1100
        assert transcript.translation.end == 4900
        assert transcript.translation.stable_id == "ENSP001"

    def test_iter_all_features_yields_one_transcript_per_group(self, tmp_path: Path) -> None:
        """Test that iter_all_features yields the same transcripts as get_all_features."""

        path = _write_gff3(tmp_path / "sample.gff3", _GFF3_TWO_TRANSCRIPTS)
        source = GFF3TranscriptSource(file_path=str(path))

        iterated = list(source.iter_all_features())

        assert len(iterated) == 2
        assert all(isinstance(t, Transcript) for t in iterated)
        assert {t.stable_id for t in iterated} == {"transcript:ENST001", "transcript:ENST002"}


class TestGetFeaturesByRegion:
    """Integration-style tests for GFF3-backed transcript building via region queries."""

    def test_get_features_by_region_uses_tabix(self, tmp_path: Path, monkeypatch) -> None:
        """Test that get_features_by_region uses tabix to query a GFF3 file."""

        path = _write_query_gff3(tmp_path / "sample.gff3.gz", "\n".join(HEADER_LINES) + "\n")
        _write_query_mock_tabix(tmp_path, monkeypatch)
        source = GFF3TranscriptSource(file_path=str(path))

        transcripts = source.get_features_by_region(seq_region_name="chr1", start=1000, end=2000)

        assert len(transcripts) == 1
        assert transcripts[0].stable_id == "transcript:ENST001"
        assert len(transcripts[0].spliced_exons) == 1

    def test_get_features_by_region_empty_region_returns_empty(self, tmp_path: Path, monkeypatch) -> None:
        """Test that get_features_by_region returns an empty list for a region with no features."""

        path = _write_query_gff3(tmp_path / "sample.gff3.gz", "\n".join(HEADER_LINES) + "\n")
        _write_query_mock_tabix(tmp_path, monkeypatch)
        source = GFF3TranscriptSource(file_path=str(path))

        transcripts = source.get_features_by_region(seq_region_name="chr3", start=1, end=10)

        assert transcripts == []


class TestIterFeaturesByRegion:
    """Integration-style tests for GFF3-backed transcript building via region queries."""

    def test_iter_features_uses_tabix(self, tmp_path: Path, monkeypatch) -> None:
        """Test that iter_features_by_region uses tabix to query a GFF3 file."""

        path = _write_query_gff3(tmp_path / "sample.gff3.gz", "\n".join(HEADER_LINES) + "\n")
        _write_query_mock_tabix(tmp_path, monkeypatch)
        source = GFF3TranscriptSource(file_path=str(path))

        transcripts = list(source.iter_features_by_region(seq_region_name="chr1", start=1000, end=2000))

        assert len(transcripts) == 1
        assert transcripts[0].stable_id == "transcript:ENST001"
        assert len(transcripts[0].spliced_exons) == 1

    def test_iter_features_empty_region_returns_empty(self, tmp_path: Path, monkeypatch) -> None:
        """Test that iter_features_by_region returns an empty iterator for a region with no features."""

        path = _write_query_gff3(tmp_path / "sample.gff3.gz", "\n".join(HEADER_LINES) + "\n")
        _write_query_mock_tabix(tmp_path, monkeypatch)
        source = GFF3TranscriptSource(file_path=str(path))

        transcripts = list(source.iter_features_by_region(seq_region_name="chr3", start=1, end=10))

        assert transcripts == []


class TestContextManager:
    """Tests for the GFF3TranscriptSource context manager."""

    def test_close_closes_attached_sequence_source(self, tmp_path: Path) -> None:
        """Test that closing the GFF3TranscriptSource also closes the attached sequence source."""

        path = _write_gff3(tmp_path / "sample.gff3")
        sequence_source = DummySequenceSource({})
        source = GFF3TranscriptSource(file_path=str(path), sequence_source=sequence_source)

        source.close()

        assert sequence_source.closed is True

    def test_context_manager_closes_on_exit(self, tmp_path: Path) -> None:
        """Test that the GFF3TranscriptSource context manager closes on exit."""

        path = _write_gff3(tmp_path / "sample.gff3")

        with GFF3TranscriptSource(file_path=str(path)) as source:
            transcripts = source.get_all_features()

        assert len(transcripts) == 1
