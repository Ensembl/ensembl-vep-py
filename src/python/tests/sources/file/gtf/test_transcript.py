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
"""Unit tests for the GTF transcript source."""

import os
from pathlib import Path

from attrs import define, field
import pytest

from ensembl.vep.common import ExonPhase, Strand
from ensembl.vep.parsers import GTFParser
from ensembl.vep.models.genic import Transcript
from ensembl.vep.sources.sequence_source import SequenceSource
from ensembl.vep.sources.file.gtf.transcript import GTFTranscriptSource


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


_GTF_ONE_TRANSCRIPT = """\
chr1\tensembl\tgene\t1000\t5000\t.\t+\t.\tgene_id "gene:1"; gene_name "BRCA1"; gene_biotype "protein_coding";
chr1\tensembl\ttranscript\t1000\t5000\t.\t+\t.\tgene_id "gene:1"; transcript_id "transcript:1"; transcript_name "BRCA1-001";
chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tgene_id "gene:1"; transcript_id "transcript:1"; exon_number "1";
chr1\tensembl\tCDS\t1100\t2000\t.\t+\t0\tgene_id "gene:1"; transcript_id "transcript:1"; protein_id "ENSP0001";
chr1\tensembl\texon\t3000\t5000\t.\t+\t.\tgene_id "gene:1"; transcript_id "transcript:1"; exon_number "2";
chr1\tensembl\tCDS\t3000\t4900\t.\t+\t0\tgene_id "gene:1"; transcript_id "transcript:1"; protein_id "ENSP0001";
"""

_GTF_TWO_TRANSCRIPTS = """\
chr1\tensembl\ttranscript\t1000\t5000\t.\t+\t.\tgene_id "gene:1"; transcript_id "transcript:1";
chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tgene_id "gene:1"; transcript_id "transcript:1";
chr1\tensembl\ttranscript\t2000\t8000\t.\t-\t.\tgene_id "gene:2"; transcript_id "transcript:2";
chr1\tensembl\texon\t2000\t3000\t.\t-\t.\tgene_id "gene:2"; transcript_id "transcript:2";
chr1\tensembl\texon\t6000\t8000\t.\t-\t.\tgene_id "gene:2"; transcript_id "transcript:2";
"""

_GTF_REVERSE_STRAND = """\
chr2\tensembl\ttranscript\t500\t1100\t.\t-\t.\tgene_id "gene:3"; transcript_id "transcript:3";
chr2\tensembl\texon\t500\t600\t.\t-\t.\tgene_id "gene:3"; transcript_id "transcript:3";
chr2\tensembl\texon\t800\t900\t.\t-\t.\tgene_id "gene:3"; transcript_id "transcript:3";
chr2\tensembl\texon\t1000\t1100\t.\t-\t.\tgene_id "gene:3"; transcript_id "transcript:3";
"""

_GTF_CDS_WITH_PHASE = """\
chr1\tensembl\ttranscript\t1000\t2000\t.\t+\t.\tgene_id "gene:1"; transcript_id "transcript:1";
chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tgene_id "gene:1"; transcript_id "transcript:1";
chr1\tensembl\tCDS\t1000\t1900\t.\t+\t1\tgene_id "gene:1"; transcript_id "transcript:1";
"""

_GTF_NO_TRANSCRIPT_ID = """\
chr1\tensembl\ttranscript\t1000\t2000\t.\t+\t.\tgene_id "gene:1";
"""

HEADER_LINES = ["#!genome-version test"]
CHR1_TRANSCRIPT_ROW = (
    'chr1\tensembl\ttranscript\t1000\t2000\t.\t+\t.\tgene_id "gene:1"; transcript_id "transcript:1";'
)
CHR1_EXON_ROW = 'chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tgene_id "gene:1"; transcript_id "transcript:1";'
CHR2_TRANSCRIPT_ROW = (
    'chr2\tensembl\ttranscript\t100\t200\t.\t+\t.\tgene_id "gene:2"; transcript_id "transcript:2";'
)


def _write_gtf(path: Path, data: str = _GTF_ONE_TRANSCRIPT) -> Path:
    """Write a GTF fixture file."""

    path.write_text(data)
    return path


def _write_query_gtf(path: Path, data: str) -> Path:
    """Write a query-mode GTF fixture with a placeholder tabix index."""

    path.write_text(data)
    path.with_suffix(path.suffix + ".tbi").write_text("")
    return path


def _write_query_mock_tabix(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Write a mock tabix executable on PATH that returns deterministic GTF records.

    GTFTranscriptSource's adaptor is constructed without access to the source's own
    tabix_path field, so the mock is put on PATH under the literal name "tabix" rather
    than passed via the tabix_path constructor argument.
    """

    tabix = tmp_path / "tabix"
    tabix.write_text(f"""#!/bin/sh
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
""")
    tabix.chmod(0o755)
    monkeypatch.setenv("PATH", f"{tmp_path}{os.pathsep}{os.environ['PATH']}")
    return tabix


class TestGetRecordId:
    """Tests for the GTF-specific record identifier resolution."""

    def test_get_record_id_uses_feature_specific_id_attribute(self) -> None:
        """Test the <feature>_id attribute is used as the record identifier."""

        records = GTFParser.parse_text(_GTF_ONE_TRANSCRIPT)
        source = GTFTranscriptSource(file_path="")

        assert source.get_record_id(records[1]) == "transcript:1"

    def test_get_record_id_synthesizes_id_when_missing(self) -> None:
        """Test a synthetic identifier is built when the feature has no matching *_id attribute."""

        record = GTFParser.parse_text(_GTF_NO_TRANSCRIPT_ID)[0]
        source = GTFTranscriptSource(file_path="")

        assert source.get_record_id(record) == "transcript:chr1:1000-2000"


class TestGetParentId:
    """Tests for the GTF-specific parent identifier resolution."""

    def test_get_parent_id_returns_parent_feature_id(self) -> None:
        """Test the parent feature type's id attribute is returned."""

        records = GTFParser.parse_text(_GTF_ONE_TRANSCRIPT)
        source = GTFTranscriptSource(file_path="")

        assert source.get_parent_id(records[1]) == "gene:1"

    def test_get_parent_id_returns_none_for_gene_records(self) -> None:
        """Test None is returned for a feature type with no configured parent."""

        records = GTFParser.parse_text(_GTF_ONE_TRANSCRIPT)
        source = GTFTranscriptSource(file_path="")

        assert source.get_parent_id(records[0]) is None


class TestGetAllFeatures:
    """Integration-style tests for GTF-backed transcript building."""

    def test_get_all_features_builds_transcript(self, tmp_path: Path) -> None:
        """Test that a transcript object is built from a GTF file."""

        path = _write_gtf(tmp_path / "sample.gtf")
        source = GTFTranscriptSource(file_path=str(path))

        transcripts = source.get_all_features()

        assert len(transcripts) == 1
        transcript = transcripts[0]
        assert isinstance(transcript, Transcript)
        assert transcript.stable_id == "transcript:1"
        assert transcript.gene is not None
        assert transcript.gene.stable_id == "gene:1"
        assert transcript.gene.name == "BRCA1"
        assert transcript.gene.biotype == "protein_coding"
        assert len(transcript.spliced_exons) == 2
        assert len(transcript.introns) == 1
        assert transcript.translation is not None
        # Offsets from the transcript's start, not genomic: CDS 1100-4900, transcript starts at 1000.
        assert transcript.translation.start == 101
        assert transcript.translation.end == 3901
        assert transcript.translation.stable_id == "ENSP0001"

    def test_get_all_features_returns_one_transcript_per_group(self, tmp_path: Path) -> None:
        """Test that multiple transcripts are built from a GTF file with multiple transcript groups."""

        path = _write_gtf(tmp_path / "sample.gtf", _GTF_TWO_TRANSCRIPTS)
        source = GTFTranscriptSource(file_path=str(path))

        transcripts = source.get_all_features()

        assert len(transcripts) == 2
        assert {t.stable_id for t in transcripts} == {"transcript:1", "transcript:2"}

    def test_transcript_without_cds_has_no_translation(self, tmp_path: Path) -> None:
        """Test that a transcript without a CDS feature has no translation object."""

        path = _write_gtf(tmp_path / "sample.gtf", _GTF_TWO_TRANSCRIPTS)
        source = GTFTranscriptSource(file_path=str(path))

        transcripts = source.get_all_features()

        assert all(t.translation is None for t in transcripts)

    def test_exon_phase_is_read_from_its_cds_row(self, tmp_path: Path) -> None:
        """GFF phase 1 (skip one base) is Ensembl phase 2 (two bases of the codon came before)."""

        path = _write_gtf(tmp_path / "sample.gtf", _GTF_CDS_WITH_PHASE)
        source = GTFTranscriptSource(file_path=str(path))

        transcript = source.get_all_features()[0]

        assert transcript.spliced_exons[0].phase == ExonPhase.TWO

    def test_reverse_strand_transcript_orders_exons_and_introns(self, tmp_path: Path) -> None:
        """Test that a reverse strand transcript orders exons and introns correctly."""

        path = _write_gtf(tmp_path / "sample.gtf", _GTF_REVERSE_STRAND)
        source = GTFTranscriptSource(file_path=str(path))

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

        path = _write_gtf(tmp_path / "sample.gtf")
        source = GTFTranscriptSource(
            file_path=str(path), sequence_source=DummySequenceSource({"chr1": 42}, {"chr1": "ACGT"})
        )

        transcript = source.get_all_features()[0]

        assert transcript.slice is not None
        assert transcript.slice.end == 42
        assert transcript.slice.seq_region_length == 42
        assert transcript.slice.seq == "ACGT"

    def test_transcript_and_exon_share_slice_instance(self, tmp_path: Path) -> None:
        """Test that the transcript and its exons share the same Slice instance."""

        path = _write_gtf(tmp_path / "sample.gtf")
        source = GTFTranscriptSource(file_path=str(path))

        transcript = source.get_all_features()[0]

        assert transcript.slice is transcript.spliced_exons[0].slice
        assert transcript.gene is not None
        assert transcript.slice is transcript.gene.slice


class TestIterAllFeatures:
    """Integration-style tests for GTF-backed transcript building via iteration."""

    def test_iter_all_features_yields_transcript(self, tmp_path: Path) -> None:
        """Test that iter_all_features yields a transcript object from a GTF file."""

        path = _write_gtf(tmp_path / "sample.gtf")
        source = GTFTranscriptSource(file_path=str(path))

        transcripts = list(source.iter_all_features())

        assert len(transcripts) == 1
        transcript = transcripts[0]
        assert transcript.stable_id == "transcript:1"
        assert transcript.gene is not None
        assert transcript.gene.stable_id == "gene:1"
        assert transcript.gene.name == "BRCA1"
        assert transcript.gene.biotype == "protein_coding"
        assert len(transcript.spliced_exons) == 2
        assert len(transcript.introns) == 1
        assert transcript.translation is not None
        # Offsets from the transcript's start, not genomic: CDS 1100-4900, transcript starts at 1000.
        assert transcript.translation.start == 101
        assert transcript.translation.end == 3901
        assert transcript.translation.stable_id == "ENSP0001"

    def test_iter_all_features_yields_one_transcript_per_group(self, tmp_path: Path) -> None:
        """Test that iter_all_features yields the same transcripts as get_all_features."""

        path = _write_gtf(tmp_path / "sample.gtf", _GTF_TWO_TRANSCRIPTS)
        source = GTFTranscriptSource(file_path=str(path))

        iterated = list(source.iter_all_features())

        assert len(iterated) == 2
        assert all(isinstance(t, Transcript) for t in iterated)
        assert {t.stable_id for t in iterated} == {"transcript:1", "transcript:2"}


class TestGetFeaturesByRegion:
    """Integration-style tests for GTF-backed transcript building via region queries."""

    def test_get_features_by_region_uses_tabix(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        """Test that get_features_by_region uses tabix to query a GTF file."""

        path = _write_query_gtf(tmp_path / "sample.gtf.gz", "\n".join(HEADER_LINES) + "\n")
        _write_query_mock_tabix(tmp_path, monkeypatch)
        source = GTFTranscriptSource(file_path=str(path))

        transcripts = source.get_features_by_region(seq_region_name="chr1", start=1000, end=2000)

        assert len(transcripts) == 1
        assert transcripts[0].stable_id == "transcript:1"
        assert len(transcripts[0].spliced_exons) == 1

    def test_get_features_by_region_empty_region_returns_empty(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Test that get_features_by_region returns an empty list for a region with no features."""

        path = _write_query_gtf(tmp_path / "sample.gtf.gz", "\n".join(HEADER_LINES) + "\n")
        _write_query_mock_tabix(tmp_path, monkeypatch)
        source = GTFTranscriptSource(file_path=str(path))

        transcripts = source.get_features_by_region(seq_region_name="chr3", start=1, end=10)

        assert transcripts == []


class TestIterFeaturesByRegion:
    """One lazy region is served by the same sequential pass as a batch of regions."""

    def test_iter_features_by_region_streams_complete_transcripts(self, tmp_path: Path) -> None:
        """A single region yields whole transcripts, both exons included, not just the overlapping one."""

        path = _write_gtf(tmp_path / "sample.gtf", _GTF_TWO_TRANSCRIPTS)
        source = GTFTranscriptSource(file_path=str(path))

        transcripts = list(source.iter_features_by_region(seq_region_name="chr1", start=7000, end=7000))

        assert [transcript.stable_id for transcript in transcripts] == ["transcript:2"]
        assert [(exon.start, exon.end) for exon in transcripts[0].spliced_exons] == [
            (6000, 8000),
            (2000, 3000),
        ]

    def test_iter_features_by_region_on_another_region_yields_nothing(self, tmp_path: Path) -> None:
        """A region the file does not carry yields nothing rather than everything."""

        path = _write_gtf(tmp_path / "sample.gtf", _GTF_TWO_TRANSCRIPTS)
        source = GTFTranscriptSource(file_path=str(path))

        assert not list(source.iter_features_by_region(seq_region_name="chr9", start=1, end=100))


class TestContextManager:
    """Tests for the GTFTranscriptSource context manager."""

    def test_close_closes_attached_sequence_source(self, tmp_path: Path) -> None:
        """Test that closing the GTFTranscriptSource also closes the attached sequence source."""

        path = _write_gtf(tmp_path / "sample.gtf")
        sequence_source = DummySequenceSource({})
        source = GTFTranscriptSource(file_path=str(path), sequence_source=sequence_source)

        source.close()

        assert sequence_source.closed is True

    def test_context_manager_closes_on_exit(self, tmp_path: Path) -> None:
        """Test that the GTFTranscriptSource context manager closes on exit."""

        path = _write_gtf(tmp_path / "sample.gtf")

        with GTFTranscriptSource(file_path=str(path)) as source:
            transcripts = source.get_all_features()

        assert len(transcripts) == 1
