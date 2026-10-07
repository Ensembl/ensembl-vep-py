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
import subprocess
from typing import Protocol

from attrs import define, field
import pytest

from ensembl.vep.common import ExonPhase, Strand
from ensembl.vep.parsers import GFF3Parser
from ensembl.vep.models.genic import Transcript
from ensembl.vep.sources.sequence_source import SequenceSource
from ensembl.vep.sources.file.gff3.transcript import GFF3TranscriptSource
from ensembl.vep.sources.file.gxf.transcript import _GXFTranscriptGroup
from ensembl.vep.sources.synonym import ChromosomeSynonyms, SequenceRegionResolver


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

_GFF3_CDS_WITH_PHASE = """\
chr1\tensembl\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript:ENST001;Parent=gene:ENSG001
chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tID=exon:ENSE001;Parent=transcript:ENST001
chr1\tensembl\tCDS\t1000\t1900\t.\t+\t1\tID=cds:P001;Parent=transcript:ENST001
"""

_GFF3_REVERSE_CDS_PHASES = """\
chr1\tensembl\tmRNA\t500\t1100\t.\t-\t.\tID=transcript:ENST002;Parent=gene:ENSG002
chr1\tensembl\texon\t1000\t1100\t.\t-\t.\tID=exon:ENSE011;Parent=transcript:ENST002
chr1\tensembl\texon\t800\t900\t.\t-\t.\tID=exon:ENSE012;Parent=transcript:ENST002
chr1\tensembl\texon\t500\t600\t.\t-\t.\tID=exon:ENSE013;Parent=transcript:ENST002
chr1\tensembl\tCDS\t1000\t1050\t.\t-\t2\tID=cds:P002;Parent=transcript:ENST002
chr1\tensembl\tCDS\t800\t900\t.\t-\t0\tID=cds:P002;Parent=transcript:ENST002
"""
"""A reverse-strand transcript whose CDS starts part-way through a codon (GFF phase 2) and whose
last exon is entirely 3' UTR."""

_GFF3_NO_ID_OR_NAME = """\
chr1\tensembl\tmRNA\t1000\t2000\t.\t+\t.\tParent=gene:ENSG001
"""

HEADER_LINES = ["##gff-version 3"]
CHR1_TRANSCRIPT_ROW = "chr1\tensembl\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript:ENST001;Parent=gene:ENSG001"
CHR1_EXON_ROW = "chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tID=exon:ENSE001;Parent=transcript:ENST001"
CHR2_TRANSCRIPT_ROW = "chr2\tensembl\tmRNA\t100\t200\t.\t+\t.\tID=transcript:ENST002;Parent=gene:ENSG002"

_GFF3_TWO_TRANSCRIPTS_SHARED_EXON_ONE_LINE = """\
chr1\tensembl\tgene\t1000\t5000\t.\t+\t.\tID=gene:ENSG001;Name=BRCA1
chr1\tensembl\tmRNA\t1000\t5000\t.\t+\t.\tID=transcript:ENST001;Parent=gene:ENSG001;biotype=protein_coding
chr1\tensembl\tmRNA\t1000\t5000\t.\t+\t.\tID=transcript:ENST002;Parent=gene:ENSG001;biotype=protein_coding
chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tID=exon:ENSE001;Parent=transcript:ENST001,transcript:ENST002
chr1\tensembl\tCDS\t1100\t2000\t.\t+\t0\tParent=transcript:ENST001,transcript:ENST002
chr1\tensembl\texon\t2500\t4000\t.\t+\t.\tID=exon:ENSE002;Parent=transcript:ENST002
chr1\tensembl\tCDS\t2500\t3900\t.\t+\t.\tParent=transcript:ENST002
chr1\tensembl\texon\t3000\t5000\t.\t+\t.\tID=exon:ENSE003;Parent=transcript:ENST001
chr1\tensembl\tCDS\t3000\t4900\t.\t+\t0\tParent=transcript:ENST001
"""


def _write_gff3(path: Path, data: str = _GFF3_ONE_TRANSCRIPT) -> Path:
    """Write a GFF3 fixture file."""

    path.write_text(data)
    return path


def _write_query_gff3(path: Path, data: str) -> Path:
    """Write a query-mode GFF3 fixture with a placeholder tabix index."""

    path.write_text(data)
    path.with_suffix(path.suffix + ".tbi").write_text("")
    return path


def _write_query_mock_tabix(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Write a mock tabix executable on PATH that returns deterministic GFF3 records.

    GFF3TranscriptSource's adaptor is constructed without access to the source's own
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
    """Tests for the GFF3-specific record identifier resolution."""

    def test_get_record_id_uses_id_attribute(self) -> None:
        """Test the ID attribute is used verbatim as the record identifier."""

        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        source = GFF3TranscriptSource(file_path="")

        assert source.get_record_id(records[1]) == "transcript:ENST001"

    def test_get_record_id_falls_back_to_name(self) -> None:
        """Test Name is used when ID is absent."""

        record = GFF3Parser.parse_text("chr1\tensembl\tgene\t1\t10\t.\t+\t.\tName=BRCA1\n")[0]
        source = GFF3TranscriptSource(file_path="")

        assert source.get_record_id(record) == "BRCA1"

    def test_get_record_id_synthesizes_id_when_missing(self) -> None:
        """Test a synthetic identifier is built when ID and Name are both absent."""

        record = GFF3Parser.parse_text(_GFF3_NO_ID_OR_NAME)[0]
        source = GFF3TranscriptSource(file_path="")

        assert source.get_record_id(record) == "mrna:chr1:1000-2000"


class TestGetParentId:
    """Tests for the GFF3-specific parent identifier resolution."""

    def test_get_parent_ids_returns_parent_attribute(self) -> None:
        """Test the Parent attribute is returned as the parent identifier."""

        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        source = GFF3TranscriptSource(file_path="")

        assert source.get_parent_ids(records[1]) == ("gene:ENSG001",)

    def test_get_parent_ids_returns_empty_tuple_when_absent(self) -> None:
        """Test None is returned when the record has no Parent attribute."""

        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        source = GFF3TranscriptSource(file_path="")

        assert not source.get_parent_ids(records[0])

    def test_get_parent_ids_shared_exon(self) -> None:
        """Test the Parent attribute is returned as the parent identifier."""

        records = GFF3Parser.parse_text(_GFF3_TWO_TRANSCRIPTS_SHARED_EXON_ONE_LINE)
        source = GFF3TranscriptSource(file_path="")

        assert source.get_parent_ids(records[3]) == ("transcript:ENST001", "transcript:ENST002")


class TestGetAllFeatures:
    """Integration-style tests for GFF3-backed transcript building."""

    def test_get_all_features_builds_transcript(self, tmp_path: Path) -> None:
        """Test that a transcript object is built from a GFF3 file."""

        path = _write_gff3(tmp_path / "sample.gff3")
        source = GFF3TranscriptSource(
            file_path=str(path), sequence_source=DummySequenceSource({"chr1": 9000, "chr2": 8000})
        )

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
        # Offsets from the transcript's start, not genomic: CDS 1100-4900, transcript starts at 1000.
        assert transcript.translation.start == 101
        assert transcript.translation.end == 3901
        assert transcript.translation.stable_id == "ENSP001"

    def test_get_all_features_excludes_shared_exon(self, tmp_path: Path) -> None:
        """Test that a transcript object misses a shared exon from a GFF3 file."""

        path = _write_gff3(tmp_path / "sample.gff3", _GFF3_TWO_TRANSCRIPTS_SHARED_EXON_ONE_LINE)
        source = GFF3TranscriptSource(
            file_path=str(path), sequence_source=DummySequenceSource({"chr1": 9000, "chr2": 8000})
        )

        transcripts = source.get_all_features()
        print([exon.stable_id for exon in transcripts[0].spliced_exons])
        print([exon.stable_id for exon in transcripts[1].spliced_exons])

        assert len(transcripts) == 2
        transcript = transcripts[0]
        assert len(transcript.spliced_exons) == 2

        transcript_2 = transcripts[1]
        assert len(transcript_2.spliced_exons) == 2

    def test_sequence_region_names_come_from_gff_directives(self, tmp_path: Path) -> None:
        """Expose declared query names without scanning feature records."""
        path = _write_gff3(
            tmp_path / "sample.gff3",
            "##gff-version 3\n"
            "##sequence-region 1 1 5000\n"
            "##sequence-region MT 1 100\n" + _GFF3_ONE_TRANSCRIPT.replace("chr1", "1"),
        )
        source = GFF3TranscriptSource(
            file_path=str(path), sequence_source=DummySequenceSource({"chr1": 9000, "chr2": 8000})
        )

        assert source.sequence_region_names == frozenset({"1", "MT"})

    def test_sequence_region_names_unions_index_and_header_names(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Combine the index's declared names with ``##sequence-region`` directives.

        The index is exact and always in sync with what ``seek`` accepts, but it only lists
        contigs with at least one row -- a header can still declare a contig the index omits.
        """
        path = _write_gff3(
            tmp_path / "sample.gff3",
            "##gff-version 3\n##sequence-region 1 1 5000\n" + _GFF3_ONE_TRANSCRIPT.replace("chr1", "1"),
        )
        monkeypatch.setattr(GFF3Parser, "index_sequence_names", property(lambda _self: ["1", "2"]))
        source = GFF3TranscriptSource(
            file_path=str(path), sequence_source=DummySequenceSource({"chr1": 9000, "chr2": 8000})
        )

        assert source.sequence_region_names == frozenset({"1", "2"})

    def test_sequence_region_names_includes_a_header_only_contig_the_index_omits(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Keep a header-declared contig the index has no rows for.

        A tabix index only lists contigs it has chunks for, so a contig declared in a
        ``##sequence-region`` line but with zero features (for example a short patch or MT) would
        otherwise disappear from the resolver's known names once the index takes part.
        """
        path = _write_gff3(
            tmp_path / "sample.gff3",
            "##gff-version 3\n##sequence-region 1 1 5000\n##sequence-region MT 1 100\n"
            + _GFF3_ONE_TRANSCRIPT.replace("chr1", "1"),
        )
        monkeypatch.setattr(GFF3Parser, "index_sequence_names", property(lambda _self: ["1"]))
        source = GFF3TranscriptSource(
            file_path=str(path), sequence_source=DummySequenceSource({"chr1": 9000, "chr2": 8000})
        )

        assert source.sequence_region_names == frozenset({"1", "MT"})

    def test_get_all_features_returns_one_transcript_per_group(self, tmp_path: Path) -> None:
        """Test that multiple transcripts are built from a GFF3 file with multiple transcript groups."""

        path = _write_gff3(tmp_path / "sample.gff3", _GFF3_TWO_TRANSCRIPTS)
        source = GFF3TranscriptSource(
            file_path=str(path), sequence_source=DummySequenceSource({"chr1": 9000, "chr2": 8000})
        )

        transcripts = source.get_all_features()

        assert len(transcripts) == 2
        assert {t.stable_id for t in transcripts} == {"transcript:ENST001", "transcript:ENST002"}

    def test_transcript_without_cds_has_no_translation(self, tmp_path: Path) -> None:
        """Test that a transcript without a CDS feature has no translation object."""

        path = _write_gff3(tmp_path / "sample.gff3", _GFF3_TWO_TRANSCRIPTS)
        source = GFF3TranscriptSource(
            file_path=str(path), sequence_source=DummySequenceSource({"chr1": 9000, "chr2": 8000})
        )

        transcripts = source.get_all_features()

        assert all(t.translation is None for t in transcripts)

    def test_exon_phase_is_read_from_its_cds_row(self, tmp_path: Path) -> None:
        """GFF phase 1 (skip one base) is Ensembl phase 2 (two bases of the codon came before)."""

        path = _write_gff3(tmp_path / "sample.gff3", _GFF3_CDS_WITH_PHASE)
        source = GFF3TranscriptSource(
            file_path=str(path), sequence_source=DummySequenceSource({"chr1": 9000, "chr2": 8000})
        )

        transcript = source.get_all_features()[0]

        assert transcript.spliced_exons[0].phase == ExonPhase.TWO
        assert transcript.cds_start_phase == 2

    def test_reverse_strand_phases_follow_transcript_order(self, tmp_path: Path) -> None:
        """Each exon takes its own CDS row's phase, and an exon with no CDS row has none."""

        path = _write_gff3(tmp_path / "sample.gff3", _GFF3_REVERSE_CDS_PHASES)
        source = GFF3TranscriptSource(
            file_path=str(path), sequence_source=DummySequenceSource({"chr1": 9000, "chr2": 8000})
        )

        transcript = source.get_all_features()[0]

        phases = [exon.phase for exon in transcript.spliced_exons]
        assert phases == [ExonPhase.ONE, ExonPhase.ZERO, ExonPhase.UNKNOWN]
        assert transcript.cds_start_phase == 1

    def test_a_transcript_without_cds_rows_leaves_phases_unset(self, tmp_path: Path) -> None:
        """No CDS row, nothing to read a phase from."""

        path = _write_gff3(tmp_path / "sample.gff3", _GFF3_REVERSE_STRAND)
        source = GFF3TranscriptSource(
            file_path=str(path), sequence_source=DummySequenceSource({"chr1": 9000, "chr2": 8000})
        )

        transcript = source.get_all_features()[0]

        assert all(exon.phase is None for exon in transcript.spliced_exons)
        assert transcript.cds_start_phase == 0

    def test_reverse_strand_transcript_orders_exons_and_introns(self, tmp_path: Path) -> None:
        """Test that a reverse strand transcript orders exons and introns correctly."""

        path = _write_gff3(tmp_path / "sample.gff3", _GFF3_REVERSE_STRAND)
        source = GFF3TranscriptSource(
            file_path=str(path), sequence_source=DummySequenceSource({"chr1": 9000, "chr2": 8000})
        )

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

        assert transcript.slice is not None
        assert transcript.slice.end == 42
        assert transcript.slice.seq_region_length == 42
        assert transcript.slice.seq == "ACGT"

    def test_transcript_and_exon_share_slice_instance(self, tmp_path: Path) -> None:
        """Test that the transcript and its exons share the same Slice instance."""

        path = _write_gff3(tmp_path / "sample.gff3")
        source = GFF3TranscriptSource(
            file_path=str(path), sequence_source=DummySequenceSource({"chr1": 9000, "chr2": 8000})
        )

        transcript = source.get_all_features()[0]

        assert transcript.slice is transcript.spliced_exons[0].slice
        assert transcript.gene is not None
        assert transcript.slice is transcript.gene.slice


class TestIterAllFeatures:
    """Integration-style tests for GFF3-backed transcript building via iteration."""

    def test_iter_all_features_yields_transcript(self, tmp_path: Path) -> None:
        """Test that iter_all_features yields a transcript object from a GFF3 file."""

        path = _write_gff3(tmp_path / "sample.gff3")
        source = GFF3TranscriptSource(
            file_path=str(path), sequence_source=DummySequenceSource({"chr1": 9000, "chr2": 8000})
        )

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
        # Offsets from the transcript's start, not genomic: CDS 1100-4900, transcript starts at 1000.
        assert transcript.translation.start == 101
        assert transcript.translation.end == 3901
        assert transcript.translation.stable_id == "ENSP001"

    def test_iter_all_features_yields_one_transcript_per_group(self, tmp_path: Path) -> None:
        """Test that iter_all_features yields the same transcripts as get_all_features."""

        path = _write_gff3(tmp_path / "sample.gff3", _GFF3_TWO_TRANSCRIPTS)
        source = GFF3TranscriptSource(
            file_path=str(path), sequence_source=DummySequenceSource({"chr1": 9000, "chr2": 8000})
        )

        iterated = list(source.iter_all_features())

        assert len(iterated) == 2
        assert all(isinstance(t, Transcript) for t in iterated)
        assert {t.stable_id for t in iterated} == {"transcript:ENST001", "transcript:ENST002"}

    def test_gene_grouping_allows_transcripts_before_their_children(self, tmp_path: Path) -> None:
        """Assemble a valid gene block whose transcript declarations precede all exons."""

        data = (
            "chr1\tensembl\tgene\t100\t500\t.\t+\t.\tID=gene:G1\n"
            "chr1\tensembl\tmRNA\t100\t300\t.\t+\t.\tID=transcript:T1;Parent=gene:G1\n"
            "chr1\tensembl\tmRNA\t200\t500\t.\t+\t.\tID=transcript:T2;Parent=gene:G1\n"
            "chr1\tensembl\texon\t100\t300\t.\t+\t.\tID=exon:E1;Parent=transcript:T1\n"
            "chr1\tensembl\texon\t200\t500\t.\t+\t.\tID=exon:E2;Parent=transcript:T2\n"
        )
        path = _write_gff3(tmp_path / "gene-grouped.gff3", data)
        source = GFF3TranscriptSource(
            file_path=str(path), sequence_source=DummySequenceSource({"chr1": 9000, "chr2": 8000})
        )

        transcripts = list(source.iter_all_features())

        assert [transcript.stable_id for transcript in transcripts] == [
            "transcript:T1",
            "transcript:T2",
        ]
        assert [transcript.spliced_exons[0].stable_id for transcript in transcripts] == [
            "exon:E1",
            "exon:E2",
        ]

    def test_interval_stream_filters_before_building_transcript_graphs(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """An unindexed scan does not construct models outside the payload intervals."""

        second_gene = (
            "chr2\tensembl\tgene\t100\t200\t.\t+\t.\tID=gene:ENSG002\n"
            "chr2\tensembl\tmRNA\t100\t200\t.\t+\t.\t"
            "ID=transcript:ENST002;Parent=gene:ENSG002\n"
            "chr2\tensembl\texon\t100\t200\t.\t+\t.\t"
            "ID=exon:ENSE002;Parent=transcript:ENST002\n"
        )
        path = _write_gff3(tmp_path / "sample.gff3", _GFF3_ONE_TRANSCRIPT + second_gene)
        source = GFF3TranscriptSource(
            file_path=str(path), sequence_source=DummySequenceSource({"chr1": 9000, "chr2": 8000})
        )
        built: list[str | None] = []
        # pylint: disable-next=protected-access
        build_transcript = GFF3TranscriptSource._build_transcript_from_record_group

        def record_build(source: GFF3TranscriptSource, group: _GXFTranscriptGroup) -> Transcript:
            built.append(source.get_record_id(group.transcript_record))
            return build_transcript(source, group)

        monkeypatch.setattr(
            GFF3TranscriptSource,
            "_build_transcript_from_record_group",
            record_build,
        )

        transcripts = list(source.iter_features_by_regions({"chr2": [(150, 150)]}))

        assert [transcript.stable_id for transcript in transcripts] == ["transcript:ENST002"]
        assert built == ["transcript:ENST002"]


class TestGetFeaturesByRegion:
    """Integration-style tests for GFF3-backed transcript building via region queries."""

    def test_get_features_by_region_uses_tabix(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        """Test that get_features_by_region uses tabix to query a GFF3 file."""

        path = _write_query_gff3(tmp_path / "sample.gff3.gz", "\n".join(HEADER_LINES) + "\n")
        _write_query_mock_tabix(tmp_path, monkeypatch)
        source = GFF3TranscriptSource(
            file_path=str(path), sequence_source=DummySequenceSource({"chr1": 9000, "chr2": 8000})
        )

        transcripts = source.get_features_by_region(seq_region_name="chr1", start=1000, end=2000)

        assert len(transcripts) == 1
        assert transcripts[0].stable_id == "transcript:ENST001"
        assert len(transcripts[0].spliced_exons) == 1

    def test_get_features_by_region_empty_region_returns_empty(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Test that get_features_by_region returns an empty list for a region with no features."""

        path = _write_query_gff3(tmp_path / "sample.gff3.gz", "\n".join(HEADER_LINES) + "\n")
        _write_query_mock_tabix(tmp_path, monkeypatch)
        source = GFF3TranscriptSource(
            file_path=str(path), sequence_source=DummySequenceSource({"chr1": 9000, "chr2": 8000})
        )

        transcripts = source.get_features_by_region(seq_region_name="chr3", start=1, end=10)

        assert transcripts == []


class TestQueryParserReuse:
    """Tests for the persistent region-query parser."""

    def test_region_queries_open_the_file_once(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        """Test repeated complete-feature queries reuse one open query parser."""

        path = _write_query_gff3(tmp_path / "sample.gff3.gz", "\n".join(HEADER_LINES) + "\n")
        _write_query_mock_tabix(tmp_path, monkeypatch)
        source = GFF3TranscriptSource(
            file_path=str(path), sequence_source=DummySequenceSource({"chr1": 9000, "chr2": 8000})
        )
        opened: list[str] = []
        original_open = GFF3Parser.open

        def counting_open(parser: GFF3Parser, file_path: str, mode: str = "auto") -> GFF3Parser:
            opened.append(mode)
            return original_open(parser, file_path, mode)

        monkeypatch.setattr(GFF3Parser, "open", counting_open)

        first = source.get_features_by_region(seq_region_name="chr1", start=1000, end=2000)
        second = source.get_features_by_region(seq_region_name="chr1", start=1000, end=2000)

        assert [transcript.stable_id for transcript in first] == ["transcript:ENST001"]
        assert [transcript.stable_id for transcript in second] == ["transcript:ENST001"]
        assert opened == ["query"]

    def test_close_releases_query_parser(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        """Test closing the source drops its query parser so the next query opens a new one."""

        path = _write_query_gff3(tmp_path / "sample.gff3.gz", "\n".join(HEADER_LINES) + "\n")
        _write_query_mock_tabix(tmp_path, monkeypatch)
        source = GFF3TranscriptSource(
            file_path=str(path), sequence_source=DummySequenceSource({"chr1": 9000, "chr2": 8000})
        )
        source.get_features_by_region(seq_region_name="chr1", start=1000, end=2000)

        source.close()

        assert source._query_parser is None  # pylint: disable=protected-access
        assert source.get_features_by_region(seq_region_name="chr1", start=1000, end=2000)


class TestIterFeaturesByRegion:
    """One lazy region is served by the same sequential pass as a batch of regions."""

    def test_iter_features_by_region_streams_complete_transcripts(self, tmp_path: Path) -> None:
        """A single region yields whole transcripts, exons included, and drops the other region's."""

        second_gene = (
            "chr2\tensembl\tgene\t100\t200\t.\t+\t.\tID=gene:ENSG002\n"
            "chr2\tensembl\tmRNA\t100\t200\t.\t+\t.\t"
            "ID=transcript:ENST002;Parent=gene:ENSG002\n"
            "chr2\tensembl\texon\t100\t200\t.\t+\t.\t"
            "ID=exon:ENSE002;Parent=transcript:ENST002\n"
        )
        path = _write_gff3(tmp_path / "sample.gff3", _GFF3_ONE_TRANSCRIPT + second_gene)
        source = GFF3TranscriptSource(
            file_path=str(path), sequence_source=DummySequenceSource({"chr1": 9000, "chr2": 8000})
        )

        transcripts = list(source.iter_features_by_region(seq_region_name="chr2", start=150, end=150))

        assert [transcript.stable_id for transcript in transcripts] == ["transcript:ENST002"]
        assert [exon.stable_id for exon in transcripts[0].spliced_exons] == ["exon:ENSE002"]

    def test_iter_features_by_region_without_coordinates_takes_the_whole_region(self, tmp_path: Path) -> None:
        """A region with no bounds means every transcript on it, as ``get_features_by_region`` does."""

        path = _write_gff3(tmp_path / "sample.gff3")
        source = GFF3TranscriptSource(
            file_path=str(path), sequence_source=DummySequenceSource({"chr1": 9000, "chr2": 8000})
        )

        transcripts = list(source.iter_features_by_region("chr1"))

        assert [transcript.stable_id for transcript in transcripts] == ["transcript:ENST001"]

    def test_iter_features_by_region_on_another_region_yields_nothing(self, tmp_path: Path) -> None:
        """A region the file does not carry yields nothing rather than everything."""

        path = _write_gff3(tmp_path / "sample.gff3")
        source = GFF3TranscriptSource(
            file_path=str(path), sequence_source=DummySequenceSource({"chr1": 9000, "chr2": 8000})
        )

        assert not list(source.iter_features_by_region(seq_region_name="chr9", start=1, end=100))

    def test_iter_features_by_region_returns_empty_list_on_bad_region(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """Test that iter_features_by_region returns an empty list and logs an error message if a
        bad region is supplied."""
        path = _write_gff3(tmp_path / "sample.gff3")
        source = GFF3TranscriptSource(
            file_path=str(path), sequence_source=DummySequenceSource({"chr1": 9000, "chr2": 8000})
        )

        with caplog.at_level("DEBUG"):
            transcript = list(source.iter_features_by_region(seq_region_name="chr1", start=-1, end=100))
            assert not transcript
            assert "Cannot parse region" in caplog.text

    def test_iter_features_by_regions_with_resolver(self, tmp_path: Path) -> None:
        """Test using a reference_region_resolver to convert contig to chr names."""
        _GFF3_CONTIG_NAMES = """\
NC_000001.11\tensembl\tgene\t1000\t5000\t.\t+\t.\tID=gene:ENSG001;Name=BRCA1;biotype=protein_coding
NC_000001.11\tensembl\tmRNA\t1000\t5000\t.\t+\t.\tID=transcript:ENST001;Parent=gene:ENSG001
NC_000001.11\tensembl\texon\t1000\t5000\t.\t+\t.\tID=exon:ENSE001;Parent=transcript:ENST001
NC_000001.11\tensembl\tCDS\t1100\t4900\t.\t+\t0\tID=cds:ENSP001;Parent=transcript:ENST001;protein_id=ENSP001
"""
        path = _write_gff3(tmp_path / "sample.gff3", data=_GFF3_CONTIG_NAMES)

        chr_synonyms = ChromosomeSynonyms(
            by_name={
                "chr1": frozenset({"NC_000001.11", "CM000663.2"}),
                **{synonym: frozenset({"chr1"}) for synonym in ("NC_000001.11", "CM000663.2")},
            }
        )
        chr_resolver = SequenceRegionResolver(
            source_names=frozenset({"chr1"}),
            synonyms=chr_synonyms,
        )
        source = GFF3TranscriptSource(
            file_path=str(path),
            sequence_source=DummySequenceSource({"NC_000001.11": 9000, "NC_000002.12": 8000}),
            reference_region_resolver=chr_resolver,
        )

        (transcript,) = list(source.iter_features_by_region(seq_region_name="chr1", start=2000, end=4000))
        assert transcript.stable_id == "transcript:ENST001"


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

        with GFF3TranscriptSource(
            file_path=str(path), sequence_source=DummySequenceSource({"chr1": 9000, "chr2": 8000})
        ) as source:
            transcripts = source.get_all_features()

        assert len(transcripts) == 1


class HtslibToolchain(Protocol):
    """Protocol for HTSlib toolchain fixtures."""

    bgzip: str
    tabix: str


def _create_tabix_gff3_cli(path: Path, htslib_toolchain: HtslibToolchain) -> Path:
    """Create and tabix-index a real GFF3 fixture with a real HTSlib toolchain.

    Only ``chr1`` is declared -- a query for any other contig exercises the unindexed-contig
    path against a genuine, in-process (noodles) tabix index rather than a mock.
    """

    unzip_path = path.with_suffix(path.suffix.replace(".gz", ""))
    unzip_path.write_text(_GFF3_ONE_TRANSCRIPT)

    subprocess.run([htslib_toolchain.bgzip, "-f", str(unzip_path)], check=True)
    subprocess.run([htslib_toolchain.tabix, "-f", "-p", "gff", str(path)], check=True)

    return path


@pytest.mark.integration
class TestUnindexedContigQuery:
    """A variant on a contig absent from the GFF index must not raise."""

    def test_get_features_by_region_on_unindexed_contig_returns_empty(
        self, tmp_path: Path, htslib_toolchain: HtslibToolchain
    ) -> None:
        """Test querying a contig the index does not declare returns no transcripts."""

        path = _create_tabix_gff3_cli(tmp_path / "sample.gff3.gz", htslib_toolchain)
        source = GFF3TranscriptSource(
            file_path=str(path), sequence_source=DummySequenceSource({"chr1": 9000, "chr2": 8000})
        )

        transcripts = source.get_features_by_region(seq_region_name="chrZ", start=1, end=100)

        assert transcripts == []

        # The source stays usable for a subsequent, known contig.
        known = source.get_features_by_region(seq_region_name="chr1", start=1000, end=5000)
        assert [transcript.stable_id for transcript in known] == ["transcript:ENST001"]

    def test_sequence_region_names_uses_the_real_index(
        self, tmp_path: Path, htslib_toolchain: HtslibToolchain
    ) -> None:
        """Test the resolver source names come from the real index, not file headers."""

        path = _create_tabix_gff3_cli(tmp_path / "sample.gff3.gz", htslib_toolchain)
        source = GFF3TranscriptSource(
            file_path=str(path), sequence_source=DummySequenceSource({"chr1": 9000, "chr2": 8000})
        )

        assert source.sequence_region_names == frozenset({"chr1"})
