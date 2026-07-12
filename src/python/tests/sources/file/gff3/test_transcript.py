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
"""Unit tests for GFF3TranscriptSource."""

import pytest

from ensembl.vep.common import Strand
from ensembl.vep.parsers import GFF3Parser
from ensembl.vep.models.genic import Transcript, Exon, Intron
from ensembl.vep.sources.file.gff3.transcript_source import (
    GFF3TranscriptSource,
    _group_gff3_records,
    _build_transcript_from_record_group,
    _build_translation,
    _GFF3TranscriptGroup,
)

_GFF3_ONE_TRANSCRIPT = """\
chr1\tensembl\tgene\t1000\t5000\t.\t+\t.\tID=gene:ENSG001;Name=BRCA1
chr1\tensembl\tmRNA\t1000\t5000\t.\t+\t.\tID=transcript:ENST001;Parent=gene:ENSG001;biotype=protein_coding
chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tID=exon:ENSE001;Parent=transcript:ENST001
chr1\tensembl\tCDS\t1100\t2000\t.\t+\t0\tParent=transcript:ENST001
chr1\tensembl\texon\t3000\t5000\t.\t+\t.\tID=exon:ENSE002;Parent=transcript:ENST001
chr1\tensembl\tCDS\t3000\t4900\t.\t+\t0\tParent=transcript:ENST001
"""

_GFF3_TWO_TRANSCRIPTS = """\
chr1\tensembl\tmRNA\t1000\t5000\t.\t+\t.\tID=transcript:ENST001;Parent=gene:ENSG001
chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tID=exon:ENSE001;Parent=transcript:ENST001
chr1\tensembl\tmRNA\t2000\t8000\t.\t-\t.\tID=transcript:ENST002;Parent=gene:ENSG002
chr1\tensembl\texon\t2000\t3000\t.\t-\t.\tID=exon:ENSE002;Parent=transcript:ENST002
chr1\tensembl\texon\t6000\t8000\t.\t-\t.\tID=exon:ENSE003;Parent=transcript:ENST002
"""

_GFF3_REVERSE_STRAND = """\
chr2\tensembl\tmRNA\t500\t900\t.\t-\t.\tID=transcript:ENST003;Parent=gene:ENSG003
chr2\tensembl\texon\t500\t600\t.\t-\t.\tID=exon:ENSE010;Parent=transcript:ENST003
chr2\tensembl\texon\t800\t900\t.\t-\t.\tID=exon:ENSE011;Parent=transcript:ENST003
"""


class TestGroupGFF3Transcripts:
    """Tests for the private _group_gff3_records helper."""

    def test_one_transcript_one_group(self) -> None:
        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        groups = _group_gff3_records(records)
        assert len(groups) == 1

    def test_exons_collected_in_group(self) -> None:
        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        groups = _group_gff3_records(records)
        assert len(groups[0].exon_records) == 2

    def test_cds_collected_in_group(self) -> None:
        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        groups = _group_gff3_records(records)
        assert len(groups[0].cds_records) == 2

    def test_gene_records_not_treated_as_transcripts(self) -> None:
        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        groups = _group_gff3_records(records)
        assert groups[0].transcript_record.feature_type == "mRNA"

    def test_gene_record_captured_on_group(self) -> None:
        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        groups = _group_gff3_records(records)
        assert groups[0].gene_record is not None
        assert groups[0].gene_record.feature_type == "gene"

    def test_gene_record_none_when_absent(self) -> None:
        records = GFF3Parser.parse_text(_GFF3_TWO_TRANSCRIPTS)
        groups = _group_gff3_records(records)
        assert all(g.gene_record is None for g in groups)

    def test_two_transcripts_two_groups(self) -> None:
        records = GFF3Parser.parse_text(_GFF3_TWO_TRANSCRIPTS)
        groups = _group_gff3_records(records)
        assert len(groups) == 2

    def test_exons_assigned_to_correct_parent(self) -> None:
        records = GFF3Parser.parse_text(_GFF3_TWO_TRANSCRIPTS)
        groups = _group_gff3_records(records)
        by_id = {g.transcript_record.attributes["ID"][0]: g for g in groups}
        assert len(by_id["transcript:ENST001"].exon_records) == 1
        assert len(by_id["transcript:ENST002"].exon_records) == 2

    def test_empty_input_returns_empty(self) -> None:
        groups = _group_gff3_records([])
        assert groups == []


class TestBuildTranscript:
    """Tests for the private _build_transcript helper."""

    def _make_group(self, gff_text: str) -> _GFF3TranscriptGroup:
        records = GFF3Parser.parse_text(gff_text)
        groups = _group_gff3_records(records)
        assert len(groups) == 1
        return groups[0]

    def test_stable_id_stripped_of_prefix(self) -> None:
        group = self._make_group(_GFF3_ONE_TRANSCRIPT)
        transcript = _build_transcript_from_record_group(group, {})
        assert transcript.stable_id == "ENST001"

    def test_coordinates(self) -> None:
        group = self._make_group(_GFF3_ONE_TRANSCRIPT)
        transcript = _build_transcript_from_record_group(group, {})
        assert transcript.start == 1000
        assert transcript.end == 5000

    def test_forward_strand(self) -> None:
        group = self._make_group(_GFF3_ONE_TRANSCRIPT)
        transcript = _build_transcript_from_record_group(group, {})
        assert transcript.strand == Strand.FORWARD

    def test_reverse_strand(self) -> None:
        group = self._make_group(_GFF3_REVERSE_STRAND)
        transcript = _build_transcript_from_record_group(group, {})
        assert transcript.strand == Strand.REVERSE

    def test_seq_region_name(self) -> None:
        group = self._make_group(_GFF3_ONE_TRANSCRIPT)
        transcript = _build_transcript_from_record_group(group, {})
        assert transcript.slice is not None
        assert transcript.slice.seq_region_name == "chr1"

    def test_spliced_exons_count(self) -> None:
        group = self._make_group(_GFF3_ONE_TRANSCRIPT)
        transcript = _build_transcript_from_record_group(group, {})
        assert len(transcript.spliced_exons) == 2

    def test_spliced_exons_sorted_by_start(self) -> None:
        group = self._make_group(_GFF3_ONE_TRANSCRIPT)
        transcript = _build_transcript_from_record_group(group, {})
        starts = [e.start for e in transcript.spliced_exons]
        assert starts == sorted(starts)

    def test_exon_stable_ids(self) -> None:
        group = self._make_group(_GFF3_ONE_TRANSCRIPT)
        transcript = _build_transcript_from_record_group(group, {})
        stable_ids = {e.stable_id for e in transcript.spliced_exons}
        assert stable_ids == {"ENSE001", "ENSE002"}

    def test_exon_coordinates(self) -> None:
        group = self._make_group(_GFF3_ONE_TRANSCRIPT)
        transcript = _build_transcript_from_record_group(group, {})
        first = transcript.spliced_exons[0]
        assert first.start == 1000
        assert first.end == 2000

    def test_transcript_is_transcript_type(self) -> None:
        group = self._make_group(_GFF3_ONE_TRANSCRIPT)
        result = _build_transcript_from_record_group(group, {})
        assert isinstance(result, Transcript)

    def test_exons_are_exon_type(self) -> None:
        group = self._make_group(_GFF3_ONE_TRANSCRIPT)
        transcript = _build_transcript_from_record_group(group, {})
        for exon in transcript.spliced_exons:
            assert isinstance(exon, Exon)


class TestBuildTranslation:
    """Tests for the private _build_translation helper."""

    def test_translation_built_from_cds(self) -> None:
        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        groups = _group_gff3_records(records)
        transcript = _build_transcript_from_record_group(groups[0], {})
        assert transcript.translation is not None
        assert transcript.translation.start == 1100
        assert transcript.translation.end == 4900

    def test_translation_none_without_cds(self) -> None:
        records = GFF3Parser.parse_text(_GFF3_TWO_TRANSCRIPTS)
        groups = _group_gff3_records(records)
        transcript = _build_transcript_from_record_group(groups[0], {})
        assert transcript.translation is None

    def test_translation_back_reference_to_transcript(self) -> None:
        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        groups = _group_gff3_records(records)
        transcript = _build_transcript_from_record_group(groups[0], {})
        assert transcript.translation.transcript is transcript

    def test_translation_stable_id_from_protein_id_attribute(self) -> None:
        gff_text = _GFF3_ONE_TRANSCRIPT.replace(
            "CDS\t1100\t2000\t.\t+\t0\tParent=transcript:ENST001",
            "CDS\t1100\t2000\t.\t+\t0\tParent=transcript:ENST001;protein_id=ENSP001",
        )
        records = GFF3Parser.parse_text(gff_text)
        groups = _group_gff3_records(records)
        transcript = _build_transcript_from_record_group(groups[0], {})
        assert transcript.translation.stable_id == "ENSP001"

    def test_empty_cds_list_returns_none(self) -> None:
        assert _build_translation([]) is None


class TestBuildGene:
    """Tests for gene construction and linking, exercised via GFF3TranscriptSource."""

    def test_gene_linked_to_transcript(self, tmp_path) -> None:
        gff_file = tmp_path / "test.gff3"
        gff_file.write_text(_GFF3_ONE_TRANSCRIPT)
        source = GFF3TranscriptSource(file_path=str(gff_file))
        transcripts = source.get_all_features()
        assert transcripts[0].gene is not None
        assert transcripts[0].gene.stable_id == "ENSG001"

    def test_gene_fields(self, tmp_path) -> None:
        gff_file = tmp_path / "test.gff3"
        gff_file.write_text(_GFF3_ONE_TRANSCRIPT)
        source = GFF3TranscriptSource(file_path=str(gff_file))
        gene = source.get_all_features()[0].gene
        assert gene.name == "BRCA1"
        assert gene.start == 1000
        assert gene.end == 5000

    def test_gene_transcripts_list_contains_transcript(self, tmp_path) -> None:
        gff_file = tmp_path / "test.gff3"
        gff_file.write_text(_GFF3_ONE_TRANSCRIPT)
        source = GFF3TranscriptSource(file_path=str(gff_file))
        transcripts = source.get_all_features()
        assert transcripts[0].gene.transcripts == transcripts

    def test_no_gene_record_leaves_gene_unset(self, tmp_path) -> None:
        gff_file = tmp_path / "test.gff3"
        gff_file.write_text(_GFF3_TWO_TRANSCRIPTS)
        source = GFF3TranscriptSource(file_path=str(gff_file))
        transcripts = source.get_all_features()
        assert all(t.gene is None for t in transcripts)

    def test_shared_gene_is_deduplicated(self, tmp_path) -> None:
        gff_text = (
            "chr1\tensembl\tgene\t1000\t9000\t.\t+\t.\tID=gene:ENSG999;Name=FOO\n"
            "chr1\tensembl\tmRNA\t1000\t5000\t.\t+\t.\tID=transcript:ENSTA;Parent=gene:ENSG999\n"
            "chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tID=exon:EA;Parent=transcript:ENSTA\n"
            "chr1\tensembl\tmRNA\t6000\t9000\t.\t+\t.\tID=transcript:ENSTB;Parent=gene:ENSG999\n"
            "chr1\tensembl\texon\t6000\t7000\t.\t+\t.\tID=exon:EB;Parent=transcript:ENSTB\n"
        )
        gff_file = tmp_path / "shared.gff3"
        gff_file.write_text(gff_text)
        source = GFF3TranscriptSource(file_path=str(gff_file))
        transcripts = source.get_all_features()
        assert len(transcripts) == 2
        assert transcripts[0].gene is transcripts[1].gene


class TestBuildIntrons:
    """Tests for intron construction, exercised via _build_transcript_from_record_group."""

    def test_two_exons_yield_one_intron(self) -> None:
        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        groups = _group_gff3_records(records)
        transcript = _build_transcript_from_record_group(groups[0], {})
        assert len(transcript.introns) == 1

    def test_intron_coordinates(self) -> None:
        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        groups = _group_gff3_records(records)
        transcript = _build_transcript_from_record_group(groups[0], {})
        intron = transcript.introns[0]
        # exon1 ends at 2000, exon2 starts at 3000
        assert intron.start == 2001
        assert intron.end == 2999

    def test_intron_prev_next_exon_forward_strand(self) -> None:
        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        groups = _group_gff3_records(records)
        transcript = _build_transcript_from_record_group(groups[0], {})
        intron = transcript.introns[0]
        assert intron.prev_exon is transcript.spliced_exons[0]
        assert intron.next_exon is transcript.spliced_exons[1]

    def test_intron_prev_next_exon_reverse_strand(self) -> None:
        records = GFF3Parser.parse_text(_GFF3_REVERSE_STRAND)
        groups = _group_gff3_records(records)
        transcript = _build_transcript_from_record_group(groups[0], {})
        intron = transcript.introns[0]
        # spliced_exons is genomic-ascending; reverse strand reads from the higher
        # genomic coordinate first, so prev_exon should be the later-starting exon
        assert intron.prev_exon is transcript.spliced_exons[1]
        assert intron.next_exon is transcript.spliced_exons[0]

    def test_single_exon_yields_no_introns(self) -> None:
        gff_text = (
            "chr1\tensembl\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript:ENST900;Parent=gene:ENSG900\n"
            "chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tID=exon:ENSE900;Parent=transcript:ENST900\n"
        )
        records = GFF3Parser.parse_text(gff_text)
        groups = _group_gff3_records(records)
        transcript = _build_transcript_from_record_group(groups[0], {})
        assert transcript.introns == []

    def test_introns_are_intron_type(self) -> None:
        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        groups = _group_gff3_records(records)
        transcript = _build_transcript_from_record_group(groups[0], {})
        for intron in transcript.introns:
            assert isinstance(intron, Intron)


class TestSliceSharing:
    """Tests for the shared, full-chromosome Slice."""

    def test_slice_start_is_one(self, tmp_path) -> None:
        gff_file = tmp_path / "test.gff3"
        gff_file.write_text(_GFF3_ONE_TRANSCRIPT)
        source = GFF3TranscriptSource(file_path=str(gff_file))
        transcript = source.get_all_features()[0]
        assert transcript.slice.start == 1

    def test_transcript_and_exon_share_slice_instance(self, tmp_path) -> None:
        gff_file = tmp_path / "test.gff3"
        gff_file.write_text(_GFF3_ONE_TRANSCRIPT)
        source = GFF3TranscriptSource(file_path=str(gff_file))
        transcript = source.get_all_features()[0]
        assert transcript.slice is transcript.spliced_exons[0].slice

    def test_transcript_and_gene_share_slice_instance(self, tmp_path) -> None:
        gff_file = tmp_path / "test.gff3"
        gff_file.write_text(_GFF3_ONE_TRANSCRIPT)
        source = GFF3TranscriptSource(file_path=str(gff_file))
        transcript = source.get_all_features()[0]
        assert transcript.slice is transcript.gene.slice

    def test_different_seq_regions_get_different_slices(self, tmp_path) -> None:
        gff_text = (
            "chr1\tensembl\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript:ENSTA;Parent=gene:ENSGA\n"
            "chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tID=exon:EA;Parent=transcript:ENSTA\n"
            "chr2\tensembl\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript:ENSTB;Parent=gene:ENSGB\n"
            "chr2\tensembl\texon\t1000\t2000\t.\t+\t.\tID=exon:EB;Parent=transcript:ENSTB\n"
        )
        gff_file = tmp_path / "two_regions.gff3"
        gff_file.write_text(gff_text)
        source = GFF3TranscriptSource(file_path=str(gff_file))
        transcripts = source.get_all_features()
        assert transcripts[0].slice is not transcripts[1].slice
        assert transcripts[0].slice.seq_region_name != transcripts[1].slice.seq_region_name


class TestFlushBufferStreaming:
    """Tests for _flush_buffer, the building block of the streaming iterator methods.

    Exercised directly (bypassing the parser) so gene-cache behavior across multiple
    flushes can be verified without depending on GFF3Parser.open(mode=...).
    """

    def test_single_flush_builds_transcript(self) -> None:
        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        source = GFF3TranscriptSource()
        transcripts = source._flush_buffer(list(records), {}, {}, {})  # pylint: disable=protected-access
        assert len(transcripts) == 1
        assert transcripts[0].stable_id == "ENST001"

    def test_gene_reused_across_separate_flushes_for_same_gene(self) -> None:
        gene_line = "chr1\tensembl\tgene\t1000\t9000\t.\t+\t.\tID=gene:ENSG999;Name=FOO\n"
        transcript_a = (
            "chr1\tensembl\tmRNA\t1000\t5000\t.\t+\t.\tID=transcript:ENSTA;Parent=gene:ENSG999\n"
            "chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tID=exon:EA;Parent=transcript:ENSTA\n"
        )
        transcript_b = (
            "chr1\tensembl\tmRNA\t6000\t9000\t.\t+\t.\tID=transcript:ENSTB;Parent=gene:ENSG999\n"
            "chr1\tensembl\texon\t6000\t7000\t.\t+\t.\tID=exon:EB;Parent=transcript:ENSTB\n"
        )

        gene_records: dict = {}
        gene_objects: dict = {}
        slice_cache: dict = {}
        source = GFF3TranscriptSource()

        # First block (gene line + transcript A) is flushed on its own, as the
        # streaming iterator would do when it sees the next top-level record.
        first_records = GFF3Parser.parse_text(gene_line + transcript_a)
        first_transcripts = source._flush_buffer(  # pylint: disable=protected-access
            list(first_records), gene_records, gene_objects, slice_cache
        )

        # Second block (transcript B alone -- its gene line was already consumed
        # and flushed above) must still resolve to the SAME Gene object.
        second_records = GFF3Parser.parse_text(transcript_b)
        gene_records["gene:ENSG999"] = GFF3Parser.parse_text(gene_line)[0]
        second_transcripts = source._flush_buffer(  # pylint: disable=protected-access
            list(second_records), gene_records, gene_objects, slice_cache
        )

        assert first_transcripts[0].gene is second_transcripts[0].gene
        assert second_transcripts[0].gene.transcripts == [first_transcripts[0], second_transcripts[0]]

    def test_records_to_features_uses_default_context(self) -> None:
        records = GFF3Parser.parse_text(_GFF3_ONE_TRANSCRIPT)
        source = GFF3TranscriptSource()

        transcripts = source._records_to_features(list(records))  # pylint: disable=protected-access

        assert len(transcripts) == 1
        assert transcripts[0].stable_id == "ENST001"

    def test_records_to_features_with_context_uses_shared_streaming_caches(self) -> None:
        gene_line = "chr1\tensembl\tgene\t1000\t9000\t.\t+\t.\tID=gene:ENSG999;Name=FOO\n"
        transcript_block = (
            "chr1\tensembl\tmRNA\t1000\t5000\t.\t+\t.\tID=transcript:ENSTA;Parent=gene:ENSG999\n"
            "chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tID=exon:EA;Parent=transcript:ENSTA\n"
        )

        gene_records = {"gene:ENSG999": GFF3Parser.parse_text(gene_line)[0]}
        gene_objects: dict = {}
        slice_cache: dict = {}
        source = GFF3TranscriptSource()

        transcripts = source._records_to_features_with_context(  # pylint: disable=protected-access
            list(GFF3Parser.parse_text(transcript_block)),
            gene_records=gene_records,
            gene_objects=gene_objects,
            slice_cache=slice_cache,
        )

        assert len(transcripts) == 1
        assert transcripts[0].gene is not None
        assert transcripts[0].gene.stable_id == "ENSG999"
        assert gene_objects["gene:ENSG999"].transcripts == [transcripts[0]]


class TestGFF3TranscriptSource:
    """Integration tests for GFF3TranscriptSource.

    get_features_by_region()/iter_features_by_region() require a tabix-indexed
    (bgzipped) file -- GFF3Parser's default "auto" mode falls back to "sequential"
    for a plain, non-indexed file, and seek() only works in "query" mode. Those
    are covered by the tabix-backed integration tests, not here.
    """

    def test_get_all_features_returns_transcripts(self, tmp_path) -> None:
        gff_file = tmp_path / "test.gff3"
        gff_file.write_text(_GFF3_TWO_TRANSCRIPTS)
        source = GFF3TranscriptSource(file_path=str(gff_file))
        features = source.get_all_features()
        assert len(features) == 2
        assert all(isinstance(f, Transcript) for f in features)

    def test_get_all_features_with_exons(self, tmp_path) -> None:
        gff_file = tmp_path / "test.gff3"
        gff_file.write_text(_GFF3_ONE_TRANSCRIPT)
        source = GFF3TranscriptSource(file_path=str(gff_file))
        transcripts = source.get_all_features()
        assert len(transcripts) == 1
        assert len(transcripts[0].spliced_exons) == 2

    def test_iter_all_features_yields_same_transcripts_as_get_all_features(self, tmp_path) -> None:
        gff_file = tmp_path / "test.gff3"
        gff_file.write_text(_GFF3_TWO_TRANSCRIPTS)
        source = GFF3TranscriptSource(file_path=str(gff_file))

        iterated = list(source.iter_all_features())

        assert len(iterated) == 2
        assert all(isinstance(f, Transcript) for f in iterated)
        assert {t.stable_id for t in iterated} == {"ENST001", "ENST002"}

    def test_close_can_be_called(self, tmp_path) -> None:
        gff_file = tmp_path / "test.gff3"
        gff_file.write_text(_GFF3_ONE_TRANSCRIPT)
        source = GFF3TranscriptSource(file_path=str(gff_file))
        source.close()

    def test_context_manager_closes_on_exit(self, tmp_path) -> None:
        gff_file = tmp_path / "test.gff3"
        gff_file.write_text(_GFF3_ONE_TRANSCRIPT)
        with GFF3TranscriptSource(file_path=str(gff_file)) as source:
            transcripts = source.get_all_features()
        assert len(transcripts) == 1

    def test_context_manager_returns_source(self, tmp_path) -> None:
        gff_file = tmp_path / "test.gff3"
        gff_file.write_text(_GFF3_ONE_TRANSCRIPT)
        with GFF3TranscriptSource(file_path=str(gff_file)) as source:
            assert isinstance(source, GFF3TranscriptSource)
