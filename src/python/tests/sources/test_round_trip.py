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
"""The overlap-parity gate: a corpus saved to disk and read back must annotate identically.

This is what makes the file-backed path trustworthy rather than merely runnable. The generated
corpus and the file source are two independent routes to the same overlaps, so each checks the
other: a corpus is written as FASTA, GFF3 and VCF, read back through the readers, and compared
against the corpus it came from -- same pairs, and the same terms from the same annotator.

It earns its place. It caught two real defects that every other test missed, because the readers'
own tests assert what the readers produce rather than whether it means the same thing: the
translation was seeded with genomic coordinates where the model wants offsets from the transcript's
start, and the transcript biotype was never read at all. Pairs matched perfectly throughout --
665 of 665 -- while 510 of those overlaps silently lost every coding consequence.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
import subprocess
from typing import Protocol

import pytest

import ensembl.vep.annotators  # noqa: F401  pylint: disable=unused-import

# ^ registers the strategies, and keeps the annotators-before-caches import order.
from ensembl.vep.annotators.registry import get_annotator
from ensembl.vep.common import Strand
from ensembl.vep.corpus import build_corpus
from ensembl.vep.corpus.exporter import get_exporter
from ensembl.vep.corpus.exporter.formats.file.fasta import FASTA_NAME
from ensembl.vep.corpus.exporter.formats.file.gff3 import GFF3_NAME, gff3_lines
from ensembl.vep.corpus.exporter.formats.file.vcf import VCF_NAME
from ensembl.vep.corpus.spec import CorpusSpec
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.sources.file.gff3.transcript import GFF3TranscriptSource
from ensembl.vep.sources.file.gtf.transcript import GTFTranscriptSource
from ensembl.vep.config import Config
from ensembl.vep.explain.selection import NamedOverlap, OverlapSelection
from ensembl.vep.sources.synonym import SequenceRegionResolver

REFERENCE_ANNOTATOR = "consequence.linear"
"""The reference strategy. Every other strategy is already checked against it by conformance."""


class _HtslibToolchain(Protocol):
    """Structural stand-in for the session ``htslib_toolchain`` fixture."""

    bgzip: str
    tabix: str


def _feature_id(overlap: TranscriptVariationAllele) -> str:
    """Return a feature's id in a form both sources agree on.

    A GFF3 names a transcript ``ID=transcript:ENST0001`` and the reader keeps that spelling, while
    the generated corpus calls the same transcript ``ENST0001``. The prefix is a spelling of the
    id, not a different feature, so it is normalised away rather than asserted on.
    """
    return overlap.feature.stable_id.split(":")[-1]


def _allele(value: str | None) -> str:
    """Return an allele with either spelling of "empty" collapsed to one.

    The corpus spells an empty allele ``""`` and the VCF reader spells it ``-``, following VEP.
    """
    return "" if value is None or value == "-" else value


def _locus(overlap: TranscriptVariationAllele) -> tuple[str | None, int, int, str, str]:
    """Return the variant's position and alleles, as the two sources should agree on them."""
    allele = overlap.overlapping_feature
    seq_region_name = allele.slice.seq_region_name if allele.slice is not None else None
    return (
        seq_region_name,
        allele.start,
        allele.end,
        _allele(allele.reference_allele),
        _allele(allele.alternate_allele),
    )


def _terms(named: Sequence[NamedOverlap]) -> dict[tuple[str, str], tuple[str, ...]]:
    """Annotate a set of overlaps and return its terms, keyed by variant name and feature.

    A fresh annotator per set: the codon cache is process-global, and a shared one could hide a
    disagreement by answering the second set from the first set's entries.
    """
    named = list(named)
    annotator = get_annotator(REFERENCE_ANNOTATOR)
    results = annotator.annotate_feature_overlaps([overlap for _name, overlap in named])
    return {
        (name, _feature_id(overlap)): tuple(consequence.term for consequence in consequences)
        for (name, overlap), consequences in zip(named, results, strict=True)
    }


@pytest.fixture(name="generated", scope="module")
def fixture_generated() -> tuple[NamedOverlap, ...]:
    """The corpus the round trip starts from: small, but every boundary overlap is in it."""
    return tuple(build_corpus(CorpusSpec(kind="synthetic.default", fuzz=0)).named_overlaps)


@pytest.fixture(name="restored", scope="module")
def fixture_restored(tmp_path_factory: pytest.TempPathFactory) -> tuple[NamedOverlap, ...]:
    """The same corpus written to disk and read back through the file readers."""
    directory: Path = tmp_path_factory.mktemp("corpus")
    get_exporter("file").export(build_corpus(CorpusSpec(kind="synthetic.default", fuzz=0)), directory)
    selection = OverlapSelection(
        config=Config.for_command(
            "explain",
            input_file=directory / VCF_NAME,
            gff=directory / GFF3_NAME,
            fasta=directory / FASTA_NAME,
        )
    )
    return tuple(selection.all_overlaps())


class TestPairParity:
    """The same overlaps come back."""

    def test_the_same_number_of_overlaps(
        self, generated: tuple[NamedOverlap, ...], restored: tuple[NamedOverlap, ...]
    ) -> None:
        """A pair lost or invented in the round trip is a defect in the readers or the writer."""
        assert len(restored) == len(generated)

    def test_the_same_variant_names(
        self, generated: tuple[NamedOverlap, ...], restored: tuple[NamedOverlap, ...]
    ) -> None:
        """The generator's name goes into the VCF ID column so that --name keeps selecting the same
        variant once the corpus has been read back."""
        assert {name for name, _overlap in restored} == {name for name, _overlap in generated}

    def test_the_same_features(
        self, generated: tuple[NamedOverlap, ...], restored: tuple[NamedOverlap, ...]
    ) -> None:
        """Every transcript survives the GFF3 round trip, prefix aside."""
        assert {_feature_id(overlap) for _name, overlap in restored} == {
            _feature_id(overlap) for _name, overlap in generated
        }

    def test_the_same_loci(
        self, generated: tuple[NamedOverlap, ...], restored: tuple[NamedOverlap, ...]
    ) -> None:
        """Positions and alleles survive the anchor base a VCF adds to an indel and the reader
        strips again -- including at position 1, where there is no preceding base to anchor to."""
        assert {_locus(overlap) for _name, overlap in restored} == {
            _locus(overlap) for _name, overlap in generated
        }


class TestTermParity:
    """The same overlaps still mean the same thing.

    Pair parity is not enough on its own: it held perfectly while the coding consequences were
    being dropped, because a transcript read back with broken coding boundaries still pairs with
    exactly the same variants.
    """

    def test_every_overlap_yields_the_same_terms(
        self, generated: tuple[NamedOverlap, ...], restored: tuple[NamedOverlap, ...]
    ) -> None:
        """Same annotator, same catalogue, same terms in the same order."""
        from_generated = _terms(generated)
        from_restored = _terms(restored)

        assert set(from_restored) == set(from_generated)
        disagreements = [
            f"{key}: restored {from_restored[key]} != generated {from_generated[key]}"
            for key in sorted(from_generated)
            if from_restored[key] != from_generated[key]
        ]
        assert not disagreements, "\n".join(disagreements[:20])

    def test_the_corpus_actually_exercises_coding_consequences(
        self, generated: tuple[NamedOverlap, ...]
    ) -> None:
        """Without this, the gate above could pass on a corpus that never coded for anything --
        which is the exact failure it exists to catch."""
        emitted = {term for terms in _terms(generated).values() for term in terms}
        assert {"missense_variant", "synonymous_variant", "frameshift_variant"} <= emitted


def _sorted_by_position(lines: list[str]) -> list[str]:
    """Return GFF3/GTF data lines sorted by (seqid, start), headers first.

    Tabix requires a coordinate-sorted file. The generators this test draws fixtures from group
    rows by feature (a gene immediately followed by its transcripts), not by position, so the
    output needs re-sorting before it can be indexed -- interleaving unrelated transcripts' rows
    this way is realistic anyway, and exercises the grouping code's tolerance for it.
    """
    header = [line for line in lines if line.startswith("#")]
    data = [line for line in lines if not line.startswith("#")]
    data.sort(key=lambda line: (line.split("\t")[0], int(line.split("\t")[3])))
    return header + data


def _fingerprint(transcript: Transcript) -> tuple:
    """Return the fields a grouping/building change could plausibly disturb.

    Ids, order (via the caller comparing sequences of these), exon spans/phase/frame, intron spans
    and their prev/next exon, coding boundaries (genomic and cDNA), translation bounds, and gene
    identity including name and biotype.
    """
    gene = transcript.gene
    return (
        transcript.stable_id,
        transcript.start,
        transcript.end,
        transcript.strand,
        transcript.biotype,
        tuple(
            (exon.stable_id, exon.start, exon.end, exon.phase, exon.frame)
            for exon in transcript.spliced_exons
        ),
        tuple(
            (
                intron.start,
                intron.end,
                intron.prev_exon.stable_id if intron.prev_exon is not None else None,
                intron.next_exon.stable_id if intron.next_exon is not None else None,
            )
            for intron in transcript.introns
        ),
        transcript.coding_start,
        transcript.coding_end,
        transcript.cdna_coding_start,
        transcript.cdna_coding_end,
        (
            (transcript.translation.start, transcript.translation.end)
            if transcript.translation is not None
            else None
        ),
        gene.stable_id if gene is not None else None,
        gene.name if gene is not None else None,
        gene.biotype if gene is not None else None,
    )


class TestIndexedMatchesSequentialGFF3:
    """The indexed query path (``get_features_by_region``) and the plain
    sequential/streaming path (``iter_all_features``) share the same grouping and model-building
    code, so they must agree on every transcript they both see.

    Sharing that code means a regression in the grouping/building logic itself would not show up
    here (both paths would be wrong the same way) -- see ``TestFixedExpectations`` below for a
    check against independently hand-computed values instead.
    """

    @pytest.fixture(scope="class", name="indexed_gff3_path")
    def fixture_indexed_gff3_path(
        self, tmp_path_factory: pytest.TempPathFactory, htslib_toolchain: _HtslibToolchain
    ) -> Path:
        """A synthetic-corpus GFF3, bgzipped and tabix-indexed with a real HTSlib toolchain."""
        directory = tmp_path_factory.mktemp("stage3b_gff3")
        corpus = build_corpus(CorpusSpec(kind="synthetic.default", fuzz=0))
        gff3_path = directory / "fixture.gff3"
        lines = _sorted_by_position(list(gff3_lines((corpus,))))
        gff3_path.write_text("\n".join(lines) + "\n")

        subprocess.run([htslib_toolchain.bgzip, "-f", str(gff3_path)], check=True)
        gz_path = gff3_path.with_name(gff3_path.name + ".gz")
        subprocess.run([htslib_toolchain.tabix, "-f", "-p", "gff", str(gz_path)], check=True)
        return gz_path

    @pytest.mark.integration
    def test_indexed_matches_sequential_over_several_windows(self, indexed_gff3_path: Path) -> None:
        """Compare ``get_features_by_region`` against ``iter_all_features`` over several
        windows: a single-transcript hit, a multi-transcript window, and an intergenic gap."""

        sequential_source = GFF3TranscriptSource(file_path=str(indexed_gff3_path))
        sequential_transcripts = list(sequential_source.iter_all_features())
        assert len(sequential_transcripts) >= 3, "the default corpus should have several transcripts"
        by_seq_region: dict[str, list[Transcript]] = {}
        for transcript in sequential_transcripts:
            assert transcript.slice is not None
            by_seq_region.setdefault(transcript.slice.seq_region_name, []).append(transcript)

        # Declaring the source names is what lets an intergenic-gap query return no records
        # directly, instead of probing a "chr"-toggled alternate spelling that does not exist in
        # this fixture's index.
        indexed_source = GFF3TranscriptSource(
            file_path=str(indexed_gff3_path),
            source_region_resolver=SequenceRegionResolver(source_names=frozenset(by_seq_region)),
        )

        windows: list[tuple[str, int, int]] = []
        for seq_region_name, transcripts in by_seq_region.items():
            # A narrow window inside the first transcript.
            first = transcripts[0]
            windows.append((seq_region_name, first.start, first.start))
            # A wide window spanning every transcript on this chromosome, likely to pull in
            # several at once (and their shared genes).
            windows.append(
                (seq_region_name, min(t.start for t in transcripts), max(t.end for t in transcripts))
            )
            # An intergenic gap just past the last transcript: no candidates on either path.
            last = transcripts[-1]
            windows.append((seq_region_name, last.end + 10_000, last.end + 10_001))

        assert windows, "the corpus should cover at least one sequence region"

        for seq_region_name, start, end in windows:
            indexed = indexed_source.get_features_by_region(
                seq_region_name=seq_region_name, start=start, end=end
            )
            expected_ids = [
                transcript.stable_id
                for transcript in by_seq_region[seq_region_name]
                if transcript.start <= end and transcript.end >= start
            ]

            assert [transcript.stable_id for transcript in indexed] == expected_ids

            sequential_by_id = {transcript.stable_id: transcript for transcript in sequential_transcripts}
            assert [_fingerprint(t) for t in indexed] == [
                _fingerprint(sequential_by_id[stable_id]) for stable_id in expected_ids
            ]

        indexed_source.close()
        sequential_source.close()


# Two overlapping transcripts on one gene, a reverse-strand transcript on another, and a
# single-exon non-coding one, spread across three separate regions so narrow-window,
# multi-transcript-window and intergenic-gap windows are all exercised.
_GTF_FIXTURE = """\
chr1\tensembl\tgene\t1000\t5000\t.\t+\t.\tgene_id "gene:1"; gene_name "GENE1"; gene_biotype "protein_coding";
chr1\tensembl\ttranscript\t1000\t5000\t.\t+\t.\tgene_id "gene:1"; transcript_id "transcript:1a";
chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tgene_id "gene:1"; transcript_id "transcript:1a";
chr1\tensembl\tCDS\t1100\t2000\t.\t+\t0\tgene_id "gene:1"; transcript_id "transcript:1a"; protein_id "P1A";
chr1\tensembl\texon\t3000\t5000\t.\t+\t.\tgene_id "gene:1"; transcript_id "transcript:1a";
chr1\tensembl\tCDS\t3000\t4900\t.\t+\t0\tgene_id "gene:1"; transcript_id "transcript:1a"; protein_id "P1A";
chr1\tensembl\ttranscript\t1200\t4800\t.\t+\t.\tgene_id "gene:1"; transcript_id "transcript:1b";
chr1\tensembl\texon\t1200\t4800\t.\t+\t.\tgene_id "gene:1"; transcript_id "transcript:1b";
chr1\tensembl\tCDS\t1300\t4700\t.\t+\t0\tgene_id "gene:1"; transcript_id "transcript:1b"; protein_id "P1B";
chr1\tensembl\tgene\t8000\t9000\t.\t-\t.\tgene_id "gene:2"; gene_name "GENE2"; gene_biotype "lncRNA";
chr1\tensembl\ttranscript\t8000\t9000\t.\t-\t.\tgene_id "gene:2"; transcript_id "transcript:2";
chr1\tensembl\texon\t8000\t8300\t.\t-\t.\tgene_id "gene:2"; transcript_id "transcript:2";
chr1\tensembl\texon\t8700\t9000\t.\t-\t.\tgene_id "gene:2"; transcript_id "transcript:2";
chr1\tensembl\tgene\t15000\t16000\t.\t+\t.\tgene_id "gene:3"; gene_name "GENE3"; gene_biotype "protein_coding";
chr1\tensembl\ttranscript\t15000\t16000\t.\t+\t.\tgene_id "gene:3"; transcript_id "transcript:3";
chr1\tensembl\texon\t15000\t16000\t.\t+\t.\tgene_id "gene:3"; transcript_id "transcript:3";
"""

_GTF_WINDOWS: tuple[tuple[str, int, int], ...] = (
    ("chr1", 1400, 1600),  # inside both transcript:1a and transcript:1b
    ("chr1", 8100, 8200),  # inside transcript:2 only
    ("chr1", 15500, 15500),  # inside transcript:3 only
    ("chr1", 6000, 7000),  # intergenic gap: no candidates
)


class TestIndexedMatchesSequentialGTF:
    """The GTF adapter shares the same grouping code as GFF3, through the same base class."""

    @pytest.fixture(scope="class", name="indexed_gtf_path")
    def fixture_indexed_gtf_path(
        self, tmp_path_factory: pytest.TempPathFactory, htslib_toolchain: _HtslibToolchain
    ) -> Path:
        """The hand-authored GTF fixture, bgzipped and tabix-indexed with a real HTSlib toolchain."""
        directory = tmp_path_factory.mktemp("stage3b_gtf")
        gtf_path = directory / "fixture.gtf"
        lines = _sorted_by_position(_GTF_FIXTURE.splitlines())
        gtf_path.write_text("\n".join(lines) + "\n")

        subprocess.run([htslib_toolchain.bgzip, "-f", str(gtf_path)], check=True)
        gz_path = gtf_path.with_name(gtf_path.name + ".gz")
        subprocess.run([htslib_toolchain.tabix, "-f", "-p", "gff", str(gz_path)], check=True)
        return gz_path

    @pytest.mark.integration
    def test_indexed_matches_sequential_over_several_windows(self, indexed_gtf_path: Path) -> None:
        """Compare the indexed and sequential paths over the fixture's three regions and a gap."""

        sequential_source = GTFTranscriptSource(file_path=str(indexed_gtf_path))
        sequential_transcripts = list(sequential_source.iter_all_features())
        sequential_by_id = {transcript.stable_id: transcript for transcript in sequential_transcripts}
        assert len(sequential_transcripts) == 4

        indexed_source = GTFTranscriptSource(
            file_path=str(indexed_gtf_path),
            source_region_resolver=SequenceRegionResolver(source_names=frozenset({"chr1"})),
        )

        for seq_region_name, start, end in _GTF_WINDOWS:
            indexed = indexed_source.get_features_by_region(
                seq_region_name=seq_region_name, start=start, end=end
            )
            expected_ids = [
                transcript.stable_id
                for transcript in sequential_transcripts
                if transcript.start <= end and transcript.end >= start
            ]

            assert [transcript.stable_id for transcript in indexed] == expected_ids
            assert [_fingerprint(t) for t in indexed] == [
                _fingerprint(sequential_by_id[stable_id]) for stable_id in expected_ids
            ]

        indexed_source.close()
        sequential_source.close()


# Four independent transcripts, each pinned to its own values by hand rather than by the reader
# under test: (A) a duplicate transcript row whose second copy has no Parent, (B) a transcript
# with no Parent at all, (C) a gene row that appears after its own transcript and exon, and (D) a
# reverse-strand transcript with a single-exon CDS.
_GFF3_FIXED_EXPECTATIONS_FIXTURE = """\
##gff-version 3
chr1\tensembl\tgene\t1000\t2000\t.\t+\t.\tID=gene:GA;Name=GENE_A;biotype=protein_coding
chr1\tensembl\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript:TA;Parent=gene:GA
chr1\tensembl\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript:TA
chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tID=exon:EA;Parent=transcript:TA
chr1\tensembl\tmRNA\t3000\t4000\t.\t+\t.\tID=transcript:TB
chr1\tensembl\texon\t3000\t4000\t.\t+\t.\tID=exon:EB;Parent=transcript:TB
chr1\tensembl\tmRNA\t5000\t6000\t.\t+\t.\tID=transcript:TC;Parent=gene:GC
chr1\tensembl\texon\t5000\t6000\t.\t+\t.\tID=exon:EC;Parent=transcript:TC
chr1\tensembl\tgene\t5000\t6000\t.\t+\t.\tID=gene:GC;Name=GENE_C;biotype=lncRNA
chr1\tensembl\tgene\t9000\t9500\t.\t-\t.\tID=gene:GD;Name=GENE_D;biotype=protein_coding
chr1\tensembl\tmRNA\t9000\t9500\t.\t-\t.\tID=transcript:TD;Parent=gene:GD
chr1\tensembl\texon\t9000\t9500\t.\t-\t.\tID=exon:ED;Parent=transcript:TD
chr1\tensembl\tCDS\t9101\t9400\t.\t-\t0\tParent=transcript:TD;protein_id=PROT_D
"""


class TestFixedExpectations:
    """Independently hand-computed values for a few grouping edge cases.

    ``TestIndexedMatchesSequentialGFF3``/``GTF`` compare two readers that share the same grouping
    and model-building code, so a bug in that shared code would not show up there -- both sides
    would agree on the same wrong answer. These values are worked out by hand instead, against
    the plain, non-streaming read path (``get_all_features``).
    """

    @pytest.fixture(scope="class", name="transcripts_by_id")
    def fixture_transcripts_by_id(self, tmp_path_factory: pytest.TempPathFactory) -> dict[str, Transcript]:
        """The fixture's transcripts, keyed by id, built through the plain (unindexed) reader."""
        directory = tmp_path_factory.mktemp("fixed_expectations")
        gff3_path = directory / "fixture.gff3"
        gff3_path.write_text(_GFF3_FIXED_EXPECTATIONS_FIXTURE)

        source = GFF3TranscriptSource(file_path=str(gff3_path))
        transcripts = source.get_all_features()
        source.close()
        return {transcript.stable_id: transcript for transcript in transcripts}

    def test_all_four_transcripts_are_built(self, transcripts_by_id: dict[str, Transcript]) -> None:
        """Nothing is dropped or duplicated."""
        assert set(transcripts_by_id) == {
            "transcript:TA",
            "transcript:TB",
            "transcript:TC",
            "transcript:TD",
        }

    def test_duplicate_transcript_row_keeps_the_earlier_gene(
        self, transcripts_by_id: dict[str, Transcript]
    ) -> None:
        """The second, Parent-less copy of transcript:TA does not clear the gene the first
        copy declared."""
        transcript = transcripts_by_id["transcript:TA"]

        assert transcript.gene is not None
        assert transcript.gene.stable_id == "gene:GA"
        assert transcript.gene.name == "GENE_A"

    def test_transcript_with_no_parent_has_no_gene(self, transcripts_by_id: dict[str, Transcript]) -> None:
        """transcript:TB never declares a Parent, so it is built with no gene at all."""
        transcript = transcripts_by_id["transcript:TB"]

        assert transcript.gene is None

    def test_gene_row_after_its_transcript_still_resolves(
        self, transcripts_by_id: dict[str, Transcript]
    ) -> None:
        """gene:GC's row comes after transcript:TC's own rows in the file; it is still attached."""
        transcript = transcripts_by_id["transcript:TC"]

        assert transcript.gene is not None
        assert transcript.gene.stable_id == "gene:GC"
        assert transcript.gene.name == "GENE_C"
        assert transcript.gene.biotype == "lncRNA"

    def test_reverse_strand_coding_transcript_matches_hand_computed_bounds(
        self, transcripts_by_id: dict[str, Transcript]
    ) -> None:
        """A single-exon reverse-strand CDS (9101-9400) on a 9000-9500 transcript.

        By hand: cds_start=9101, cds_end=9400 (min/max of the CDS row); reverse-strand
        translation offsets are measured from the transcript's end, so translation.start =
        9500 - 9400 + 1 = 101 and translation.end = 9500 - 9101 + 1 = 400. With a single exon,
        the cDNA and genomic coding boundaries fall out the same way: coding_start = 9400,
        coding_end = 9101, cdna_coding_start = 101, cdna_coding_end = 400.
        """
        transcript = transcripts_by_id["transcript:TD"]

        assert (transcript.start, transcript.end, transcript.strand) == (9000, 9500, Strand.REVERSE)
        assert transcript.translation is not None
        assert (transcript.translation.start, transcript.translation.end) == (101, 400)
        assert (transcript.coding_start, transcript.coding_end) == (9400, 9101)
        assert (transcript.cdna_coding_start, transcript.cdna_coding_end) == (101, 400)
