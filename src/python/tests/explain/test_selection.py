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
"""Tests for naming, pairing, finding and choosing the overlaps ``vep explain`` selects from."""

from __future__ import annotations

import subprocess
from collections.abc import Iterator
from pathlib import Path
from typing import Protocol

import pytest

import ensembl.vep.annotators  # noqa: F401  pylint: disable=unused-import

# ^ keeps the annotators-before-caches import order, as the sibling runner tests do.
from ensembl.vep.common import Strand
from ensembl.vep.config import Config
from ensembl.vep.config.sections.common.source import SourceConfig
from ensembl.vep.explain.selection import (
    FeatureSelector,
    NamedOverlap,
    OverlapSelection,
    VariantSelector,
)
from ensembl.vep.models.genic.transcript import Transcript
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.models.variation.variation import Variation
from ensembl.vep.sources.file.gff3.transcript import GFF3TranscriptSource
from ensembl.vep.sources.file.vcf.variation import PicklableVCFVariationSource, VCFVariationSource
from tests.models.conftest import DummySequenceSource

DEFAULT_FLANK = SourceConfig().upstream_distance
"""The flank a run uses when none was configured, which these fixtures are built around."""

GFF3 = """\
chr1\tensembl\tgene\t1000\t5000\t.\t+\t.\tID=gene:G1;Name=GENE1;biotype=protein_coding
chr1\tensembl\tmRNA\t1000\t5000\t.\t+\t.\tID=transcript:ENST001;Parent=gene:G1
chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tID=exon:E1;Parent=transcript:ENST001
chr1\tensembl\tCDS\t1100\t2000\t.\t+\t0\tID=cds:P1;Parent=transcript:ENST001;protein_id=P1
chr1\tensembl\texon\t3000\t5000\t.\t+\t.\tID=exon:E2;Parent=transcript:ENST001
chr1\tensembl\tCDS\t3000\t4900\t.\t+\t0\tID=cds:P1;Parent=transcript:ENST001;protein_id=P1
chr1\tensembl\tgene\t20000\t21000\t.\t+\t.\tID=gene:G2;Name=GENE2;biotype=protein_coding
chr1\tensembl\tmRNA\t20000\t21000\t.\t+\t.\tID=transcript:ENST002;Parent=gene:G2
chr1\tensembl\texon\t20000\t21000\t.\t+\t.\tID=exon:E3;Parent=transcript:ENST002
chr2\tensembl\tgene\t50\t900\t.\t-\t.\tID=gene:G3;Name=GENE3;biotype=protein_coding
chr2\tensembl\tmRNA\t50\t900\t.\t-\t.\tID=transcript:ENST003;Parent=gene:G3
chr2\tensembl\texon\t50\t900\t.\t-\t.\tID=exon:E4;Parent=transcript:ENST003
"""
"""Two transcripts far apart on chr1 and one on chr2, so both the flank bound and the sequence
region split are exercised by the same fixture."""

VCF = """\
##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t1500\trs_inside\tA\tC\t.\tPASS\t.
chr1\t9000\trs_flank\tA\tC\t.\tPASS\t.
chr1\t30000\trs_far\tA\tC\t.\tPASS\t.
chr1\t1600\t.\tA\tC\t.\tPASS\t.
chr1\t1700\trs_multi\tA\tC,G\t.\tPASS\t.
chr2\t100\trs_chr2\tC\tCAT\t.\tPASS\t.
"""
"""One variant inside a transcript, one only within its flank, one beyond every transcript, one
with no ID, one multi-allelic, and one on the other sequence region."""


class _Sources:
    """A selection over the fixtures, with its two readers supplied rather than built.

    Assigning ``sources`` is the seam that keeps these tests off the filesystem for sequence: the
    readers share one in-memory sequence source, which is all the pairing needs.
    """

    def __init__(self, variants: Path, transcripts: Path, *, flank: int = DEFAULT_FLANK) -> None:
        sequence_source = DummySequenceSource(data={"chr1": "N" * 25000, "chr2": "N" * 2000})
        self.variation_source = PicklableVCFVariationSource(
            file_path=str(variants), sequence_source=sequence_source
        )
        self.transcript_source = GFF3TranscriptSource(
            file_path=str(transcripts), sequence_source=sequence_source
        )
        self.selection = OverlapSelection(
            config=Config.for_command("explain", input_file=variants, gff=transcripts, distance=flank)
        )
        self.selection.sources = (self.variation_source, self.transcript_source)

    def paired(self, *, flank: int = DEFAULT_FLANK) -> list[NamedOverlap]:
        """Every overlap the two files pair into."""
        if flank != DEFAULT_FLANK:
            narrowed = OverlapSelection(
                config=Config.for_command(
                    "explain",
                    input_file=self.variation_source.file_path,
                    gff=self.transcript_source.file_path,
                    distance=flank,
                )
            )
            narrowed.sources = (self.variation_source, self.transcript_source)
            return list(narrowed.all_overlaps())
        return list(self.selection.all_overlaps())

    def found(self, **selectors: object) -> list[NamedOverlap]:
        """The overlaps of the first variant a name or a variant selector identifies."""
        variant = selectors.pop("variant", None)
        allele_num = selectors.pop("allele_num", None)
        selector = (
            VariantSelector.parse(str(variant), allele_num=allele_num)  # type: ignore[arg-type]
            if variant is not None
            else None
        )
        return list(
            self.selection.candidate_overlaps(
                name=selectors.pop("name", None),  # type: ignore[arg-type]
                selector=selector,
            )
        )


def _write(directory: Path) -> tuple[Path, Path]:
    """Write the VCF and GFF3 fixtures into a directory, returning their paths."""
    variants = directory / "variants.vcf"
    variants.write_text(VCF, encoding="utf-8")
    transcripts = directory / "transcripts.gff3"
    transcripts.write_text(GFF3, encoding="utf-8")
    return variants, transcripts


@pytest.fixture(name="sources")
def fixture_sources(tmp_path: Path) -> _Sources:
    """Readers over the fixtures, with sequence supplied in memory."""
    return _Sources(*_write(tmp_path))


def _refuse(*_arguments: object, **_keywords: object) -> None:
    """Stand in for a reader call a targeted lookup must not make."""
    raise AssertionError("read more than the one variant's overlaps")


def _unindexed(*_arguments: object, **_keywords: object) -> None:
    """Stand in for a region query against a file with no index."""
    raise OSError("no tabix index found")


class _HtslibToolchain(Protocol):
    """Structural stand-in for the session ``htslib_toolchain`` fixture."""

    bgzip: str
    tabix: str


class TestPairing:
    """What the two files pair into."""

    def test_it_pairs_variants_with_transcripts_in_range(self, sources: _Sources) -> None:
        """Both the variant inside a transcript and the one only within its flank pair with it."""
        by_name = dict(sources.paired())
        assert by_name["rs_inside"].feature.stable_id == "transcript:ENST001"
        assert by_name["rs_flank"].feature.stable_id == "transcript:ENST001"

    def test_a_variant_beyond_every_transcript_pairs_with_nothing(self, sources: _Sources) -> None:
        """rs_far sits past ENST001's flank and short of ENST002's, so it yields no overlap."""
        assert "rs_far" not in {name for name, _overlap in sources.paired()}

    def test_regions_do_not_cross(self, sources: _Sources) -> None:
        """The chr2 variant pairs only with the chr2 transcript."""
        assert dict(sources.paired())["rs_chr2"].feature.stable_id == "transcript:ENST003"

    def test_a_variant_without_an_id_is_named_by_locus(self, sources: _Sources) -> None:
        """The locus name is the same shape VariantSelector.parse reads."""
        assert "chr1:1600:A>C" in {name for name, _overlap in sources.paired()}

    def test_a_multi_allelic_variant_names_each_allele(self, sources: _Sources) -> None:
        """One VCF ID covers both alleles, so the alternate is appended to keep names unique."""
        assert {"rs_multi:C", "rs_multi:G"} <= {name for name, _overlap in sources.paired()}

    def test_a_name_and_a_feature_together_identify_one_overlap(self, sources: _Sources) -> None:
        """The same uniqueness the generated source promises."""
        pairs = [(name, overlap.feature.stable_id) for name, overlap in sources.paired()]
        assert len(pairs) == len(set(pairs))

    def test_the_overlaps_are_transcript_variation_alleles(self, sources: _Sources) -> None:
        """What a FactContext is built from."""
        assert all(isinstance(overlap, TranscriptVariationAllele) for _name, overlap in sources.paired())

    def test_the_flank_is_honoured(self, sources: _Sources) -> None:
        """A smaller flank drops the variant that only reached its transcript through the flank."""
        assert "rs_flank" in {name for name, _overlap in sources.paired()}
        assert "rs_flank" not in {name for name, _overlap in sources.paired(flank=10)}

    def test_pairing_is_in_input_order(self, sources: _Sources) -> None:
        """``--index`` counts these, so the order is part of the contract."""
        assert [name for name, _overlap in sources.paired()][:2] == ["rs_inside", "rs_flank"]


class TestFinding:
    """Selecting one variant's overlaps without pairing every variant in the file."""

    def test_an_snv_is_found_by_locus(self, sources: _Sources) -> None:
        """The selector shape is the same one the pairing names its variants with."""
        assert [name for name, _overlap in sources.found(variant="chr1:1500:A>C")] == ["rs_inside"]

    def test_an_indel_is_found_by_locus(self, sources: _Sources) -> None:
        """The VCF reader spells an empty allele ``-`` and the corpus spells it ``""``. Finding
        must not care which, or every indel read from a file would be unselectable."""
        assert [name for name, _overlap in sources.found(variant="chr2:101:->AT")] == ["rs_chr2"]

    def test_a_name_still_finds_it(self, sources: _Sources) -> None:
        """A VCF ID is a name like any other."""
        found = sources.found(name="rs_inside")
        assert [name for name, _overlap in found] == ["rs_inside"]
        assert found[0][1].feature.stable_id == "transcript:ENST001"

    def test_an_id_selects_the_first_alternate_allele(self, sources: _Sources) -> None:
        """CHR:POS:ID names a record, and its first alternate allele is the default."""
        name, overlap = sources.found(variant="chr1:1700:rs_multi")[0]
        assert name == "rs_multi:C"
        assert overlap.overlapping_feature.alternate_allele == "C"

    def test_an_allele_number_selects_a_later_alternate_allele(self, sources: _Sources) -> None:
        """Alternate alleles are counted from 1, as VEP's ALLELE_NUM counts them."""
        assert sources.found(variant="chr1:1700:rs_multi", allele_num=2)[0][0] == "rs_multi:G"

    def test_an_allele_number_past_the_last_allele_finds_nothing(self, sources: _Sources) -> None:
        """rs_multi has two alternate alleles, so a third matches nothing."""
        assert not sources.found(variant="chr1:1700:rs_multi", allele_num=3)

    def test_an_id_at_another_position_finds_nothing(self, sources: _Sources) -> None:
        """The position is part of the selector, not a hint."""
        assert not sources.found(variant="chr1:1701:rs_multi")

    def test_an_id_matches_the_input_position_of_an_indel(self, sources: _Sources) -> None:
        """An ID selector uses the VCF POS, which for an insertion is the anchor base before the
        minimised start a REF>ALT selector uses."""
        assert sources.found(variant="chr2:100:rs_chr2")[0][0] == "rs_chr2"

    def test_finding_never_pairs_the_whole_file(
        self, sources: _Sources, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Not every transcript is built to explain one variant."""
        monkeypatch.setattr(GFF3TranscriptSource, "get_all_features", _refuse)

        name, overlap = sources.found(variant="chr1:1500:A>C")[0]

        assert name == "rs_inside"
        assert overlap.feature.stable_id == "transcript:ENST001"

    def test_an_unindexed_input_is_read_only_up_to_the_first_match(
        self, sources: _Sources, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """rs_inside is the first record, so nothing after it is read."""
        read: list[Variation] = []
        read_all = VCFVariationSource.iter_all_features

        def counting(reader: VCFVariationSource) -> Iterator[Variation]:
            for variation in read_all(reader):
                read.append(variation)
                yield variation

        monkeypatch.setattr(VCFVariationSource, "iter_all_features", counting)

        assert sources.found(variant="chr1:1500:rs_inside")[0][0] == "rs_inside"
        assert len(read) == 1

    def test_an_indexed_input_is_queried_by_region(
        self, sources: _Sources, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A selector's position is enough to ask the index, so the input is never read in order."""
        records = list(sources.variation_source.iter_all_features())
        queried: list[tuple[str | None, int | None, int | None]] = []

        def query(
            _reader: VCFVariationSource,
            *,
            seq_region_name: str | None = None,
            start: int | None = None,
            end: int | None = None,
        ) -> list[Variation]:
            queried.append((seq_region_name, start, end))
            return records

        monkeypatch.setattr(VCFVariationSource, "get_features_by_region", query)
        monkeypatch.setattr(VCFVariationSource, "iter_all_features", _refuse)

        assert sources.found(variant="chr2:101:->AT")[0][0] == "rs_chr2"
        assert queried == [("chr2", 100, 101)]

    def test_an_unindexed_transcript_file_is_streamed_for_the_window(
        self, sources: _Sources, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Without an index the transcripts are streamed, and only the variant's window is kept."""
        monkeypatch.setattr(GFF3TranscriptSource, "get_features_by_region", _unindexed)

        assert sources.found(variant="chr1:1500:A>C")[0][1].feature.stable_id == "transcript:ENST001"

    def test_an_indexed_transcript_file_is_queried_for_the_flank_window(
        self, sources: _Sources, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """The query covers the variant and its flank, clamped to the start of the region."""
        transcripts = sources.transcript_source.get_all_features()
        queried: list[tuple[str, int, int]] = []

        def query(
            _reader: GFF3TranscriptSource, *, seq_region_name: str, start: int, end: int
        ) -> list[Transcript]:
            queried.append((seq_region_name, start, end))
            return transcripts

        monkeypatch.setattr(GFF3TranscriptSource, "get_features_by_region", query)
        monkeypatch.setattr(GFF3TranscriptSource, "iter_features_by_region", _refuse)

        assert sources.found(variant="chr1:1500:A>C")[0][1].feature.stable_id == "transcript:ENST001"
        assert queried == [("chr1", 1, 1500 + DEFAULT_FLANK)]

    def test_a_name_not_in_the_file_finds_nothing(self, sources: _Sources) -> None:
        """Nothing stops the read early, and the miss is reported by the caller, not raised here."""
        assert not sources.found(name="rs_absent")

    @pytest.mark.integration
    def test_indexed_files_find_the_same_overlap(
        self, tmp_path: Path, htslib_toolchain: _HtslibToolchain
    ) -> None:
        """Real bgzipped, tabix-indexed inputs answer through the region queries."""
        variants = tmp_path / "variants.vcf"
        variant_lines = VCF.splitlines()
        body = sorted(variant_lines[2:], key=lambda line: (line.split("\t")[0], int(line.split("\t")[1])))
        variants.write_text("\n".join([*variant_lines[:2], *body]) + "\n", encoding="utf-8")
        transcripts = tmp_path / "transcripts.gff3"
        transcripts.write_text(GFF3, encoding="utf-8")
        for path, preset in ((variants, "vcf"), (transcripts, "gff")):
            subprocess.run([htslib_toolchain.bgzip, "-f", str(path)], check=True)
            subprocess.run([htslib_toolchain.tabix, "-f", "-p", preset, f"{path}.gz"], check=True)

        sources = _Sources(Path(f"{variants}.gz"), Path(f"{transcripts}.gz"))

        name, overlap = sources.found(variant="chr1:1700:rs_multi", allele_num=2)[0]

        assert name == "rs_multi:G"
        assert overlap.feature.stable_id == "transcript:ENST001"


class TestFeatureSelector:
    """Naming a transcript with or without the type prefix a GFF3 gives it."""

    def test_either_spelling_names_the_same_feature(self) -> None:
        """The GFF3 id is ``transcript:ENST001``; a caller typing ``ENST001`` means the same one."""
        feature = Transcript(start=1, end=10, strand=Strand.FORWARD, stable_id="transcript:ENST001")

        assert FeatureSelector.parse("ENST001").matches(feature)
        assert FeatureSelector.parse("transcript:ENST001").matches(feature)

    def test_another_feature_is_not_named(self) -> None:
        """A prefix match is not a match: the id after the prefix has to be the same one."""
        feature = Transcript(start=1, end=10, strand=Strand.FORWARD, stable_id="transcript:ENST001")

        assert not FeatureSelector.parse("ENST002").matches(feature)
        assert not FeatureSelector.parse("ENST00").matches(feature)

    def test_an_empty_selector_is_rejected(self) -> None:
        """A selector has to name something, so the empty string is the caller's mistake."""
        with pytest.raises(ValueError, match="feature stable id"):
            FeatureSelector.parse("")


class TestVariantSelector:
    """Reading the ``CHR:POS:REF>ALT`` and ``CHR:POS:ID`` selectors."""

    def test_a_selector_reads_a_dash_as_empty(self) -> None:
        """VEP spells an empty allele as a dash, so the selector does too."""
        insertion = VariantSelector.parse("chr1:10:->acg")
        deletion = VariantSelector.parse("chr1:10:A>-")
        assert (insertion.seq_region_name, insertion.position) == ("chr1", 10)
        assert (insertion.reference_allele, insertion.alternate_allele) == ("", "ACG")
        assert (deletion.reference_allele, deletion.alternate_allele) == ("A", "")

    def test_a_selector_without_alleles_is_an_identifier(self) -> None:
        """CHR:POS:ID names a record by the input's own identifier, and its first allele by default."""
        selector = VariantSelector.parse("22:15528194:rs123")
        assert selector.identifier == "rs123"
        assert selector.position == 15528194
        assert selector.allele_num == 1
        assert selector.reference_allele is None

    def test_an_identifier_may_contain_a_colon(self) -> None:
        """Only the first two colons separate parts, so an ID keeps its own."""
        assert VariantSelector.parse("chr1:10:chr1:10:A:T").identifier == "chr1:10:A:T"

    @pytest.mark.parametrize(
        ("variant", "allele_num", "message"),
        [
            pytest.param("chr1:10:A>T", 2, "names its allele", id="alleles-already-named"),
            pytest.param("chr1:10:rs1", 0, "from 1", id="counted-from-one"),
            pytest.param("chr1:10:", None, "CHR:POS:ID", id="empty-detail"),
            pytest.param("chr1:x:A>T", None, "integer position", id="non-numeric-position"),
            pytest.param("chr1:10", None, "CHR:POS:ID", id="missing-detail"),
        ],
    )
    def test_an_unusable_selector_is_rejected(
        self, variant: str, allele_num: int | None, message: str
    ) -> None:
        """A selector's shape is fixed, so anything else is reported rather than guessed at."""
        with pytest.raises(ValueError, match=message):
            VariantSelector.parse(variant, allele_num=allele_num)

    def test_the_query_span_covers_an_anchor_base(self) -> None:
        """A VCF anchors an indel one base before its minimised start."""
        assert VariantSelector.parse("chr1:101:->AT").query_span == (100, 101)
        assert VariantSelector.parse("chr1:1:->AT").query_span == (1, 1)


class TestFileSources:
    """Reaching real files through the run path's own builders."""

    @staticmethod
    def _selection(**options: object) -> OverlapSelection:
        """Return a selection over the given options, building its own readers."""
        return OverlapSelection(config=Config.for_command("explain", **options))

    def test_it_builds_both_sources_from_the_options(self, tmp_path: Path) -> None:
        """The seam explain reads its input through."""
        variants, transcripts = _write(tmp_path)
        fasta = tmp_path / "ref.fa"
        fasta.write_text(">chr1\n" + "N" * 25000 + "\n>chr2\n" + "N" * 2000 + "\n", encoding="utf-8")

        selection = self._selection(input_file=variants, gff=transcripts, fasta=fasta)

        paired = selection.all_overlaps()
        assert {name for name, _overlap in paired} >= {"rs_inside", "rs_chr2"}

    def test_a_missing_transcript_file_is_rejected_before_reading(self, tmp_path: Path) -> None:
        """Failing at build time says which option was missing; failing later would not."""
        variants, _transcripts = _write(tmp_path)

        with pytest.raises(ValueError, match="--gff or --gtf"):
            list(self._selection(input_file=variants).all_overlaps())

    def test_no_reference_is_allowed_but_reported(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """A run refuses a missing FASTA; explaining goes on without it and says what is lost.

        Only the build is exercised: reading an overlap without a reference cannot resolve a
        sequence region's length, which is exactly what the warning is about.
        """
        variants, transcripts = _write(tmp_path)
        selection = self._selection(input_file=variants, gff=transcripts)

        with caplog.at_level("WARNING"):
            selection._file_sources()  # pylint: disable=protected-access

        assert any("no --fasta given" in record.getMessage() for record in caplog.records)

    def test_the_readers_are_built_once(self, tmp_path: Path) -> None:
        """Both files are opened once per run, however many times the selection is asked."""
        variants, transcripts = _write(tmp_path)
        fasta = tmp_path / "ref.fa"
        fasta.write_text(">chr1\n" + "N" * 25000 + "\n>chr2\n" + "N" * 2000 + "\n", encoding="utf-8")
        selection = self._selection(input_file=variants, gff=transcripts, fasta=fasta)

        list(selection.all_overlaps())
        built = selection.sources
        list(selection.all_overlaps())

        assert built is not None
        assert selection.sources is built
