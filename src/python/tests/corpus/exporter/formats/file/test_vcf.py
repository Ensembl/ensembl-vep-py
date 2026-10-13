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
"""Tests for VCF corpus export."""

from __future__ import annotations

from pathlib import Path

import pytest

from ensembl.vep.common import CoordSystemLevel, Strand
from ensembl.vep.corpus import Corpus
from ensembl.vep.corpus.exporter.formats.file.vcf import VCF_NAME, distinct_variants, vcf_alleles
from ensembl.vep.corpus.exporter.manifest import CorpusManifest
from ensembl.vep.models.sequence import CoordSystem, Slice
from ensembl.vep.models.variation import VariationAllele
from ensembl.vep.parsers.vcf import VCFParser


@pytest.fixture(name="records", scope="module")
def fixture_records(saved: tuple[Path, CorpusManifest]) -> list:
    """Parse exported VCF records once through the production parser."""
    directory, _manifest = saved
    return VCFParser.parse_file(directory / VCF_NAME, limit=None)


def _retrim(position: int, reference: str, alternate: str) -> tuple[int, str, str]:
    """Independently undo the anchor padding used by VCF indels."""
    if len(reference) == len(alternate) == 1:
        return position, reference, alternate
    if reference[0] == alternate[0]:
        return position + 1, reference[1:], alternate[1:]
    if reference[-1] == alternate[-1]:
        return position, reference[:-1], alternate[:-1]
    return position, reference, alternate


def _slice_with(sequence: str) -> Slice:
    """Build a slice carrying the supplied known reference sequence."""
    sequence_region = Slice(
        coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
        seq_region_name="chrTest",
        start=1,
        end=len(sequence),
        strand=Strand.FORWARD,
    )
    with sequence_region.unsafe() as unsafe_slice:
        setattr(unsafe_slice, "_seq", sequence)
    return sequence_region


class TestVcf:
    """Variant records exported as VCF."""

    def test_writes_each_distinct_variant_once(self, source: Corpus, records: list) -> None:
        """A variant overlapping several transcripts remains one VCF record."""
        expected = {name for name, _allele in distinct_variants(source)}
        written = [record.id[0] for record in records]
        assert len(written) == len(expected)
        assert set(written) == expected

    def test_records_retrim_to_the_model_alleles(self, source: Corpus, records: list) -> None:
        """VCF anchor padding reverses exactly to the model's trimmed alleles."""
        alleles_by_name = dict(distinct_variants(source))
        for record in records:
            allele = alleles_by_name[record.id[0]]
            assert _retrim(record.pos, record.reference, record.alt[0]) == (
                allele.start,
                allele.reference_allele or "",
                allele.alternate_allele or "",
            )

    def test_records_are_position_sorted_within_each_sequence_region(self, records: list) -> None:
        """The order meets the requirement of indexed VCF consumers."""
        for chromosome in {record.chrom for record in records}:
            positions = [record.pos for record in records if record.chrom == chromosome]
            assert positions == sorted(positions)


class TestVcfAlleles:
    """The padding rule for individual VCF alleles."""

    def test_snv_needs_no_padding(self) -> None:
        """A single-base substitution retains its own coordinate and alleles."""
        sequence_region = _slice_with("ACGTACGT")
        allele = VariationAllele(
            start=3,
            end=3,
            strand=Strand.FORWARD,
            slice=sequence_region,
            reference_allele="G",
            alternate_allele="T",
        )
        assert vcf_alleles(allele, "ACGTACGT") == (3, "G", "T")

    def test_indels_are_anchored_before_the_change(self) -> None:
        """Deletions and insertions use the preceding reference base as their anchor."""
        sequence_region = _slice_with("ACGTACGT")
        deletion = VariationAllele(
            start=3,
            end=4,
            strand=Strand.FORWARD,
            slice=sequence_region,
            reference_allele="GT",
            alternate_allele="",
        )
        insertion = VariationAllele(
            start=3,
            end=2,
            strand=Strand.FORWARD,
            slice=sequence_region,
            reference_allele="",
            alternate_allele="GG",
        )
        assert vcf_alleles(deletion, "ACGTACGT") == (2, "CGT", "C")
        assert vcf_alleles(insertion, "ACGTACGT") == (2, "C", "CGG")

    def test_position_one_indels_are_anchored_after_the_change(self) -> None:
        """The first base uses right-padding because no preceding anchor exists."""
        sequence_region = _slice_with("ACGTACGT")
        deletion = VariationAllele(
            start=1,
            end=1,
            strand=Strand.FORWARD,
            slice=sequence_region,
            reference_allele="A",
            alternate_allele="",
        )
        insertion = VariationAllele(
            start=1,
            end=0,
            strand=Strand.FORWARD,
            slice=sequence_region,
            reference_allele="",
            alternate_allele="TT",
        )
        assert vcf_alleles(deletion, "ACGTACGT") == (1, "AC", "C")
        assert vcf_alleles(insertion, "ACGTACGT") == (1, "A", "TTA")
