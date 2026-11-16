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
"""Tests for how VCF spells alleles."""

from __future__ import annotations

from ensembl.vep.common.vcf import trim_right, vcf_record_alleles

SEQUENCE = "ACGTACGTAC"


def _base_at(position: int) -> str:
    """Return the base at a 1-based position of ``SEQUENCE``."""
    return SEQUENCE[position - 1]


class TestTrimRight:
    """The trailing-base trim, as ensembl-variation does it."""

    def test_shared_trailing_bases_are_removed(self) -> None:
        """Bases every allele ends with go, down to one base."""
        assert trim_right(["ACG", "AG"]) == ["AC", "A"]

    def test_it_never_trims_an_allele_away(self) -> None:
        """Trimming stops at a one-base allele."""
        assert trim_right(["AG", "G"]) == ["AG", "G"]

    def test_different_last_bases_stop_it(self) -> None:
        """Nothing shared, nothing trimmed."""
        assert trim_right(["ACG", "ACT"]) == ["ACG", "ACT"]

    def test_a_single_allele_is_returned_as_is(self) -> None:
        """There is nothing to compare against."""
        assert trim_right(["ACG"]) == ["ACG"]


class TestVcfRecordAlleles:
    """The padding rules."""

    def test_same_length_alleles_are_written_as_they_are(self) -> None:
        """SNVs and MNVs need no padding."""
        assert vcf_record_alleles(3, ["G", "T"], _base_at) == (3, ["G", "T"])
        assert vcf_record_alleles(3, ["GT", "CA"], _base_at) == (3, ["GT", "CA"])

    def test_a_deletion_is_padded_with_the_base_before(self) -> None:
        """The anchor moves the position back by one."""
        assert vcf_record_alleles(3, ["GT", "-"], _base_at) == (2, ["CGT", "C"])

    def test_an_insertion_is_padded_with_the_base_before_its_position(self) -> None:
        """An insertion before position 3 sits after base 2."""
        assert vcf_record_alleles(3, ["", "GG"], _base_at) == (2, ["C", "CGG"])

    def test_at_position_one_the_following_base_pads(self) -> None:
        """The spec's rule at the start of a sequence; the one deliberate difference from Perl."""
        assert vcf_record_alleles(1, ["A", ""], _base_at) == (1, ["AC", "C"])
        assert vcf_record_alleles(1, ["", "TT"], _base_at) == (1, ["A", "TTA"])

    def test_alleles_sharing_a_first_base_are_not_padded(self) -> None:
        """After trimming, ACG>AG is AC>A, which already starts with a shared base."""
        assert vcf_record_alleles(5, ["ACG", "AG"], _base_at) == (5, ["AC", "A"])

    def test_multiple_alternates_are_padded_together(self) -> None:
        """One anchor for the whole record."""
        assert vcf_record_alleles(3, ["G", "T", "GA"], _base_at) == (2, ["CG", "CT", "CGA"])

    def test_a_symbolic_allele_is_bracketed_and_not_trimmed(self) -> None:
        """Symbolic alleles are named, never padded with a base."""
        assert vcf_record_alleles(3, ["G", "DEL"], _base_at) == (2, ["CG", "<DEL>"])

    def test_the_base_lookup_is_only_used_when_padding(self) -> None:
        """An SNV never reads the sequence."""

        def refuse(_position: int) -> str:
            raise AssertionError("no padding base should be needed")

        assert vcf_record_alleles(3, ["G", "T"], refuse) == (3, ["G", "T"])
