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
"""VCF-back Variation feature source"""

import logging
import re
from attrs import define, field

from ensembl.vep.common import Locus, Strand
from ensembl.vep.models.sequence import Slice
from ensembl.vep.models.variation.variation import Variation
from ensembl.vep.models.variation.structural_variation import StructuralVariation, StructuralVariationAllele
from ensembl.vep.parsers.vcf import VCFParser
from ensembl.vep.sources.file.file_feature_source import FileFeatureSource
from ensembl.vep.ensembl_vep_core import VCFRecord  # noqa: F401,F403 pylint: disable=E0401,E0611

logger = logging.getLogger(__name__)

BREAKEND_REGEX = re.compile(
    r"""\[[^*]+\[ # one or more non-asterisk characters inside opening square brackets
    | \][^*]+\]   # OR one or more non-asterisk characters inside closing square brackets
    | \.[A-Z]+    # OR a dot followed by base characters
    | [A-Z]+\.    # OR base characters followed by a dot""",
    re.X,
)
SYMBOLIC_ALLELES = ("DEL", "INS", "DUP", "INV", "CNV")
SVTYPES = ("BND", *SYMBOLIC_ALLELES)


@define(kw_only=True, auto_attribs=True)
class VCFVariationSource(FileFeatureSource[VCFParser, Variation]):
    """Variation source backed by a VCF file.

    Attributes:
        adaptor: The parser adaptor used to stream VCF records.
    """

    adaptor: VCFParser = field(init=False, factory=VCFParser)

    def _records_to_features(self, records: list[VCFRecord]) -> list[Variation]:
        """Convert VCFRecords into Variation or StructuralVariation objects.

        Args:
            records (list[VCFRecord]): VCF records to convert into variation features.

        Returns:
            list[Variation]: Variation and StructuralVariation objects created from the VCF records.
        """

        features: list[Variation] = []
        for record in records:
            if self._is_structural_variation_record(record):
                features.extend(self._records_to_structural_variation([record]))
            else:
                features.extend(self._records_to_variation([record]))

        return features

    def _records_to_variation(self, records: list[VCFRecord]) -> list[Variation]:
        """Convert VCFRecords into Variation objects.

        Args:
            records (list[VCFRecord]): VCF records to convert into Variation objects.

        Returns:
            list[Variation]: Variation objects created from the VCF records.
        """

        variations: list[Variation] = []
        for record in records:
            uploaded_allele_string = "/".join([record.reference, *record.alt])
            uploaded_end = record.pos + len(record.reference) - 1
            ref, alts, start, end, _ = self._remove_anchor_bases(
                record.reference,
                record.alt,
                start=record.pos,
                end=uploaded_end,
            )
            assert isinstance(alts, list)
            allele_string = "/".join([ref, *alts])
            variations.append(
                Variation(
                    stable_id=self._stable_id_from_vcf_record(record),
                    start=start,
                    end=end,
                    strand=Strand.FORWARD,
                    slice=self.slice_cache.get_slice(
                        record.chrom,
                        self.sequence_source,
                    ),
                    allele_string=allele_string,
                    uploaded_locus=Locus(
                        start=record.pos,
                        end=uploaded_end,
                        allele_string=uploaded_allele_string,
                    ),
                )
            )
        return variations

    def _records_to_structural_variation(self, records: list[VCFRecord]) -> list[StructuralVariation]:
        """Convert VCFRecords into StructuralVariation objects.

        Args:
            records (list[VCFRecord]): VCF records to convert into StructuralVariation objects.

        Returns:
            list[StructuralVariation]: StructuralVariation objects created from the VCF records.
        """

        structural_variations: list[StructuralVariation] = []
        for record in records:
            seq_region_slice = self.slice_cache.get_slice(
                record.chrom,
                self.sequence_source,
            )

            variation_alleles = []
            for i in range(len(record.alt)):
                variation_alleles.append(
                    self._create_sv_allele_from_record(
                        record, allele_index=i, seq_region_slice=seq_region_slice
                    )
                )
            start = record.pos
            end = max(var.end for var in variation_alleles)

            structural_variation = StructuralVariation(
                stable_id=self._stable_id_from_vcf_record(record),
                start=start,
                end=end,
                strand=Strand.FORWARD,
                slice=seq_region_slice,
                allele_string="/".join([record.reference, *record.alt]),
                uploaded_locus=Locus(
                    start=record.pos,
                    end=end,
                    allele_string="/".join([record.reference, *record.alt]),
                ),
                copy_number=self._int_info_value(record, "CN"),
            )

            structural_variation.variation_alleles = variation_alleles

            structural_variations.append(structural_variation)

        return structural_variations

    def _create_sv_allele_from_record(
        self, record: VCFRecord, *, allele_index: int, seq_region_slice: Slice | None
    ) -> StructuralVariationAllele:
        """Create a StructuralVariationAllele object from a VCFRecord.

        Args:
            record (VCFRecord): VCF record to convert into StructuralVariation objects.
            allele_index (int): index of the alt allele to create object for.
            slice (Slice | None): Slice object to add to StructuralVariationAllele.

        Returns:
            StructuralVariationAllele: StructuralVariationAllele object created from the VCF record.
        """
        end, svlen = self._end_and_svlen_for_sva_from_vcf_record(record, allele_index)
        ref, alt, start, end, _ = self._remove_anchor_bases_symbolic_allele(
            record.reference, record.alt[allele_index], start=record.pos, end=end
        )

        return StructuralVariationAllele(
            start=start,
            end=end,
            svlen=svlen,
            strand=Strand.FORWARD,
            slice=seq_region_slice,
            reference_allele=ref,
            alternate_allele=alt,
        )

    def _is_structural_variation_record(self, record: VCFRecord) -> bool:
        """Return whether the VCF record should be represented as a StructuralVariation.

        Args:
            record (VCFRecord): VCF record to inspect for structural-variation hints.

        Returns:
            bool: True if the record should be converted to a StructuralVariation, otherwise False.
        """

        return not self._is_bases_only(record) and (
            self._has_svtype(record) or self._is_symbolic(record) or self._is_breakend(record)
        )

    @staticmethod
    def _is_bases_only(record: VCFRecord) -> bool:
        """Is the VCF record a SNV or short INDEL

        Args:
            record (VCFRecord): VCF record to inspect for structural-variation hints.

        Returns:
            bool: True if the record only contains nucleotide bases in REF and ALT, otherwise False.
        """

        return all(re.fullmatch(r"^[ACGT]+$", allele) for allele in [record.reference, *record.alt])

    @staticmethod
    def _has_svtype(record: VCFRecord) -> bool:
        """Does the record have the SVTYPE info key

        Args:
            record (VCFRecord): VCF record to inspect for structural-variation hints.

        Returns:
            bool: True if the record has a valid value for SVTYPE, otherwise False.
        """

        svtype = record.info.get("SVTYPE", "").upper()
        return svtype in SVTYPES

    @staticmethod
    def _is_symbolic(record: VCFRecord) -> bool:
        """Do any of the ALT alleles match a symbolic allele

        Args:
            record (VCFRecord): VCF record to inspect for structural-variation hints.

        Returns:
            bool: True if ANY of ALT alleles match a symbolic allele, otherwise False.
        """

        # Remove opening angle bracket for comparison
        alts = [allele.replace("<", "") for allele in record.alt]
        return any(allele[0:3] in SYMBOLIC_ALLELES for allele in alts)

    @staticmethod
    def _is_breakend(record: VCFRecord) -> bool:
        """Do any of the ALT alleles match a breakend allele

        Args:
            record (VCFRecord): VCF record to inspect for structural-variation hints.

        Returns:
            bool: True if any of ALT alleles match breakend notation, otherwise False.
        """

        return any(re.search(BREAKEND_REGEX, allele) for allele in record.alt)

    @staticmethod
    def _is_symbolic_insertion(alt_allele_string: str) -> bool:
        """Is the symbolic allele an insertion

        Args:
            alt_allele_string (str): Alternate allele string

        Returns:
            bool: True if the alt string matches an insertion allele.
        """

        return alt_allele_string in ("INS", "<INS>")

    @staticmethod
    def _stable_id_from_vcf_record(record: VCFRecord) -> str | None:
        """Return a stable identifier from the VCF ID column, when present.

        Args:
            record (VCFRecord): VCF record containing ID values.

        Returns:
            str | None: Semicolon-delimited stable identifier, or None when the record has no IDs.
        """

        return ";".join(record.id) or None

    def _end_and_svlen_for_sva_from_vcf_record(
        self, record: VCFRecord, allele_index: int
    ) -> tuple[int | None, int | None]:
        """Return the VCF record end coordinate and length for an allele.

        Args:
            record (VCFRecord): VCF record.
            allele_index: index for the allele.

        Returns:
            tuple: (End coordinate, SV length).
        """

        if allele_index >= len(record.alt):
            logger.warning(f"Allele index, {allele_index}, is too large, {';'.join(record.id)}")
            return (None, None)

        lengths = record.info.get("SVLEN")
        if lengths is not None:
            svlen = lengths.split(",")[allele_index]
            svlen = None if svlen == "." else abs(int(svlen))
            if svlen is None:
                end = record.pos + len(record.reference) - 1
            else:
                end = record.pos + svlen
        else:
            end = self._int_info_value(record, "END")
            if end is None:
                svlen = None
            else:
                svlen = abs(end - record.pos)

        # For insertions, end is not POS + SVLEN
        if self._is_symbolic_insertion(record.alt[allele_index]):
            end = record.pos + len(record.reference) - 1

        if end is None or svlen is None:
            logger.warning(f"Could not determine END/SVLEN for {record.id}")

        return end, svlen

    @staticmethod
    def _int_info_value(record: VCFRecord, key: str) -> int | None:
        """Return an INFO value parsed as an integer.

        Args:
            record (VCFRecord): VCF record containing INFO values.
            key (str): INFO field key to parse.

        Returns:
            int | None: Parsed integer value, or None when the key is absent or not an integer.
        """

        value = record.info.get(key)
        if value is None:
            return None

        try:
            return int(value)
        except (TypeError, ValueError):
            return None

    @staticmethod
    def _remove_anchor_bases(
        ref: str,
        alt: str | list[str],
        *,
        start: int = 0,
        end: int | None = None,
        strand: Strand = Strand.FORWARD,
    ) -> tuple[str, str | list[str], int, int, bool]:
        """Remove common leading anchor bases from non-symbolic alleles.

        For multi-allelic variants, only bases shared by the reference and every alternate allele are
        removed so that all alleles retain a common locus.

        Args:
            ref: Reference allele sequence.
            alt: One alternate allele or a list of alternate alleles.
            start: Start coordinate.
            end: End coordinate. Defaults to ``start + len(ref) - 1``.
            strand: Strand on which the alleles are represented.

        Returns:
            Tuple of (new_ref, new_alt, new_start, new_end, changed). Empty alleles are represented by
            ``-`` and the shape of ``new_alt`` matches the input.
        """

        end = end if end is not None else start + len(ref) - 1

        # Alternate alleles can be a single string or a list of strings. Normalize to a list for processing.
        if isinstance(alt, str):
            alt_is_string = True
            alt_alleles = [alt]
        else:
            alt_is_string = False
            alt_alleles = alt.copy()
        if not alt_alleles:
            return (ref, alt, start, end, False)

        if strand == Strand.UNKNOWN:
            logger.warning(
                "Strand is unknown, cannot remove anchor bases with confidence. "
                "Returning unmodified alleles."
            )
            return (ref, alt, start, end, False)

        alleles = [ref, *alt_alleles]
        changed = False

        while all(alleles) and len({allele[0] for allele in alleles}) == 1:
            alleles = [allele[1:] for allele in alleles]
            if strand == Strand.REVERSE:
                end -= 1
            else:
                start += 1
            changed = True

        new_ref = alleles[0] or "-"
        new_alts = [allele or "-" for allele in alleles[1:]]
        new_alt: str | list[str] = new_alts[0] if alt_is_string else new_alts
        return (new_ref, new_alt, start, end, changed)

    @staticmethod
    def _remove_anchor_bases_symbolic_allele(
        ref: str,
        alt: str,
        *,
        start: int = 0,
        end: int | None,
        strand: Strand = Strand.FORWARD,
    ) -> tuple[str, str, int, int, bool]:
        """Remove the reference anchor bases from a symbolic allele.

        Converts symbolic allele in the form A/<DEL> to -/<DEL>, adjusting coordinates accordingly.

        Args:
            ref: Reference allele sequence.
            alt: Alternate allele sequence.
            start: Start coordinate (default 0).
            end: End coordinate (default: start + len(ref) - 1).
            strand: 1 for forward, -1 for reverse strand.

        Returns:
            Tuple of (new_ref, new_alt, new_start, new_end, changed).
        """

        end = end or (start + (len(ref) - 1))

        changed = False

        if strand == Strand.UNKNOWN:
            logger.warning(
                "Strand is unknown, cannot trim alleles with confidence. Returning untrimmed alleles."
            )
            return (ref, alt, start, end, changed)

        if ref == "-":
            logger.warning("Allele is already trimmed. Returning untrimmed alleles.")
            return (ref, alt, start, end, changed)

        if strand == Strand.REVERSE:
            end -= len(ref)
        else:
            start += len(ref)
        ref = "-"
        changed = True

        return (ref, alt, start, end, changed)
