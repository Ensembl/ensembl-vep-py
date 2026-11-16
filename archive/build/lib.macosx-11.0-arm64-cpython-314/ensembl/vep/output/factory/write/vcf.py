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
"""VCF output: the input line back, with a ``CSQ`` INFO field of pipe-separated consequences.

A port of Perl VEP's ``OutputFactory/VCF.pm``. When the variant came from a VCF, its original line is
written back with any previous ``CSQ`` replaced; otherwise a record is built from the variation.
Each (allele, feature) becomes one comma-separated chunk whose fields follow the ``Format`` the
header declares -- which is exactly what ``output/parser/vcf.py`` reads back.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Iterator, Sequence

from attrs import define, field

from ensembl.vep.common.methods import reverse_complement
from ensembl.vep.common.vcf import EMPTY_ALLELE, vcf_record_alleles
from ensembl.vep.models.variation import Variation
from ensembl.vep.output.factory.base import AnnotatedVariant
from ensembl.vep.output.factory.write.base import OutputFactory

__all__ = ["DEFAULT_INFO_FIELD", "VCF_COLS", "VCFOutputFactory", "VCFOutputWriter"]

VCF_COLS = (
    "Allele",
    "Consequence",
    "IMPACT",
    "SYMBOL",
    "Gene",
    "Feature_type",
    "Feature",
    "BIOTYPE",
    "EXON",
    "INTRON",
    "HGVSc",
    "HGVSp",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "Existing_variation",
)
"""The ``CSQ`` fields VEP writes by default, in its order."""

DEFAULT_INFO_FIELD = "CSQ"
FILEFORMAT = "##fileformat=VCFv4.1"
COLUMN_HEADING = "#" + "\t".join(("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))
_INFO_COLUMN = 7


@define(kw_only=True, auto_attribs=True)
class VCFOutputFactory(OutputFactory):
    """Writes annotated variations as VCF.

    Attributes:
        info_field (str): The INFO key consequences are written under.
        keep_csq (bool): Keep a ``CSQ`` already on the input line, and its header, instead of
            replacing them.
    """

    info_field: str = field(default=DEFAULT_INFO_FIELD)
    keep_csq: bool = field(default=False)

    def fields(self) -> list[str]:
        """Return the fields each chunk carries: VEP's defaults, then what the flags add."""
        fields = list(VCF_COLS)
        for name in self.flag_fields():
            if name not in fields:
                fields.append(name)
        return fields

    def headers(self, input_headers: Sequence[str] = ()) -> list[str]:
        """Return the header lines, carrying the input file's own when there are any.

        Args:
            input_headers (Sequence[str]): The input VCF's header lines, column heading last.

        Returns:
            list[str]: The meta lines, then VEP's own, then the column heading.
        """
        headers: list[str] = []
        if input_headers:
            if not any(line.lower().startswith("##fileformat=vcfv4") for line in input_headers):
                headers.append(FILEFORMAT)
            info_declaration = re.compile(rf"ID={re.escape(self.info_field)},", re.IGNORECASE)
            for line in input_headers:
                if self.keep_csq or not (info_declaration.search(line) or line.startswith("##VEP")):
                    headers.append(line)
            column_heading = headers.pop()
        else:
            headers.append(FILEFORMAT)
            column_heading = COLUMN_HEADING

        headers.append(f'##VEP="v{self.version}" time="{self.timestamp()}"')
        headers.append(
            f"##INFO=<ID={self.info_field},Number=.,Type=String,"
            f'Description="Consequence annotations from Ensembl VEP. Format: {"|".join(self.fields())}">'
        )
        if self.command_line:
            headers.append(f"##VEP-command-line='{self.command_line}'")
        headers.append(column_heading)
        return headers

    def info_chunk(self, output_hash: dict[str, object], strand: int = 1) -> str:
        """Encode one output hash as a pipe-separated chunk.

        The escaping is VEP's: the characters VCF reserves inside INFO -- ``,`` ``;`` ``|`` and
        whitespace -- are replaced so a chunk can never split, and ``-`` means empty in every field
        but ``Allele``. A zero is kept rather than blanked.
        """
        chunk = []
        for column in self.fields():
            value = output_hash.get(column)
            if value is None or value == "":
                chunk.append("")
                continue
            if not isinstance(value, (list, tuple)) and str(value) == "0":
                chunk.append("0")
                continue

            data = "&".join(str(item) for item in value) if isinstance(value, (list, tuple)) else str(value)
            if column == "Allele":
                if strand < 0 and data != EMPTY_ALLELE:
                    data = reverse_complement(data)
            elif data == EMPTY_ALLELE:
                data = ""
            if column == "SYMBOL":
                data = re.sub(r"\s", "%20", data)
            data = data.replace(",", "&").replace(";", "%3B")
            data = re.sub(r"\s+", "_", data)
            chunk.append(data.replace("|", "&"))
        return "|".join(chunk)

    def vcf_record(self, variation: Variation) -> list[str]:
        """Build the eight VCF columns for a variation that did not come from a VCF line."""
        chromosome = variation.slice
        sequence_holder: dict[str, str] = {}

        def base_at(position: int) -> str:
            if chromosome is None:
                return "N"
            if "seq" not in sequence_holder:
                sequence_holder["seq"] = chromosome.seq or ""
            sequence = sequence_holder["seq"]
            return sequence[position - 1] if 0 < position <= len(sequence) else "N"

        alleles = (variation.allele_string or "").split("/")
        position, spelled = vcf_record_alleles(variation.start, alleles, base_at)
        return [
            chromosome.seq_region_name if chromosome is not None else ".",
            str(position),
            variation.stable_id or ".",
            spelled[0] if spelled else ".",
            ",".join(spelled[1:]) or ".",
            ".",
            ".",
            ".",
        ]

    def line(self, variation: Variation, record_columns: Sequence[str] | None = None) -> str:
        """Return one variation's output line.

        Args:
            variation (Variation): The annotated variation.
            record_columns (Sequence[str] | None): Its original VCF columns, written back with the
                consequences added; or None to build a record from the variation.

        Returns:
            str: The tab-separated line.
        """
        if record_columns is not None:
            columns = self.prepare_input_columns(record_columns)
            info = "" if columns[_INFO_COLUMN] in {"", "."} else columns[_INFO_COLUMN]
        else:
            columns = self.vcf_record(variation)
            info = ""

        strand = int(variation.strand.value)
        chunks = [
            self.info_chunk(output_hash, strand)
            for output_hash in self.output_hashes(variation, record_columns)
        ]
        if chunks:
            info = f"{info};" if info else ""
            info += f"{self.info_field}={','.join(chunks)}"
        columns[_INFO_COLUMN] = info or "."
        return "\t".join(columns)

    def prepare_input_columns(self, record_columns: Sequence[str]) -> list[str]:
        """Copy input columns and remove an annotation field scheduled for replacement."""
        columns = [str(column) for column in record_columns]
        while len(columns) <= _INFO_COLUMN:
            columns.append(".")
        info = "" if columns[_INFO_COLUMN] in {"", "."} else columns[_INFO_COLUMN]
        key = re.escape(self.info_field)
        if not self.keep_csq and re.search(rf"(^|;){key}=", info):
            info = re.sub(rf"(^|;){key}=\S+?(;\S|$)", r"\2", info, count=1)
            info = info.removeprefix(";")
        columns[_INFO_COLUMN] = info or "."
        return columns

    def lines(
        self, variations: Iterable[Variation], records: Iterable[Sequence[str] | None] | None = None
    ) -> Iterator[str]:
        """Yield one line per variation, pairing each with its original record when given."""
        if records is None:
            for variation in variations:
                yield self.line(variation)
            return
        for variation, record_columns in zip(variations, records, strict=True):
            yield self.line(variation, record_columns)


@define(kw_only=True, auto_attribs=True)
class VCFOutputWriter:
    """Write annotated-variant snapshots as VCF records.

    Attributes:
        factory (VCFOutputFactory): The VCF formatter that supplies headers and CSQ encoding.
    """

    factory: VCFOutputFactory = field(factory=VCFOutputFactory)

    def headers(self, input_headers: Sequence[str] = ()) -> list[str]:
        """Return VCF headers for this writer."""
        return self.factory.headers(input_headers)

    def line(self, variant: AnnotatedVariant) -> str:
        """Return one annotated VCF record."""
        columns = self._columns(variant)
        info = "" if columns[_INFO_COLUMN] in {"", "."} else columns[_INFO_COLUMN]
        chunks = [self.factory.info_chunk(annotation, variant.strand) for annotation in variant.annotations]
        if chunks:
            info = f"{info};" if info else ""
            info += f"{self.factory.info_field}={','.join(chunks)}"
        columns[_INFO_COLUMN] = info or "."
        return "\t".join(columns)

    def lines(self, variants: Iterable[AnnotatedVariant]) -> Iterator[str]:
        """Yield one VCF record per annotated variant."""
        for variant in variants:
            yield self.line(variant)

    def _columns(self, variant: AnnotatedVariant) -> list[str]:
        """Return input columns ready for a replacement annotation field."""
        if variant.input_columns:
            return self.factory.prepare_input_columns(variant.input_columns)
        raise ValueError("VCF output requires the original input columns")
