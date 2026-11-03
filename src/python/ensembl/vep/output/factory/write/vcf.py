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

A port of Perl VEP's ``OutputFactory/VCF.pm``. A snapshot's original VCF line is written back with
any previous ``CSQ`` replaced. Each selected (allele, feature) record becomes one comma-separated
chunk whose fields follow the ``Format`` the header declares -- which is exactly what
``output/parser/vcf.py`` reads back.

A snapshot that did not come from a VCF line carries no columns, and VCF cannot invent them: an
indel's spelling needs reference bases, which a detached snapshot deliberately does not hold. So
whoever builds such a snapshot supplies the columns, with the padding rule
``corpus/exporter/formats/file/vcf.py`` uses; writing one without them is refused rather than
guessed at.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from datetime import datetime

from attrs import define, field

from ensembl.vep import __version__
from ensembl.vep.common.constants import EMPTY_ALLELE
from ensembl.vep.common.methods import reverse_complement
from ensembl.vep.output.factory.base import AnnotatedVariant, AnnotatedVariantFactory

__all__ = ["DEFAULT_INFO_FIELD", "VCF_COLS", "VCFOutputWriter"]

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
class VCFOutputWriter:
    """VCF serializer for annotated-variant snapshots.

    Attributes:
        records (AnnotatedVariantFactory): The factory whose flags the declared ``Format`` must
            describe, so a reader of the header finds the fields the records actually carry.
        info_field (str): The INFO key consequences are written under.
        keep_csq (bool): Keep a ``CSQ`` already on the input line, and its header, instead of
            replacing them.
        command_line (str): The command written into the header.
        time (str | None): The timestamp written into the header, or None for now.
    """

    records: AnnotatedVariantFactory = field(factory=AnnotatedVariantFactory)
    info_field: str = field(default=DEFAULT_INFO_FIELD)
    keep_csq: bool = field(default=False)
    command_line: str = field(default="")
    time: str | None = field(default=None)

    @property
    def version(self) -> str:
        """Return the version written into the header."""
        return __version__

    def timestamp(self) -> str:
        """Return the header timestamp, formatted as VEP formats it."""
        return self.time if self.time is not None else datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def fields(self) -> list[str]:
        """Return the fields each chunk carries: VEP's defaults, then what the flags add."""
        fields = list(VCF_COLS)
        for name in self.records.flag_fields():
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

    def line(self, variant: AnnotatedVariant) -> str:
        """Return one annotated VCF record."""
        columns = self._columns(variant)
        info = "" if columns[_INFO_COLUMN] in {"", "."} else columns[_INFO_COLUMN]
        chunks = [self.info_chunk(annotation, variant.strand) for annotation in variant.annotations]
        if chunks:
            info = f"{info};" if info else ""
            info += f"{self.info_field}={','.join(chunks)}"
        columns[_INFO_COLUMN] = info or "."
        return "\t".join(columns)

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

    def _columns(self, variant: AnnotatedVariant) -> list[str]:
        """Return input columns ready for a replacement annotation field."""
        if variant.input_columns:
            return self.prepare_input_columns(variant.input_columns)
        raise ValueError("VCF output requires the original input columns")
