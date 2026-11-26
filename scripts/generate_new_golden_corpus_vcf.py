# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Build golden-corpus VCF input from a tab-separated test-case table.

The table is deliberately transcript-oriented: genomic allele coordinates may be
written relative to a transcript, its CDS, or one of its exons or introns.  This
module resolves those expressions using the same file-backed ``Transcript``
objects that annotation uses, then outputs each allele as a VCF record.
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

from ensembl.vep.annotators.consequence.rules.catalogue import (
    catalogue_digest,
    load_consequence_definitions,
)
from ensembl.vep.common.methods import reverse_complement
from ensembl.vep.corpus.golden.terms import GoldenTerms
from ensembl.vep.corpus.spec import CorpusSpec
from ensembl.vep.models.genic import Transcript
from ensembl.vep.sources.file.fasta.sequence import FastaSequenceSource
from ensembl.vep.sources.file.gff3.transcript import GFF3TranscriptSource

__all__ = ["GoldenVcfInputError", "REQUIRED_COLUMNS", "write_vcf"]

logger = logging.getLogger(__name__)

REQUIRED_COLUMNS = (
    "VariantID",
    "TranscriptID",
    "GeneID",
    "GeneName",
    "Chr",
    "Start",
    "End",
    "Strand",
    "AlleleStart",
    "AlleleEnd",
    "RefAllele",
    "AltAllele",
    "ExpectedConsequences",
    "Comment",
    "TranscriptContext",
    "ScenarioGroup",
    "GRCh37TranscriptID",
)
"""The required TSV columns, in their required order."""

_COORDINATE_EXPRESSION = re.compile(
    r"^\s*(?P<base>[A-Za-z_][A-Za-z0-9_]*|\d+)(?P<offsets>(?:\s*[+-]\s*\d+)*)\s*$"
)
_EXON_OR_INTRON = re.compile(r"^(?P<kind>exon|intron)_(?P<index>[1-9]\d*)_(?P<edge>start|end)$")
_OFFSET = re.compile(r"(?P<operator>[+-])\s*(?P<value>\d+)")

_DEBUG = False


class GoldenVcfInputError(ValueError):
    """A test-case table cannot be represented as a valid golden VCF."""


@dataclass(frozen=True, kw_only=True)
class _Row:
    """One validated TSV row together with its physical line number."""

    line_number: int
    values: dict[str, str]

    def value(self, column: str) -> str:
        """Return one stripped required value, rejecting an empty value."""

        value = self.values[column].strip()
        if not value:
            raise GoldenVcfInputError(f"row {self.line_number}: {column} must not be empty")
        return value


def _read_rows(path: Path) -> list[_Row]:
    """Read and validate the test-case TSV before producing any output."""

    try:
        handle = path.open(newline="", encoding="utf-8")
    except OSError as error:
        raise GoldenVcfInputError(f"cannot read test-case TSV {path}: {error}") from error

    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise GoldenVcfInputError("test-case TSV has no header")
        fieldnames = reader.fieldnames
        if fieldnames != list(REQUIRED_COLUMNS):
            expected = ", ".join(REQUIRED_COLUMNS)
            actual = ", ".join(fieldnames)
            raise GoldenVcfInputError(f"unexpected TSV columns; expected [{expected}], got [{actual}]")

        rows: list[_Row] = []
        for line_number, values in enumerate(reader, start=2):
            if None in values:
                raise GoldenVcfInputError(f"row {line_number}: too many tab-separated fields")
            assert values is not None
            if any(value is None for value in values.values()):
                raise GoldenVcfInputError(f"row {line_number}: too few tab-separated fields")
            rows.append(
                _Row(line_number=line_number, values={key: value or "" for key, value in values.items()})
            )

    if not rows:
        raise GoldenVcfInputError("test-case TSV has no data rows")
    return sorted(rows, key=_input_row_sort_key)


def _input_row_sort_key(row: _Row) -> tuple[str, int, int]:
    """Return the required chromosome and numeric transcript-start sort key."""

    return (
        row.value("Chr"),
        _integer(row.value("Start"), row=row, column="Start"),
        row.line_number,
    )


def _requested_intervals(rows: Iterable[_Row]) -> dict[str, list[tuple[int, int]]]:
    """Return the transcript spans whose GFF3 models the test-case rows need."""

    intervals: dict[str, list[tuple[int, int]]] = {}
    for row in rows:
        chromosome = row.value("Chr")
        start = _integer(row.value("Start"), row=row, column="Start")
        end = _integer(row.value("End"), row=row, column="End")
        if end < start:
            raise GoldenVcfInputError(f"row {row.line_number}: End must not be less than Start")
        intervals.setdefault(chromosome, []).append((start, end))
    return intervals


def _integer(value: str, *, row: _Row, column: str) -> int:
    """Parse one positive genomic coordinate from the table."""

    try:
        coordinate = int(value)
    except ValueError as error:
        raise GoldenVcfInputError(
            f"row {row.line_number}: {column} must be an integer, got {value!r}"
        ) from error
    if coordinate < 1:
        raise GoldenVcfInputError(f"row {row.line_number}: {column} must be at least 1, got {coordinate}")
    return coordinate


def _strand(value: str, *, row: _Row) -> int:
    """Parse a table strand into the corresponding ``Strand.value``."""

    symbols = {"+": 1, "+1": 1, "1": 1, "-": -1, "-1": -1}
    strand = symbols.get(value)
    if strand is None:
        raise GoldenVcfInputError(
            f"row {row.line_number}: Strand must be one of '+', '1', '-', or '-1', got {value!r}"
        )
    return strand


def _coordinate_base(name: str, transcript: Transcript, *, row: _Row, column: str) -> int:
    """Resolve the non-arithmetic portion of a transcript coordinate expression."""

    if name.startswith("t_"):
        if name.endswith("start"):
            return transcript.seq_region_start
        else:
            return transcript.seq_region_end
    elif name.startswith("cds_"):
        if name.endswith("start"):
            return transcript.cds.start
        else:
            return transcript.cds.end

    match = _EXON_OR_INTRON.fullmatch(name)
    if match is None:
        choices = "t_start, t_end, cds_start, cds_end, exon_<n>_<start|end>, or intron_<n>_<start|end>"
        raise GoldenVcfInputError(
            f"row {row.line_number}: {column} has unknown coordinate {name!r}; use {choices}"
        )

    features = transcript.spliced_exons if match["kind"] == "exon" else transcript.introns
    index = int(match["index"]) - 1
    if index >= len(features):
        raise GoldenVcfInputError(
            f"row {row.line_number}: {column} requests {match['kind']} {index + 1}, but transcript "
            f"{transcript.stable_id!r} has {len(features)} {match['kind']}s"
        )
    return getattr(features[index], match["edge"])


def _resolve_coordinate(expression: str, transcript: Transcript, *, row: _Row, column: str) -> int:
    """Resolve an integer or a simple transcript-relative coordinate expression."""

    match = _COORDINATE_EXPRESSION.fullmatch(expression)
    if match is None:
        raise GoldenVcfInputError(
            f"row {row.line_number}: {column} has invalid coordinate expression {expression!r}"
        )
    base = match["base"]
    value = int(base) if base.isdecimal() else _coordinate_base(base, transcript, row=row, column=column)
    for offset in _OFFSET.finditer(match["offsets"]):
        adjustment = int(offset["value"])
        value += adjustment if offset["operator"] == "+" else -adjustment
    if value < 1:
        raise GoldenVcfInputError(
            f"row {row.line_number}: {column} resolves to {value}, but genomic coordinates start at 1"
        )
    return value


def _transcript_for_row(row: _Row, by_id: dict[str, list[Transcript]]) -> Transcript:
    """Return the one GFF3 transcript identified by all coordinate metadata in a row."""

    transcript_id = row.value("TranscriptID")
    start = _integer(row.value("Start"), row=row, column="Start")
    end = _integer(row.value("End"), row=row, column="End")
    if end < start:
        raise GoldenVcfInputError(f"row {row.line_number}: End must not be less than Start")
    chromosome = row.value("Chr")
    strand = _strand(row.value("Strand"), row=row)
    candidates = [
        transcript
        for transcript in by_id.get(transcript_id, [])
        if transcript.start == start
        and transcript.end == end
        and transcript.strand.value == strand
        and transcript.slice is not None
        and transcript.slice.seq_region_name == chromosome
    ]
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise GoldenVcfInputError(
            f"row {row.line_number}: GFF3 has no transcript {transcript_id!r} on {chromosome}:"
            f"{start}-{end} with strand {strand}"
        )
    raise GoldenVcfInputError(
        f"row {row.line_number}: GFF3 has multiple transcripts matching {transcript_id!r} on "
        f"{chromosome}:{start}-{end}"
    )


def _reference_sequence(
    sequence_source: FastaSequenceSource,
    chromosome: str,
    start: int,
    end: int,
    *,
    row: _Row,
) -> str:
    """Fetch one validated, forward-strand reference span."""

    if chromosome not in sequence_source.sequence_region_names:
        raise GoldenVcfInputError(f"row {row.line_number}: chromosome {chromosome!r} is not in the FASTA")
    length = sequence_source.get_seq_region_length(chromosome)
    assert length is not None
    if end > length:
        raise GoldenVcfInputError(
            f"row {row.line_number}: coordinate {end} is beyond FASTA chromosome "
            f"{chromosome!r} length {length}"
        )
    sequence = sequence_source.get_sequence_by_start_end_strand(chromosome, start, end)
    if sequence is None or len(sequence) != end - start + 1:
        raise GoldenVcfInputError(
            f"row {row.line_number}: unable to fetch {chromosome}:{start}-{end} from the FASTA"
        )
    return sequence


def _vcf_alleles(
    row: _Row,
    transcript: Transcript,
    sequence_source: FastaSequenceSource,
) -> tuple[int, str, str]:
    """Resolve one row's genomic allele and write it in VCF's anchored convention."""

    chromosome = row.value("Chr")
    allele_start = _resolve_coordinate(row.value("AlleleStart"), transcript, row=row, column="AlleleStart")
    allele_end = _resolve_coordinate(row.value("AlleleEnd"), transcript, row=row, column="AlleleEnd")
    reference = row.value("RefAllele")
    alternate = row.value("AltAllele")

    insertion = allele_start == allele_end + 1
    if insertion:
        anchor = _reference_sequence(sequence_source, chromosome, allele_end, allele_end, row=row)
        if reference not in {"NA", "-"}:
            raise GoldenVcfInputError(
                f"row {row.line_number}: an insertion (AlleleStart = AlleleEnd + 1) must use "
                "RefAllele 'NA' or '-'"
            )
        if alternate in {"NA", "-"}:
            raise GoldenVcfInputError(
                f"row {row.line_number}: an insertion must provide a non-empty AltAllele"
            )
        return allele_end, anchor, anchor + alternate

    if allele_end < allele_start:
        raise GoldenVcfInputError(
            f"row {row.line_number}: AlleleEnd must be at least AlleleStart, except for an insertion"
        )
    if reference == "-":
        raise GoldenVcfInputError(f"row {row.line_number}: RefAllele '-' is only valid for an insertion")
    fetched_reference = _reference_sequence(sequence_source, chromosome, allele_start, allele_end, row=row)
    if reference == "NA":
        reference = fetched_reference
    if alternate == "NA":
        alternate = reverse_complement(reference)

    if alternate != "-":
        return allele_start, reference, alternate
    if allele_start == 1:
        raise GoldenVcfInputError(
            f"row {row.line_number}: cannot left-anchor a deletion starting at position 1"
        )
    anchor = _reference_sequence(sequence_source, chromosome, allele_start - 1, allele_start - 1, row=row)
    return allele_start - 1, anchor + reference, anchor


def _record_id(row: _Row) -> str:
    """Return the VCF identifier supplied by the test-case schema."""

    return row.value("VariantID")


def _consequence_ranks() -> dict[str, int]:
    """Return the consequence catalogue's rank for each term."""

    return {definition.term: definition.rank for definition in load_consequence_definitions()}


def _expected_terms(row: _Row, consequence_ranks: dict[str, int]) -> tuple[str, ...]:
    """Return expected terms in catalogue-rank order, preserving unknown-label order."""

    terms = tuple(term.strip() for term in row.value("ExpectedConsequences").split(","))
    if any(not term for term in terms):
        raise GoldenVcfInputError(
            f"row {row.line_number}: ExpectedConsequences must be a comma-separated list of terms"
        )
    return tuple(sorted(terms, key=lambda term: consequence_ranks.get(term, float("inf"))))


def _header(contigs: Iterable[tuple[str, int]]) -> list[str]:
    """Return the VCF header for exactly the chromosomes used by the TSV."""

    lines = ["##fileformat=VCFv4.2", "##source=ensembl_vep_golden"]
    lines.extend(f"##contig=<ID={chromosome},length={length}>" for chromosome, length in contigs)
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    return lines


def _has_tabix_index(path: Path) -> bool:
    """Return whether a GFF3 has an index for complete indexed transcript queries."""

    return Path(f"{path}.tbi").is_file() or Path(f"{path}.csi").is_file()


def write_vcf(
    test_cases: str | Path,
    fasta: str | Path,
    gff3: str | Path,
    output_vcf: str | Path,
    output_json: str | Path,
) -> int:
    """Generate VCF and golden-terms files from transcript-relative test cases.

    All rows are resolved before either output is opened, so malformed input
    never leaves a partial VCF or terms file behind.  The return value is the
    number of VCF records written.
    """

    rows = _read_rows(Path(test_cases))
    fasta_path = Path(fasta)
    gff3_path = Path(gff3)
    output_vcf_path = Path(output_vcf)
    output_json_path = Path(output_json)
    if not fasta_path.is_file():
        raise GoldenVcfInputError(f"FASTA file does not exist: {fasta_path}")
    if not gff3_path.is_file():
        raise GoldenVcfInputError(f"GFF3 file does not exist: {gff3_path}")

    sequence_source = FastaSequenceSource(file_path=fasta_path)
    transcript_source = GFF3TranscriptSource(file_path=str(gff3_path), sequence_source=sequence_source)
    try:
        consequence_ranks = _consequence_ranks()
        intervals_by_chromosome = _requested_intervals(rows)
        transcripts_by_id: dict[str, list[Transcript]] = {}
        seen_transcripts: set[tuple[str, str, int, int, int]] = set()

        def add_transcript(transcript: Transcript) -> None:
            """Store one model once even when overlapping queries return it repeatedly."""

            chromosome = transcript.slice.seq_region_name if transcript.slice is not None else ""
            key = (
                transcript.stable_id,
                chromosome,
                transcript.start,
                transcript.end,
                transcript.strand.value,
            )
            if key not in seen_transcripts:
                seen_transcripts.add(key)
                transcripts_by_id.setdefault(transcript.stable_id, []).append(transcript)

        if _has_tabix_index(gff3_path):
            for chromosome, intervals in intervals_by_chromosome.items():
                for start, end in sorted(set(intervals)):
                    for transcript in transcript_source.get_features_by_region(
                        seq_region_name=chromosome, start=start, end=end
                    ):
                        add_transcript(transcript)
        else:
            for transcript in transcript_source.iter_features_by_regions(intervals_by_chromosome):
                add_transcript(transcript)

        contigs: dict[str, int] = {}
        records: list[tuple[str, int, int, str]] = []
        expected_by_overlap: dict[str, tuple[str, ...]] = {}
        used_variant_ids: set[str] = set()
        for row in rows:
            transcript = _transcript_for_row(row, transcripts_by_id)
            variant_id = _record_id(row)
            if variant_id in used_variant_ids:
                raise GoldenVcfInputError(
                    f"row {row.line_number}: VariantID {variant_id!r} is repeated; VCF IDs must be unique"
                )
            used_variant_ids.add(variant_id)
            overlap_key = f"{transcript.stable_id}::{variant_id}"
            expected_by_overlap[overlap_key] = _expected_terms(row, consequence_ranks)

            chromosome = row.value("Chr")
            if chromosome not in contigs:
                if chromosome not in sequence_source.sequence_region_names:
                    raise GoldenVcfInputError(
                        f"row {row.line_number}: chromosome {chromosome!r} is not in the FASTA"
                    )
                length = sequence_source.get_seq_region_length(chromosome)
                assert length is not None
                contigs[chromosome] = length
            position, reference, alternate = _vcf_alleles(row, transcript, sequence_source)
            line = "\t".join((chromosome, str(position), variant_id, reference, alternate, ".", ".", "."))
            records.append((chromosome, position, row.line_number, line))
    finally:
        transcript_source.close()

    try:
        records.sort(key=lambda record: record[:3])
        output_vcf_path.write_text(
            "\n".join([*_header(contigs.items()), *(record[3] for record in records)]) + "\n",
            encoding="utf-8",
        )
        GoldenTerms(
            spec=CorpusSpec(kind="golden.file", parameters={"dataset": "default"}),
            annotator="consequence.linear",
            catalogue=catalogue_digest(),
            terms=expected_by_overlap,
        ).write(output_json_path)
    except OSError as error:
        raise GoldenVcfInputError(f"cannot write output files: {error}") from error
    return len(records)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the golden VCF generator command line."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("test_cases", type=Path, help="tab-separated golden test-case table")
    parser.add_argument("fasta", type=Path, help="reference FASTA file")
    parser.add_argument("gff3", type=Path, help="transcript annotation GFF3 file")
    parser.add_argument("output_vcf", type=Path, help="VCF file to write")
    parser.add_argument("output_json", type=Path, help="golden terms JSON file to write")
    parser.add_argument("--debug", action="store_true", help="Output debugging info")
    arguments = parser.parse_args(argv)
    # global _DEBUG = True if arguments.debug else False

    try:
        count = write_vcf(
            arguments.test_cases,
            arguments.fasta,
            arguments.gff3,
            arguments.output_vcf,
            arguments.output_json,
        )
    except GoldenVcfInputError as error:
        parser.error(str(error))
    print(f"wrote {count} VCF records to {arguments.output_vcf} and terms to {arguments.output_json}")
    return 0


if __name__ == "__main__":  # pragma: no cover - exercised through the thin script wrapper
    raise SystemExit(main())
