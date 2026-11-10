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
"""Text and HTML summaries of explained overlaps.

The summary describes the variant, transcript, affected positions, codons, and full-sequence
changes before presenting how the annotator reached its result.
"""

from __future__ import annotations

from html import escape

from ensembl.vep.common import Strand
from ensembl.vep.common.constants import EMPTY_ALLELE, UNKNOWN_VALUE
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.models.overlap.transcript_variation_allele import Position
from ensembl.vep.common.methods import format_coords

__all__ = ["format_summary", "summary_html", "summary_rows"]

_CONTEXT_LENGTH = 12
"""Bases or residues shown either side of the first difference in a full CDS or peptide."""


def summary_rows(overlap: TranscriptVariationAllele) -> list[tuple[str, str]]:
    """Return applicable overlap details as labelled rows.

    Args:
        overlap (TranscriptVariationAllele): The overlap to describe.

    Returns:
        list[tuple[str, str]]: Label-value pairs in reading order.
    """
    rows = [("variant", _variant_text(overlap)), ("transcript", _transcript_text(overlap))]
    rows.append(("location", _location_text(overlap)))
    rows.append(("coding", _coding_text(overlap)))

    positions = _positions_text(overlap)
    if positions:
        rows.append(("positions", positions))

    if overlap.reference_codon is not None:
        rows.append(("codons", f"{overlap.reference_codon or '-'}/{overlap.alternate_codon or '-'}"))
        rows.append(("amino acids", f"{overlap.reference_peptide or '-'}/{overlap.alternate_peptide or '-'}"))

    if overlap.reference_cds is not None and overlap.alternate_cds is not None:
        rows.append(("full CDS", _sequence_change(overlap.reference_cds, overlap.alternate_cds, "bp")))
    if overlap.reference_cds_peptide is not None and overlap.alternate_cds_peptide is not None:
        peptide_change = _sequence_change(overlap.reference_cds_peptide, overlap.alternate_cds_peptide, "aa")
        rows.append(("full peptide", peptide_change))
    return rows


def format_summary(overlap: TranscriptVariationAllele) -> str:
    """Return the summary as an aligned text block headed ``summary:``."""
    rows = summary_rows(overlap)
    width = max(len(label) for label, _value in rows)
    lines = ["summary:"]
    for label, value in rows:
        lines.append(f"  {label:<{width}}  {value}")
    return "\n".join(lines)


def summary_html(overlap: TranscriptVariationAllele) -> str:
    """Return the summary as an HTML table for an explanation page."""
    cells = []
    for label, value in summary_rows(overlap):
        cells.append(f"<tr><th>{escape(label)}</th><td><code>{escape(value)}</code></td></tr>")
    return f'<table class="vep-table">{"".join(cells)}</table>'


def _variant_text(overlap: TranscriptVariationAllele) -> str:
    """Return the allele's locus and alleles, with the record's id and the allele's number."""
    allele = overlap.overlapping_feature
    region = allele.slice.seq_region_name if allele.slice is not None else UNKNOWN_VALUE
    text = f"{region}:{format_coords(allele.start, allele.end)} "
    text += f"{allele.reference_allele or EMPTY_ALLELE}>{allele.alternate_allele or EMPTY_ALLELE}"

    variation = allele.variation
    if variation is None:
        return text
    details = []
    if variation.stable_id:
        details.append(f"id {variation.stable_id}")
    alleles = list(variation.variation_alleles or ())
    for number, candidate in enumerate(alleles, start=1):
        if candidate is allele:
            details.append(f"alternate allele {number} of {len(alleles)}")
    if variation.uploaded_locus is not None and variation.uploaded_locus.allele_string:
        details.append(f"input {variation.uploaded_locus.start} {variation.uploaded_locus.allele_string}")
    return f"{text} ({', '.join(details)})" if details else text


def _transcript_text(overlap: TranscriptVariationAllele) -> str:
    """Return the transcript's id, span, strand, biotype and gene."""
    transcript = overlap.feature
    strand = "+" if transcript.strand == Strand.FORWARD else "-"
    details = [f"{transcript.start}-{transcript.end}", f"strand {strand}"]
    gene = transcript.gene
    biotype = transcript.biotype or (gene.biotype if gene is not None else None)
    if biotype:
        details.append(biotype)
    if gene is not None:
        gene_text = gene.stable_id
        if gene.symbol:
            gene_text += f" {gene.symbol}"
        details.append(f"gene {gene_text}")
    return f"{transcript.stable_id} ({', '.join(details)})"


def _location_text(overlap: TranscriptVariationAllele) -> str:
    """Return which exons and introns the variant touches, or how far outside it lies."""
    transcript = overlap.feature
    allele = overlap.overlapping_feature
    variant_min, variant_max = min(allele.start, allele.end), max(allele.start, allele.end)

    touched = []
    for label, features in (("exon", transcript.spliced_exons), ("intron", transcript.introns)):
        numbers = []
        for number, feature in enumerate(features, start=1):
            if variant_min <= feature.end and variant_max >= feature.start:
                numbers.append(number)
        if numbers:
            touched.append(f"{label} {_number_range(numbers)}/{len(features)}")
    if touched:
        return ", ".join(touched)

    if variant_max < transcript.start:
        distance = transcript.start - variant_max
    else:
        distance = variant_min - transcript.end
    upstream = (variant_max < transcript.start) == (transcript.strand == Strand.FORWARD)
    return f"{distance} bp {'upstream' if upstream else 'downstream'} of the transcript"


def _number_range(numbers: list[int]) -> str:
    """Return ``[3]`` as ``3`` and ``[3, 4]`` as ``3-4``."""
    return str(numbers[0]) if len(numbers) == 1 else f"{numbers[0]}-{numbers[-1]}"


def _coding_text(overlap: TranscriptVariationAllele) -> str:
    """Return the coding span and length, and the padding of an incomplete 5' end."""
    transcript = overlap.feature
    if transcript.coding_start is None or transcript.coding_end is None:
        return "non-coding"

    low, high = sorted((transcript.coding_start, transcript.coding_end))
    text = f"genomic {low}-{high}"
    if transcript.cdna_coding_start is not None and transcript.cdna_coding_end is not None:
        text += f", cDNA {transcript.cdna_coding_start}-{transcript.cdna_coding_end}"
    if transcript.cds_length is not None:
        text += f", {transcript.cds_length} bp"
    phase = transcript.cds_start_phase
    if phase:
        text += f" including {phase} N (5' incomplete, start phase {phase})"
    return text


def _positions_text(overlap: TranscriptVariationAllele) -> str:
    """Return the cDNA, CDS and protein positions VEP would write, where there are any."""
    positions: list[tuple[str, Position | None]] = [
        ("cDNA", overlap.cdna_position),
        ("CDS", overlap.cds_position),
        ("protein", overlap.protein_position),
    ]
    texts = []
    for label, position in positions:
        if position is not None:
            texts.append(f"{label} {format_coords(*position)}")
    return ", ".join(texts)


def _sequence_change(reference: str, alternate: str, unit: str) -> str:
    """Return sequence lengths and the first difference, or report no change."""
    if reference == alternate:
        return f"unchanged, {len(reference)} {unit}"

    first_difference = 0
    for reference_symbol, alternate_symbol in zip(reference, alternate):
        if reference_symbol != alternate_symbol:
            break
        first_difference += 1

    window_start = max(0, first_difference - _CONTEXT_LENGTH)
    window_end = first_difference + _CONTEXT_LENGTH + 1
    reference_window = _marked_window(reference, window_start, first_difference, window_end)
    alternate_window = _marked_window(alternate, window_start, first_difference, window_end)
    return (
        f"{len(reference)} -> {len(alternate)} {unit}, first difference at {first_difference + 1}: "
        f"{reference_window} -> {alternate_window}"
    )


def _marked_window(sequence: str, window_start: int, marked: int, window_end: int) -> str:
    """Return ``sequence[window_start:window_end]``, the symbol at ``marked`` in brackets."""
    before = sequence[window_start:marked]
    symbol = sequence[marked : marked + 1] or "-"
    after = sequence[marked + 1 : window_end]
    leading = "..." if window_start > 0 else ""
    trailing = "..." if window_end < len(sequence) else ""
    return f"{leading}{before}[{symbol}]{after}{trailing}"
