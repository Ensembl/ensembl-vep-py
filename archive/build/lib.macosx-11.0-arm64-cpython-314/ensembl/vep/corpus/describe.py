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
"""Describing what a corpus holds, before any annotator has said anything about it.

Shape only: chromosomes, transcripts, variants by class, how the variants pair into overlaps, and
the overlap listing that gives ``vep explain`` its ``--name``/``--index``/``--variant`` selectors.
What an annotator makes of the corpus is a run summary, and lives with the stats.
"""

from __future__ import annotations

from collections.abc import Sequence

from attrs import frozen

from ensembl.vep.common import Strand
from ensembl.vep.corpus.base import Corpus, as_corpus
from ensembl.vep.corpus.store import distinct_variants
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.models.sequence import Slice
from ensembl.vep.models.variation import VariationAllele

__all__ = [
    "CHROMOSOME_HEADERS",
    "CLASS_HEADERS",
    "OVERLAP_HEADERS",
    "TRANSCRIPT_HEADERS",
    "ChromosomeShape",
    "CorpusDescription",
    "OverlapEntry",
    "TranscriptShape",
    "describe_variant",
    "variant_class",
]

CHROMOSOME_HEADERS = ("Chromosome", "Length", "Transcripts", "Variants", "Overlaps")
TRANSCRIPT_HEADERS = (
    "Transcript",
    "Chromosome",
    "Strand",
    "Start",
    "End",
    "Exons",
    "Coding",
    "Biotype",
    "Overlaps",
)
CLASS_HEADERS = ("Class", "Variants")
OVERLAP_HEADERS = ("Index", "Name", "Feature", "Locus", "Alleles")

EMPTY_ALLELE = "-"
_ALLELE_DISPLAY_LENGTH = 12
"""Longest allele shown in full. A whole-transcript deletion is thousands of bases, and nobody
selects one of those by typing it out; ``--name`` is the selector for it."""


def _abbreviate(allele: str) -> str:
    """Return an allele short enough for a table cell, keeping its length visible."""
    if len(allele) <= _ALLELE_DISPLAY_LENGTH:
        return allele
    return f"{allele[:_ALLELE_DISPLAY_LENGTH]}...({len(allele)}bp)"


def variant_class(allele: VariationAllele) -> str:
    """Return a variant's class from its trimmed alleles alone."""
    reference_allele = (allele.reference_allele or "").replace(EMPTY_ALLELE, "")
    alternate_allele = (allele.alternate_allele or "").replace(EMPTY_ALLELE, "")
    if not alternate_allele:
        return "deletion"
    if not reference_allele:
        return "insertion"
    if len(reference_allele) == 1 and len(alternate_allele) == 1:
        return "SNV"
    return "substitution"


def describe_variant(allele: VariationAllele) -> tuple[str, str]:
    """Return a variant's locus and alleles, spelled the way ``--variant`` accepts them.

    A long allele is abbreviated with its length, so the listing stays a table; the full selector
    for such a variant is its name.
    """
    seq_region_name = allele.slice.seq_region_name if allele.slice is not None else "?"
    reference_allele = _abbreviate(allele.reference_allele or EMPTY_ALLELE)
    alternate_allele = _abbreviate(allele.alternate_allele or EMPTY_ALLELE)
    return f"{seq_region_name}:{allele.start}", f"{reference_allele}>{alternate_allele}"


@frozen(kw_only=True, auto_attribs=True)
class ChromosomeShape:
    """One chromosome's totals.

    Attributes:
        name (str): Sequence region name.
        length (int): Length in bases.
        transcripts (int): Transcripts on it.
        variants (int): Distinct variants on it.
        overlaps (int): Overlaps those variants pair into.
    """

    name: str
    length: int
    transcripts: int
    variants: int
    overlaps: int


@frozen(kw_only=True, auto_attribs=True)
class TranscriptShape:
    """One transcript's shape and how many overlaps reach it.

    Attributes:
        stable_id (str): The transcript's stable id.
        chromosome (str): Sequence region name.
        strand (str): ``+`` or ``-``.
        start (int): Genomic start.
        end (int): Genomic end.
        exons (int): Exon count.
        coding (bool): Whether it carries a translation.
        biotype (str): Its own biotype, or its gene's, or empty.
        overlaps (int): Overlaps pairing a variant with it.
    """

    stable_id: str
    chromosome: str
    strand: str
    start: int
    end: int
    exons: int
    coding: bool
    biotype: str
    overlaps: int


@frozen(kw_only=True, auto_attribs=True)
class OverlapEntry:
    """One selectable overlap, as ``vep explain`` addresses it.

    Attributes:
        index (int): Position in the source's order, for ``--index``.
        name (str): The variant's name, for ``--name``.
        feature (str): The transcript's stable id, for ``--transcript``.
        locus (str): ``CHR:POS``, the first part of a ``--variant`` selector.
        alleles (str): ``REF>ALT``, the rest of it, with ``-`` for an empty allele.
    """

    index: int
    name: str
    feature: str
    locus: str
    alleles: str


@frozen(kw_only=True, auto_attribs=True)
class CorpusDescription:
    """A corpus's shape, as tables.

    Attributes:
        kind (str): Which corpus kind this describes.
        chromosomes (tuple[ChromosomeShape, ...]): Per-chromosome totals.
        transcripts (tuple[TranscriptShape, ...]): Every transcript.
        variant_classes (tuple[tuple[str, int], ...]): Class -> distinct variants, largest first.
        batch_size (int): The chunk size ``transcripts_per_chunk`` was measured at.
        transcripts_per_chunk (float): Mean distinct transcripts per chunk of overlaps -- the
            number that says whether a batch shares work between overlaps or not.
        overlaps (tuple[OverlapEntry, ...]): The selectable overlaps, in source order.
    """

    kind: str
    chromosomes: tuple[ChromosomeShape, ...]
    transcripts: tuple[TranscriptShape, ...]
    variant_classes: tuple[tuple[str, int], ...]
    batch_size: int
    transcripts_per_chunk: float
    overlaps: tuple[OverlapEntry, ...]

    @classmethod
    def build(cls, corpus: Corpus, *, select: str | None = None, batch_size: int = 256) -> CorpusDescription:
        """Describe a source's corpora.

        Args:
            corpus (Corpus): The corpus to describe.
            select (str | None): Keep only the overlaps of the variant with this name.
            batch_size (int): Overlaps per chunk, for the sharing measure.

        Returns:
            CorpusDescription: The description.
        """
        corpus = as_corpus(corpus)
        class_totals: dict[str, int] = {}
        for _name, allele in distinct_variants(corpus):
            kind = variant_class(allele)
            class_totals[kind] = class_totals.get(kind, 0) + 1
        classes = sorted(class_totals.items(), key=lambda entry: (-entry[1], entry[0]))

        named_overlaps = list(corpus.named_overlaps)
        entries = []
        for index, (name, overlap) in enumerate(named_overlaps):
            if select is not None and name != select:
                continue
            locus, alleles = describe_variant(overlap.overlapping_feature)
            entries.append(
                OverlapEntry(
                    index=index, name=name, feature=overlap.feature.stable_id, locus=locus, alleles=alleles
                )
            )

        return cls(
            kind=corpus.spec.kind,
            chromosomes=tuple(
                _chromosome_shape(corpus, sequence_region) for sequence_region in corpus.sequence_regions
            ),
            transcripts=tuple(_transcript_shapes(corpus)),
            variant_classes=tuple(classes),
            batch_size=batch_size,
            transcripts_per_chunk=_transcripts_per_chunk(
                [overlap for _name, overlap in named_overlaps], batch_size
            ),
            overlaps=tuple(entries),
        )

    def chromosome_rows(self) -> list[tuple[object, ...]]:
        """Return the chromosome table, in ``CHROMOSOME_HEADERS`` order."""
        rows: list[tuple[object, ...]] = []
        for shape in self.chromosomes:
            rows.append((shape.name, shape.length, shape.transcripts, shape.variants, shape.overlaps))
        return rows

    def transcript_rows(self) -> list[tuple[object, ...]]:
        """Return the transcript table, in ``TRANSCRIPT_HEADERS`` order."""
        rows: list[tuple[object, ...]] = []
        for shape in self.transcripts:
            rows.append(
                (
                    shape.stable_id,
                    shape.chromosome,
                    shape.strand,
                    shape.start,
                    shape.end,
                    shape.exons,
                    "yes" if shape.coding else "no",
                    shape.biotype,
                    shape.overlaps,
                )
            )
        return rows

    def class_rows(self) -> list[tuple[object, ...]]:
        """Return the variant-class table, in ``CLASS_HEADERS`` order."""
        return list(self.variant_classes)

    def overlap_rows(self) -> list[tuple[object, ...]]:
        """Return the overlap listing, in ``OVERLAP_HEADERS`` order."""
        rows: list[tuple[object, ...]] = []
        for entry in self.overlaps:
            rows.append((entry.index, entry.name, entry.feature, entry.locus, entry.alleles))
        return rows


def _chromosome_shape(corpus: Corpus, sequence_region: Slice) -> ChromosomeShape:
    """Return one corpus chromosome's totals."""
    transcripts = [transcript for transcript in corpus.transcripts if transcript.slice is sequence_region]
    variants = [allele for _name, allele in corpus.named_variants if allele.slice is sequence_region]
    overlaps = [
        overlap for _name, overlap in corpus.named_overlaps if overlap.feature.slice is sequence_region
    ]
    return ChromosomeShape(
        name=sequence_region.seq_region_name,
        length=sequence_region.end,
        transcripts=len(transcripts),
        variants=len(variants),
        overlaps=len(overlaps),
    )


def _transcript_shapes(corpus: Corpus) -> list[TranscriptShape]:
    """Return every transcript's shape, with its overlap count."""
    overlaps_by_transcript: dict[str, int] = {}
    for _name, overlap in corpus.named_overlaps:
        stable_id = overlap.feature.stable_id
        overlaps_by_transcript[stable_id] = overlaps_by_transcript.get(stable_id, 0) + 1

    shapes = []
    for transcript in corpus.transcripts:
        if transcript.slice is None:
            continue
        biotype = transcript.biotype
        if biotype is None and transcript.gene is not None:
            biotype = transcript.gene.biotype
        shapes.append(
            TranscriptShape(
                stable_id=transcript.stable_id,
                chromosome=transcript.slice.seq_region_name,
                strand="+" if transcript.strand == Strand.FORWARD else "-",
                start=transcript.start,
                end=transcript.end,
                exons=len(transcript.spliced_exons),
                coding=transcript.translation is not None,
                biotype=biotype or "",
                overlaps=overlaps_by_transcript.get(transcript.stable_id, 0),
            )
        )
    return shapes


def _transcripts_per_chunk(overlaps: Sequence[TranscriptVariationAllele], batch_size: int) -> float:
    """Return the mean number of distinct transcripts per chunk of ``batch_size`` overlaps."""
    if not overlaps or batch_size <= 0:
        return 0.0
    distinct_counts = []
    for chunk_start in range(0, len(overlaps), batch_size):
        chunk = overlaps[chunk_start : chunk_start + batch_size]
        distinct_counts.append(len({overlap.feature.stable_id for overlap in chunk}))
    return sum(distinct_counts) / len(distinct_counts)
