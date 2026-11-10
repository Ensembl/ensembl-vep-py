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
"""Format-neutral annotated variants and their factory.

The factory is a port of the format-independent half of Perl VEP's ``OutputFactory.pm``. A variation
first becomes a hash of what is known from its position alone (``variation_hash``). That hash is
copied and extended once per (allele, feature) overlap (``allele_hash``). The copies are then
filtered, picked or squashed as the options ask (``output_hashes``), and ``build_annotated_variant``
hands them to ``AnnotatedVariant.build``, which detaches the variant's own fields. The snapshot is
the only thing that crosses a worker boundary and the only thing an output format sees.

Only fields this model can actually supply are ever populated. A column VEP fills from data not
carried here yet -- HGVS, SIFT, frequencies, co-located variants -- stays empty rather than being
approximated, which is VEP's own convention for "not available". The *column layout* for a given
set of flags is VEP's, so a downstream tool reading the ``Format`` declaration finds what it
expects.

Consequences reach an allele through ``VariationAllele.transcript_variation_alleles``. Assigning
that list invalidates the allele's aggregate, so assign each overlap's consequences first and
attach the overlaps afterwards.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Iterator, Sequence
from functools import cache
from typing import Any

from attrs import define, field, frozen

from ensembl.vep.annotators.consequence.rules.catalogue import load_consequence_definitions
from ensembl.vep.common.constants import EMPTY_ALLELE, UNKNOWN_VALUE
from ensembl.vep.common.methods import format_coords
from ensembl.vep.models.overlap import OverlapConsequence, TranscriptVariationAllele
from ensembl.vep.models.variation import Variation, VariationAllele
from ensembl.vep.stats import RunSummary, RunSummaryCollector

__all__ = [
    "DEFAULT_PICK_ORDER",
    "DISTANCE_TERMS",
    "FLAG_FIELDS",
    "AlleleOverlap",
    "AnnotatedVariant",
    "AnnotatedVariantFactory",
    "AnnotatedVariants",
    "AnnotationRecord",
]

AnnotationRecord = dict[str, object]
"""One VEP output hash, copied from a selected allele-feature annotation."""


DEFAULT_PICK_ORDER = (
    "mane_select",
    "mane_plus_clinical",
    "canonical",
    "appris",
    "tsl",
    "biotype",
    "ccds",
    "rank",
    "length",
    "ensembl",
    "refseq",
)
"""The order ``--pick`` weighs a transcript by, as VEP defaults it. A criterion this model has no
data for ties for every candidate and so defers to the next one -- exactly VEP's behaviour when
two candidates are equal on a criterion."""

DISTANCE_TERMS = frozenset({"upstream_gene_variant", "downstream_gene_variant"})
"""The terms whose output carries a ``DISTANCE`` to the transcript."""

FLAG_FIELDS = (
    ("allele_number", ("ALLELE_NUM",)),
    ("show_ref_allele", ("REF_ALLELE",)),
    ("flag_pick", ("PICK",)),
    ("flag_pick_allele", ("PICK",)),
    ("variant_class", ("VARIANT_CLASS",)),
    ("symbol", ("SYMBOL", "SYMBOL_SOURCE", "HGNC_ID")),
    ("biotype", ("BIOTYPE",)),
    ("numbers", ("EXON", "INTRON")),
)
"""Flag -> the output fields it adds, in VEP's order, for the flags this factory supports.

The whole field list is declared even where only part of it can be filled (``SYMBOL_SOURCE``,
``HGNC_ID``), so the column layout matches VEP's for the same flags."""

_MISSING_RANK = 1000
_UNKNOWN_PICK_VALUE = 1
_TSL_PATTERN = re.compile(r"tsl(\d+)")
_APPRIS_PATTERN = re.compile(r"([A-Za-z]).+(\d+)")
_VERSIONED = re.compile(r"\.\d+$")


def _as_iterator(variants: Iterable[AnnotatedVariant]) -> Iterator[AnnotatedVariant]:
    """Return a single-pass iterator over annotated variants."""
    return iter(variants)


@cache
def _intergenic_consequence() -> OverlapConsequence:
    """Return ``intergenic_variant`` as the catalogue defines it: the term for an allele with no overlap."""
    for definition in load_consequence_definitions():
        if definition.term == "intergenic_variant":
            return OverlapConsequence(
                term=definition.term,
                accession=definition.accession,
                rank=definition.rank,
                tier=definition.tier,
                impact=definition.impact,
            )
    raise LookupError("the catalogue defines no intergenic_variant term")


@frozen(kw_only=True, auto_attribs=True)
class AnnotatedVariant:
    """Format-neutral snapshot of one annotated variant.

    Attributes:
        id (str | None): The uploaded or stable variant identifier.
        seq_region_name (str): The sequence-region name, as the input spelled it.
        start (int): The inclusive genomic start.
        end (int): The inclusive genomic end.
        strand (int): The genomic strand.
        allele_string (str): The reference and alternate alleles separated by ``/``.
        input_columns (tuple[str, ...]): The original input columns, when present.
        annotations (tuple[AnnotationRecord, ...]): Selected VEP output hashes.
        most_severe_consequence (str | None): The most severe emitted consequence term.
        variant_class (str | None): The Sequence Ontology class term.
        assembly_name (str | None): The assembly name, when known.
        skipped (bool): Whether annotation was skipped for this input variant.
    """

    id: str | None
    seq_region_name: str
    start: int
    end: int
    strand: int
    allele_string: str
    input_columns: tuple[str, ...] = ()
    annotations: tuple[AnnotationRecord, ...] = ()
    most_severe_consequence: str | None = None
    variant_class: str | None = None
    assembly_name: str | None = None
    skipped: bool = False

    @classmethod
    def build(
        cls,
        variation: Variation,
        record_columns: Sequence[str] | None = None,
        *,
        annotations: tuple[AnnotationRecord, ...] = (),
        assembly_name: str | None = None,
        skipped: bool = False,
    ) -> AnnotatedVariant:
        """Return the snapshot of one variation, carrying records chosen elsewhere.

        Which records a run writes is an output-format question, answered by the output factory;
        everything here is the variant's own.

        Args:
            variation (Variation): The variation whose annotation is complete.
            record_columns (Sequence[str] | None): Original input columns, when available.
            annotations (tuple[AnnotationRecord, ...]): The selected output hashes.
            assembly_name (str | None): Assembly name to expose in structured output.
            skipped (bool): Whether sources could not annotate this variation.

        Returns:
            AnnotatedVariant: The detached snapshot.
        """
        uploaded_locus = variation.uploaded_locus
        uploaded_name = uploaded_locus.seq_region_name if uploaded_locus is not None else None
        sequence_region_name = uploaded_name or UNKNOWN_VALUE
        class_term = variation.class_so_term
        most_severe = None if skipped else variation.most_severe_consequence
        most_severe_term = most_severe.term if most_severe is not None else None
        if most_severe_term is None and annotations:
            # A term the model reports only inside the records, such as one a flag squashed.
            first_terms = annotations[0].get("Consequence")
            if isinstance(first_terms, (list, tuple)) and first_terms:
                most_severe_term = str(first_terms[0])

        return cls(
            id=variation.stable_id,
            seq_region_name=sequence_region_name,
            start=variation.start,
            end=variation.end,
            strand=int(variation.strand.value),
            allele_string=variation.allele_string or (class_term.value if class_term is not None else ""),
            input_columns=tuple(str(column) for column in record_columns or ()),
            annotations=annotations,
            most_severe_consequence=most_severe_term,
            variant_class=class_term.value if class_term is not None else None,
            assembly_name=assembly_name,
            skipped=skipped,
        )


@define(kw_only=True, auto_attribs=True)
class AnnotatedVariants(Iterator[AnnotatedVariant]):
    """Closable annotation iterator with its completed run summary.

    Attributes:
        variants (Iterable[AnnotatedVariant]): Ordered annotated variants.
        collector (RunSummaryCollector): Summary state updated as variants become available.
    """

    _iterator: Iterator[AnnotatedVariant] = field(alias="variants", converter=_as_iterator)
    _collector: RunSummaryCollector = field(alias="collector", repr=False)
    _completed: bool = field(default=False, init=False)
    _closed: bool = field(default=False, init=False)
    _variant_count: int = field(default=0, init=False)

    def __iter__(self) -> AnnotatedVariants:
        """Return this single-pass annotation stream."""
        return self

    def __next__(self) -> AnnotatedVariant:
        """Return the next annotated variant in input order."""
        if self._closed:
            raise StopIteration
        try:
            variant = next(self._iterator)
        except StopIteration:
            self._completed = True
            self.close()
            raise
        self._variant_count += 1
        return variant

    def __enter__(self) -> AnnotatedVariants:
        """Return this stream for bounded consumption."""
        return self

    def __exit__(self, *_error: Any) -> None:
        """Close worker resources when consumption stops."""
        self.close()

    @property
    def completed(self) -> bool:
        """Return whether every input variant has been consumed."""
        return self._completed

    @property
    def variant_count(self) -> int:
        """Return the number of variants yielded to the consumer."""
        return self._variant_count

    @property
    def collector(self) -> RunSummaryCollector:
        """Return the collector after complete consumption.

        Raises:
            RuntimeError: If the stream has not been exhausted.
        """
        if not self._completed:
            raise RuntimeError("the run summary is available after the annotation stream is exhausted")
        return self._collector

    @property
    def summary(self) -> RunSummary:
        """Return the final summary after complete consumption.

        Raises:
            RuntimeError: If the stream has not been exhausted.
        """
        if not self._completed:
            raise RuntimeError("the run summary is available after the annotation stream is exhausted")
        return self._collector.finish()

    def close(self) -> None:
        """Release resources held by the underlying annotation iterator."""
        if self._closed:
            return
        close = getattr(self._iterator, "close", None)
        if close is not None:
            close()
        self._closed = True


@frozen(kw_only=True, auto_attribs=True)
class AlleleOverlap:
    """One allele of a variation against one feature, with the consequences it was assigned.

    Attributes:
        allele_number (int): The allele's 1-based position among the variation's alternates.
        allele (VariationAllele): The allele.
        overlap (TranscriptVariationAllele | None): The feature overlap, or None for an allele that
            overlaps nothing -- which VEP reports as intergenic.
        consequences (tuple[OverlapConsequence, ...]): The consequences, most severe first.
    """

    allele_number: int
    allele: VariationAllele
    overlap: TranscriptVariationAllele | None
    consequences: tuple[OverlapConsequence, ...]

    @property
    def allele_seq(self) -> str:
        """Return the allele as VEP writes it, using ``-`` when empty."""
        return self.allele.alternate_allele or EMPTY_ALLELE


@define(kw_only=True, auto_attribs=True)
class AnnotatedVariantFactory:  # pylint: disable=too-many-instance-attributes
    """Factory for detached annotated-variant snapshots.

    Attributes:
        most_severe (bool): One record per variation, carrying only its most severe term.
        summary (bool): One record per variation, carrying every distinct term, most severe first.
        pick (bool): Keep only the single best (allele, feature) per variation.
        pick_allele (bool): Keep the best feature per allele.
        flag_pick (bool): Keep everything, and mark the single best with ``PICK``.
        flag_pick_allele (bool): Keep everything, and mark the best per allele with ``PICK``.
        pick_order (tuple[str, ...]): The criteria the picking flags weigh, in order.
        no_intergenic (bool): Write nothing for an allele that overlaps no feature.
        variant_class (bool): Add ``VARIANT_CLASS``.
        symbol (bool): Add the gene symbol.
        biotype (bool): Add the transcript biotype.
        numbers (bool): Add exon and intron numbers.
        allele_number (bool): Add ``ALLELE_NUM``.
        show_ref_allele (bool): Add ``REF_ALLELE``.
        transcript_version (bool): Append the version to transcript identifiers.
        gene_version (bool): Append the version to gene identifiers.
    """

    most_severe: bool = field(default=False)
    summary: bool = field(default=False)
    pick: bool = field(default=False)
    pick_allele: bool = field(default=False)
    flag_pick: bool = field(default=False)
    flag_pick_allele: bool = field(default=False)
    pick_order: tuple[str, ...] = field(default=DEFAULT_PICK_ORDER)
    no_intergenic: bool = field(default=False)
    variant_class: bool = field(default=False)
    symbol: bool = field(default=False)
    biotype: bool = field(default=False)
    numbers: bool = field(default=False)
    allele_number: bool = field(default=False)
    show_ref_allele: bool = field(default=False)
    transcript_version: bool = field(default=False)
    gene_version: bool = field(default=False)

    def flag_fields(self) -> list[str]:
        """Return the output fields the enabled flags add, in order and without repeats."""
        fields: list[str] = []
        for flag, names in FLAG_FIELDS:
            if not getattr(self, flag):
                continue
            for name in names:
                if name not in fields:
                    fields.append(name)
        return fields

    def output_hashes(
        self, variation: Variation, record_columns: Sequence[str] | None = None
    ) -> list[dict[str, object]]:
        """Return the output records of one variation, filtered and picked as the flags ask.

        Args:
            variation (Variation): The variation, its alleles' overlaps carrying consequences.
            record_columns (Sequence[str] | None): The original VCF columns, if it came from one.

        Returns:
            list[dict[str, object]]: One hash per written (allele, feature), or a single squashed
                hash under ``most_severe``/``summary``. Empty under ``no_intergenic`` for a
                variation that overlaps nothing.
        """
        variation_fields = self.variation_hash(variation, record_columns)
        overlaps = self.allele_overlaps(variation)
        if not overlaps:
            return []
        if self.summary or self.most_severe:
            return [self._summary_only(variation_fields, overlaps)]

        selected, picked = self._filter(overlaps)
        hashes = []
        for entry in selected:
            output_hash = dict(variation_fields)
            self.allele_hash(output_hash, entry)
            if any(entry is chosen for chosen in picked):
                output_hash["PICK"] = 1
            hashes.append(output_hash)
        return hashes

    def build_annotated_variant(
        self,
        variation: Variation,
        record_columns: Sequence[str] | None = None,
        *,
        assembly_name: str | None = None,
        skipped: bool = False,
    ) -> AnnotatedVariant:
        """Detach one annotated variation, carrying the records the flags select.

        Args:
            variation (Variation): The variation whose annotations are complete.
            record_columns (Sequence[str] | None): Original input columns, when available.
            assembly_name (str | None): Assembly name to expose in structured output.
            skipped (bool): Whether sources could not annotate this variation.

        Returns:
            AnnotatedVariant: A serialisable snapshot for stream and file factories.
        """
        annotations = () if skipped else tuple(self.output_hashes(variation, record_columns))
        return AnnotatedVariant.build(
            variation,
            record_columns,
            annotations=annotations,
            assembly_name=assembly_name,
            skipped=skipped,
        )

    def variation_hash(
        self, variation: Variation, record_columns: Sequence[str] | None = None
    ) -> dict[str, object]:
        """Return what is known about a variation from its position alone."""
        uploaded_locus = variation.uploaded_locus
        uploaded_name = uploaded_locus.seq_region_name if uploaded_locus is not None else None
        chromosome = uploaded_name or UNKNOWN_VALUE
        class_term = variation.class_so_term
        alleles = variation.allele_string or (class_term.value if class_term is not None else "")
        positional_name = f"{chromosome}_{variation.start}_{alleles}"

        name = variation.stable_id
        if name and name not in {".", positional_name}:
            uploaded = name
        elif record_columns is not None and len(record_columns) >= 5:
            reference, alternates = record_columns[3], record_columns[4].replace(",", "/")
            uploaded = "_".join(
                (
                    record_columns[0],
                    str(int(record_columns[1])),
                    f"{reference}/{alternates}" if alternates else reference,
                )
            )
        else:
            uploaded = positional_name

        fields: dict[str, object] = {
            "Uploaded_variation": uploaded,
            "Location": f"{chromosome}:{format_coords(variation.start, variation.end)}",
        }
        if self.variant_class and class_term is not None:
            fields["VARIANT_CLASS"] = class_term.value
        return fields

    def allele_overlaps(self, variation: Variation) -> list[AlleleOverlap]:
        """Return every (allele, feature) of a variation, grouped by feature.

        Grouped by feature and then in allele order, as VEP groups them. An allele with no overlap
        yields one intergenic entry, unless ``no_intergenic`` is set.
        """
        order: list[str] = []
        by_feature: dict[str, list[AlleleOverlap]] = {}
        intergenic: list[AlleleOverlap] = []
        for index, allele in enumerate(variation.variation_alleles):
            overlaps = allele.transcript_variation_alleles
            if not overlaps:
                if not self.no_intergenic:
                    intergenic.append(
                        AlleleOverlap(
                            allele_number=index + 1,
                            allele=allele,
                            overlap=None,
                            consequences=(_intergenic_consequence(),),
                        )
                    )
                continue
            for overlap in overlaps:
                key = overlap.feature.stable_id
                if key not in by_feature:
                    by_feature[key] = []
                    order.append(key)
                consequences = tuple(sorted(overlap.overlap_consequences, key=lambda item: item.rank))
                by_feature[key].append(
                    AlleleOverlap(
                        allele_number=index + 1, allele=allele, overlap=overlap, consequences=consequences
                    )
                )
        grouped = [entry for key in order for entry in by_feature[key]]
        return grouped + intergenic

    def allele_hash(self, output_hash: dict[str, object], entry: AlleleOverlap) -> None:
        """Extend a copy of the variation's hash with one (allele, feature)'s fields, in place."""
        output_hash["Allele"] = entry.allele_seq
        output_hash["Consequence"] = [consequence.term for consequence in entry.consequences]
        if entry.consequences and entry.consequences[0].impact is not None:
            output_hash["IMPACT"] = entry.consequences[0].impact.value
        if self.allele_number:
            output_hash["ALLELE_NUM"] = entry.allele_number
        if self.show_ref_allele:
            output_hash["REF_ALLELE"] = entry.allele.reference_allele or EMPTY_ALLELE
        if entry.overlap is not None:
            self._transcript_fields(output_hash, entry)

    def _transcript_fields(self, output_hash: dict[str, object], entry: AlleleOverlap) -> None:
        """Add what a transcript overlap contributes."""
        overlap = entry.overlap
        assert overlap is not None
        transcript = overlap.feature
        gene = transcript.gene

        output_hash["Feature_type"] = "Transcript"
        feature = transcript.stable_id
        if self.transcript_version and transcript.version and not _VERSIONED.search(feature):
            feature = f"{feature}.{transcript.version}"
        output_hash["Feature"] = feature
        if gene is not None:
            gene_id = gene.stable_id
            if self.gene_version and gene.version and not _VERSIONED.search(gene_id):
                gene_id = f"{gene_id}.{gene.version}"
            output_hash["Gene"] = gene_id
        output_hash["STRAND"] = int(transcript.strand.value)

        if self.numbers:
            exon_number = getattr(overlap, "exon_number", None)
            if exon_number:
                output_hash["EXON"] = f"{exon_number}/{len(transcript.spliced_exons)}"
            intron_number = getattr(overlap, "intron_number", None)
            if intron_number:
                output_hash["INTRON"] = f"{intron_number}/{len(transcript.introns)}"

        if DISTANCE_TERMS.intersection(consequence.term for consequence in entry.consequences):
            distance = getattr(overlap, "distance_to_transcript", None)
            if distance is not None:
                output_hash["DISTANCE"] = distance

        if self.symbol and gene is not None and gene.symbol not in (None, EMPTY_ALLELE):
            output_hash["SYMBOL"] = gene.symbol
        if self.biotype:
            # This model may carry the biotype on the gene only; VEP reads the transcript's.
            biotype = transcript.biotype or (gene.biotype if gene is not None else None)
            if biotype:
                output_hash["BIOTYPE"] = biotype

        self._position_fields(output_hash, overlap)

    @staticmethod
    def _position_fields(output_hash: dict[str, object], overlap: TranscriptVariationAllele) -> None:
        """Add cDNA/CDS/protein positions, amino acids and codons, where the overlap has them.

        A position is written only when the variant touches that sequence at all, as VEP writes it,
        and an end off the sequence prints as ``?``. Codons come already cased the way VEP writes
        them, the allele's own bases upper case.
        """
        for column, position in (
            ("cDNA_position", overlap.cdna_position),
            ("CDS_position", overlap.cds_position),
            ("Protein_position", overlap.protein_position),
        ):
            if position is not None:
                output_hash[column] = format_coords(*position)

        if overlap.feature.translation is None:
            return
        reference_peptide, alternate_peptide = overlap.reference_peptide, overlap.alternate_peptide
        if reference_peptide is not None or alternate_peptide is not None:
            if reference_peptide == alternate_peptide:
                output_hash["Amino_acids"] = reference_peptide
            else:
                output_hash["Amino_acids"] = (
                    f"{reference_peptide or EMPTY_ALLELE}/{alternate_peptide or EMPTY_ALLELE}"
                )
        reference_codon, alternate_codon = overlap.reference_codon, overlap.alternate_codon
        if reference_codon is not None or alternate_codon is not None:
            output_hash["Codons"] = f"{reference_codon or EMPTY_ALLELE}/{alternate_codon or EMPTY_ALLELE}"

    def _summary_only(
        self, variation_fields: dict[str, object], overlaps: Sequence[AlleleOverlap]
    ) -> dict[str, object]:
        """Squash every overlap's consequences into one record: all distinct terms, or the worst."""
        output_hash = dict(variation_fields)
        consequences = sorted(
            (consequence for entry in overlaps for consequence in entry.consequences),
            key=lambda item: item.rank,
        )
        if not consequences:
            return output_hash
        if self.summary:
            terms: list[str] = []
            for consequence in consequences:
                if consequence.term not in terms:
                    terms.append(consequence.term)
            output_hash["Consequence"] = terms
        else:
            output_hash["Consequence"] = [consequences[0].term]
        return output_hash

    def _filter(self, overlaps: list[AlleleOverlap]) -> tuple[list[AlleleOverlap], list[AlleleOverlap]]:
        """Apply the picking flags, returning what to write and what to mark ``PICK``."""
        if self.pick:
            return [self.pick_worst(overlaps)], []
        if self.pick_allele:
            return [self.pick_worst(group) for group in _by_allele(overlaps)], []
        if self.flag_pick:
            return overlaps, [self.pick_worst(overlaps)]
        if self.flag_pick_allele:
            return overlaps, [self.pick_worst(group) for group in _by_allele(overlaps)]
        return overlaps, []

    def pick_worst(self, overlaps: Sequence[AlleleOverlap]) -> AlleleOverlap:
        """Return the single overlap the pick order prefers.

        A port of VEP's ``pick_worst_VariationFeatureOverlapAllele``: narrow the candidates one
        criterion at a time, keeping every candidate tied on it, until one remains; the first of
        any tie that survives every criterion wins.

        Raises:
            ValueError: If there is nothing to pick from.
        """
        if not overlaps:
            raise ValueError("nothing to pick from")
        candidates = list(overlaps)
        for criterion in self.pick_order:
            if len(candidates) == 1:
                break
            scores = [_pick_value(entry, criterion) for entry in candidates]
            best = min(scores)
            candidates = [entry for entry, score in zip(candidates, scores, strict=True) if score == best]
        return candidates[0]


def _by_allele(overlaps: Sequence[AlleleOverlap]) -> list[list[AlleleOverlap]]:
    """Group overlaps by allele, in allele order."""
    groups: dict[str, list[AlleleOverlap]] = {}
    for entry in overlaps:
        groups.setdefault(entry.allele_seq, []).append(entry)
    return list(groups.values())


def _attribute_values(entry: AlleleOverlap, name: str) -> list[str]:
    """Return the values of every transcript attribute with one name."""
    if entry.overlap is None:
        return []
    values = []
    for attribute in entry.overlap.feature.attributes or []:
        if attribute.name == name:
            values.append(str(attribute.value))
    return values


# pylint: disable-next=too-many-return-statements
def _pick_value(entry: AlleleOverlap, criterion: str) -> int:
    """Return an overlap's score on one pick criterion; lower is preferred."""
    if criterion == "rank":
        return entry.consequences[0].rank if entry.consequences else _MISSING_RANK
    if entry.overlap is None:
        return _UNKNOWN_PICK_VALUE if criterion not in {"tsl", "appris"} else 100
    transcript = entry.overlap.feature
    if criterion == "mane_select":
        return 0 if _attribute_values(entry, "MANE_Select") else 1
    if criterion == "mane_plus_clinical":
        return 0 if _attribute_values(entry, "MANE_Plus_Clinical") else 1
    if criterion == "biotype":
        biotype = transcript.biotype or (transcript.gene.biotype if transcript.gene is not None else None)
        return 0 if biotype == "protein_coding" else 1
    if criterion == "length":
        length = transcript.cds_length if transcript.translation is not None else transcript.length
        return -(length or 0)
    if criterion == "tsl":
        for value in _attribute_values(entry, "TSL"):
            match = _TSL_PATTERN.search(value)
            if match and int(match.group(1)):
                return int(match.group(1))
        return 100
    if criterion == "appris":
        for value in _attribute_values(entry, "appris"):
            match = _APPRIS_PATTERN.search(value)
            if match:
                grade = int(match.group(2)) + (10 if match.group(1)[0] == "a" else 0)
                if grade:
                    return grade
        return 100
    # canonical, ccds, ensembl, refseq: nothing in this model distinguishes candidates.
    return _UNKNOWN_PICK_VALUE
