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
"""Naming and finding the overlaps ``vep explain`` can select from.

Explaining one derivation needs one ``(variant, transcript)`` pair out of a real input, which is a
different job from annotating a whole file: there is no buffering, no output, and reading the whole
input to pair every variant would be waste when the caller named one. What it must not be is a
second opinion about which pairs exist -- so the pairing here is the same
``TranscriptOverlapFinder`` over the same sources the run path builds, and only the *search* is
explain's own.
"""

from __future__ import annotations

import logging
from collections.abc import Iterator, Sequence
from contextlib import closing
from typing import Generator

from attrs import define

from ensembl.vep.buffer import InputBuffer
from ensembl.vep.common import DEFAULT_DISTANCE
from ensembl.vep.config.sections.common.source import SourceConfig
from ensembl.vep.matcher import OverlapMatcher
from ensembl.vep.models.genic.transcript import Transcript
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.models.variation.variation import Variation
from ensembl.vep.models.variation.variation_allele import VariationAllele
from ensembl.vep.overlaps import TranscriptOverlapFinder
from ensembl.vep.sources.build import (
    TranscriptSource,
    build_annotation_sources,
    build_chromosome_synonyms,
    build_input_source,
    build_sequence_source,
    transcript_source_of,
)
from ensembl.vep.sources.file.vcf.variation import VCFVariationSource

__all__ = [
    "EMPTY_ALLELE",
    "NamedOverlap",
    "VariantSelector",
    "build_explain_sources",
    "find_named_overlaps",
    "named_overlaps",
    "names_feature",
    "variant_name",
]

logger = logging.getLogger(__name__)

NamedOverlap = tuple[str, TranscriptVariationAllele]
"""One overlap and the name of the variant it came from."""

EMPTY_ALLELE = "-"
"""How an empty allele is spelled in a ``CHR:POS:REF>ALT`` selector, as VEP spells it."""


def _normalised_allele(value: str | None) -> str:
    """Return an allele in the single spelling a selector comparison can rely on.

    An empty allele is spelled ``""`` by the generated corpus and ``"-"`` by the VCF reader, which
    follows VEP. Selecting an overlap must work the same way whichever source produced it, so both
    collapse to ``""`` here. Without this, a selector matched SNVs from a file but never indels.

    Args:
        value (str | None): The allele as the model stores it.

    Returns:
        str: The allele upper-cased, with either spelling of "empty" as an empty string.
    """
    if value is None or value == EMPTY_ALLELE:
        return ""
    return value.upper()


@define(frozen=True, kw_only=True, auto_attribs=True)
class VariantSelector:
    """A variant named as ``CHR:POS:REF>ALT`` or as ``CHR:POS:ID``.

    The two forms locate a variant differently. ``REF>ALT`` names one allele by its minimised
    coordinates, as the model stores them. ``ID`` names a record by the identifier and position the
    input gives it -- a VCF's ID and POS columns -- and ``allele_num`` picks one of its alternate
    alleles.

    Attributes:
        seq_region_name (str): Sequence region name.
        position (int): The allele's start for ``REF>ALT``, the input's own position for ``ID``.
        reference_allele (str | None): Upper-cased reference allele, ``""`` when empty; None for ``ID``.
        alternate_allele (str | None): Upper-cased alternate allele, ``""`` when empty; None for ``ID``.
        identifier (str | None): The variant's identifier in the input; None for ``REF>ALT``.
        allele_num (int): Which alternate allele an ``ID`` selects, counting from 1.
    """

    seq_region_name: str
    position: int
    reference_allele: str | None = None
    alternate_allele: str | None = None
    identifier: str | None = None
    allele_num: int = 1

    @classmethod
    def parse(cls, selector: str, *, allele_num: int | None = None) -> VariantSelector:
        """Parse a selector, with ``-`` for an empty allele.

        Args:
            selector (str): ``CHR:POS:REF>ALT``, e.g. ``chr1:1052:->ACG`` for an insertion before 1052,
                or ``CHR:POS:ID``, e.g. ``22:15528194:rs123``.
            allele_num (int | None): The alternate allele an ``ID`` selects, counting from 1. None means
                the first.

        Returns:
            VariantSelector: The parsed selector.

        Raises:
            ValueError: If the selector does not have either shape, or ``allele_num`` is given for a
                ``REF>ALT`` selector, which already names its allele, or is below 1.
        """
        parts = selector.split(":", 2)
        if len(parts) != 3 or not parts[2]:
            raise ValueError(f"expected a variant as CHR:POS:REF>ALT or CHR:POS:ID, got {selector!r}")
        seq_region_name, position_text, detail = parts
        try:
            position = int(position_text)
        except ValueError as error:
            raise ValueError(f"expected an integer position in {selector!r}") from error

        if ">" not in detail:
            if allele_num is not None and allele_num < 1:
                raise ValueError(f"allele numbers count alternate alleles from 1, got {allele_num}")
            return cls(
                seq_region_name=seq_region_name,
                position=position,
                identifier=detail,
                allele_num=allele_num or 1,
            )

        if allele_num is not None:
            raise ValueError(
                f"an allele number selects from a CHR:POS:ID variant; {selector!r} names its allele"
            )
        reference_allele, alternate_allele = detail.split(">", 1)
        return cls(
            seq_region_name=seq_region_name,
            position=position,
            reference_allele=_normalised_allele(reference_allele),
            alternate_allele=_normalised_allele(alternate_allele),
        )

    @property
    def query_span(self) -> tuple[int, int]:
        """The 1-based span a region query must cover to find this variant's input record.

        One base before the position as well, because a VCF anchors an indel on the preceding base
        and the minimised start of an insertion lies after it.
        """
        return max(1, self.position - 1), self.position

    def matches(self, allele: VariationAllele) -> bool:
        """Return whether an allele is the one this selector names."""
        seq_region_name = allele.slice.seq_region_name if allele.slice is not None else None
        if seq_region_name != self.seq_region_name:
            return False
        if self.identifier is None:
            return (
                allele.start == self.position
                and _normalised_allele(allele.reference_allele) == self.reference_allele
                and _normalised_allele(allele.alternate_allele) == self.alternate_allele
            )

        variation = allele.variation
        if variation is None or not variation.stable_id:
            return False
        if self.identifier not in variation.stable_id.split(";"):
            return False
        uploaded_start = (
            variation.uploaded_locus.start if variation.uploaded_locus is not None else variation.start
        )
        if uploaded_start != self.position:
            return False
        alleles = variation.variation_alleles or ()
        return len(alleles) >= self.allele_num and alleles[self.allele_num - 1] is allele


def names_feature(stable_id: str | None, feature: str) -> bool:
    """Return whether a feature selector names a stable id, with or without its type prefix.

    A GFF3 spells a transcript's id ``transcript:ENST00000775241`` and the reader keeps it, while a
    caller naturally types ``ENST00000775241``; both name the same feature.

    Args:
        stable_id (str | None): The feature's stable id.
        feature (str): The name the caller typed.

    Returns:
        bool: Whether they name the same feature.
    """
    if stable_id is None:
        return False
    return feature in (stable_id, stable_id.split(":", 1)[-1])


def _allele_locus_name(allele: VariationAllele) -> str:
    """Return an allele's ``CHR:POS:REF>ALT`` name.

    Deliberately the exact shape ``VariantSelector.parse`` reads, so a name printed by one command
    can be passed straight back to another as ``--variant``.

    Args:
        allele (VariationAllele): The allele to name.

    Returns:
        str: The selector-shaped name.
    """
    seq_region_name = allele.slice.seq_region_name if allele.slice is not None else "?"
    reference_allele = allele.reference_allele or EMPTY_ALLELE
    alternate_allele = allele.alternate_allele or EMPTY_ALLELE
    return f"{seq_region_name}:{allele.start}:{reference_allele}>{alternate_allele}"


def variant_name(variation: Variation, allele: VariationAllele, *, several: bool) -> str:
    """Return the name an allele's overlaps are yielded under.

    A VCF ID names the *variant*, so a multi-allelic record would give two alleles the same name and
    the ``(name, feature)`` pair would stop being unique. The alternate allele is appended in that
    case, and only in that case, so the common single-allele record keeps its plain ID.

    Args:
        variation (Variation): The variant the allele belongs to.
        allele (VariationAllele): The allele being named.
        several (bool): Whether the variant carries more than one alternate allele.

    Returns:
        str: The name to yield this allele's overlaps under.
    """
    if not variation.stable_id:
        return _allele_locus_name(allele)
    if not several:
        return variation.stable_id
    return f"{variation.stable_id}:{allele.alternate_allele or EMPTY_ALLELE}"


def build_explain_sources(config: SourceConfig) -> tuple[VCFVariationSource, TranscriptSource]:
    """Return the input source and the transcript source the options name.

    The same builders the run path uses, so a variant explained here is read, named and resolved
    exactly as a run would read, name and resolve it.

    Args:
        config (SourceConfig): The configured sources.

    Returns:
        tuple[VCFVariationSource, TranscriptSource]: The input source and the transcript source.

    Raises:
        ValueError: If the transcript files are ambiguous or absent, or no input file was given.
    """

    synonyms = build_chromosome_synonyms(config)
    if config.fasta is None:
        # A run refuses this; explaining is allowed to go on without it, and says what is lost. The
        # sequence facts then have no bases and a slice cannot resolve its region's length.
        logger.warning(
            "no --fasta given: slices cannot resolve a sequence-region length and sequence facts "
            "have no bases to read"
        )
        sequence_source = None
    else:
        sequence_source = build_sequence_source(config)

    sources = build_annotation_sources(config, sequence_source=sequence_source, synonyms=synonyms)
    return (
        build_input_source(config, sequence_source=sequence_source, synonyms=synonyms),
        transcript_source_of(sources),
    )


def named_overlaps(
    variation_source: VCFVariationSource,
    transcript_source: TranscriptSource,
    *,
    upstream_distance: int = DEFAULT_DISTANCE,
    downstream_distance: int = DEFAULT_DISTANCE,
) -> Iterator[NamedOverlap]:
    """Yield every overlap the two files pair into, in input order.

    Pairing every variant means reading both files whole: an un-indexed GFF3 or GTF cannot serve a
    region query at all, and a caller asking for the listing has asked for all of them. Selecting
    one variant does not go through here.

    Args:
        variation_source (VCFVariationSource): Reader for the variants.
        transcript_source (TranscriptSource): Reader for the transcripts.
        upstream_distance (int): Bases before a transcript's 5' end that still pair with it.
        downstream_distance (int): Bases past a transcript's 3' end that still pair with it.

    Yields:
        NamedOverlap: Each overlap and the name of its variant.
    """

    transcripts_by_region: dict[str, list[Transcript]] = {}
    for transcript in transcript_source.get_all_features():
        region = transcript.slice.seq_region_name if transcript.slice is not None else None
        if region is not None:
            transcripts_by_region.setdefault(region, []).append(transcript)

    for variation in variation_source.get_all_features():
        alleles = list(variation.variation_alleles or ())
        for allele in alleles:
            region = allele.slice.seq_region_name if allele.slice is not None else None
            if region is None:
                logger.warning("skipping an allele with no sequence region: %r", allele)
                continue
            name = variant_name(variation, allele, several=len(alleles) > 1)
            for overlap in _pair(
                allele,
                transcripts_by_region.get(region, ()),
                upstream_distance=upstream_distance,
                downstream_distance=downstream_distance,
            ):
                yield name, overlap


def find_named_overlaps(
    variation_source: VCFVariationSource,
    transcript_source: TranscriptSource,
    *,
    upstream_distance: int = DEFAULT_DISTANCE,
    downstream_distance: int = DEFAULT_DISTANCE,
    name: str | None = None,
    selector: VariantSelector | None = None,
) -> Iterator[NamedOverlap]:
    """Yield the overlaps of the first variant a name or selector identifies, pairing no other.

    Both files are queried by region where they are indexed. Where they are not, the variants are
    read in order only as far as the first match, and just that variant's window of transcripts is
    streamed.

    Args:
        variation_source (VCFVariationSource): Reader for the variants.
        transcript_source (TranscriptSource): Reader for the transcripts.
        upstream_distance (int): Bases before a transcript's 5' end that still pair with it.
        downstream_distance (int): Bases past a transcript's 3' end that still pair with it.
        name (str | None): The variant name, as ``named_overlaps`` would yield it.
        selector (VariantSelector | None): The variant's locus and alleles, or its identifier.

    Yields:
        NamedOverlap: Each overlap of the matching alleles.
    """

    for allele_name, allele in _first_matching_alleles(variation_source, name=name, selector=selector):
        assert allele.slice is not None
        region = allele.slice.seq_region_name
        low, high = min(allele.start, allele.end), max(allele.start, allele.end)
        # The window has to cover whichever direction applies, which depends on each transcript's
        # strand; the pairing below is what decides, so the wider distance is fetched.
        widest = max(upstream_distance, downstream_distance)
        nearby = _transcripts_near(transcript_source, region, low - widest, high + widest)
        for overlap in _pair(
            allele,
            nearby,
            upstream_distance=upstream_distance,
            downstream_distance=downstream_distance,
        ):
            yield allele_name, overlap


def _pair(
    allele: VariationAllele,
    transcripts: Sequence[Transcript],
    *,
    upstream_distance: int,
    downstream_distance: int,
) -> Iterator[TranscriptVariationAllele]:
    """Pair one allele with the transcripts it meets, by the run path's own rule."""

    finder = TranscriptOverlapFinder(
        candidates=lambda _buffer: transcripts,
        matcher=OverlapMatcher(),
        upstream_distance=upstream_distance,
        downstream_distance=downstream_distance,
    )
    yield from finder.find(InputBuffer(objects=(allele,)))


def _first_matching_alleles(
    variation_source: VCFVariationSource,
    *,
    name: str | None,
    selector: VariantSelector | None,
) -> list[tuple[str, VariationAllele]]:
    """Return the named alleles of the first variant with any allele matching, or none."""

    with closing(_candidate_variations(variation_source, selector)) as variations:
        for variation in variations:
            alleles = list(variation.variation_alleles or ())
            matched = []
            for allele in alleles:
                if allele.slice is None:
                    continue
                allele_name = variant_name(variation, allele, several=len(alleles) > 1)
                if name is not None and allele_name != name:
                    continue
                if selector is not None and not selector.matches(allele):
                    continue
                matched.append((allele_name, allele))
            if matched:
                return matched
    return []


def _candidate_variations(
    variation_source: VCFVariationSource,
    selector: VariantSelector | None,
) -> Generator[Variation, None, None]:
    """Yield the variants that could match, from a region query when the input has an index."""

    if selector is None:
        yield from variation_source.iter_all_features()
        return

    start, end = selector.query_span
    try:
        found = variation_source.get_features_by_region(
            seq_region_name=selector.seq_region_name, start=start, end=end
        )
    except OSError:
        logger.info("the variant input has no region index; reading it in order to the first match")
        yield from variation_source.iter_all_features()
        return
    yield from found


def _transcripts_near(
    transcript_source: TranscriptSource,
    region: str,
    start: int,
    end: int,
) -> list[Transcript]:
    """Return the transcripts meeting a window, from a region query when the file has an index."""

    start = max(1, start)
    try:
        return transcript_source.get_features_by_region(seq_region_name=region, start=start, end=end)
    except OSError:
        logger.info("the transcript file has no region index; streaming it for %s:%d-%d", region, start, end)
        return list(transcript_source.iter_features_by_region(seq_region_name=region, start=start, end=end))
