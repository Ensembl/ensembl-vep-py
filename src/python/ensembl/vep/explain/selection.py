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
"""Overlap selection for consequence explanations.

Explaining one derivation needs one ``(variant, transcript)`` pair out of a real input, which is a
different job from annotating a whole file: there is no buffering, no output, and reading the whole
input to pair every variant would be waste when the caller named one. What it must not be is a
second opinion about which pairs exist -- so the pairing here is built by the transcript track's
own ``build_transcript_overlap_finder`` over the same sources the run path builds, and only the
*search* is explain's own.

The two kinds of source answer the same question differently: a file is searched for the named
variant alone, pairing no more of it than that takes, while a generated corpus is walked, there
being nothing to query. ``OverlapSelection`` holds that choice for one run, reading each source
lazily and at most once, so selecting by name costs one pass and selecting nothing costs none.
"""

from __future__ import annotations

import logging
from collections.abc import Iterator, Sequence
from contextlib import closing
from pathlib import Path
from typing import Generator

from attrs import define, field

from ensembl.vep.buffer import FixedSizeInputBuffer
from ensembl.vep.common.constants import EMPTY_ALLELE
from ensembl.vep.config import Config
from ensembl.vep.corpus.base import Corpus
from ensembl.vep.corpus.registry import build_corpus
from ensembl.vep.corpus.spec import CorpusSpec
from ensembl.vep.models.genic.transcript import Transcript
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.models.variation.variation import Variation
from ensembl.vep.models.variation.variation_allele import VariationAllele
from ensembl.vep.sources.build import (
    InputSource,
    TranscriptSource,
    build_chromosome_synonyms,
    build_input_source,
    build_sequence_source,
    build_transcript_source,
)
from ensembl.vep.workflow.tracks import build_transcript_overlap_finder, transcript_candidate_window

__all__ = ["FeatureSelector", "NamedOverlap", "OverlapSelection", "VariantSelector"]

logger = logging.getLogger(__name__)

NamedOverlap = tuple[str, TranscriptVariationAllele]
"""One overlap and the name of the variant it came from."""


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
    """Variant selector written as ``CHR:POS:REF>ALT`` or ``CHR:POS:ID``.

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
        """Return the 1-based query span covering this variant's input record.

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


@define(frozen=True, kw_only=True, auto_attribs=True)
class FeatureSelector:
    """Feature selector using a stable identifier and optional GFF3 type prefix.

    A GFF3 spells a transcript's id ``transcript:ENST00000775241`` and the reader keeps it, while a
    caller naturally types ``ENST00000775241``; both name the same feature.

    Attributes:
        stable_id (str): The stable id the caller typed, either spelling.
    """

    stable_id: str

    @classmethod
    def parse(cls, selector: str) -> FeatureSelector:
        """Parse a feature selector.

        Args:
            selector (str): A stable id, e.g. ``ENST00000775241`` or ``transcript:ENST00000775241``.

        Returns:
            FeatureSelector: The parsed selector.

        Raises:
            ValueError: If the selector is empty.
        """
        if not selector:
            raise ValueError("expected a feature stable id, got an empty selector")
        return cls(stable_id=selector)

    def matches(self, feature: Transcript) -> bool:
        """Return whether a feature is the one this selector names.

        A prefix match is not a match: the id after the prefix has to be the same one.
        """
        return self.stable_id in (feature.stable_id, feature.stable_id.split(":", 1)[-1])


@define(kw_only=True, auto_attribs=True)
class OverlapSelection:
    """Selectable overlaps and the sources used to find them.

    Attributes:
        config (Config): The run's configuration; its ``source``, ``corpus`` and ``explain``
            sections say where the overlaps come from and which one is wanted.
        sources (tuple[InputSource, TranscriptSource] | None): The two readers, built from
            the configuration on first use. Assignable so a caller can supply readers of its own
            -- a test over in-memory sequence, say -- instead of real files.
    """

    config: Config
    sources: tuple[InputSource, TranscriptSource] | None = field(init=False, default=None)
    _built_corpus: Corpus | None = field(init=False, default=None)

    def close(self) -> None:
        """Release the generated corpus, if one was built."""
        if self._built_corpus is not None:
            self._built_corpus.close()
            self._built_corpus = None

    def select(self) -> NamedOverlap:
        """Return the single overlap the selectors identify.

        A name or a variant identifies the variant, not the overlap, so one that pairs with several
        features needs ``--transcript`` as well. Saying which stable ids matched is the whole value
        of the error, so the message names them.

        Returns:
            NamedOverlap: The selected overlap and the name it came under.

        Raises:
            LookupError: If nothing matches, if several do and no feature was named, or if the
                index is past the end.
            ValueError: If the variant selector is malformed, or an allele number was given
                without a ``CHR:POS:ID`` variant.
        """
        explain = self.config.explain
        if explain.index is not None:
            for position, named_overlap in enumerate(self.all_overlaps()):
                if position == explain.index:
                    return named_overlap
            raise LookupError(f"index {explain.index} is past the end of {self.description()}")

        if explain.allele_num is not None and explain.variant is None:
            raise ValueError("an allele number selects from a CHR:POS:ID variant; give the variant too")
        selector = (
            VariantSelector.parse(explain.variant, allele_num=explain.allele_num)
            if explain.variant is not None
            else None
        )
        selector_text = explain.name if explain.name is not None else explain.variant
        feature = FeatureSelector.parse(explain.transcript) if explain.transcript is not None else None

        matches: list[NamedOverlap] = []
        for candidate_name, overlap in self.candidate_overlaps(explain.name, selector):
            if explain.name is not None and candidate_name != explain.name:
                continue
            if selector is not None and not selector.matches(overlap.overlapping_feature):
                continue
            if feature is not None and not feature.matches(overlap.feature):
                continue
            matches.append((candidate_name, overlap))
            if feature is not None:
                # A named feature makes the pair unique, so nothing later can make this ambiguous.
                break

        if not matches:
            suffix = f" for feature {explain.transcript!r}" if feature is not None else ""
            raise LookupError(f"no overlap matching {selector_text!r} in {self.description()}{suffix}")
        if len(matches) > 1:
            stable_ids = sorted({overlap.feature.stable_id for _name, overlap in matches})
            raise LookupError(
                f"{len(matches)} overlaps match {selector_text!r} (features: {', '.join(stable_ids)}); "
                "disambiguate by naming one"
            )
        return matches[0]

    def candidate_overlaps(
        self, name: str | None = None, selector: VariantSelector | None = None
    ) -> Iterator[NamedOverlap]:
        """Yield the overlaps that could match, pairing no more of the input than it takes.

        A real input is searched for the named variant alone; a generated corpus is walked, there
        being nothing to query, and is returned unfiltered for the caller to sift.

        Args:
            name (str | None): The variant name to look for.
            selector (VariantSelector | None): The variant's locus and alleles, or its identifier.

        Yields:
            NamedOverlap: Every overlap that still has to be filtered by the caller.
        """
        if self.config.source.input_file is None:
            yield from self._corpus().named_overlaps
            return
        for allele_name, allele in self._get_candidate_variation_alleles(name=name, selector=selector):
            for overlap in self._pair(allele, self._get_candidate_features(allele)):
                yield allele_name, overlap

    def all_overlaps(self) -> Iterator[NamedOverlap]:
        """Yield every selectable overlap, in the order ``--index`` counts.

        Pairing every variant means reading both files whole: an un-indexed GFF3 or GTF cannot
        serve a region query at all, and a caller asking for the listing has asked for all of
        them. The transcripts are grouped by sequence region once, rather than queried per
        variant, which is what keeps this one pass over each file.
        """
        if self.config.source.input_file is None:
            yield from self._corpus().named_overlaps
            return

        transcripts_by_region: dict[str, list[Transcript]] = {}
        for transcript in self._get_candidate_features():
            region = transcript.slice.seq_region_name if transcript.slice is not None else None
            if region is not None:
                transcripts_by_region.setdefault(region, []).append(transcript)

        for name, allele in self._get_candidate_variation_alleles():
            assert allele.slice is not None
            nearby = transcripts_by_region.get(allele.slice.seq_region_name, ())
            for overlap in self._pair(allele, nearby):
                yield name, overlap

    def _get_candidate_variation_alleles(
        self, *, name: str | None = None, selector: VariantSelector | None = None
    ) -> Iterator[tuple[str, VariationAllele]]:
        """Yield named alleles from the input: every one, or the ones a name or selector identifies.

        Filtered, this stops at the first variation with any allele matching, which is what makes
        selecting by name cost a partial pass rather than a whole file. The input is queried by
        region where a selector gives a position and the file has an index.

        Args:
            name (str | None): The variant name to look for, as ``display_name`` spells it.
            selector (VariantSelector | None): The variant's locus and alleles, or its identifier.

        Yields:
            tuple[str, VariationAllele]: Each allele and the name it is known by.
        """
        variation_source, _transcript_source = self._file_sources()
        wanted = name is not None or selector is not None

        def variations() -> Generator[Variation, None, None]:
            """Yield the variations that could match, from a region query where there is an index."""
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

        with closing(variations()) as candidates:
            for variation in candidates:
                matched = []
                for allele in variation.variation_alleles or ():
                    if allele.slice is None:
                        logger.warning("skipping an allele with no sequence region: %r", allele)
                        continue
                    allele_name = allele.display_name
                    if name is not None and allele_name != name:
                        continue
                    if selector is not None and not selector.matches(allele):
                        continue
                    matched.append((allele_name, allele))
                yield from matched
                if matched and wanted:
                    return

    def _get_candidate_features(self, allele: VariationAllele | None = None) -> Sequence[Transcript]:
        """Return the transcripts one allele could pair with, or every transcript in the file.

        One allele's window is a region query where the file has an index, and a stream of just
        that window where it has not. The window is the transcript track's own, so a transcript
        fetched here is one a run would have fetched for the same allele.

        Args:
            allele (VariationAllele | None): The allele to fetch a window around, or None for all.

        Returns:
            Sequence[Transcript]: The candidate transcripts.
        """
        _variation_source, transcript_source = self._file_sources()
        if allele is None:
            return list(transcript_source.get_all_features())

        assert allele.slice is not None
        region = allele.slice.seq_region_name
        widest = transcript_candidate_window(self.config)
        start = max(1, min(allele.start, allele.end) - widest)
        end = max(allele.start, allele.end) + widest
        try:
            return transcript_source.get_features_by_region(seq_region_name=region, start=start, end=end)
        except OSError:
            logger.info(
                "the transcript file has no region index; streaming it for %s:%d-%d", region, start, end
            )
            return list(
                transcript_source.iter_features_by_region(seq_region_name=region, start=start, end=end)
            )

    def _pair(
        self, allele: VariationAllele, transcripts: Sequence[Transcript]
    ) -> Iterator[TranscriptVariationAllele]:
        """Pair one allele with the transcripts it meets, by the run path's own rule.

        The candidates are the transcripts already found here, so the track's finder pairs them
        without any of the caching a run wants around it.
        """
        finder = build_transcript_overlap_finder(self.config, lambda _buffer: transcripts)
        yield from finder.find(FixedSizeInputBuffer(objects=(allele,)))

    def _file_sources(self) -> tuple[InputSource, TranscriptSource]:
        """Build the input and transcript sources the options name, once per run.

        The same builders the run path uses, so a variant explained here is read, named and
        resolved exactly as a run would read, name and resolve it.

        Returns:
            tuple[InputSource, TranscriptSource]: The input source and transcript source.

        Raises:
            ValueError: If the transcript files are ambiguous or absent, or no input file was given.
        """
        if self.sources is not None:
            return self.sources

        config = self.config.source
        synonyms = build_chromosome_synonyms(config)
        if config.fasta is None:
            # A run refuses this; explaining is allowed to go on without it, and says what is
            # lost. The sequence facts then have no bases and a slice cannot resolve its region's
            # length.
            logger.warning(
                "no --fasta given: slices cannot resolve a sequence-region length and sequence facts "
                "have no bases to read"
            )
            sequence_source = None
        else:
            sequence_source = build_sequence_source(config)

        self.sources = (
            build_input_source(config, sequence_source=sequence_source, synonyms=synonyms),
            build_transcript_source(config, sequence_source=sequence_source, synonyms=synonyms),
        )
        return self.sources

    def _corpus(self) -> Corpus:
        """Return the generated corpus the options describe."""
        if self._built_corpus is None:
            self._built_corpus = build_corpus(CorpusSpec.from_options(self.config.corpus))
        return self._built_corpus

    def description(self) -> str:
        """Describe where the overlaps come from, for an error or a log line."""
        source = self.config.source
        if source.input_file is None:
            return self._corpus().description
        features = source.gff if source.gff is not None else source.gtf
        feature_name = "no transcripts" if features is None else Path(features).name
        distance = ",".join(str(value) for value in source.distance)
        return f"{Path(source.input_file).name} against {feature_name} (distance {distance})"

    def listing(self) -> str:
        """Return every selectable overlap, one per line, as index, name and feature."""
        lines = []
        for index, (name, overlap) in enumerate(self.all_overlaps()):
            lines.append(f"[{index:>5}] {name:<42} {overlap.feature.stable_id}")
        return "\n".join(lines)
