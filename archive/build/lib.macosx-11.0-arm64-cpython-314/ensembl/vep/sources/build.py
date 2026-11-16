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
"""Building a run's sources from its configuration.

Two kinds of source, named apart on purpose: the **input source** holds the variants being
annotated, and the **annotation sources** hold what they are annotated against. Perl VEP draws the
same line -- its input is a Parser and its annotation sources come from
``AnnotationSourceAdaptor::get_all``, which reads the configuration and returns whichever set of
sources it names (a cache, a database, a GFF3/GTF file, custom files).

These are functions rather than a class holding the set: a run's sources are built once, wired to
each other, and then used directly by whoever asked for them. The one thing they must share is the
reference sequence and the synonym map, which is why those are built first and passed in.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence

from ensembl.vep.config.sections.common.source import SourceConfig
from ensembl.vep.sources.feature_source import FeatureSource
from ensembl.vep.sources.file.fasta.sequence import FastaSequenceSource
from ensembl.vep.sources.file.gff3.transcript import GFF3TranscriptSource
from ensembl.vep.sources.file.gtf.transcript import GTFTranscriptSource
from ensembl.vep.sources.file.vcf.variation import VCFVariationSource
from ensembl.vep.sources.synonym import ChromosomeSynonyms, SequenceRegionResolver

__all__ = [
    "TranscriptSource",
    "build_annotation_sources",
    "build_chromosome_synonyms",
    "build_input_source",
    "build_sequence_source",
    "transcript_source_of",
]

logger = logging.getLogger(__name__)

TranscriptSource = GFF3TranscriptSource | GTFTranscriptSource
"""The file-backed transcript annotation sources."""


def build_chromosome_synonyms(config: SourceConfig) -> ChromosomeSynonyms:
    """Return the run's sequence-region synonyms.

    Args:
        config (SourceConfig): The configured sources.

    Returns:
        ChromosomeSynonyms: The map ``--synonyms`` names, or an empty one when it was not given.
    """

    if config.synonyms is None:
        return ChromosomeSynonyms()
    return ChromosomeSynonyms.from_file(config.synonyms)


def build_sequence_source(config: SourceConfig) -> FastaSequenceSource:
    """Return the reference sequence source every other source shares.

    Args:
        config (SourceConfig): The configured sources.

    Returns:
        FastaSequenceSource: The reference source.

    Raises:
        ValueError: If no FASTA was configured.
    """

    if config.fasta is None:
        raise ValueError("annotating a variant file needs a reference: give --fasta")
    return FastaSequenceSource(file_path=config.fasta)


def build_annotation_sources(
    config: SourceConfig,
    *,
    sequence_source: FastaSequenceSource | None,
    synonyms: ChromosomeSynonyms,
) -> list[FeatureSource]:
    """Return the sources a run annotates against.

    Args:
        config (SourceConfig): The configured sources.
        sequence_source (FastaSequenceSource | None): The shared reference source. None is for a
            caller that can work without reference sequence -- ``vep explain`` allows it, a run
            does not -- and leaves every source with no sequence and no reference names.
        synonyms (ChromosomeSynonyms): The run's sequence-region synonyms.

    Returns:
        list[FeatureSource]: One source per configured annotation source.

    Raises:
        ValueError: If the transcript files are ambiguous or absent.
        NotImplementedError: If a configured source has no reader yet.
    """

    if config.cache:
        raise NotImplementedError("annotating from a VEP cache (--cache, --dir_cache) is not implemented yet")
    if config.database:
        raise NotImplementedError("annotating from an Ensembl database is not implemented yet")
    if config.custom:
        raise NotImplementedError("custom annotation sources are not implemented yet")

    return [
        _build_transcript_source(config, sequence_source=sequence_source, synonyms=synonyms),
    ]


def build_input_source(
    config: SourceConfig,
    *,
    sequence_source: FastaSequenceSource | None,
    synonyms: ChromosomeSynonyms,
) -> VCFVariationSource:
    """Return the source holding the variants being annotated.

    Args:
        config (SourceConfig): The configured sources.
        sequence_source (FastaSequenceSource | None): The shared reference source, or None.
        synonyms (ChromosomeSynonyms): The run's sequence-region synonyms.

    Returns:
        VCFVariationSource: The input source.

    Raises:
        ValueError: If no input file was configured.
    """

    if config.input_file is None:
        raise ValueError("annotating needs variants: give --input_file")
    return VCFVariationSource(
        file_path=str(config.input_file),
        sequence_source=sequence_source,
        sequence_region_resolver=SequenceRegionResolver(
            source_names=_reference_names(sequence_source),
            synonyms=synonyms,
        ),
    )


def transcript_source_of(sources: Sequence[FeatureSource]) -> TranscriptSource:
    """Return the transcript source among a run's annotation sources.

    Consequence calling is defined against transcripts, so the runner asks for that one source by
    what it is rather than by where it came from. Once another kind of annotation source exists,
    this is the one place that has to learn how to choose.

    Args:
        sources (Sequence[FeatureSource]): The run's annotation sources.

    Returns:
        TranscriptSource: The transcript source.

    Raises:
        ValueError: If the configured sources hold no transcripts.
    """

    for source in sources:
        if isinstance(source, (GFF3TranscriptSource, GTFTranscriptSource)):
            return source
    raise ValueError("no transcript source was configured")


def _build_transcript_source(
    config: SourceConfig,
    *,
    sequence_source: FastaSequenceSource | None,
    synonyms: ChromosomeSynonyms,
) -> TranscriptSource:
    """Build the configured transcript source, wired to resolve names in both directions.

    Raises:
        ValueError: If the transcript files are ambiguous or absent.
    """

    if config.gff is not None and config.gtf is not None:
        raise ValueError("give --gff or --gtf, not both")
    if config.gff is None and config.gtf is None:
        raise ValueError("annotating a variant file needs transcripts: give --gff or --gtf")

    # Two directions, one map: records read from the file get reference names, and a region query
    # built from a reference name goes out in the file's own spelling.
    reference_region_resolver = SequenceRegionResolver(
        source_names=_reference_names(sequence_source),
        synonyms=synonyms,
    )
    source: TranscriptSource
    if config.gff is not None:
        source = GFF3TranscriptSource(
            file_path=str(config.gff),
            sequence_source=sequence_source,
            reference_region_resolver=reference_region_resolver,
        )
    else:
        assert config.gtf is not None
        source = GTFTranscriptSource(
            file_path=str(config.gtf),
            sequence_source=sequence_source,
            reference_region_resolver=reference_region_resolver,
        )
    source.source_region_resolver = SequenceRegionResolver(
        source_names=source.sequence_region_names,
        synonyms=synonyms,
    )
    return source


def _reference_names(sequence_source: FastaSequenceSource | None) -> frozenset[str]:
    """Return the FASTA's declared names, or none when there is no FASTA or no readable index."""

    if sequence_source is None:
        return frozenset()
    try:
        return sequence_source.sequence_region_names
    except OSError:
        return frozenset()
