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
"""Source construction from run configuration.

Two kinds of source, named apart on purpose: the **input source** holds the variants being
annotated, and the **annotation sources** hold what they are annotated against. Perl VEP draws the
same line -- its input is a Parser and its annotation sources come from
``AnnotationSourceAdaptor::get_all``, which reads the configuration and returns whichever set of
sources it names (a cache, a database, a GFF3/GTF file, custom files).

These are functions rather than a class holding the set: a run's sources are built once, wired to
each other, and then used directly by whoever asked for them. The one thing they must share is the
reference sequence and the synonym map, which is why those are built first and passed in.

There is one builder per kind of annotation source rather than one call returning the whole set.
Each annotation track builds the source it annotates against, so nothing has to pick its own source
out of a list of everything configured.

Which input format to read is decided here, by name, against the closed set in
``common/constants.py``. One conditional in one place, rather than a registry each source
registers with: the configuration package offers the same names as ``--format`` choices and cannot
import the sources, which import it.
"""

from __future__ import annotations

import logging
from pathlib import Path

from ensembl.vep.common.constants import COMPRESSION_SUFFIXES, INPUT_FORMAT_EXTENSIONS
from ensembl.vep.config.sections.common.source import SourceConfig
from ensembl.vep.sources.file.fasta.sequence import FastaSequenceSource
from ensembl.vep.sources.file.gff3.transcript import GFF3TranscriptSource
from ensembl.vep.sources.file.gtf.transcript import GTFTranscriptSource
from ensembl.vep.sources.file.vcf.variation import PicklableVCFVariationSource
from ensembl.vep.sources.sequence_source import SequenceSource
from ensembl.vep.sources.synonym import ChromosomeSynonyms, SequenceRegionResolver

__all__ = [
    "DEFAULT_INPUT_FORMAT",
    "InputSource",
    "TranscriptSource",
    "build_chromosome_synonyms",
    "build_input_source",
    "build_sequence_source",
    "build_transcript_source",
    "detect_input_format",
]

logger = logging.getLogger(__name__)

DEFAULT_INPUT_FORMAT = "vcf"
"""What an input file whose name names no format is read as."""

TranscriptSource = GFF3TranscriptSource | GTFTranscriptSource
"""The file-backed transcript annotation sources."""

InputSource = PicklableVCFVariationSource
"""The sources a run can read its variants from.

What a run asks of one: ``iter_all_records`` and ``records_to_features`` either side of a process
boundary, ``region_of`` and ``record_columns`` to place and echo a record, and
``header_lines`` and ``record_count`` for what output and progress need before reading."""


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


def _reference_names(sequence_source: SequenceSource | None) -> frozenset[str]:
    """Return the names the reference sequence declares, or none without one.

    The reference and not the input, because this is the spelling every model slice is built in,
    and a slice is what sequence is fetched with: the sequence source is asked for a region by the
    name it declares, with nothing mapping the name on its way in. Every other side maps into or
    out of it -- a GXF record's own name is resolved in, a region query back out -- and the input's
    own spelling is kept on the variation's ``uploaded_locus``, which is what output reports.

    Asked of the source rather than of a FASTA: which kind of reference a run has is
    ``build_sequence_source``'s business, and a cache or a database answers this as well as a
    file does. A source that lists no regions, and one whose index cannot be read, both leave
    every uploaded name resolving to itself.
    """

    if sequence_source is None:
        return frozenset()
    try:
        return sequence_source.sequence_region_names
    except OSError:
        return frozenset()


def _region_resolver(
    sequence_source: SequenceSource | None,
    synonyms: ChromosomeSynonyms,
    source_names: frozenset[str] | None = None,
) -> SequenceRegionResolver:
    """Return a resolver naming regions the way one side of a run spells them.

    Args:
        sequence_source (SequenceSource | None): The shared reference source, or None for a
            caller that can work without reference sequence -- ``vep explain`` allows it, a run
            does not -- which leaves the resolver with no names.
        synonyms (ChromosomeSynonyms): The run's sequence-region synonyms.
        source_names (frozenset[str] | None): The names to resolve into. None asks for the
            reference's own names, which is what a record read from a file is given.

    Returns:
        SequenceRegionResolver: The resolver for those names.
    """

    names = _reference_names(sequence_source) if source_names is None else source_names
    return SequenceRegionResolver(source_names=names, synonyms=synonyms)


def _format_for_path(input_file: Path | str) -> str | None:
    """Return the format owning an input path's extension, or None when no format claims it."""

    name = Path(input_file).name.lower()
    for suffix in COMPRESSION_SUFFIXES:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    for format_name, extensions in sorted(INPUT_FORMAT_EXTENSIONS.items()):
        if any(name.endswith(extension) for extension in extensions):
            return format_name
    return None


def detect_input_format(input_file: Path | str, configured_format: str | None = None) -> str:
    """Return the format an input file is read as.

    A configured format wins, since only the caller knows what an unqualified file holds. An
    unrecognised name is refused rather than guessed at, but an unrecognised *extension* falls
    back to ``DEFAULT_INPUT_FORMAT``: VEP has always read ``--input_file`` of any name as VCF, and
    ``--format`` is how another format is named.

    Args:
        input_file (Path | str): The configured input file.
        configured_format (str | None): What ``--format`` named, or None.

    Returns:
        str: A registered input format name.

    Raises:
        ValueError: If a configured format names no registered input source.
    """

    if configured_format is not None:
        name = configured_format.strip().lower()
        if name not in INPUT_FORMAT_EXTENSIONS:
            supported = ", ".join(sorted(INPUT_FORMAT_EXTENSIONS)) or "<none>"
            raise ValueError(f"no reader for input format {configured_format!r} yet; supported: {supported}")
        return name

    detected = _format_for_path(input_file)
    if detected is not None:
        return detected
    logger.debug("%s names no input format; reading it as %s", input_file, DEFAULT_INPUT_FORMAT)
    return DEFAULT_INPUT_FORMAT


def build_input_source(
    config: SourceConfig,
    *,
    sequence_source: SequenceSource | None,
    synonyms: ChromosomeSynonyms,
) -> InputSource:
    """Return the source holding the variants being annotated.

    Args:
        config (SourceConfig): The configured sources.
        sequence_source (SequenceSource | None): The shared reference source, or None.
        synonyms (ChromosomeSynonyms): The run's sequence-region synonyms.

    Returns:
        InputSource: The input source the configured format reads with.

    Raises:
        ValueError: If no input file was configured, or its format has no reader.
    """

    if config.input_file is None:
        raise ValueError("annotating needs variants: give --input_file")

    source_format = detect_input_format(config.input_file, config.format)
    # source_format can't be anything other than vcf at the moment
    # keep as a guard
    if source_format != "vcf":
        raise ValueError(f"no reader for input format {source_format!r} yet")
    return PicklableVCFVariationSource(
        file_path=str(config.input_file),
        sequence_source=sequence_source,
        sequence_region_resolver=_region_resolver(sequence_source, synonyms),
    )


def build_transcript_source(
    config: SourceConfig,
    *,
    sequence_source: SequenceSource | None,
    synonyms: ChromosomeSynonyms,
) -> TranscriptSource:
    """Return the transcript source a run annotates against, wired to resolve names both ways.

    Args:
        config (SourceConfig): The configured sources.
        sequence_source (SequenceSource | None): The shared reference source. None is for a
            caller that can work without reference sequence -- ``vep explain`` allows it, a run
            does not -- and leaves the source with no sequence and no reference names.
        synonyms (ChromosomeSynonyms): The run's sequence-region synonyms.

    Returns:
        TranscriptSource: The transcript source.

    Raises:
        ValueError: If the transcript files are ambiguous or absent.
        NotImplementedError: If the configured transcript provider has no reader yet.
    """

    if config.cache:
        raise NotImplementedError("annotating from a VEP cache (--cache, --dir_cache) is not implemented yet")
    if config.database:
        raise NotImplementedError("annotating from an Ensembl database is not implemented yet")
    if config.gff is not None and config.gtf is not None:
        raise ValueError("give --gff or --gtf, not both")
    if config.gff is None and config.gtf is None:
        raise ValueError("annotating a variant file needs transcripts: give --gff or --gtf")

    # Two directions, one map: records read from the file get reference names, and a region query
    # built from a reference name goes out in the file's own spelling.
    reference_region_resolver = _region_resolver(sequence_source, synonyms)
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
    source.source_region_resolver = _region_resolver(sequence_source, synonyms, source.sequence_region_names)
    return source
