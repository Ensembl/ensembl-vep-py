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
"""Load a corpus from FASTA, GFF3 or GTF, and VCF files."""

from __future__ import annotations

from pathlib import Path

from attrs import define

from ensembl.vep.caches import WeakSliceCache
from ensembl.vep.common import Strand
from ensembl.vep.corpus.loaded.base import LoadedCorpus, collect_sequence_regions
from ensembl.vep.corpus.overlap import build_transcript_overlaps
from ensembl.vep.corpus.registry import register_corpus
from ensembl.vep.corpus.spec import CorpusSpec
from ensembl.vep.models.sequence import Slice
from ensembl.vep.models.variation import Variation, VariationAllele
from ensembl.vep.sources.file.fasta.sequence import FastaSequenceSource
from ensembl.vep.sources.file.gff3.transcript import GFF3TranscriptSource
from ensembl.vep.sources.file.gtf.transcript import GTFTranscriptSource
from ensembl.vep.sources.file.vcf.variation import VCFVariationSource
from ensembl.vep.sources.sequence_source import SequenceSource
from ensembl.vep.sources.synonym import ChromosomeSynonyms, SequenceRegionResolver

__all__ = ["load_file_corpus"]


@define(kw_only=True)
class _ResolvingSliceCache:
    """Resolve uploaded names before delegating to the shared slice cache."""

    cache: WeakSliceCache
    resolver: SequenceRegionResolver

    def get_slice(
        self,
        seq_region_name: str,
        sequence_source: SequenceSource | None = None,
        *,
        start: int = 1,
        end: int | None = None,
        strand: Strand = Strand.FORWARD,
    ) -> Slice:
        """Return a slice using the reference source's canonical region name."""
        resolved_name = self.resolver.resolve(seq_region_name)
        return self.cache.get_slice(
            resolved_name,
            sequence_source,
            start=start,
            end=end,
            strand=strand,
        )

    def clear(self) -> None:
        """Clear the shared cache."""
        self.cache.clear()


def _required_path(specification: CorpusSpec, name: str) -> Path:
    """Return a required path parameter."""
    value = specification.parameter(name)
    if not isinstance(value, str) or not value:
        raise ValueError(f"{specification.kind} requires a non-empty {name!r} path")
    path = Path(value)
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def _named_variants(variations: list[Variation]) -> tuple[tuple[str, VariationAllele], ...]:
    """Name every loaded allele uniquely within its variation."""
    named_variants = []
    used_names: set[str] = set()
    for variation_index, variation in enumerate(variations, start=1):
        base_name = variation.stable_id or f"variation-{variation_index}"
        allele_count = len(variation.variation_alleles)
        for allele_index, allele in enumerate(variation.variation_alleles, start=1):
            name = base_name if allele_count == 1 else f"{base_name}:{allele_index}"
            candidate = name
            suffix = 2
            while candidate in used_names:
                candidate = f"{name}:{suffix}"
                suffix += 1
            name = candidate
            used_names.add(name)
            named_variants.append((name, allele))
    return tuple(named_variants)


def _share_sequence_regions(variations: list[Variation], canonical_by_name: dict[str, Slice]) -> None:
    """Attach loaded variations to the transcript source's sequence-region objects."""
    for variation in variations:
        sequence_region = variation.slice
        if sequence_region is None:
            continue
        canonical = canonical_by_name.setdefault(sequence_region.seq_region_name, sequence_region)
        variation.slice = canonical
        for allele in variation.variation_alleles:
            allele.slice = canonical


@register_corpus("loaded.files")
def load_file_corpus(specification: CorpusSpec) -> LoadedCorpus:
    """Load a corpus from FASTA, VCF, and GFF3 or GTF files.

    ``specification`` must provide ``fasta`` and ``vcf`` parameters plus exactly one of ``gff3``
    or ``gtf``. Its optional ``window`` parameter controls overlap construction.

    Args:
        specification (CorpusSpec): File paths and loader parameters.

    Returns:
        LoadedCorpus: Materialised corpus owning its sequence source.

    Raises:
        FileNotFoundError: If a required input path does not exist.
        TypeError: If a path or overlap window has an invalid type.
        ValueError: If the transcript format selection or loaded memberships are invalid.
    """
    fasta_path = _required_path(specification, "fasta")
    vcf_path = _required_path(specification, "vcf")
    gff3_value = specification.parameter("gff3")
    gtf_value = specification.parameter("gtf")
    if bool(gff3_value) == bool(gtf_value):
        raise ValueError("loaded.files requires exactly one of 'gff3' or 'gtf'")
    transcript_path = _required_path(specification, "gff3" if gff3_value else "gtf")

    sequence_source = FastaSequenceSource(file_path=fasta_path)
    try:
        slice_cache = WeakSliceCache()
        reference_region_resolver = SequenceRegionResolver(
            source_names=sequence_source.sequence_region_names,
            synonyms=ChromosomeSynonyms(),
        )
        transcript_source_class = GFF3TranscriptSource if gff3_value else GTFTranscriptSource
        transcript_source = transcript_source_class(
            file_path=str(transcript_path),
            sequence_source=sequence_source,
            slice_cache=slice_cache,
            reference_region_resolver=reference_region_resolver,
        )
        variation_source = VCFVariationSource(
            file_path=str(vcf_path),
            sequence_source=sequence_source,
            slice_cache=_ResolvingSliceCache(
                cache=slice_cache,
                resolver=reference_region_resolver,
            ),
        )
        transcripts = tuple(transcript_source.iter_all_features())
        variations = list(variation_source.iter_all_features())
        canonical_by_name = {
            transcript.slice.seq_region_name: transcript.slice
            for transcript in transcripts
            if transcript.slice is not None
        }
        _share_sequence_regions(variations, canonical_by_name)
        named_variants = _named_variants(variations)
        window = specification.parameter("window", 5000)
        if not isinstance(window, int) or isinstance(window, bool):
            raise TypeError("loaded corpus overlap window must be an integer")
        named_overlaps = build_transcript_overlaps(named_variants, transcripts, window=window)
        sequence_regions = collect_sequence_regions(
            [*transcripts, *(allele for _name, allele in named_variants)]
        )
        return LoadedCorpus(
            spec=specification,
            sequence_regions=sequence_regions,
            named_variants=named_variants,
            transcripts=transcripts,
            named_overlaps=named_overlaps,
            resources=(sequence_source,),
        )
    except Exception:
        sequence_source.close()
        raise
