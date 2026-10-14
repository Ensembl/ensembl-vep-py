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
"""VCF rendering for corpora."""

from __future__ import annotations

from collections.abc import Iterator, Sequence

from ensembl.vep.common.vcf import vcf_record_alleles
from ensembl.vep.corpus.base import Corpus
from ensembl.vep.models.variation import VariationAllele

__all__ = ["VCF_NAME", "vcf_alleles", "vcf_lines"]

VCF_NAME = "corpus.vcf"

_SOURCE = "ensembl_vep_synthetic"


def distinct_variants(corpus: Corpus) -> list[tuple[str, VariationAllele]]:
    """Return each variant once, under its name, in generated order."""
    return list(corpus.named_variants)


def vcf_alleles(allele: VariationAllele, sequence: str) -> tuple[int, str, str]:
    """Return the ``(POS, REF, ALT)`` a VCF record spells an allele with."""
    position, (reference_allele, alternate_allele) = vcf_record_alleles(
        allele.start,
        (allele.reference_allele or "", allele.alternate_allele or ""),
        lambda base_position: sequence[base_position - 1],
    )
    return position, reference_allele, alternate_allele


def vcf_lines(corpora: Sequence[Corpus]) -> Iterator[str]:
    """Yield the VCF records for corpora."""
    yield "##fileformat=VCFv4.2"
    yield f"##source={_SOURCE}"
    for corpus in corpora:
        for sequence_region in corpus.sequence_regions:
            yield f"##contig=<ID={sequence_region.seq_region_name},length={sequence_region.end}>"
    yield "\t".join(("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))
    for corpus in corpora:
        for sequence_region in corpus.sequence_regions:
            sequence = sequence_region.seq or ""
            records = []
            for name, allele in distinct_variants(corpus):
                if allele.slice is not sequence_region:
                    continue
                position, reference_allele, alternate_allele = vcf_alleles(allele, sequence)
                records.append((position, name, reference_allele, alternate_allele))
            records.sort(key=lambda record: (record[0], record[1]))
            for position, name, reference_allele, alternate_allele in records:
                yield "\t".join(
                    (
                        sequence_region.seq_region_name,
                        str(position),
                        name,
                        reference_allele,
                        alternate_allele,
                        ".",
                        ".",
                        ".",
                    )
                )
