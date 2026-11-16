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
"""Build deterministic synthetic corpora."""

from __future__ import annotations

from ensembl.vep.corpus.combine import CombinedCorpus, combine_corpora
from ensembl.vep.corpus.registry import build_corpus, register_corpus
from ensembl.vep.corpus.spec import CorpusSpec
from ensembl.vep.corpus.synthetic.base import SyntheticCorpus
from ensembl.vep.corpus.synthetic.components.sequence import build_transcript
from ensembl.vep.corpus.synthetic.components.variation import (
    boundary_variant_alleles,
    build_transcript_variation_alleles,
    fuzz_variant_alleles,
)
from ensembl.vep.corpus.synthetic.default import build_synthetic_genome, default_corpus
from ensembl.vep.corpus.synthetic.heavy import (
    build_heavy_locus_genome,
    build_heavy_locus_transcripts,
    heavy_locus_corpus,
    heavy_locus_variant_alleles,
)

__all__ = [
    "SyntheticCorpus",
    "boundary_variant_alleles",
    "build_heavy_locus_genome",
    "build_heavy_locus_transcripts",
    "build_synthetic_genome",
    "build_transcript",
    "build_transcript_variation_alleles",
    "default_corpus",
    "fuzz_variant_alleles",
    "heavy_locus_corpus",
    "heavy_locus_variant_alleles",
]


@register_corpus("synthetic.all", aliases=("both",))
def _all_corpora(specification: CorpusSpec) -> CombinedCorpus:
    """Build the default and heavy synthetic corpora as one corpus."""
    default_spec = CorpusSpec(kind="synthetic.default", parameters=specification.parameters)
    heavy_spec = CorpusSpec(kind="synthetic.heavy", parameters=specification.parameters)
    return combine_corpora(
        specification,
        (build_corpus(default_spec), build_corpus(heavy_spec)),
    )
