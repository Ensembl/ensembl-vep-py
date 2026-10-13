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
from ensembl.vep.corpus.synthetic import default as _default  # noqa: F401
from ensembl.vep.corpus.synthetic import heavy as _heavy  # noqa: F401

__all__ = ["SyntheticCorpus"]


@register_corpus("synthetic.all", aliases=("both",))
def _all_corpora(specification: CorpusSpec) -> CombinedCorpus:
    """Build the default and heavy synthetic corpora as one corpus."""
    default_spec = CorpusSpec(kind="synthetic.default", parameters=specification.parameters)
    heavy_spec = CorpusSpec(kind="synthetic.heavy", parameters=specification.parameters)
    return combine_corpora(
        specification,
        (build_corpus(default_spec), build_corpus(heavy_spec)),
    )
