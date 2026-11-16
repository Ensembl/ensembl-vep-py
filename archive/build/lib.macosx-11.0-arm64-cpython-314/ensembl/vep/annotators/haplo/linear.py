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
"""Haplotype-aware linear annotator scaffold.

Keeps the ``haplo.`` namespace live and registered so a real haplotype-aware annotator (one that
evaluates consequences over a phased combination of alleles rather than one variant at a time)
can be dropped in later without touching ``consequence/strategies/``. Until then this is a plain
alias for the consequence annotator's LinearFactsAnnotator -- a deliberate cross-kind dependency
on a placeholder, not a shape the two annotator kinds are meant to keep.
"""

from __future__ import annotations

from ensembl.vep.annotators.consequence.strategies.linear import (
    LinearFactsAnnotator as ConsequenceLinearFactsAnnotator,
)
from ensembl.vep.annotators.registry import register_annotator

__all__ = ["LinearFactsAnnotator"]


@register_annotator("haplo.linear")
class LinearFactsAnnotator(ConsequenceLinearFactsAnnotator):
    """Alias for the consequence linear annotator until haplotype-specific rules exist.

    Attributes:
        name (str): Registry name, "haplo.linear".
    """

    name = "haplo.linear"
