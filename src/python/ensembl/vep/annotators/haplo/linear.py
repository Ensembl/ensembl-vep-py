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

This placeholder keeps the ``haplo.`` registry namespace available until annotation evaluates
phased allele combinations. It currently delegates to the consequence linear annotator.
"""

from __future__ import annotations

from ensembl.vep.annotators.consequence.strategies.linear import (
    LinearFactsAnnotator as ConsequenceLinearFactsAnnotator,
)
from ensembl.vep.annotators.registry import register_annotator

__all__ = ["LinearFactsAnnotator"]


@register_annotator("haplo.linear")
class LinearFactsAnnotator(ConsequenceLinearFactsAnnotator):
    """Provide the consequence linear annotator under the haplotype namespace.

    Attributes:
        name (str): Registry name, "haplo.linear".
    """

    name = "haplo.linear"
