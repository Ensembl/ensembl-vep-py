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

from __future__ import annotations

from attrs import define, field

from ensembl.vep.models import Feature

@define(kw_only=True, auto_attribs=True)
class OverlapConsequence:
    """
    A predicted consequence of an overlap between a feature and a base feature. 
    
    This is typically used to represent the consequence of a variant on a transcript, but it can be used 
    for any type of feature overlap.

    Attributes:
        term: The SO term for the consequence, e.g. "missense_variant".
        accession: The SO accession for the consequence, e.g. "SO:0001583".
        is_default: Whether this consequence is the default consequence for the overlapping feature.
            Only one consequence can be default, for transcript variation allele it is generally intergenic_variant.
        rank: The rank of the consequence, with 1 being the most severe.
        tier: The tier of the consequence. If tier 1 consequence are assigned, then no tier 2 consequences will be
            checked or assigned.
        impact: The predicted impact of the consequence, e.g. "HIGH", "MODERATE", "LOW", "MODIFIER".
        feature_type: The type of the overlapped feature that is affected by this consequence type, e.g. Transcript.
        overlapping_feature_type: The type of the overlapping feature that causes this consequence type, e.g. VariationAllele.
    """
    term: str
    accession: str
    is_default: bool = field(default=False, init=False)
    rank: int = field(default=0, init=False)
    tier: int = field(default=0, init=False)
    impact: str | None = field(default=None, init=False)
    feature_type: type[Feature] | None = field(default=None, init=False)
    overlapping_feature_type: type[Feature] | None = field(default=None, init=False)