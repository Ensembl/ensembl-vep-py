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
"""TranscriptStructuralVariationAllele class."""

from __future__ import annotations

from attr import define, field

from ensembl.vep.models.overlap.transcript_variation_allele import TranscriptVariationAllele


@define(kw_only=True, auto_attribs=True)

__all__ = ["TranscriptStructuralVariationAllele"]


class TranscriptStructuralVariationAllele(TranscriptVariationAllele):
    """
    Represents the overlap between a transcript and a structural variation allele.

    Attributes:
        exon_numbers (list[int]):
            The exon numbers of the transcript that are overlapped by the structural variation allele.
        intron_numbers (list[int]):
            The intron numbers of the transcript that are overlapped by the structural variation allele.
    """

    exon_numbers: list[int] = field(factory=list, init=False)
    intron_numbers: list[int] = field(factory=list, init=False)
    # more ...
