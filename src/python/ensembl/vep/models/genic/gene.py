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
"""Gene class."""

from __future__ import annotations

from typing import TYPE_CHECKING
from attrs import define, field

from ensembl.vep.common import ExternalReference
from ensembl.vep.models import Feature

if TYPE_CHECKING:
    from ensembl.vep.models.genic import Transcript


@define(kw_only=True, auto_attribs=True)

__all__ = ["Gene"]


class Gene(Feature):
    """
    Represents a gene.

    Attributes:
        stable_id (str): Ensembl stable gene identifier (e.g. "ENSG00000116016").
        version (int | None): Version of the stable identifier.
        name (str | None): Gene name from known nomenclature such as HGNC(e.g. "BRCA2").
        symbol (str | None): Gene symbol or short name (e.g. "BRCA2").
        alternative_names (list[ExternalReference] | None): Alternative names for the gene.
        transcripts (list[Transcript] | None): List of transcripts associated with the gene.
        biotype (str | None): Biotype of the gene (e.g. "protein_coding").
    """

    stable_id: str
    version: int | None = field(default=None, init=False)
    name: str | None = field(default=None, init=False)
    symbol: str | None = field(default=None, init=False)
    alternative_names: list[ExternalReference] | None = field(factory=list, init=False)
    transcripts: list[Transcript] | None = field(factory=list, init=False)
    biotype: str | None = field(default=None, init=False)
