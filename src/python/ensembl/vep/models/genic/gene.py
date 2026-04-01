from __future__ import annotations

from attrs import define, field
from typing import TYPE_CHECKING

from ensembl.vep.common import ExternalReference
from ensembl.vep.models import Feature
if TYPE_CHECKING:
    from ensembl.vep.models.genic import Transcript

@define(kw_only=True, auto_attribs=True)
class Gene(Feature):
    """
    Represents a gene.

    Attributes:
        stable_id:  Ensembl stable transcript identifier (e.g. "ENST00000366667").
        gene_id:    Ensembl stable gene identifier (e.g. "ENSG00000116016").
        biotype:    Transcript biotype (e.g. "protein_coding").
        slice:      Genomic location of the transcript.
        version:    Version of the stable identifier.
    """
    stable_id: str
    version: int | None = field(default=None, init=False)
    name: str | None = field(default=None, init=False)
    symbol: str | None = field(default=None, init=False)
    alternative_names: list[ExternalReference] | None = field(default=[], init=False)
    transcripts: list[Transcript] | None = field(default=[], init=False)
    biotype: str | None = field(default=None, init=False)