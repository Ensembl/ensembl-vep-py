from attrs import define, field
from typing import TYPE_CHECKING

from ensembl.vep.common import ExternalReference
if TYPE_CHECKING:
    from ensembl.vep.models.genic import Transcript

@define(kw_only=True, auto_attribs=True)
class Gene:
    """
    Represents a gene in the Ensembl database.

    Attributes:
        stable_id:  Ensembl stable transcript identifier (e.g. "ENST00000366667").
        gene_id:    Ensembl stable gene identifier (e.g. "ENSG00000116016").
        biotype:    Transcript biotype (e.g. "protein_coding").
        slice:      Genomic location of the transcript.
        version:    Version of the stable identifier.
    """
    stable_id: str
    version: int | None = None
    name: str | None = None
    symbol: str | None = None
    alternative_names: list[ExternalReference] | None = None
    transcripts: list[Transcript] | None = None
    biotype: str | None = None