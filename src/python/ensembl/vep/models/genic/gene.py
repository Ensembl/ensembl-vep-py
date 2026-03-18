from __future__ import annotations

from dataclasses import dataclass, field 
from enum import Enum, auto
from typing import Optional
from ..sequence import Slice 
from ..genic import Transcript

@dataclass
class Gene:
    """
    Lightweight reference to an Ensembl transcript (Core API object).

    For full transcript data (exons, CDS, UTRs …) the actual Core API
    transcript object should be fetched via the transcript adaptor.

    Attributes:
        stable_id:  Ensembl stable transcript identifier (e.g. "ENST00000366667").
        gene_id:    Ensembl stable gene identifier (e.g. "ENSG00000116016").
        biotype:    Transcript biotype (e.g. "protein_coding").
        slice:      Genomic location of the transcript.
        version:    Version of the stable identifier.
    """
    stable_id: str
    name: str
    symbol: str
    alternative_symbols: list[str]
    alternative_names: list[str]
    version: int
    slice: Slice
    transcripts: list[Transcript]
    biotype: str