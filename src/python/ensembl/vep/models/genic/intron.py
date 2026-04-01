from __future__ import annotations

from attrs import define, field
from typing import TYPE_CHECKING

from ensembl.vep.models import Feature
if TYPE_CHECKING:
    from ensembl.vep.models.genic import Exon

@define(kw_only=True, auto_attribs=True)
class Intron(Feature):
    """
    Represents an intron.
    """
    next_exon: Exon | None = field(default=None, init=False)
    prev_exon: Exon | None = field(default=None, init=False)