from __future__ import annotations

from attrs import define, field
from typing import TYPE_CHECKING

from ensembl.vep.models import Feature
if TYPE_CHECKING:
    from ensembl.vep.models.genic import Transcript

@define(kw_only=True, auto_attribs=True)
class UTR(Feature):
    """
    Represents a UTR.

    Attributes:
        transcript: The transcript to which this UTR belongs.
    """
    transcript: Transcript | None = field(default=None, init=False)