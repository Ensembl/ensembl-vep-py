from __future__ import annotations

from attrs import define, field
from typing import TYPE_CHECKING
from functools import cached_property
from threading import Lock

from ensembl.vep.models import Feature
if TYPE_CHECKING:
    from ensembl.vep.models.genic import Transcript, Translation

@define(kw_only=True, auto_attribs=True)
class cDNA(Feature):
    """
    Represents a CDS.

    Attributes:
        transcript: The transcript to which this cDNA belongs.
    """
    transcript: Transcript | None = field(default=None, init=False)