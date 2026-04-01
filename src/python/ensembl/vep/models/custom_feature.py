from __future__ import annotations

from attrs import define

from ensembl.vep.models import Feature


@define(kw_only=True, auto_attribs=True)
class CustomFeature(Feature):
    pass