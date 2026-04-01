from attrs import define

from ensembl.vep.models import Feature

@define(kw_only=True, auto_attribs=True)
class Motif(Feature):
    pass