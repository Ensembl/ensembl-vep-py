from attrs import define

from ensembl.vep.models.variation import Variation

@define(kw_only=True, auto_attribs=True)
class StructuralVariation(Variation):
    pass