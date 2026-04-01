from attrs import define

from ensembl.vep.models.overlap import MultiFeatureOverlap

@define(kw_only=True, auto_attribs=True)
class TranscriptHaplotype(MultiFeatureOverlap):
    pass