@dataclass
class TranscriptVariation:
    """
    All allele-level consequences of a VariationFeature overlapping a Transcript.

    Groups all TranscriptVariationAllele objects that share the same VF and
    Transcript, providing a convenient per-transcript view.

    Attributes:
        variation_feature: The source VF.
        transcript:        The overlapped transcript.
        alleles:           All TVA objects for each alt (and optionally ref) VFA.
    """
    variation_feature: VariationFeature
    transcript: Transcript
    alleles: list[TranscriptVariationAllele] = field(default_factory=list)

    def get_all_consequences(self) -> list[OverlapConsequence]:
        """Collect all unique consequences across all alleles."""
        seen: set[str] = set()
        out: list[OverlapConsequence] = []
        for tva in self.alleles:
            for c in tva.consequences:
                if c.so_term.value not in seen:
                    seen.add(c.so_term.value)
                    out.append(c)
        return out
    
@dataclass
class RegulatoryFeatureVariation:
    """
    All allele-level consequences of a VF overlapping a RegulatoryFeature.

    Attributes:
        variation_feature:  The source VF.
        regulatory_feature: The overlapped regulatory feature.
        alleles:            All RFVA objects.
    """
    variation_feature: VariationFeature
    regulatory_feature: RegulatoryFeature
    alleles: list[RegulatoryFeatureVariationAllele] = field(default_factory=list)