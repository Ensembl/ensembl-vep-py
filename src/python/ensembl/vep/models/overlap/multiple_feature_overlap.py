@dataclass
class BaseMultiFeatureOverlap:
    """
    Holds ONE overlapped genomic feature alongside MULTIPLE VFA objects.

    This is the generic replacement for the legacy TranscriptHaplotypeContainer
    (THA) and the foundation for haplotype-aware consequence computation and
    complex SV representation.

    Unlike BaseFeatureOverlap (BFO) which captures a 1-to-1 feature overlap,
    BMFO captures a 1-to-many relationship (one feature : many VFAs).

    Primary use-cases:
    1. Haplotype-aware consequence calling – multiple phased VFAs on a single
       transcript.  Use TranscriptMultiVariation (subclass) for this.
    2. Complex structural variation – multiple events (BND, DEL within INV …)
       each represented as a separate VFA; BMFO aggregates them for a unified
       consequence call on the overlapped feature.

    A new set of predicates operating on BMFO is required to compute
    combined / joint consequences from the set of VFAs.

    Attributes:
        overlapped_feature: The genomic feature being overlapped (Transcript,
                            RegulatoryFeature, …).  Typed as ``object`` here;
                            subclasses should narrow the type.
        variation_feature_alleles: The VFAs overlapping the feature.
        consequences:       Combined consequences derived from all VFAs acting
                            together on the overlapped feature.
    """
    overlapped_feature: object  # Transcript | RegulatoryFeature | …
    variation_feature_alleles: list[VariationFeatureAllele] = field(default_factory=list)
    consequences: list[OverlapConsequence] = field(default_factory=list)

    def add_vfa(self, vfa: VariationFeatureAllele) -> None:
        """Add a VFA to this multi-overlap group."""
        self.variation_feature_alleles.append(vfa)

    def get_vfa_count(self) -> int:
        return len(self.variation_feature_alleles)
    
@dataclass
class TranscriptMultiVariation(BaseMultiFeatureOverlap):
    """
    Haplotype-aware consequence container for a single transcript.

    Replaces TranscriptHaplotypeContainer (THA) from the legacy Perl API.

    Holds a Transcript and the set of phased VFAs that co-occur on the same
    haplotype, enabling compound-heterozygous and *cis*-combination effects
    to be modelled correctly.

    Attributes:
        overlapped_feature: The Transcript (narrowed type from BMFO).
        haplotype_id:       Optional identifier for the haplotype (e.g. from a
                            VCF PS tag or derived phasing algorithm).
        is_reference_haplotype: True if all VFAs are ref alleles (i.e. the
                            wild-type haplotype is being modelled for comparison).
    """
    # overlapped_feature is Transcript but kept as `object` in base for genericity
    haplotype_id: Optional[str] = None
    is_reference_haplotype: bool = False

    @property
    def transcript(self) -> Transcript:
        """Typed convenience accessor for the overlapped transcript."""
        assert isinstance(self.overlapped_feature, Transcript), (
            "overlapped_feature must be a Transcript for TranscriptMultiVariation"
        )
        return self.overlapped_feature  # type: ignore[return-value]