class OverlapConsequence:
    """
    A predicted consequence of a VariationFeatureAllele overlapping a genomic
    feature (transcript, regulatory element, etc.).

    In the redesigned model consequences are computed on VFA objects rather
    than on VariationFeature, removing the need for the legacy OverlapAllele
    intermediary objects.

    Attributes:
        so_term:         Sequence Ontology consequence term.
        so_accession:    SO accession (e.g. "SO:0001583").
        impact:          Impact rating: HIGH / MODERATE / LOW / MODIFIER.
        display_term:    Human-readable label.
        description:     Longer description of the consequence.
        is_default:      True if this is the default/primary consequence.
        rank:            Severity rank (lower = more severe).
    """
    # so_term: ConsequenceTerm
    # so_accession: str
    # impact: str  # "HIGH" | "MODERATE" | "LOW" | "MODIFIER"
    # display_term: str
    # description: Optional[str] = None
    # is_default: bool = False
    # rank: int = 0
    pass