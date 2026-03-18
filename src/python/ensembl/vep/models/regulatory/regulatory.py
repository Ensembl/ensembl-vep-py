@dataclass
class RegulatoryFeature:
    """
    Lightweight reference to an Ensembl regulatory feature (Regulation API).

    Attributes:
        stable_id:       Ensembl stable regulatory feature identifier.
        feature_type:    Regulatory feature type (e.g. "promoter", "enhancer").
        slice:           Genomic location.
    """
    stable_id: str
    feature_type: str
    slice: Optional[Slice] = None