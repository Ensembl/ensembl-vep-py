from attrs import frozen

@frozen
class VariationLocus:
    """Container for a variation's genomic position and alleles.
    
    Args:
        start (int): Start position of the variation.
        end (int): End position of the variation.
        allele_string (str): String representation of the alleles.
    """
    start: int
    end: int
    allele_string: str


@frozen
class ExternalDB:
    """Container for external database information.
    
    Args:
        id (str): Unique identifier for the external database.
        name (str): Name of the external database.
        description (str | None): Description of the external database.
        url (str | None): URL where this external database can be accessed.
        release (str | None): Release version of the external database.
    """
    id: str
    name: str
    description: str | None = None
    url: str | None = None
    release: str | None = None


@frozen
class ExternalReference:
    """Container for external reference information.
    
    Args:
        accession_id (str): Accession ID of the external reference.
        name (str): Name of the external reference.
        description (str | None): Description of the external reference.
        assignment_method (str | None): Method used to assign this reference, e.g. - DIRECT.
        url (str | None): URL where this external reference can be accessed.
        source (ExternalDB): Source database for this reference.
    """
    accession_id: str
    name: str
    description: str | None = None
    assignment_method: str | None = None
    url: str | None = None
    source: ExternalDB | None = None