from attrs import frozen

@frozen
class Locus:
    """Container for a genomic position and (optionally) alleles.
    
    Args:
        start (int): Start position of the locus.
        end (int): End position of the locus.
        allele_string (str | None): String representation of the alleles. 
            It could be in the format "{reference_allele}/{alternate_allele}", or just "allele" 
            where allele could be either reference or alternate.
    """
    start: int
    end: int
    allele_string: str | None = None


@frozen
class Attribute:
    """Container for a attribute of a feature.

    Args:
        name (str): Name of the attribute.
        value (str): Value of the attribute.
        code (str | None): Code of the attribute.
        description (str | None): Description of the attribute.
    """
    name: str
    value: str
    code: str | None = None
    description: str | None = None


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


@frozen
class PredictionResult:
    """Container for a prediction result of a program.
    
    Args:
        plugin_name (str): Name of the plugin that generated this prediction.
        consequence (str): Consequence term predicted by the plugin.
        score (float | None): Score associated with the prediction, if applicable.
        description (str | None): Description of the prediction.
    """
    tool: str
    # need something more generic and can be used for both plugins and custom annotation
    #     plugin gets attached to FeatureOverlap and works on FeatureOverlap after overlap computation is done
    #     custom annotation after doing overlap computation gets attached to VFA
    # score: int
    # result: str
    # classification: float | None = None
    version: str | None = None
    modes: list[str] = []
    reference_data: list[ExternalDB] = []