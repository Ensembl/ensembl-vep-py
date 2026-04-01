from enum import Enum

class Strand(Enum):
    """Genomic strand."""
    FORWARD = 1
    REVERSE = -1
    UNKNOWN = 0

class CoordSystemLevel(Enum):
    """Coordinate system level."""
    CHROMOSOME = 1
    TOPLEVEL = 2
    CONTIG = 3

class ExonPhase(Enum):
    """Phase of an exon."""
    ZERO = 0
    ONE = 1
    TWO = 2
    UNKNOWN = -1

class ClassSOTerm(Enum):
    """Sequence Ontology term describing the class of a variation."""
    SEQUENCE_ALTERATION = "sequence_alteration"
    SNV = "SNV"
    INSERTION = "insertion"
    DELETION = "deletion"
    INDEL = "indel"
    SUBSTITUTION = "substitution"
    INVERSION = "inversion"
    COPY_NUMBER_VARIATION = "copy_number_variation"
    TANDEM_REPEAT = "tandem_repeat"