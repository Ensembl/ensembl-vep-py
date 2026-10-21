# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Common enums."""

from enum import Enum, IntEnum


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
    """Sequence Ontology term describing the class of a variation.

    The first block is what ``infer_variant_so_term`` can derive from an allele string alone --
    the classes a short variant can have. The second block is structural-variant only which may
    may need INFO/SVTYPE to be derived.
    """

    SEQUENCE_ALTERATION = "sequence_alteration"
    SNV = "SNV"
    INSERTION = "insertion"
    DELETION = "deletion"
    INDEL = "indel"
    SUBSTITUTION = "substitution"
    INVERSION = "inversion"
    COPY_NUMBER_VARIATION = "copy_number_variation"
    TANDEM_REPEAT = "tandem_repeat"

    # Structural-variant classes. Never produced by infer_variant_so_term.
    DUPLICATION = "duplication"
    TANDEM_DUPLICATION = "tandem_duplication"
    COPY_NUMBER_GAIN = "copy_number_gain"
    COPY_NUMBER_LOSS = "copy_number_loss"
    MOBILE_ELEMENT_INSERTION = "mobile_element_insertion"
    ALU_INSERTION = "Alu_insertion"
    HERV_INSERTION = "HERV_insertion"
    LINE1_INSERTION = "LINE1_insertion"
    SVA_INSERTION = "SVA_insertion"
    MOBILE_ELEMENT_DELETION = "mobile_element_deletion"
    ALU_DELETION = "Alu_deletion"
    HERV_DELETION = "HERV_deletion"
    LINE1_DELETION = "LINE1_deletion"
    SVA_DELETION = "SVA_deletion"
    NOVEL_SEQUENCE_INSERTION = "novel_sequence_insertion"
    TRANSLOCATION = "translocation"
    INTERCHROMOSOMAL_BREAKPOINT = "interchromosomal_breakpoint"
    INTRACHROMOSOMAL_BREAKPOINT = "intrachromosomal_breakpoint"
    CHROMOSOME_BREAKPOINT = "chromosome_breakpoint"
    LOSS_OF_HETEROZYGOSITY = "loss_of_heterozygosity"
    PROBE = "probe"


class ConsequenceImpact(Enum):
    """Impact of a consequence."""

    HIGH = "HIGH"
    MODERATE = "MODERATE"
    LOW = "LOW"
    MODIFIER = "MODIFIER"


class SourceType(Enum):
    """Type of annotation/input source."""

    FILE = "file"
    DB = "db"
    CACHE = "cache"


class StartMethod(str, Enum):
    """Multiprocessing worker start methods."""

    AUTO = "auto"
    SPAWN = "spawn"
    FORK = "fork"
    FORKSERVER = "forkserver"


class Readthrough(IntEnum):
    """How far translation may run when a variant destroys the stop codon.

    Opt-in, and ``NONE`` by default. An ``IntEnum`` because the levels are ordered rather
    than independent, each one extends the previous.
    """

    NONE = 0
    """Stop at the (possibly shifted) CDS end. Translation never enters the UTR."""

    UTR = 1
    """Continue into the 3' UTR, looking for the next in-frame stop. Stays entirely within the
    transcript's own spliced sequence."""

    BEYOND_TRANSCRIPT = 2
    """As ``UTR``, and if that still yields no stop, continue into raw genomic sequence past the
    transcript's own boundary, read directly from the attached slice (bounded by 
    ``MAX_READTHROUGH_BASES``). Trusts that the attached slice already covers the flank."""


class CostTier(IntEnum):
    """The declared cost tiers a fact can be registered with, cheapest first.

    An ``IntEnum`` rather than a plain ``int`` alias, so a tier is both a real integer (comparable,
    summable, sortable) and a self-naming, enumerable value: ``CostTier.MAPPING.name`` instead of
    a hand-maintained ``{20: "MAPPING"}`` lookup, and ``list(CostTier)`` for anything that wants
    every tier without duplicating this declaration.

    The numbers are relative in magnitude and not absolute. The declared relative magnitudes do
    not match measured wall-clock time. Treat the tier ordering as load-bearing and the numeric
    ratios as illustrative only.
    """

    FREE = 0
    """Arithmetic on the variant alone."""
    BOUNDS = 1
    """Comparing against a single bound (e.g. transcript start/end)."""
    INTERVAL = 5
    """Scanning exon/intron interval lists."""
    MAPPING = 20
    """Genomic -> CDS coordinate mapping."""
    SEQUENCE = 60
    """Extracting sequence, splicing in the alt allele, translating."""
    TRANSCRIPT_SEQUENCE = 200
    """Rebuilding or translating the whole CDS or cDNA, not just the touched codon window."""

    @classmethod
    def label(cls, cost: int) -> str:
        """Return the name of the tier ``cost`` belongs to, or the number itself.

        A declared cost is a member, so this is ``cost.name``. The measured-cost override
        (``data/fact_costs.json``, opt-in) sets arbitrary profiled integers that match no tier, so
        the numeric fallback is the normal case there rather than an error.
        """
        try:
            return cls(cost).name
        except ValueError:
            return str(cost)
