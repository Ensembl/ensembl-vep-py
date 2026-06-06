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

from enum import Enum
from ensembl.vep.ensembl_vep_core import Strand, ExonPhase # pylint: disable=E0611,W0611


class CoordSystemLevel(Enum):
    """Coordinate system level."""

    CHROMOSOME = 1
    TOPLEVEL = 2
    CONTIG = 3


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


class ConsequenceImpact(Enum):
    """Impact of a consequence."""

    HIGH = "HIGH"
    MODERATE = "MODERATE"
    LOW = "LOW"
    MODIFIER = "MODIFIER"
