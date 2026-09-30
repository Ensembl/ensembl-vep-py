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
"""Common constants."""

from ensembl.vep.common.enums import ExonPhase, Strand, ClassSOTerm

CLASS_SO_TERM_ABBRV_TO_ENUM: dict[str, ClassSOTerm] = {
    "INS": ClassSOTerm.INSERTION,
    "INS_ME": ClassSOTerm.MOBILE_ELEMENT_INSERTION,
    "INS_ALU": ClassSOTerm.ALU_INSERTION,
    "INS_HERV": ClassSOTerm.HERV_INSERTION,
    "INS_LINE1": ClassSOTerm.LINE1_INSERTION,
    "INS_SVA": ClassSOTerm.SVA_INSERTION,
    "DEL": ClassSOTerm.DELETION,
    "DEL_ME": ClassSOTerm.MOBILE_ELEMENT_DELETION,
    "DEL_ALU": ClassSOTerm.ALU_DELETION,
    "DEL_HERV": ClassSOTerm.HERV_DELETION,
    "DEL_LINE1": ClassSOTerm.LINE1_DELETION,
    "DEL_SVA": ClassSOTerm.SVA_DELETION,
    "TREP": ClassSOTerm.TANDEM_REPEAT,
    "TDUP": ClassSOTerm.TANDEM_DUPLICATION,
    "DUP": ClassSOTerm.DUPLICATION,
    "CNV": ClassSOTerm.COPY_NUMBER_VARIATION,
    "INV": ClassSOTerm.INVERSION,
    "BND": ClassSOTerm.CHROMOSOME_BREAKPOINT,
}

# for . (non-strandedness) we treat it as forward strand
STRAND_SYMBOL_TO_ENUM = {
    "+": Strand.FORWARD,
    "-": Strand.REVERSE,
    ".": Strand.FORWARD,
    "?": Strand.UNKNOWN,
}

GFF_TO_ENSEMBL_PHASE = {
    0: ExonPhase.ZERO,
    2: ExonPhase.ONE,
    1: ExonPhase.TWO,
}
"""A CDS row's GFF phase (bases to skip to reach a codon) as an Ensembl exon phase (bases of the
first codon on what came before), converted as ``Bio::EnsEMBL::VEP::AnnotationSource::File::
BaseGXF::_convert_phase`` does."""
