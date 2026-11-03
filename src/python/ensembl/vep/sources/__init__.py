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
"""File-backed and database-backed feature sources."""

from __future__ import annotations

from ensembl.vep.sources.feature_source import FeatureSource
from ensembl.vep.sources.file.fasta.sequence import FastaSequenceSource, UNKNOWN_REGION_LENGTH
from ensembl.vep.sources.file.file_feature_source import FileFeatureSource
from ensembl.vep.sources.file.gff3.transcript import GFF3TranscriptSource
from ensembl.vep.sources.file.gtf.transcript import GTFTranscriptSource
from ensembl.vep.sources.file.vcf.variation import PicklableVCFVariationSource, VCFVariationSource
from ensembl.vep.sources.sequence_source import SequenceSource
from .synonym import ChromosomeSynonyms, SequenceRegionResolver

__all__ = [
    "ChromosomeSynonyms",
    "FeatureSource",
    "FileFeatureSource",
    "FastaSequenceSource",
    "GFF3TranscriptSource",
    "GTFTranscriptSource",
    "PicklableVCFVariationSource",
    "SequenceSource",
    "SequenceRegionResolver",
    "UNKNOWN_REGION_LENGTH",
    "VCFVariationSource",
]
