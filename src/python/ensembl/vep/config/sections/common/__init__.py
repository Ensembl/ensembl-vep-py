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
"""Options more than one command group reads.

A section shared across groups lives here so that one spelling, one default and one environment
variable serve every command that reads it; a section only one group's verbs share stays in that
group's own package. The synthetic corpus is here for that reason: explain, stats,
corpus, catalogue and bench all read it.
"""

from ensembl.vep.config.sections.common.annotation import AnnotationConfig
from ensembl.vep.config.sections.common.corpus import CorpusConfig
from ensembl.vep.config.sections.common.execution import ExecutionConfig
from ensembl.vep.config.sections.common.logging import LoggingConfig
from ensembl.vep.config.sections.common.source import SourceConfig

__all__ = [
    "AnnotationConfig",
    "CorpusConfig",
    "ExecutionConfig",
    "LoggingConfig",
    "SourceConfig",
]
