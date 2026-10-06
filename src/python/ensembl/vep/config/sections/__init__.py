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
"""The configuration sections."""

from ensembl.vep.config.sections.annotation import AnnotationConfig
from ensembl.vep.config.sections.base import ConfigSection
from ensembl.vep.config.sections.bench import (
    BenchConfig,
    BenchMatrixConfig,
    BenchOracleConfig,
    BenchSweepConfig,
)
from ensembl.vep.config.sections.catalogue import (
    CatalogueAddConfig,
    CatalogueConfig,
    CatalogueCostsConfig,
    CatalogueRemoveConfig,
    CatalogueShowConfig,
    CatalogueValidateConfig,
    CatalogueWriteConfig,
)
from ensembl.vep.config.sections.corpus import (
    CorpusConfig,
    CorpusDescribeConfig,
    CorpusGoldenConfig,
    CorpusSaveConfig,
)
from ensembl.vep.config.sections.explain import ExplainConfig
from ensembl.vep.config.sections.execution import ExecutionConfig
from ensembl.vep.config.sections.logging import LoggingConfig
from ensembl.vep.config.sections.run import RunConfig
from ensembl.vep.config.sections.source import SourceConfig
from ensembl.vep.config.sections.stats import StatsConfig

__all__ = [
    "AnnotationConfig",
    "BenchConfig",
    "BenchMatrixConfig",
    "BenchOracleConfig",
    "BenchSweepConfig",
    "ConfigSection",
    "CatalogueAddConfig",
    "CatalogueConfig",
    "CatalogueCostsConfig",
    "CatalogueRemoveConfig",
    "CatalogueShowConfig",
    "CatalogueValidateConfig",
    "CatalogueWriteConfig",
    "CorpusConfig",
    "CorpusDescribeConfig",
    "CorpusGoldenConfig",
    "CorpusSaveConfig",
    "ExplainConfig",
    "ExecutionConfig",
    "LoggingConfig",
    "RunConfig",
    "SourceConfig",
    "StatsConfig",
]
