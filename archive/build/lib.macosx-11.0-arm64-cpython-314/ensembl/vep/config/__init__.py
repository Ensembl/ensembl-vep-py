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
"""Configuration declarations, sections and loaders."""

from ensembl.vep.config.base import Config, load
from ensembl.vep.config.option_spec import (
    OptionSpec,
    option,
    section_options,
    unwrap_optional_type,
)
from ensembl.vep.config.sections import (
    AnnotationConfig,
    BenchConfig,
    CatalogueAddConfig,
    CatalogueConfig,
    CatalogueCostsConfig,
    CatalogueRemoveConfig,
    CatalogueShowConfig,
    CatalogueValidateConfig,
    CatalogueWriteConfig,
    ConfigSection,
    CorpusConfig,
    ExecutionConfig,
    ExplainConfig,
    LoggingConfig,
    RunConfig,
    SourceConfig,
    StatsConfig,
)
from ensembl.vep.config.sources import ConfigLayer, ConfigSource

__all__ = [
    "AnnotationConfig",
    "BenchConfig",
    "CatalogueAddConfig",
    "CatalogueConfig",
    "CatalogueCostsConfig",
    "CatalogueRemoveConfig",
    "CatalogueShowConfig",
    "CatalogueValidateConfig",
    "CatalogueWriteConfig",
    "Config",
    "ConfigLayer",
    "ConfigSection",
    "ConfigSource",
    "CorpusConfig",
    "ExecutionConfig",
    "ExplainConfig",
    "LoggingConfig",
    "OptionSpec",
    "RunConfig",
    "SourceConfig",
    "StatsConfig",
    "load",
    "option",
    "section_options",
    "unwrap_optional_type",
]
