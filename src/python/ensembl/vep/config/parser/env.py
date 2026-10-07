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
"""Reading configuration from the environment."""

from __future__ import annotations

import os

from ensembl.vep.config.option_spec import section_options
from ensembl.vep.config.sections.base import ConfigSection
from ensembl.vep.config.sources import ConfigLayer, coerce_layer

__all__ = ["ENV_PREFIX", "default_env_name", "read_environment"]

ENV_PREFIX = "ENSEMBL_VEP"


def default_env_name(section: str, field_name: str) -> str:
    """Return the environment variable an option is read from, e.g. ENSEMBL_VEP_CORPUS_FUZZ."""
    return f"{ENV_PREFIX}_{section}_{field_name}".upper()


def read_environment(
    section_types: tuple[type[ConfigSection], ...], environ: dict[str, str] | None = None
) -> ConfigLayer:
    """Read whatever the environment sets, for the given sections.

    Args:
        section_types (tuple[type[ConfigSection], ...]): The sections in play.
        environ (dict[str, str] | None): Environment to read, or None for ``os.environ``.

    Returns:
        ConfigLayer: The options the environment sets, coerced to their declared types.

    Raises:
        ValueError: If a value is not of its declared type or not one of its choices.
    """
    source = os.environ if environ is None else environ
    layer: ConfigLayer = {}
    for section_type in section_types:
        for field_name, (_field, spec) in section_options(section_type).items():
            name = spec.env or default_env_name(section_type.section_name, field_name)
            if name in source:
                layer.setdefault(section_type.section_name, {})[field_name] = source[name]
    return coerce_layer(layer, section_types)
