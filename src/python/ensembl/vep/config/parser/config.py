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
"""Reading configuration from a file.

TOML is the native format and its tables map onto section names. A file that is not valid TOML is
read as legacy VEP config -- one ``option value`` per line -- so a config file from an existing VEP
install works unchanged. Legacy keys are flat, and are resolved to sections by looking up which one
declares each option.
"""

from __future__ import annotations

import tomllib
from pathlib import Path

from ensembl.vep.config.sections.base import ConfigSection
from ensembl.vep.config.parser.common import resolve_option
from ensembl.vep.config.sources import ConfigLayer, coerce_layer

__all__ = ["read_config_file", "read_legacy_config_file"]


def read_legacy_config_file(path: Path, section_types: tuple[type[ConfigSection], ...]) -> ConfigLayer:
    """Read a legacy VEP config file: one ``option value`` per line, ``#`` comments.

    Args:
        path (Path): The file to read.
        section_types (tuple[type[ConfigSection], ...]): The sections in play.

    Returns:
        ConfigLayer: The options the file sets, grouped by section and coerced. Unrecognised keys
            are ignored: a legacy config carries options this implementation has not grown yet.
    """
    layer: ConfigLayer = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.split("#", 1)[0].strip()
        if not stripped:
            continue
        key, _separator, raw_value = stripped.partition(" ")
        # Legacy VEP writes "option value", but a hand-written file may use "option = value".
        raw_value = raw_value.strip().removeprefix("=").strip()
        resolved = resolve_option(key, raw_value or "1", section_types)
        if resolved is None:
            continue
        section, field_name, value = resolved
        layer.setdefault(section, {})[field_name] = value
    return coerce_layer(layer, section_types)


def read_config_file(path: Path, section_types: tuple[type[ConfigSection], ...]) -> ConfigLayer:
    """Read a configuration file, as TOML if it parses and as legacy config otherwise.

    Args:
        path (Path): The file to read.
        section_types (tuple[type[ConfigSection], ...]): The sections in play.

    Returns:
        ConfigLayer: The options the file sets, coerced to their declared types.

    Raises:
        FileNotFoundError: If ``path`` does not exist.
    """
    text = path.read_text(encoding="utf-8")
    try:
        parsed = tomllib.loads(text)
    except tomllib.TOMLDecodeError:
        return read_legacy_config_file(path, section_types)
    tables: ConfigLayer = {}
    for section, values in parsed.items():
        if not isinstance(values, dict):
            continue
        for key, value in values.items():
            resolved = resolve_option(str(key), value, section_types)
            if resolved is None:
                continue
            resolved_section, field_name, resolved_value = resolved
            if resolved_section != str(section):
                continue
            tables.setdefault(resolved_section, {})[field_name] = resolved_value
    if not tables:
        return read_legacy_config_file(path, section_types)
    return coerce_layer(tables, section_types)
