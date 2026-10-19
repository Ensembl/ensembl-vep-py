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

TOML is the native format and its tables map onto section names. A file with no TOML table is read
as a Perl VEP config file, the way Perl VEP's own ``read_config_from_file`` reads it, so a config
file from an existing VEP install works unchanged. Legacy keys are flat, and are resolved to
sections by looking up which one declares each option.
"""

from __future__ import annotations

import logging
import re
import tomllib
from pathlib import Path
from ensembl.vep.config.option_spec import section_options
from ensembl.vep.config.parser.common import resolve_option
from ensembl.vep.config.sections.base import ConfigSection
from ensembl.vep.config.sources import ConfigLayer, coerce_layer

__all__ = ["read_config_file", "read_legacy_config_file"]

logger = logging.getLogger(__name__)

_TABLE_HEADER = re.compile(
    r"""^\s*\[\[?  # Optional whitespace and one or two opening square brackets
    [^\[\]]+       # One or more non-square bracket characters
    \]\]?          # One or two closing square brackets
    \s*(\#.*)?$    # Optional whitespace and optional comment until end of line/string""",
    flags=(re.VERBOSE | re.MULTILINE),
)
"""A TOML ``[table]`` header line, which no legacy config line can look like."""

_QUOTES = "\"'"


def _split_legacy_line(line: str) -> list[str]:
    """Split a legacy config line into words on whitespace, dropping quotes and keeping quoted spaces."""
    words: list[str] = []
    current: list[str] = []
    in_word = False
    open_quote = ""
    for character in line:
        if open_quote:
            if character == open_quote:
                open_quote = ""
            else:
                current.append(character)
        elif character in _QUOTES:
            open_quote = character
            in_word = True
        elif character.isspace():
            if in_word:
                words.append("".join(current))
                current = []
                in_word = False
        else:
            current.append(character)
            in_word = True
    if in_word:
        words.append("".join(current))
    return words


def _accumulates(section_type: type[ConfigSection], field_name: str) -> bool:
    """Return whether repeated lines add to a field rather than replace it, as Perl's ALLOW_MULTIPLE.

    Only an option declared ``repeatable`` qualifies -- Perl's ``custom`` and ``plugin``. Every
    other tuple option reads one comma-separated value, exactly as it does on the command line.
    """
    _section_field, spec = section_options(section_type)[field_name]
    return spec.repeatable


def read_legacy_config_file(path: Path, section_types: tuple[type[ConfigSection], ...]) -> ConfigLayer:
    """Read a Perl VEP config file as Perl VEP's ``read_config_from_file`` does.

    A line is a comment only if it starts with ``#``. The key and value are separated by any run of
    whitespace, a leading dash on the key is dropped, and quotes are removed from the value with the
    spaces inside them kept. The first occurrence of a repeated key wins, except for an option
    declared ``repeatable`` (Perl's ``custom`` and ``plugin``), which accumulates one item per line.
    ``option = value`` is accepted as well.

    Args:
        path (Path): The file to read.
        section_types (tuple[type[ConfigSection], ...]): The sections in play.

    Returns:
        ConfigLayer: The options the file sets, grouped by section and coerced. Unrecognised keys
            are ignored: a legacy config carries options this implementation has not grown yet.

    Raises:
        ValueError: If a value is of the wrong type or not one of the option's choices.
    """
    by_name = {section_type.section_name: section_type for section_type in section_types}
    layer: ConfigLayer = {}
    accumulated: dict[tuple[str, str], list[object]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        words = _split_legacy_line(line)
        key = words[0]
        values = words[1:]
        # Perl VEP writes "option value", but a hand-written file may use "option = value".
        if values and values[0] == "=":
            values = values[1:]
        elif values and values[0].startswith("="):
            values[0] = values[0][1:]
        first = values[0] if values else "1"
        resolved = resolve_option(key, first, section_types)
        if resolved is None:
            continue
        section, field_name, value = resolved
        if _accumulates(by_name[section], field_name):
            if values:
                # A constant flag resolves to its own value; anything else keeps the whole line.
                accumulated.setdefault((section, field_name), []).append(
                    " ".join(values) if value is first else value
                )
        else:
            layer.setdefault(section, {}).setdefault(field_name, value)
    for (section, field_name), tuple_values in accumulated.items():
        layer.setdefault(section, {})[field_name] = tuple_values
    return coerce_layer(layer, section_types)


def _read_toml_tables(
    parsed: dict[str, object], section_types: tuple[type[ConfigSection], ...], path: Path
) -> ConfigLayer:
    """Collect the options a parsed TOML document sets, warning about keys no section declares."""
    section_names = {section_type.section_name for section_type in section_types}
    tables: ConfigLayer = {}
    for section, values in parsed.items():
        if not isinstance(values, dict):
            continue
        for key, value in values.items():
            resolved = resolve_option(str(key), value, section_types)
            if resolved is None or resolved[0] != section:
                if section in section_names:
                    logger.warning("%s: ignoring unknown option '%s' in table [%s]", path, key, section)
                continue
            resolved_section, field_name, resolved_value = resolved
            tables.setdefault(resolved_section, {})[field_name] = resolved_value
    return coerce_layer(tables, section_types)


def read_config_file(path: Path, section_types: tuple[type[ConfigSection], ...]) -> ConfigLayer:
    """Read a configuration file, as TOML if it is TOML-shaped and as legacy config otherwise.

    A file is TOML if it has a ``[table]`` header or parses as TOML with at least one table. A
    syntax error in such a file is an error, not a cue to reread it as legacy config.

    Args:
        path (Path): The file to read.
        section_types (tuple[type[ConfigSection], ...]): The sections in play.

    Returns:
        ConfigLayer: The options the file sets, coerced to their declared types.

    Raises:
        FileNotFoundError: If ``path`` does not exist.
        ValueError: If a TOML file has a syntax error, or a value is invalid for its option.
    """
    text = path.read_text(encoding="utf-8")
    has_header = _TABLE_HEADER.search(text) is not None
    try:
        parsed = tomllib.loads(text)
    except tomllib.TOMLDecodeError as invalid:
        if has_header:
            raise ValueError(f"{path}: invalid TOML: {invalid}") from invalid
        return read_legacy_config_file(path, section_types)
    if not has_header and not any(isinstance(value, dict) for value in parsed.values()):
        return read_legacy_config_file(path, section_types)
    return _read_toml_tables(parsed, section_types, path)
