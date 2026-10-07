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
"""Configuration layers and how they combine."""

from __future__ import annotations

import logging
from collections.abc import Callable, Iterable
from enum import Enum
from pathlib import Path
from typing import Any, get_args, get_origin

from ensembl.vep.config.option_spec import OptionSpec, section_options, unwrap_optional_type
from ensembl.vep.config.sections.base import ConfigSection

__all__ = [
    "ConfigLayer",
    "ConfigSource",
    "coerce_layer",
    "coerce_value",
    "merge_layers",
]

logger = logging.getLogger(__name__)

ConfigLayer = dict[str, dict[str, object]]
"""One layer of configuration: section name -> field name -> value, carrying only what it sets."""

TRUE_VALUES = frozenset({"1", "true", "yes", "on"})
FALSE_VALUES = frozenset({"0", "false", "no", "off"})


class ConfigSource(Enum):
    """Where a configuration value came from."""

    DEFAULT = "default"
    FILE = "file"
    ENVIRONMENT = "environment"
    COMMAND_LINE = "command line"


def _to_bool(raw: object) -> bool:
    """Read a boolean from a value that may have arrived as text.

    Raises:
        ValueError: If the text is not a recognised truth value.
    """
    if isinstance(raw, bool):
        return raw
    text = str(raw).strip().lower()
    if text in TRUE_VALUES:
        return True
    if text in FALSE_VALUES:
        return False
    raise ValueError(f"expected a boolean, got '{raw}'")


def _to_tuple(raw: object, type_annotation: Any) -> tuple:
    """Read a tuple from a comma-separated string or an existing sequence."""
    item_types: list[Any] = []
    for argument in get_args(type_annotation):
        if argument is Ellipsis:
            continue
        item_types.append(argument)
    item_type = item_types[0] if item_types else str
    pieces: list[object] = []
    if isinstance(raw, str):
        for piece in raw.split(","):
            stripped = piece.strip()
            if stripped:
                pieces.append(stripped)
    elif isinstance(raw, Iterable):
        pieces = list(raw)
    else:
        # One value for a field that takes several: a TOML "distance = 5000", or an int from a
        # Python caller. Perl VEP's own comma options accept the single form the same way.
        pieces = [raw]
    converted: list[Any] = []
    for value in pieces:
        converted.append(item_type(value))
    return tuple(converted)


def _to_enum(raw: object, enum_type: type[Enum]) -> Enum:
    """Read an enum member from itself or from its name, in any case.

    Raises:
        ValueError: If the text names no member.
    """
    if isinstance(raw, enum_type):
        return raw
    try:
        return enum_type[str(raw).upper()]
    except KeyError:
        accepted = ", ".join(member.name.lower() for member in enum_type)
        raise ValueError(f"{str(raw)!r} is not one of {accepted}") from None


def _converter_for(type_annotation: Any) -> Callable[[object], object]:
    """Return the converter that reads one declared type.

    Anything unrecognised is passed through untouched: coercion is best-effort, not a validator.
    """
    if type_annotation is bool:
        return _to_bool
    if type_annotation is Path:
        return lambda raw: Path(str(raw))
    if isinstance(type_annotation, type) and issubclass(type_annotation, Enum):
        return lambda raw: _to_enum(raw, type_annotation)
    if get_origin(type_annotation) is tuple:
        return lambda raw: _to_tuple(raw, type_annotation)
    if type_annotation in (int, float, str):
        return type_annotation
    return lambda raw: raw


def coerce_value(raw: object, type_annotation: Any) -> object:
    """Coerce a value read from a file or the environment to a field's declared type.

    Command-line values have already been converted by argparse. These have not: TOML gives real
    ints and bools but strings for paths, and the environment gives strings for everything.

    Args:
        raw (object): The value as read.
        type_annotation (Any): The field's resolved type.

    Returns:
        object: The coerced value.

    Raises:
        ValueError: If a string cannot be read as the declared type.
    """
    type_annotation, optional = unwrap_optional_type(type_annotation)
    if optional and (raw is None or raw == ""):
        return None
    return _converter_for(type_annotation)(raw)


def _check_choices(section: str, field_name: str, value: object, spec: OptionSpec) -> None:
    """Check a coerced value against the option's declared choices, if it has any.

    Raises:
        ValueError: If the value, or any member of a tuple value, is not a permitted choice.
    """
    choices = spec.resolved_choices()
    if choices is None or value is None or isinstance(value, Enum):
        return
    members = value if isinstance(value, tuple) else (value,)
    for member in members:
        if member not in choices:
            accepted = ", ".join(str(choice) for choice in choices)
            raise ValueError(f"{section}.{field_name}: {member!r} is not one of {accepted}")


def _coerce_repeatable(raw: object, type_annotation: Any) -> tuple:
    """Coerce a repeatable option's value without splitting an item on its commas.

    One ``custom`` spec is itself comma-separated, so a lone string is one item, and a list or
    tuple is one item per element -- the TOML spelling of giving the option twice.
    """
    item_type, _optional = unwrap_optional_type(type_annotation)
    arguments = [argument for argument in get_args(item_type) if argument is not Ellipsis]
    converter = _converter_for(arguments[0] if arguments else str)
    items = raw if isinstance(raw, (list, tuple)) else [raw]
    return tuple(converter(item) for item in items)


def coerce_layer(layer: ConfigLayer, section_types: tuple[type[ConfigSection], ...]) -> ConfigLayer:
    """Coerce every value in a layer to its declared type and check it against its choices.

    Unknown sections are dropped silently, since a layer may cover sections the command does not
    read. An unknown field in a known section is dropped with a warning.

    Raises:
        ValueError: If a value cannot be read as its declared type or is not a permitted choice. The
            message names the section and field.
    """
    by_name = {section_type.section_name: section_type for section_type in section_types}
    coerced: ConfigLayer = {}
    for section, values in layer.items():
        section_type = by_name.get(section)
        if section_type is None:
            continue
        options = section_options(section_type)
        for field_name, raw in values.items():
            if field_name not in options:
                logger.warning("Ignoring unknown option '%s' in section '%s'", field_name, section)
                continue
            section_field, spec = options[field_name]
            try:
                if spec.repeatable:
                    value: object = _coerce_repeatable(raw, section_field.type)
                else:
                    value = coerce_value(raw, section_field.type)
            except ValueError as invalid:
                raise ValueError(f"{section}.{field_name}: {invalid}") from invalid
            _check_choices(section, field_name, value, spec)
            coerced.setdefault(section, {})[field_name] = value
    return coerced


def merge_layers(*layers: ConfigLayer) -> ConfigLayer:
    """Merge configuration layers in increasing precedence order."""
    merged: ConfigLayer = {}
    for layer in layers:
        for section, values in layer.items():
            merged.setdefault(section, {}).update(values)
    return merged
