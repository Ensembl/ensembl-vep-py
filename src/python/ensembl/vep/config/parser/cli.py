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
"""Build an argparse surface from the section declarations.

Each action uses ``default=argparse.SUPPRESS`` so omitted flags stay absent from the parsed
namespace. That keeps them distinct from values read from a file or the environment.
"""

from __future__ import annotations

import argparse
from enum import Enum
from pathlib import Path
from typing import Any, get_args, get_origin

from ensembl.vep.config.option_spec import OptionSpec, section_options, unwrap_optional_type
from ensembl.vep.config.sections.base import ConfigSection
from ensembl.vep.config.parser.common import flag_for
from ensembl.vep.config.sources import ConfigLayer

__all__ = [
    "add_section_arguments",
    "build_parser",
    "check_no_flag_collisions",
    "parse_args",
]


def _comma_separated(item_type: Any) -> Any:
    """Return a converter parsing a comma-separated string into a tuple."""

    def convert(value: str) -> tuple:
        if not value:
            return ()
        pieces: list[Any] = []
        for piece in value.split(","):
            pieces.append(item_type(piece.strip()))
        return tuple(pieces)

    convert.__name__ = f"comma_separated_{getattr(item_type, '__name__', 'value')}"
    return convert


def _enum_converter(enum_type: type[Enum]) -> Any:
    """Return a converter mapping a lowercased member name to an enum member."""

    def convert(value: str) -> Enum:
        try:
            return enum_type[value.upper()]
        except KeyError as error:
            names: list[str] = []
            for member in enum_type:
                names.append(member.name.lower())
            names_text = ", ".join(names)
            raise argparse.ArgumentTypeError(f"expected one of {names_text}, got '{value}'") from error

    convert.__name__ = enum_type.__name__.lower()
    return convert


def _add_argument(parser: argparse._ActionsContainer, section_field: Any, spec: OptionSpec) -> None:
    """Add one declared option to a parser or argument group."""
    flag = flag_for(section_field.name, spec)
    flags = [f"--{flag}"]
    if spec.short:
        flags.insert(0, f"-{spec.short}")

    type_annotation, _optional = unwrap_optional_type(section_field.type)
    keywords: dict[str, Any] = {
        "dest": section_field.name,
        "default": argparse.SUPPRESS,
        "help": argparse.SUPPRESS if spec.hidden else spec.help,
    }

    if type_annotation is bool:
        # A negative-sense flag (legacy VEP's --no_stats) has no positive counterpart, and
        # BooleanOptionalAction refuses to derive one for a name that already starts with "no".
        keywords["action"] = "store_true" if flag.startswith("no_") else argparse.BooleanOptionalAction
    elif isinstance(type_annotation, type) and issubclass(type_annotation, Enum):
        keywords["type"] = _enum_converter(type_annotation)
        keywords["choices"] = list(type_annotation)
        keywords["metavar"] = "{" + ",".join(member.name.lower() for member in type_annotation) + "}"
    elif get_origin(type_annotation) is tuple:
        item_types: list[Any] = []
        for argument in get_args(type_annotation):
            if argument is Ellipsis:
                continue
            item_types.append(argument)
        keywords["type"] = _comma_separated(item_types[0] if item_types else str)
    elif type_annotation is Path:
        keywords["type"] = Path
    else:
        keywords["type"] = type_annotation

    choices = spec.resolved_choices()
    if choices is not None and "choices" not in keywords:
        keywords["choices"] = list(choices)
    if spec.metavar is not None:
        keywords["metavar"] = spec.metavar

    parser.add_argument(*flags, **keywords)


def add_section_arguments(parser: argparse.ArgumentParser, section_type: type[ConfigSection]) -> None:
    """Add every option a section declares to a parser, as its own argument group."""
    group = parser.add_argument_group(f"{section_type.section_name} options")
    for section_field, spec in section_options(section_type).values():
        _add_argument(group, section_field, spec)


def check_no_flag_collisions(section_types: tuple[type[ConfigSection], ...]) -> None:
    """Raise if two of the given sections declare the same flag.

    ``output_file`` means the annotation output to ``run``, the report to ``stats`` and the
    explanation to ``explain``; that only matters when one command reads two such sections.
    Reported here because argparse's own message names neither section.

    Raises:
        ValueError: If a flag is declared by more than one of them.
    """
    owners: dict[str, list[str]] = {}
    for section_type in section_types:
        for field_name, (_field, spec) in section_options(section_type).items():
            owners.setdefault(flag_for(field_name, spec), []).append(section_type.section_name)
    collisions = {flag: sections for flag, sections in owners.items() if len(sections) > 1}
    if collisions:
        details: list[str] = []
        for flag, sections in sorted(collisions.items()):
            details.append(f"--{flag} in {', '.join(sections)}")
        detail = "; ".join(details)
        raise ValueError(f"Sections declare colliding flags: {detail}")


def build_parser(
    prog: str, description: str, section_types: tuple[type[ConfigSection], ...]
) -> argparse.ArgumentParser:
    """Return a parser exposing exactly the options the given sections declare.

    Raises:
        ValueError: If two of the sections declare the same flag.
    """
    check_no_flag_collisions(section_types)
    parser = argparse.ArgumentParser(
        prog=prog, description=description, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    for section_type in section_types:
        add_section_arguments(parser, section_type)
    return parser


def parse_args(
    parser: argparse.ArgumentParser, argv: list[str], section_types: tuple[type[ConfigSection], ...]
) -> ConfigLayer:
    """Parse command-line arguments into a sparse section -> field -> value mapping.

    Returns:
        ConfigLayer: Only the options actually present on the command line, so a caller can
            distinguish "not given" from "given the default value".
    """
    given = vars(parser.parse_args(argv))
    layer: ConfigLayer = {}
    for section_type in section_types:
        values: dict[str, object] = {}
        for name in section_options(section_type):
            if name in given:
                values[name] = given[name]
        if values:
            layer[section_type.section_name] = values
    return layer
