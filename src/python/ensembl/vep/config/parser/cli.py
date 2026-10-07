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
from ensembl.vep.config.parser.common import flag_for, option_spellings
from ensembl.vep.config.sources import ConfigLayer

__all__ = [
    "add_section_arguments",
    "destination",
    "build_parser",
    "check_no_flag_collisions",
    "parse_args",
]


class BooleanSwitchAction(argparse.Action):
    """Set a boolean to True or False according to which spelling was given.

    One action carries both senses so ``--help`` lists them together, ``--zone_map, --no_zone_map``.
    """

    def __init__(
        self,
        option_strings: list[str],
        dest: str,
        *,
        false_flags: tuple[str, ...] = (),
        default: Any = argparse.SUPPRESS,
        help: str | None = None,  # pylint: disable=redefined-builtin
    ) -> None:
        super().__init__(option_strings=option_strings, dest=dest, nargs=0, default=default, help=help)
        self.false_flags = frozenset(false_flags)

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        values: Any,
        option_string: str | None = None,
    ) -> None:
        setattr(namespace, self.dest, option_string not in self.false_flags)


class AppendItemAction(argparse.Action):
    """Collect each occurrence of a repeatable option into a tuple, keeping every value whole."""

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        values: Any,
        option_string: str | None = None,
    ) -> None:
        collected = getattr(namespace, self.dest, ())
        setattr(namespace, self.dest, (*collected, values))


def opposite_spelling(spelling: str) -> str:
    """Return the spelling that sets a boolean the other way: ``zone_map`` <-> ``no_zone_map``."""
    if spelling.startswith("no_"):
        return spelling.removeprefix("no_")
    return f"no_{spelling}"


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


def destination(section_name: str, field_name: str) -> str:
    """Return the parsed-namespace key one section's option is read back from.

    Namespaced by section because two sections of one command may declare the same field name --
    ``format`` is the input format to ``source`` and the report format to ``stats`` -- and a bare
    field name would land one value in both.
    """

    return f"{section_name}__{field_name}"


def _add_argument(
    parser: argparse._ActionsContainer, section_name: str, section_field: Any, spec: OptionSpec
) -> None:
    """Add one declared option to a parser or argument group."""
    flag = flag_for(section_field.name, spec)
    flags = [f"--{flag}"]
    for synonym in spec.synonyms:
        flags.append(f"--{synonym}")
    if spec.short:
        flags.insert(0, f"-{spec.short}")

    type_annotation, _optional = unwrap_optional_type(section_field.type)
    if type_annotation is bool:
        _add_boolean_arguments(parser, section_name, section_field, spec, flags)
        _add_constant_arguments(parser, section_name, section_field, spec, flag)
        return

    keywords: dict[str, Any] = {
        "dest": destination(section_name, section_field.name),
        "default": argparse.SUPPRESS,
        "help": argparse.SUPPRESS if spec.hidden else spec.help,
    }

    if isinstance(type_annotation, type) and issubclass(type_annotation, Enum):
        keywords["type"] = _enum_converter(type_annotation)
        keywords["choices"] = list(type_annotation)
        keywords["metavar"] = "{" + ",".join(member.name.lower() for member in type_annotation) + "}"
    elif get_origin(type_annotation) is tuple:
        item_types: list[Any] = []
        for argument in get_args(type_annotation):
            if argument is Ellipsis:
                continue
            item_types.append(argument)
        item_type = item_types[0] if item_types else str
        if spec.repeatable:
            keywords["action"] = AppendItemAction
            keywords["type"] = item_type
        else:
            keywords["type"] = _comma_separated(item_type)
    elif type_annotation is Path:
        keywords["type"] = Path
    else:
        keywords["type"] = type_annotation

    choices = spec.resolved_choices()
    if choices is not None and "choices" not in keywords:
        keywords["choices"] = list(choices)
    if spec.metavar is not None:
        keywords["metavar"] = spec.metavar
    elif (
        "metavar" not in keywords
        and "choices" not in keywords
        and keywords.get("action") in (None, AppendItemAction)
    ):
        # argparse would name the value after the namespaced dest otherwise.
        keywords["metavar"] = flag.upper()

    parser.add_argument(*flags, **keywords)
    _add_constant_arguments(parser, section_name, section_field, spec, flag)


def _add_boolean_arguments(
    parser: argparse._ActionsContainer,
    section_name: str,
    section_field: Any,
    spec: OptionSpec,
    flags: list[str],
) -> None:
    """Add the spellings of a boolean option, in Perl VEP's underscore style.

    A positive option (``zone_map``) takes ``--zone_map`` and ``--no_zone_map``, plus the hyphenated
    ``--no-zone_map`` as a hidden alias. A negative-sense option (``no_stats``) keeps ``--no_stats``
    and gains a hidden ``--stats`` that turns it back off, named in the help of the first.
    """
    dest = destination(section_name, section_field.name)
    flag = flag_for(section_field.name, spec)
    help_text = argparse.SUPPRESS if spec.hidden else spec.help
    opposites = [
        f"--{opposite_spelling(spelling.lstrip('-'))}" for spelling in flags if spelling.startswith("--")
    ]

    if flag.startswith("no_"):
        if not spec.hidden:
            help_text = f"{spec.help} ({', '.join(opposites)} turns it back on)"
        parser.add_argument(*flags, dest=dest, action="store_true", default=argparse.SUPPRESS, help=help_text)
        parser.add_argument(
            *opposites, dest=dest, action="store_false", default=argparse.SUPPRESS, help=argparse.SUPPRESS
        )
        return

    parser.add_argument(
        *flags,
        *opposites,
        dest=dest,
        action=BooleanSwitchAction,
        false_flags=tuple(opposites),
        default=argparse.SUPPRESS,
        help=help_text,
    )
    hyphenated = [opposite.replace("no_", "no-", 1) for opposite in opposites]
    parser.add_argument(
        *hyphenated, dest=dest, action="store_false", default=argparse.SUPPRESS, help=argparse.SUPPRESS
    )


def _add_constant_arguments(
    parser: argparse._ActionsContainer,
    section_name: str,
    section_field: Any,
    spec: OptionSpec,
    flag: str,
) -> None:
    """Add the flags that set one option to a fixed value, as Perl VEP's ``--vcf`` does."""

    for constant_flag, constant in spec.constants.items():
        parser.add_argument(
            f"--{constant_flag}",
            dest=destination(section_name, section_field.name),
            action="store_const",
            const=constant,
            default=argparse.SUPPRESS,
            help=argparse.SUPPRESS if spec.hidden else f"shorthand for --{flag} {constant}",
        )


def add_section_arguments(parser: argparse.ArgumentParser, section_type: type[ConfigSection]) -> None:
    """Add every option a section declares to a parser, as its own argument group."""
    group = parser.add_argument_group(f"{section_type.section_name} options")
    for section_field, spec in section_options(section_type).values():
        _add_argument(group, section_type.section_name, section_field, spec)


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
        for field_name, (section_field, spec) in section_options(section_type).items():
            for spelling in option_spellings(field_name, spec):
                if spelling == field_name and spelling != flag_for(field_name, spec):
                    # The field name is not a flag when the option renames itself.
                    continue
                owners.setdefault(spelling, []).append(section_type.section_name)
            if unwrap_optional_type(section_field.type)[0] is bool:
                for spelling in (flag_for(field_name, spec), *spec.synonyms):
                    owners.setdefault(opposite_spelling(spelling), []).append(section_type.section_name)
            for constant_flag in spec.constants:
                owners.setdefault(constant_flag, []).append(section_type.section_name)
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
            key = destination(section_type.section_name, name)
            if key in given:
                values[name] = given[key]
        if values:
            layer[section_type.section_name] = values
    return layer
