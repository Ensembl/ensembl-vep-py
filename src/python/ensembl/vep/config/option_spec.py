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
"""Option declarations for configuration sections."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from types import NoneType, UnionType
from typing import TYPE_CHECKING, Any, Union, get_args, get_origin

import attrs
from attrs import Factory, define, fields

if TYPE_CHECKING:
    from ensembl.vep.config.sections.base import ConfigSection

__all__ = ["OptionSpec", "option", "section_options", "unwrap_optional_type"]

METADATA_KEY = "vep.option"
"""Key the ``option()`` payload is stored under in a field's metadata."""


@define(frozen=True, kw_only=True, auto_attribs=True)
class OptionSpec:
    """How one option is exposed.

    Attributes:
        help (str): One-line description, used as the argparse help text.
        flag (str | None): Long flag spelling without the dashes, or None for the field name.
        short (str | None): Single-character short flag, without the dash.
        synonyms (tuple[str, ...]): Further long spellings that set the same option, for the names
            Perl VEP uses. They work on the command line and in a config file; the environment
            variable keeps the canonical name.
        constants (Mapping[str, object]): Flags that set this option to a fixed value, such as
            Perl VEP's ``--vcf`` for ``--output_format vcf``, as flag name -> value.
        choices (Sequence[object] | Callable[[], Sequence[object]] | None): The permitted values,
            or a callable returning them for a vocabulary another package owns.
        metavar (str | None): argparse metavar override.
        env (str | None): Environment variable name, or None to derive one.
        hidden (bool): True to keep the option out of ``--help``.
    """

    help: str
    flag: str | None = None
    short: str | None = None
    synonyms: tuple[str, ...] = ()
    constants: Mapping[str, object] = Factory(dict)
    choices: Sequence[object] | Callable[[], Sequence[object]] | None = None
    metavar: str | None = None
    env: str | None = None
    hidden: bool = False

    def resolved_choices(self) -> tuple[object, ...] | None:
        """Return the permitted values, calling ``choices`` if needed."""
        if self.choices is None:
            return None
        if callable(self.choices):
            choices = self.choices()
        else:
            choices = self.choices
        return tuple(choices)


def option(
    help: str,  # pylint: disable=redefined-builtin
    *,
    flag: str | None = None,
    short: str | None = None,
    synonyms: tuple[str, ...] = (),
    constants: Mapping[str, object] | None = None,
    choices: Sequence[object] | Callable[[], Sequence[object]] | None = None,
    metavar: str | None = None,
    env: str | None = None,
    hidden: bool = False,
) -> dict[str, OptionSpec]:
    """Return ``attrs`` field metadata for one option.

    Args:
        help (str): One-line description, used as the argparse help text.
        flag (str | None): Long flag spelling without the dashes, or None for the field name.
        short (str | None): Single-character short flag, without the dash.
        synonyms (tuple[str, ...]): Further long spellings that set the same option.
        constants (Mapping[str, object] | None): Flags that set this option to a fixed value, as
            flag name -> value.
        choices (Sequence[object] | Callable[[], Sequence[object]] | None): The permitted values,
            or a callable returning them.
        metavar (str | None): argparse metavar override.
        env (str | None): Environment variable name, or None to derive one.
        hidden (bool): True to keep the option out of ``--help``.

    Returns:
        dict[str, OptionSpec]: The metadata to pass as ``field(metadata=...)``.
    """
    spec = OptionSpec(
        help=help,
        flag=flag,
        short=short,
        synonyms=synonyms,
        constants={} if constants is None else dict(constants),
        choices=choices,
        metavar=metavar,
        env=env,
        hidden=hidden,
    )
    return {METADATA_KEY: spec}


def section_options(section_type: type[ConfigSection]) -> dict[str, tuple[attrs.Attribute, OptionSpec]]:
    """Return the declared options of a section, keyed by field name.

    Args:
        section_type (type[ConfigSection]): The section class to inspect.

    Returns:
        dict[str, tuple[attrs.Attribute, OptionSpec]]: Field name -> (field, spec), in declaration
            order, for the fields carrying ``option()`` metadata.
    """
    attrs.resolve_types(section_type)
    specs: dict[str, tuple[attrs.Attribute, OptionSpec]] = {}
    for section_field in fields(section_type):
        spec = section_field.metadata.get(METADATA_KEY)
        if spec is not None:
            specs[section_field.name] = (section_field, spec)
    return specs


def unwrap_optional_type(type_annotation: Any) -> tuple[Any, bool]:
    """Split ``X | None`` into ``X`` and a nullability flag.

    Args:
        type_annotation (Any): The resolved field type.

    Returns:
        tuple[Any, bool]: The non-None member and whether None was a member. A union of more than
            one non-None type is returned unchanged and reported as not optional.
    """
    origin = get_origin(type_annotation)
    if origin not in (Union, UnionType):
        return type_annotation, False
    members: list[Any] = []
    for member in get_args(type_annotation):
        if member is not NoneType:
            members.append(member)
    if len(members) != 1:
        return type_annotation, False
    return members[0], True
