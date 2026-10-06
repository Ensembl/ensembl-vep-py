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
"""Helpers shared by the command-line, environment and file readers."""

from __future__ import annotations

from ensembl.vep.config.option_spec import OptionSpec, section_options
from ensembl.vep.config.sections.base import ConfigSection

__all__ = ["flag_for", "option_spellings", "resolve_option", "resolve_option_key"]


def flag_for(field_name: str, spec: OptionSpec) -> str:
    """Return the long flag spelling for one option, without the leading dashes.

    Flags keep underscores, so ``--input_file`` matches legacy VEP and a legacy config file needs
    no key translation.
    """
    return spec.flag or field_name


def option_spellings(field_name: str, spec: OptionSpec) -> tuple[str, ...]:
    """Return every long spelling that sets one option: its field name, flag and synonyms."""

    spellings = [field_name, flag_for(field_name, spec)]
    spellings.extend(spec.synonyms)
    unique: list[str] = []
    for spelling in spellings:
        if spelling not in unique:
            unique.append(spelling)
    return tuple(unique)


def _normalised(key: str) -> str:
    """Return a flat key as a field name: no dashes, underscores throughout."""

    return key.strip().lstrip("-").replace("-", "_")


def resolve_option_key(key: str, section_types: tuple[type[ConfigSection], ...]) -> tuple[str, str] | None:
    """Resolve a flat option key to the section and field that declare it.

    Args:
        key (str): The key as written, with or without leading dashes.
        section_types (tuple[type[ConfigSection], ...]): The sections in play.

    Returns:
        tuple[str, str] | None: (section name, field name), or None if no section declares it.
    """
    resolved = resolve_option(key, None, section_types)
    if resolved is None:
        return None
    section_name, field_name, _value = resolved
    return section_name, field_name


def resolve_option(
    key: str, value: object, section_types: tuple[type[ConfigSection], ...]
) -> tuple[str, str, object] | None:
    """Resolve a flat key and its value to the section, field and value they set.

    A key naming one of an option's constants (``--vcf`` for ``--output_format vcf``) sets that
    constant, whatever value came with it: the flag carries the value, so a config file writing
    ``vcf 1`` means the same as ``vcf``.

    Args:
        key (str): The key as written, with or without leading dashes.
        value (object): The value given with it, ignored for a constant flag.
        section_types (tuple[type[ConfigSection], ...]): The sections in play.

    Returns:
        tuple[str, str, object] | None: (section name, field name, value), or None if no section
            declares the key.
    """
    normalised = _normalised(key)
    for section_type in section_types:
        for field_name, (_field, spec) in section_options(section_type).items():
            if normalised in option_spellings(field_name, spec):
                return section_type.section_name, field_name, value
            for constant_flag, constant in spec.constants.items():
                if normalised == _normalised(constant_flag):
                    return section_type.section_name, field_name, constant
    return None
