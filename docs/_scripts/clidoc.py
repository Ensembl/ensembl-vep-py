# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Generate the CLI reference pages from the registered runners and their config sections.

The CLI counterpart of ``sphinx-apidoc``: run it before ``sphinx-build`` (``make clidoc``) and it
writes one Markdown page per command group, a page for the options several groups share, and an
index, into the output directory (``docs/reference/cli`` by default).

Every option row comes from the same ``argparse`` parser ``vep`` builds at run time, so flags,
synonyms and metavars match ``vep <command> -h``. Defaults, environment variable names and
config-file keys come from the ``ConfigSection`` fields and their ``OptionSpec`` metadata.

Usage: ``python docs/_scripts/clidoc.py [OUTPUT_DIR]``
"""

from __future__ import annotations

import argparse
import os
import sys
from collections.abc import Iterable
from enum import Enum
from pathlib import Path

import attrs

from ensembl.vep.cli import COMMAND_SECTIONS, CONFIG_FLAG, DEFAULT_COMMAND, OTHER_SECTION
from ensembl.vep.config import Config
from ensembl.vep.config.option_spec import OptionSpec, section_options
from ensembl.vep.config.parser import build_parser
from ensembl.vep.config.parser.cli import destination
from ensembl.vep.config.parser.env import default_env_name
from ensembl.vep.config.sections.base import ConfigSection
from ensembl.vep.runner import get_runner, runner_default, runner_groups

CONFIG_HELP = "read options from a TOML or Perl VEP config file"
"""Help for ``--config``, which ``ensembl.vep.cli.main`` adds by hand rather than via a section."""

USAGE_WIDTH = "100"
"""Column width argparse wraps the usage lines to."""

ENTRYPOINTS_MARKERS = "<!-- entrypoints-table:start -->\n<!-- entrypoints-table:end -->"
"""Markers the ``ensembl.utils.docs`` entry-points extension fills in at build time."""


# --------------------------------------------------------------------------------------------------
# Small formatting helpers
# --------------------------------------------------------------------------------------------------


def _invocation(name: str) -> str:
    """Return a runner name as typed on the command line: ``run.default`` -> ``vep run default``."""
    return f"vep {name.replace('.', ' ', 1)}"


def _label(*parts: str) -> str:
    """Return a Sphinx cross-reference label built from the given parts."""
    return "cli-" + "-".join(part.replace(".", "-").replace("_", "-") for part in parts)


def _cell(text: str) -> str:
    """Escape text for a Markdown table cell."""
    return text.replace("|", "\\|").replace("\n", " ")


def _title(section_name: str) -> str:
    """Return the heading of a section's options: ``bench_endtoend`` -> ``Bench endtoend options``."""
    return f"{section_name.replace('_', ' ').capitalize()} options"


def _summary(section_type: type[ConfigSection]) -> str:
    """Return the first paragraph of a section's docstring, joined onto one line."""
    doc = (section_type.__doc__ or "").strip()
    paragraph = doc.split("\n\n", 1)[0]
    return " ".join(line.strip() for line in paragraph.splitlines())


def _value(value: object) -> str:
    """Render one default value the way it is typed on the command line."""
    if isinstance(value, Enum):
        return value.name.lower()
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (tuple, list, frozenset, set)):
        return ",".join(_value(item) for item in value)
    return str(value)


def _default(section_field: attrs.Attribute) -> str:
    """Return a field's default as a table cell, calling ``Factory`` defaults to get the value."""
    default = section_field.default
    if default is attrs.NOTHING:
        return "*required*"
    if isinstance(default, attrs.Factory):  # type: ignore[arg-type]
        if default.takes_self:
            return "*computed*"
        default = default.factory()
    if default is None or default == () or default == "":
        return "–"
    return f"`{_cell(_value(default))}`"


def _usage(prog: str, description: str, section_types: tuple[type[ConfigSection], ...]) -> str:
    """Return the usage line ``vep`` prints for a command, ``--config`` included."""
    parser = _command_parser(prog, description, section_types)
    return parser.format_usage().strip()


def _command_parser(
    prog: str, description: str, section_types: tuple[type[ConfigSection], ...]
) -> argparse.ArgumentParser:
    """Return the parser ``ensembl.vep.cli.main`` builds for one command."""
    parser = build_parser(prog, description, section_types)
    parser.add_argument(CONFIG_FLAG, metavar="FILE", default=None, help=CONFIG_HELP)
    return parser


# --------------------------------------------------------------------------------------------------
# Option tables
# --------------------------------------------------------------------------------------------------


def _section_rows(section_type: type[ConfigSection]) -> list[str]:
    """Return the Markdown table rows for every visible option a section declares."""
    parser = build_parser("vep", "", (section_type,))
    # The parser's own help formatter renders each flag list exactly as ``-h`` shows it.
    formatter = parser._get_formatter()  # pylint: disable=protected-access
    actions_by_dest: dict[str, list[argparse.Action]] = {}
    for action in parser._actions:  # pylint: disable=protected-access
        actions_by_dest.setdefault(action.dest, []).append(action)

    rows: list[str] = []
    for field_name, (section_field, spec) in section_options(section_type).items():
        if spec.hidden:
            continue
        actions = actions_by_dest.get(destination(section_type.section_name, field_name), [])
        main_actions = [
            action
            for action in actions
            if action.help != argparse.SUPPRESS and not isinstance(action, argparse._StoreConstAction)
        ]
        if not main_actions:
            continue
        flags = ", ".join(
            formatter._format_action_invocation(action)  # pylint: disable=protected-access
            for action in main_actions
        )
        env = spec.env or default_env_name(section_type.section_name, field_name)
        rows.append(
            f"| `{_cell(flags)}` | {_default(section_field)} | `{env}` | {_cell(_description(spec))} |"
        )
    return rows


def _description(spec: OptionSpec) -> str:
    """Return an option's help, capitalised, with its shorthand flags appended."""
    text = spec.help[:1].upper() + spec.help[1:]
    if not text.endswith("."):
        text += "."
    for constant_flag, constant in spec.constants.items():
        text += f" `--{constant_flag}` is shorthand for `{_value(constant)}`."
    return text


def _section_block(section_type: type[ConfigSection], heading: str, used_by: Iterable[str] = ()) -> list[str]:
    """Return a section's labelled heading, summary, TOML hint and option table."""
    name = section_type.section_name
    lines = [f"({_label('section', name)})=", f"{heading} {_title(name)}", ""]
    summary = _summary(section_type)
    if summary:
        lines += [summary, ""]
    used = sorted(used_by)
    if used:
        lines += ["Used by " + ", ".join(f"{{ref}}`{_invocation(n)} <{_label(n)}>`" for n in used) + ".", ""]
    lines += [f"In a TOML config file these options go under the `[{name}]` table.", ""]
    rows = _section_rows(section_type)
    if not rows:
        return [*lines, "_This section declares no command-line options._", ""]
    return [
        *lines,
        "| Option | Default | Environment variable | Description |",
        "| ------ | ------- | -------------------- | ----------- |",
        *rows,
        "",
    ]


def _section_links(section_names: Iterable[str]) -> str:
    """Return inline cross-references to the given sections' option tables."""
    return ", ".join(f"{{ref}}`{_title(name)} <{_label('section', name)}>`" for name in section_names)


# --------------------------------------------------------------------------------------------------
# Pages
# --------------------------------------------------------------------------------------------------


def _ordered_groups() -> list[str]:
    """Return the command groups in the order ``vep -h`` lists them."""
    groups = runner_groups()
    ordered: list[str] = []
    for _heading, section_groups in COMMAND_SECTIONS:
        ordered.extend(group for group in section_groups if group in groups)
    ordered.extend(group for group in groups if group not in ordered)
    return ordered


def _heading_of(group: str) -> str:
    """Return the ``vep -h`` heading a group is listed under."""
    for heading, section_groups in COMMAND_SECTIONS:
        if group in section_groups:
            return heading
    return OTHER_SECTION


def _common_sections() -> dict[str, set[str]]:
    """Return the sections read by commands from more than one group, with the commands using them."""
    used_by: dict[str, set[str]] = {}
    for names in runner_groups().values():
        for name in names:
            for section_name in get_runner(name).sections:
                used_by.setdefault(section_name, set()).add(name)
    groups_of = {
        section: {name.split(".", 1)[0] for name in commands} for section, commands in used_by.items()
    }
    return {section: used_by[section] for section, groups in groups_of.items() if len(groups) > 1}


def _group_page(group: str, common: dict[str, set[str]]) -> str:
    """Return the Markdown page documenting every command of one group."""
    names = runner_groups()[group]
    lines = [f"({_label(group)})=", f"# vep {group}", ""]
    lines += [f"Commands to {_heading_of(group)}.", ""]
    default = runner_default(group)
    if default is not None and len(names) > 1:
        lines += [f"Without a verb, `vep {group}` runs `{_invocation(default)}`.", ""]
    if group == DEFAULT_COMMAND.split(".", 1)[0]:
        lines += [f"Without any command, `vep` runs `{_invocation(DEFAULT_COMMAND)}`.", ""]

    # Sections shared by several commands of this group (but no other group) are documented once.
    section_users: dict[str, list[str]] = {}
    for name in names:
        for section_name in get_runner(name).sections:
            if section_name not in common:
                section_users.setdefault(section_name, []).append(name)
    shared = [section for section, users in section_users.items() if len(users) > 1]

    for name in names:
        runner = get_runner(name)
        section_types = Config.section_types_for(runner.sections)
        heading = "#" if len(names) == 1 else "##"
        if len(names) > 1:
            lines += [f"({_label(name)})=", f"## {_invocation(name)}", ""]
        elif _label(name) != _label(group):
            # A single-command group: the page title is the command, so label it as one too.
            lines.insert(0, f"({_label(name)})=")
        lines += [runner.description[:1].upper() + runner.description[1:] + ".", ""]
        lines += ["```text", _usage(_invocation(name), runner.description, section_types), "```", ""]
        common_used = [section for section in runner.sections if section in common]
        shared_used = [section for section in runner.sections if section in shared]
        if common_used:
            lines += [f"Common options: {_section_links(common_used)}.", ""]
        if shared_used:
            lines += [f"Shared options: {_section_links(shared_used)}.", ""]
        lines += [f"Global options: {{ref}}`--config FILE <{_label('global')}>`.", ""]
        own = [
            section_type
            for section_type in section_types
            if section_type.section_name not in common and section_type.section_name not in shared
        ]
        for section_type in own:
            lines += _section_block(section_type, heading + "#")

    if shared:
        lines += [f"## Options shared by `vep {group}` commands", ""]
        for section_name in shared:
            section_type = Config.section_types[section_name]
            lines += _section_block(section_type, "###", section_users[section_name])
    return "\n".join(lines).rstrip() + "\n"


def _common_page(common: dict[str, set[str]]) -> str:
    """Return the Markdown page documenting the sections several command groups read."""
    lines = [
        f"({_label('common')})=",
        "# Common options",
        "",
        "Option sections read by commands from more than one group. Each command's page links to",
        "the ones it accepts.",
        "",
    ]
    ordered = sorted(common, key=lambda section: (-len(common[section]), section))
    for section_name in ordered:
        lines += _section_block(Config.section_types[section_name], "##", common[section_name])
    return "\n".join(lines).rstrip() + "\n"


def _index_page(groups: list[str]) -> str:
    """Return the CLI reference landing page: entry points, option sources, command listing."""
    lines = [
        f"({_label('reference')})=",
        "# CLI reference",
        "",
        "Every command, verb and option of the `vep` command line. These pages are generated by",
        "`make clidoc` from the option declarations, so they always match `vep <command> -h`.",
        "",
        "## Entry points",
        "",
        ENTRYPOINTS_MARKERS,
        "",
        "## Setting options",
        "",
        "Each option can be set in three places. When the same option is set more than once, the",
        "later source in this list wins:",
        "",
        "1. a config file passed with `--config FILE`, either TOML with one `[section]` table per",
        "   option section, or a Perl VEP config file with one `option value` per line;",
        "2. its environment variable, listed in each option table;",
        "3. its command-line flag.",
        "",
        "Options left unset take the default shown in the tables.",
        "",
        f"({_label('global')})=",
        "## Global options",
        "",
        "| Option | Default | Description |",
        "| ------ | ------- | ----------- |",
        f"| `{CONFIG_FLAG} FILE` | – | {_cell(CONFIG_HELP[:1].upper() + CONFIG_HELP[1:])}. May appear anywhere"
        " on the command line. |",
        "| `-h, --help` | – | Show the help for `vep`, a group (`vep <group> -h`) or a command. |",
        "",
        "## Commands",
        "",
        "| Command | Description |",
        "| ------- | ----------- |",
    ]
    for group in groups:
        for name in runner_groups()[group]:
            description = get_runner(name).description
            lines.append(f"| {{ref}}`{_invocation(name)} <{_label(name)}>` | {_cell(description)} |")
    lines += [
        "",
        "```{toctree}",
        ":hidden:",
        "",
        *groups,
        "common",
        "```",
    ]
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> None:
    """Write the CLI reference pages."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument(
        "output_dir", nargs="?", type=Path, default=Path("docs/reference/cli"), help="where to write"
    )
    args = parser.parse_args(argv)
    os.environ["COLUMNS"] = USAGE_WIDTH  # argparse wraps usage lines to the terminal width

    output_dir: Path = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    groups = _ordered_groups()
    common = _common_sections()

    pages = {"index.md": _index_page(groups), "common.md": _common_page(common)}
    for group in groups:
        pages[f"{group}.md"] = _group_page(group, common)
    for filename, content in pages.items():
        (output_dir / filename).write_text(content, encoding="utf-8")
    print(f"clidoc: wrote {len(pages)} pages to {output_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
