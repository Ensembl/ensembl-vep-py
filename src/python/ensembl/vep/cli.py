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
"""The ``vep`` command line.

The first argument selects the command; when it is not a registered one the default is used, so
``vep --input_file x.vcf`` still parses. ``--config FILE`` may stand anywhere and is taken off
first, so ``vep --config f.toml corpus describe`` reaches ``corpus describe``. ``vep -h`` lists
every command and ``vep <group> -h`` the verbs of one group, laid out the way ``git -h`` lays out
its commands.
"""

from __future__ import annotations

import sys
from pathlib import Path

from ensembl.vep.common import setup_logging
from ensembl.vep.config import Config, load
from ensembl.vep.config.parser import build_parser
from ensembl.vep.runner import ExitCode, available_runners, get_runner, runner_default, runner_groups

# Importing ensembl.vep.runner above is what registers the built-in runners.

__all__ = ["main", "resolve_command"]

DEFAULT_COMMAND = "run.default"
CONFIG_FLAG = "--config"
HELP_FLAGS = frozenset({"-h", "--help"})

COMMAND_SECTIONS: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("annotate variants", ("run",)),
    ("explain and summarise annotation", ("explain", "stats")),
    ("build and check the synthetic corpus", ("corpus",)),
    ("measure annotation speed", ("bench",)),
    ("maintain the consequence catalogue", ("catalogue",)),
)
"""Headings for the help listing, each over the command groups it covers, in display order."""

OTHER_SECTION = "other commands"
"""The heading for a registered group no entry in ``COMMAND_SECTIONS`` covers."""


def _is_verb_group(token: str) -> bool:
    """Return whether a word names a group of verbs rather than a command of its own."""
    return token in runner_groups() and token not in available_runners()


def _section_of(group: str) -> str:
    """Return the help heading a command group is listed under."""
    for heading, groups in COMMAND_SECTIONS:
        if group in groups:
            return heading
    return OTHER_SECTION


def _command_listing(entries: list[tuple[str, str, str]]) -> list[str]:
    """Return ``(heading, name, description)`` entries as headed, aligned blocks.

    Headings keep the order they first appear in; names are indented three spaces and padded to
    the longest, as ``git -h`` lays out its commands.
    """
    width = max((len(name) for _heading, name, _description in entries), default=0) + 3
    by_heading: dict[str, list[str]] = {}
    for heading, name, description in entries:
        by_heading.setdefault(heading, []).append(f"   {name:<{width}}{description}")

    lines: list[str] = []
    for heading, commands in by_heading.items():
        lines.extend(["", heading, *commands])
    return lines


def _usage() -> str:
    """Return the top-level help, listing every command under the heading of what it is for."""
    groups = runner_groups()
    ordered_groups: list[str] = []
    for _heading, section_groups in COMMAND_SECTIONS:
        ordered_groups.extend(group for group in section_groups if group in groups)
    ordered_groups.extend(group for group in groups if group not in ordered_groups)

    entries = []
    for group in ordered_groups:
        for name in groups[group]:
            entries.append((_section_of(group), name.replace(".", " ", 1), get_runner(name).description))

    default_invocation = DEFAULT_COMMAND.replace(".", " ", 1)
    lines = [
        "usage: vep [-h | --help] [--config <file>] <command> [<verb>] [<options>]",
        "",
        "These are the vep commands, grouped by what they are for:",
        *_command_listing(entries),
        "",
        f"Without a command vep runs '{default_invocation}', so 'vep --input_file x.vcf' annotates x.vcf.",
        "See 'vep <command> -h' for a command's options, or 'vep <group> -h' for a group's verbs.",
    ]
    return "\n".join(lines)


def _group_usage(group: str) -> str:
    """Return the help for one group of verbs, laid out like the top-level help."""
    entries = []
    for name in runner_groups()[group]:
        entries.append((_section_of(group), name.split(".", 1)[1], get_runner(name).description))

    lines = [
        f"usage: vep {group} <verb> [<options>]",
        "",
        f"These are the 'vep {group}' verbs:",
        *_command_listing(entries),
        "",
    ]
    default = runner_default(group)
    if default is not None:
        lines.append(f"Without a verb 'vep {group}' runs 'vep {default.replace('.', ' ', 1)}'.")
    lines.append(f"See 'vep {group} <verb> -h' for a verb's options.")
    return "\n".join(lines)


def _help_requested(argv: list[str]) -> str | None:
    """Return the help page a leading ``-h`` asks for, or None when it belongs to a command.

    ``vep -h`` and ``vep <group> -h`` are the only two. ``vep --input_file x.vcf -h`` and
    ``vep run default -h`` still reach the command's own options.
    """
    if argv and argv[0] in HELP_FLAGS:
        return _usage()
    if len(argv) >= 2 and argv[1] in HELP_FLAGS and _is_verb_group(argv[0]):
        return _group_usage(argv[0])
    return None


def _unknown_command(token: str) -> str:
    """Return the error and help page for a leading word that is not a command."""
    if token in runner_groups():
        return f"vep {token}: needs a verb.\n\n{_group_usage(token)}"
    return f"vep: '{token}' is not a vep command. See 'vep --help'.\n\n{_usage()}"


def resolve_command(argv: list[str]) -> tuple[str, list[str]]:
    """Split arguments into a command name and the arguments that belong to it.

    A bare group name resolves to that group's single runner if it has one, otherwise stays
    unresolved so the caller can list its verbs.

    Only a leading flag, or no arguments at all, means the default command. ``vep`` takes no
    positional arguments, so a leading bare word can only be a command, and one that is not
    registered is an error rather than something to pass to the default.

    Args:
        argv (list[str]): Arguments without the program name.

    Returns:
        tuple[str, list[str]]: The command name and the remaining arguments. The name is empty when
            the leading word is not a command, or names a group without one of its verbs.
    """
    registered = set(available_runners())
    if not argv or argv[0].startswith("-"):
        return DEFAULT_COMMAND, argv
    if argv[0] in registered:
        return argv[0], argv[1:]

    groups = runner_groups()
    if argv[0] in groups:
        names = groups[argv[0]]
        if len(argv) > 1 and f"{argv[0]}.{argv[1]}" in registered:
            return f"{argv[0]}.{argv[1]}", argv[2:]
        default = runner_default(argv[0])
        if default is not None:
            return default, argv[1:]
        if len(names) == 1:
            return names[0], argv[1:]
    return "", argv


def _take_config_flag(argv: list[str]) -> tuple[list[str], str | None]:
    """Remove ``--config FILE`` or ``--config=FILE`` from the arguments, if present.

    Handled here rather than as a declared option: it is read before the command is resolved, so
    it may stand before or after it, and belongs to no one section. The last occurrence wins.

    Args:
        argv (list[str]): Every argument, without the program name.

    Returns:
        tuple[list[str], str | None]: The arguments without the flag, and the file it named.

    Raises:
        ValueError: If the flag has no file after it.
    """
    remaining_args: list[str] = []
    config_file: str | None = None
    index = 0
    while index < len(argv):
        argument = argv[index]
        if argument == CONFIG_FLAG:
            if index + 1 >= len(argv):
                raise ValueError(f"argument {CONFIG_FLAG}: expected one argument")
            config_file = argv[index + 1]
            index += 2
            continue
        if argument.startswith(f"{CONFIG_FLAG}="):
            config_file = argument.split("=", 1)[1]
            if not config_file:
                raise ValueError(f"argument {CONFIG_FLAG}: expected one argument")
            index += 1
            continue
        remaining_args.append(argument)
        index += 1
    return remaining_args, config_file


def main(argv: list[str] | None = None) -> ExitCode:
    """Run a command and return its exit code.

    Args:
        argv (list[str] | None): Arguments without the program name, or None for ``sys.argv[1:]``.

    Returns:
        ExitCode: The runner's exit code, ``OK`` after printing a help page, or ``USAGE`` if no
            command could be resolved, ``--config`` is malformed, or the configuration is invalid
            or unreadable.
    """
    try:
        arguments, config_file = _take_config_flag(list(sys.argv[1:] if argv is None else argv))
    except ValueError as invalid:
        print(f"vep: error: {invalid}", file=sys.stderr)
        return ExitCode.USAGE

    help_page = _help_requested(arguments)
    if help_page is not None:
        print(help_page)
        return ExitCode.OK

    command, remaining_args = resolve_command(arguments)
    if not command:
        print(_unknown_command(arguments[0]), file=sys.stderr)
        return ExitCode.USAGE

    runner_type = get_runner(command)
    section_types = Config.section_types_for(runner_type.sections)
    parser = build_parser(f"vep {command.replace('.', ' ', 1)}", runner_type.description, section_types)
    parser.add_argument(
        CONFIG_FLAG, metavar="FILE", default=None, help="read options from a TOML or Perl VEP config file"
    )

    try:
        config = load(
            remaining_args,
            command=command,
            section_names=runner_type.sections,
            parser=parser,
            config_file=Path(config_file) if config_file else None,
        )
    except OSError as unreadable:
        if config_file is None:
            # Nothing else in load() reads a file, so an OSError here is not ours to explain away.
            raise
        parser.print_usage(sys.stderr)
        reason = unreadable.strerror or unreadable
        print(f"{parser.prog}: error: cannot read config file '{config_file}': {reason}", file=sys.stderr)
        return ExitCode.USAGE
    except ValueError as invalid:
        parser.print_usage(sys.stderr)
        print(f"{parser.prog}: error: {invalid}", file=sys.stderr)
        return ExitCode.USAGE
    setup_logging(config.logging.effective_level().upper(), config.logging.log_file)
    return runner_type(config).run()


if __name__ == "__main__":
    raise SystemExit(main())
