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
"""The runner contract."""

from __future__ import annotations

from enum import IntEnum
from typing import ClassVar, Protocol, runtime_checkable

from ensembl.vep.config import Config

__all__ = ["ExitCode", "Runner"]


class ExitCode(IntEnum):
    """Process exit codes returned by runners."""

    OK = 0
    """The command completed successfully."""
    FAILURE = 1
    """The command ran and reported a failure."""
    USAGE = 2
    """The command was invoked incorrectly."""


@runtime_checkable
class Runner(Protocol):
    """One command's execution.

    Attributes:
        name (ClassVar[str]): The registry name, following ``<domain>.<verb>``.
        sections (ClassVar[tuple[str, ...]]): The configuration sections this runner reads, named
            by their attribute on ``Config``. The command line is derived from exactly these,
            so a runner sees only the options it declares an interest in.
        description (ClassVar[str]): One-line description, used as the subcommand's help text.
    """

    name: ClassVar[str]
    sections: ClassVar[tuple[str, ...]]
    description: ClassVar[str]

    def __init__(self, config: Config) -> None:
        """Initialize the runner.

        Args:
            config (Config): The populated configuration.
        """
        raise NotImplementedError

    def run(self) -> ExitCode:
        """Execute the command.

        Returns:
            ExitCode: The process exit code.
        """
        raise NotImplementedError
