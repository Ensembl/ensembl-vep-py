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
"""Logging options, read by every runner."""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

from attrs import define, field

from ensembl.vep.config.option_spec import option
from ensembl.vep.config.sections.base import ConfigSection

__all__ = ["LoggingConfig"]

LOG_LEVELS = ("debug", "info", "warning", "error", "critical")


@define(frozen=True, kw_only=True, auto_attribs=True)
class LoggingConfig(ConfigSection):
    """How much the run says, and where it says it.

    Attributes:
        level (str): Root level for the ``ensembl.vep`` logger hierarchy.
        log_file (Path | None): File to append log records to, or None for stream output only.
        verbose (bool): Shorthand for ``level="debug"``. Also what turns derivation tracing on,
            since ensembl.vep.trace keys off the strategy logger being at DEBUG.
        quiet (bool): Shorthand for ``level="error"``.
    """

    section_name: ClassVar[str] = "logging"

    level: str = field(default="info", metadata=option("logging level", choices=LOG_LEVELS))
    log_file: Path | None = field(default=None, metadata=option("write log records to this file"))
    verbose: bool = field(default=False, metadata=option("log at debug level", short="v"))
    quiet: bool = field(default=False, metadata=option("log errors only"))

    def effective_level(self) -> str:
        """Returns the level to apply once ``verbose``/``quiet`` are taken into account.

        Returns:
            str: One of ``LOG_LEVELS``. ``quiet`` wins over ``verbose`` when both are set, so the
                quieter of two conflicting requests is honoured.
        """
        if self.quiet:
            return "error"
        if self.verbose:
            return "debug"
        return self.level
