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
"""The haplotype-aware annotation runner."""

from __future__ import annotations

import logging
from typing import ClassVar

from ensembl.vep.config import Config
from ensembl.vep.runner.base import ExitCode
from ensembl.vep.runner.registry import register_runner

__all__ = ["Runner"]

logger = logging.getLogger(__name__)


@register_runner("run.haplo")
class Runner:
    """Placeholder runner for the haplotype workflow command."""

    name: ClassVar[str] = "run.haplo"
    sections: ClassVar[tuple[str, ...]] = (
        "logging",
        "source",
        "annotation",
        "execution",
        "run",
    )
    description: ClassVar[str] = "annotate phased haplotypes"

    def __init__(self, config: Config) -> None:
        """Store the resolved configuration."""
        self.config = config

    def run(self) -> ExitCode:
        """Report that the haplotype input pipeline is not implemented yet."""
        logger.error("haplotype annotation is not available yet")
        return ExitCode.USAGE
