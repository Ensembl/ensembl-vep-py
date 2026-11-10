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
"""Consequence catalogue term addition."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import ClassVar

from ensembl.vep.annotators.consequence.rules.editing import apply_add
from ensembl.vep.annotators.consequence.rules.validation import ValidationReport
from ensembl.vep.config import Config
from ensembl.vep.runner.base import ExitCode
from ensembl.vep.runner.catalogue.common import (
    destination_for,
    log_diff,
    log_validation_report,
    read_payload,
    write_payload,
)
from ensembl.vep.runner.registry import register_runner

__all__ = ["Runner"]

logger = logging.getLogger(__name__)


@register_runner("catalogue.add")
class Runner:
    """Runner for adding and validating a consequence term."""

    name: ClassVar[str] = "catalogue.add"
    sections: ClassVar[tuple[str, ...]] = ("logging", "catalogue", "catalogue_write", "catalogue_add")
    description: ClassVar[str] = "add a term to the consequence catalogue"

    def __init__(self, config: Config) -> None:
        """Hold the configuration."""
        self.config = config

    def run(self) -> ExitCode:
        """Add the term, or report why it was refused."""
        entry = self.config.catalogue_add
        if entry.term is None:
            logger.error("nothing to add: give --term")
            return ExitCode.USAGE

        baseline, payload = read_payload(self.config)
        try:
            destination = destination_for(self.config) if not self.config.catalogue_write.dry_run else None
            result = apply_add(
                payload,
                term=entry.term,
                accession=entry.accession,
                impact=entry.impact,
                rank=entry.rank,
                tier=entry.tier,
                so_parent=entry.so_parent,
                subsumes=entry.subsumes,
                no_derive_subsumes=entry.no_derive_subsumes,
                when=entry.when,
                allow_missing_when=entry.allow_missing_when,
                notes=entry.notes,
            )
        except ValueError as error:
            logger.error("%s", error)
            return ExitCode.USAGE

        if result is None:
            return ExitCode.FAILURE
        updated, messages = result
        for message in messages:
            logger.info("%s", message)

        report = ValidationReport.for_catalogue_entries(updated["consequences"], baseline["consequences"])
        log_validation_report(report)
        if report.has_errors:
            logger.error("refusing to write: validation reported errors")
            return ExitCode.FAILURE

        if destination is None:
            log_diff(baseline, updated, self.config.catalogue.catalogue_file or Path("catalogue.json"))
            logger.info("dry run: nothing written")
            return ExitCode.OK

        log_diff(baseline, updated, destination)
        write_payload(updated, destination)
        return ExitCode.OK
