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
"""Shared helpers for the catalogue runners."""

from __future__ import annotations

import copy
import difflib
import logging
from pathlib import Path
from typing import Any

from ensembl.vep.annotators.consequence.rules.catalogue import catalogue_path, dump_catalogue, load_catalogue
from ensembl.vep.annotators.consequence.rules.validation import ValidationReport
from ensembl.vep.config import Config

__all__ = ["destination_for", "log_diff", "log_validation_report", "read_payload", "write_payload"]

logger = logging.getLogger(__name__)


def read_payload(config: Config) -> tuple[dict[str, Any], dict[str, Any]]:
    """Return the catalogue payload and an untouched copy of it.

    The copy is the baseline a change is validated against, so a pre-existing finding is not
    reported as one this change introduced.
    """
    payload = load_catalogue(config.catalogue.catalogue_file)
    return copy.deepcopy(payload), payload


def destination_for(config: Config) -> Path:
    """Return where a change should be written.

    Raises:
        ValueError: If neither a path nor permission to overwrite the packaged catalogue was given.
            An installed ``vep`` must not silently rewrite its own wheel.
    """
    if config.catalogue.catalogue_file is not None:
        return config.catalogue.catalogue_file
    if config.catalogue_write.catalogue_overwrite:
        return catalogue_path()
    raise ValueError(
        "refusing to write: give --catalogue PATH, or --catalogue_overwrite to replace the "
        "packaged catalogue in place"
    )


def log_validation_report(report: ValidationReport) -> None:
    """Log a validation report, errors first."""
    for error in report.errors:
        logger.error("ERROR: %s", error)
    for finding in report.unreachable_rules:
        logger.error("ERROR: unreachable rule introduced by this change: %s", finding.term)
    for warning in report.warnings:
        logger.warning("warning: %s", warning)
    if not report.has_errors:
        logger.info("validation passed")


def log_diff(original: dict[str, Any], updated: dict[str, Any], destination: Path) -> list[str]:
    """Log a unified diff of a catalogue change, and return its lines."""
    lines = list(
        difflib.unified_diff(
            dump_catalogue(original).splitlines(keepends=True),
            dump_catalogue(updated).splitlines(keepends=True),
            fromfile=f"a/{destination.name}",
            tofile=f"b/{destination.name}",
        )
    )
    for line in lines:
        logger.info("%s", line.rstrip())
    return lines


def write_payload(payload: dict[str, Any], destination: Path) -> None:
    """Write a catalogue payload, in the checked-in formatting."""
    destination.write_text(dump_catalogue(payload), encoding="utf-8")
    logger.info("wrote %s", destination)
