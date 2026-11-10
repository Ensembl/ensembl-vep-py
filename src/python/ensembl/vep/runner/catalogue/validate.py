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
"""Consequence rule-set validation."""

from __future__ import annotations

import logging
from typing import ClassVar

from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.annotators.consequence.facts.base import verify_requires
from ensembl.vep.annotators.consequence.rules.rulesets import FactRuleSet, build_fact_ruleset
from ensembl.vep.annotators.consequence.rules.validation import ValidationReport
from ensembl.vep.config import Config
from ensembl.vep.corpus.registry import CorpusOptionsError, build_corpus_from_options
from ensembl.vep.runner.base import ExitCode
from ensembl.vep.runner.catalogue.common import log_validation_report
from ensembl.vep.runner.registry import register_runner

__all__ = ["Runner"]

logger = logging.getLogger(__name__)


@register_runner("catalogue.validate")
class Runner:
    """Runner for configured catalogue validation tiers.

    A non-zero exit is a real failure: the checked-in catalogue passes today.
    """

    name: ClassVar[str] = "catalogue.validate"
    sections: ClassVar[tuple[str, ...]] = ("logging", "catalogue", "catalogue_validate", "corpus")
    description: ClassVar[str] = "validate the consequence rule set"

    def __init__(self, config: Config) -> None:
        """Hold the configuration."""
        self.config = config

    def run(self) -> ExitCode:
        """Report every finding, and fail if any is build-breaking."""
        ruleset = build_fact_ruleset()
        report = ValidationReport.for_ruleset(ruleset)
        logger.info("checked %d rules", len(ruleset.rules))
        log_validation_report(report)
        exit_code = ExitCode.FAILURE if report.has_errors else ExitCode.OK

        if self.config.catalogue_validate.verify_requires:
            try:
                if self._report_requires_problems():
                    exit_code = ExitCode.FAILURE
            except CorpusOptionsError as unusable:
                logger.error("%s", unusable)
                return ExitCode.USAGE

        if self.config.catalogue_validate.solver == "z3" and self._report_solver_findings(ruleset):
            exit_code = ExitCode.FAILURE

        return exit_code

    def _report_requires_problems(self) -> bool:
        """Check every declared ``requires`` against real contexts, and return whether any failed."""
        with build_corpus_from_options(self.config.corpus) as corpus:
            contexts = FactContext.for_overlaps([overlap for _name, overlap in corpus.named_overlaps])
            problems = verify_requires(contexts)
            description = corpus.description
        logger.info("verify_requires over %s: %d problem(s)", description, len(problems))
        for problem in problems:
            logger.error("ERROR: %s", problem)
        return bool(problems)

    def _report_solver_findings(self, ruleset: FactRuleSet) -> bool:
        """Run the whole-ruleset SAT queries, and return whether any rule was unreachable."""
        # pylint: disable=import-outside-toplevel
        from ensembl.vep.annotators.consequence.rules.validation import SolverChecks

        try:
            findings = SolverChecks.findings(ruleset)
        except RuntimeError as unavailable:
            logger.warning("solver skipped: %s", unavailable)
            return False

        unreachable = False
        for finding in findings:
            if finding.startswith("UNREACHABLE"):
                logger.error("ERROR: %s", finding)
                unreachable = True
            else:
                logger.info("%s", finding)
        return unreachable
