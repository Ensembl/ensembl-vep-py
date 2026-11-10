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
"""Golden corpus term generation and validation."""

from __future__ import annotations

import logging
from typing import Any, ClassVar

from attrs import asdict

from ensembl.vep.annotators.factory import build_annotator
from ensembl.vep.config import Config, ConfigSource
from ensembl.vep.corpus.golden import GoldenReport, GoldenTerms, record_golden
from ensembl.vep.corpus.registry import build_corpus
from ensembl.vep.corpus.spec import CorpusSpec
from ensembl.vep.runner.base import ExitCode
from ensembl.vep.runner.registry import register_runner

__all__ = ["Runner"]

logger = logging.getLogger(__name__)


@register_runner("corpus.golden")
class Runner:
    """Runner for generating or validating golden corpus terms.

    Checking reads the corpus spec from the golden file itself, so ``--golden_file`` alone is a
    complete check. The annotator comes from the file too, unless ``--annotator`` is given
    explicitly -- checking another strategy against one annotator's golden terms is a conformance
    check against a fixed expectation, which is the point.

    Attributes:
        name (ClassVar[str]): The registry name.
        sections (ClassVar[tuple[str, ...]]): The configuration sections read.
        description (ClassVar[str]): One-line description for the subcommand list.
    """

    name: ClassVar[str] = "corpus.golden"
    sections: ClassVar[tuple[str, ...]] = ("logging", "corpus", "annotation", "corpus_golden")
    description: ClassVar[str] = "record or check the terms an annotator emits for a corpus"

    def __init__(self, config: Config) -> None:
        """Hold the configuration."""
        self.config = config

    def run(self) -> ExitCode:
        """Record or check, as asked."""
        golden = self.config.corpus_golden
        if golden.golden_file is None:
            logger.error("give --golden_file PATH to say which golden terms to check or record")
            return ExitCode.USAGE

        if golden.update_golden:
            return self._update()
        return self._check()

    def _check(self) -> ExitCode:
        """Check current annotations against an existing golden terms file."""
        golden = self.config.corpus_golden
        assert golden.golden_file is not None
        if not golden.golden_file.exists():
            logger.error("%s does not exist; record it first with --update_golden", golden.golden_file)
            return ExitCode.USAGE

        try:
            expected = GoldenTerms.load(golden.golden_file)
        except (ValueError, KeyError) as unreadable:
            logger.error("cannot read %s: %s", golden.golden_file, unreadable)
            return ExitCode.FAILURE

        annotator_name = expected.annotator
        if self.config.source_of("annotation", "annotator") is not ConfigSource.DEFAULT:
            annotator_name = self.config.annotation.annotator
        try:
            actual = self._record(expected.spec, annotator_name)
        except (KeyError, TypeError, ValueError, OSError) as unusable:
            logger.error("%s", unusable)
            return ExitCode.USAGE
        report = GoldenReport.compare(expected, actual)
        emit = logger.info if report.ok else logger.error
        for line in report.lines():
            emit("%s", line)
        return ExitCode.OK if report.ok else ExitCode.FAILURE

    def _update(self) -> ExitCode:
        """Record current annotations as the requested golden terms file."""
        golden_file = self.config.corpus_golden.golden_file
        assert golden_file is not None
        try:
            spec = CorpusSpec.from_options(self.config.corpus)
            recorded = self._record(spec, self.config.annotation.annotator)
        except (KeyError, TypeError, ValueError, OSError) as unusable:
            logger.error("%s", unusable)
            return ExitCode.USAGE
        recorded.write(golden_file)
        logger.info(
            "recorded %d overlaps from %s over %s to %s",
            len(recorded.terms),
            recorded.annotator,
            spec.kind,
            golden_file,
        )
        return ExitCode.OK

    def _record(self, spec: CorpusSpec, annotator_name: str) -> GoldenTerms:
        """Annotate a corpus with the named annotator, configured caches applied."""
        options: dict[str, Any] = asdict(self.config.annotation)
        options["annotator"] = annotator_name
        annotator = build_annotator(
            **options,
            upstream_distance=self.config.source.upstream_distance,
            downstream_distance=self.config.source.downstream_distance,
        )
        with build_corpus(spec) as corpus:
            return record_golden(corpus, annotator, annotator_name)
