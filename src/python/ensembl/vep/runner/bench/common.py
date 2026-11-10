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
"""Shared benchmark corpora, gates, and result handling."""

from __future__ import annotations

import logging
from collections.abc import Iterator, Sequence
from contextlib import ExitStack, contextmanager

from ensembl.vep.annotators.factory import build_annotator
from ensembl.vep.annotators.registry import available_annotators
from ensembl.vep.bench.baseline import (
    BenchResults,
    RunInfo,
    compare_results,
    git_state,
    machine_info,
    utc_timestamp,
)
from ensembl.vep.bench.harness import (
    REFERENCE_STRATEGY,
    TimingPolicy,
    conformance_gate,
    consequence_strategies,
)
from ensembl.vep.config import Config
from ensembl.vep.config.report import resolve_report_format
from ensembl.vep.corpus.registry import (
    CorpusOptionsError,
    build_corpus,
    corpus_spec_from_options,
)
from ensembl.vep.corpus.spec import CorpusSpec
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.runner.base import ExitCode

__all__ = [
    "corpora",
    "gate",
    "policy",
    "record_results",
    "resolve_strategies",
    "run_info",
    "table_format",
]

logger = logging.getLogger(__name__)


def resolve_strategies(requested: Sequence[str], default: Sequence[str]) -> list[str]:
    """Return the requested strategies, or the registered defaults when none were asked for.

    Raises:
        ValueError: If a requested strategy is not registered.
    """
    registered = set(available_annotators())
    if not requested:
        return [name for name in default if name in registered]
    unknown = [name for name in requested if name not in registered]
    if unknown:
        raise ValueError(f"not registered: {', '.join(unknown)}")
    return list(requested)


def table_format(config: Config) -> str:
    """Return the table format the bench report asks for.

    Raises:
        ValueError: If a diagram format was asked for; a timing table has no diagram form.
    """

    return resolve_report_format(
        config.bench.output_format, graph=False, table_default="text", graph_formats=()
    )


@contextmanager
def corpora(config: Config) -> Iterator[list[tuple[str, list[TranscriptVariationAllele]]]]:
    """Keep each configured corpus alive while yielding its labelled overlaps."""
    spec = corpus_spec_from_options(config.corpus)
    with ExitStack() as resources:
        try:
            labelled = []
            for kind in spec.kinds:
                selected = CorpusSpec(kind=kind, parameters=spec.parameters)
                corpus = resources.enter_context(build_corpus(selected))
                labelled.append((kind, [overlap for _name, overlap in corpus.named_overlaps]))
        except (KeyError, TypeError, ValueError, OSError) as error:
            if isinstance(error, CorpusOptionsError):
                raise
            message = error.args[0] if isinstance(error, KeyError) and error.args else str(error)
            raise CorpusOptionsError(str(message)) from error
        yield labelled


def policy(config: Config) -> TimingPolicy:
    """Return the configured timing policy."""
    return TimingPolicy(repeats=config.bench.repeats, warmup=config.bench.warmup)


def gate(labelled: Sequence[tuple[str, Sequence[TranscriptVariationAllele]]]) -> bool:
    """Check every registered consequence strategy against the reference over every corpus, and log it."""
    reference = build_annotator(annotator=REFERENCE_STRATEGY)
    candidates = {
        name: build_annotator(annotator=name)
        for name in consequence_strategies()
        if name != REFERENCE_STRATEGY
    }
    for label, overlaps in labelled:
        report = conformance_gate(overlaps, reference=reference, candidates=candidates)
        if not report.ok:
            suffix = "+" if report.truncated else ""
            logger.error(
                "conformance gate failed on %s (%d%s mismatches):", label, len(report.mismatches), suffix
            )
            for mismatch in report.mismatches:
                logger.error("  %s", mismatch)
            return False
        logger.info(
            "conformance gate passed on %s: %d overlaps, %d strategies",
            label,
            report.overlaps,
            len(report.candidates),
        )
    return True


def run_info(config: Config, verb: str) -> RunInfo:
    """Return what produced this run."""
    spec = CorpusSpec.from_options(config.corpus)
    return RunInfo(
        verb=verb,
        corpus=spec.as_dict(),
        corpus_parameters=spec.parameter_hash(),
        git=git_state(),
        machine=machine_info(),
        timestamp=utc_timestamp(),
        repeats=config.bench.repeats,
        warmup=config.bench.warmup,
    )


def record_results(config: Config, results: BenchResults) -> ExitCode:
    """Write the results if asked, and compare them against a baseline if one was given."""
    bench = config.bench
    if bench.results_file is not None:
        results.write(bench.results_file)
        logger.info("wrote results to %s", bench.results_file)
    if bench.baseline_file is None:
        return ExitCode.OK
    if not bench.baseline_file.exists():
        logger.error("baseline %s does not exist; write one first with --results_file", bench.baseline_file)
        return ExitCode.USAGE

    comparison = compare_results(
        BenchResults.load(bench.baseline_file), results, threshold=bench.regression_threshold
    )
    if comparison.corpus_changed:
        logger.error("the baseline was measured over a different corpus; the comparison would mean nothing")
        return ExitCode.USAGE
    for change in comparison.cost_changes:
        logger.warning("cost changed -- a behaviour change, not a timing one: %s", change)
    for row in comparison.unmatched:
        logger.warning("no baseline row for %s", row)
    for regression in comparison.regressions:
        logger.error("regression beyond %.0f%%: %s", bench.regression_threshold * 100, regression.line())
    if comparison.ok:
        logger.info(
            "no regression beyond %.0f%% against %s", bench.regression_threshold * 100, bench.baseline_file
        )
    return ExitCode.OK if comparison.ok else ExitCode.FAILURE
