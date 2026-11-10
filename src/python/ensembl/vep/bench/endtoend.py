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
"""Timing whole runs of the command line, one configuration at a time.

``bench matrix`` times an annotator over overlaps that are already built. This times what a user
runs: reading the VCF, loading transcripts from the indexed GFF3, fetching sequence, annotating,
and writing the output. The difference between the two is the point -- a strategy that wins several
times over in isolation moves the end-to-end number by whatever share of the run annotation is.

Each run is a fresh process, so no configuration is timed against another's warmed module state,
and the process cost a real run pays is counted rather than amortised away. The child is given this
tree on ``PYTHONPATH``, so the benchmark measures the working copy it was started from and never an
``ensembl-vep`` installed elsewhere.
"""

from __future__ import annotations

import logging
import os
import statistics
import subprocess
import sys
import time
from pathlib import Path

from attrs import field, frozen

import ensembl.vep
from ensembl.vep.bench.suites import Configuration
from ensembl.vep.bench.workspace import Workspace

__all__ = [
    "MISMATCH_LIMIT",
    "EndToEndRow",
    "RunOutcome",
    "consequence_sets",
    "output_mismatches",
    "repository_path",
    "run_command",
    "run_once",
    "time_runs",
]

logger = logging.getLogger(__name__)

MISMATCH_LIMIT = 20
"""Output mismatches collected before the gate stops, so one broken configuration cannot flood it."""

_CSQ_PREFIX = "CSQ="
"""The INFO key the consequence terms are written under."""

_RSS_UNIT_BYTES = 1 if sys.platform == "darwin" else 1024
"""``ru_maxrss`` is bytes on macOS and kibibytes on Linux."""


@frozen(kw_only=True, auto_attribs=True)
class RunOutcome:
    """What one run of the command line cost.

    Attributes:
        seconds (float): Wall clock of the whole process.
        peak_rss_bytes (int): The high-water mark of the largest process in the tree -- not the sum
            across workers, which is what a sampling profiler would report instead.
    """

    seconds: float
    peak_rss_bytes: int


@frozen(kw_only=True, auto_attribs=True)
class EndToEndRow:
    """One configuration's measurements over one corpus.

    Attributes:
        corpus (str): The corpus kind annotated.
        configuration (str): The configuration's label, as ``bench matrix`` spells it.
        processes (int): Worker processes the run was given.
        variants (int): Variants in the input.
        seconds (float): Median wall clock over the timed repeats.
        ms_per_1k (float): Median wall clock per thousand variants.
        variants_per_second (float): Throughput at the median.
        peak_rss_bytes (int): Largest peak seen over the timed repeats.
    """

    corpus: str
    configuration: str
    processes: int
    variants: int
    seconds: float
    ms_per_1k: float
    variants_per_second: float
    peak_rss_bytes: int = field(default=0)

    def to_json(self) -> dict[str, object]:
        """Return the row as plain JSON values."""
        return {
            "corpus": self.corpus,
            "configuration": self.configuration,
            "processes": self.processes,
            "variants": self.variants,
            "seconds": self.seconds,
            "ms_per_1k": self.ms_per_1k,
            "variants_per_second": self.variants_per_second,
            "peak_rss_bytes": self.peak_rss_bytes,
        }


def repository_path() -> Path:
    """Return the directory holding the ``ensembl`` package this process imported.

    What makes the benchmark measure the working copy it runs from: the child gets this on
    ``PYTHONPATH``, which takes precedence over anything installed into the interpreter.
    """
    return Path(ensembl.vep.__file__).resolve().parents[2]


def child_environment() -> dict[str, str]:
    """Return the environment for a run, with this tree ahead of any installed copy."""
    environment = dict(os.environ)
    existing = environment.get("PYTHONPATH")
    tree = str(repository_path())
    environment["PYTHONPATH"] = f"{tree}{os.pathsep}{existing}" if existing else tree
    return environment


def run_command(
    configuration: Configuration, *, processes: int, workspace: Workspace, output: Path
) -> list[str]:
    """Return the command line for one run of one configuration.

    Args:
        configuration (Configuration): The strategy and caches to run with.
        processes (int): Worker processes for the run.
        workspace (Workspace): The corpus the run annotates.
        output (Path): Where the run writes its annotated VCF.

    Returns:
        list[str]: The argv, invoking this tree's CLI module rather than an installed entry point.
    """
    caches = configuration.caches
    return [
        sys.executable, "-m", "ensembl.vep.cli", "run", "default",
        "--input_file", str(workspace.vcf_path),
        "--gff", str(workspace.indexed_gff3_path),
        "--fasta", str(workspace.fasta_path),
        "--output_file", str(output),
        "--force_overwrite",
        "--no_stats",
        "--no_progress",
        "--annotator", configuration.strategy,
        "--batch_size", str(configuration.batch_size),
        "--processes", str(processes),
        "--zone_map" if "zonemap" in caches else "--no-zone_map",
        "--codon_cache" if "codon" in caches else "--no-codon_cache",
        "--mapper_cache" if "mapper" in caches else "--no-mapper_cache",
    ]  # fmt: skip


def run_once(command: list[str], *, log_path: Path) -> RunOutcome:
    """Run one command to completion, timing it and reading its resource usage.

    Stderr goes to a file rather than a pipe: the run is waited for with ``os.wait4``, which cannot
    drain a pipe, and a chatty run would otherwise fill it and deadlock.

    Args:
        command (list[str]): The argv to run.
        log_path (Path): Where the run's stderr is written.

    Returns:
        RunOutcome: The wall clock and peak resident size.

    Raises:
        subprocess.CalledProcessError: If the run fails, carrying its log as output.
    """
    started = time.perf_counter()
    with log_path.open("wb") as log:
        process = subprocess.Popen(  # pylint: disable=consider-using-with
            command, env=child_environment(), stdout=subprocess.DEVNULL, stderr=log
        )
        _pid, status, usage = os.wait4(process.pid, 0)
    elapsed = time.perf_counter() - started

    # The child is already reaped, so tell Popen rather than letting it wait for a process that
    # no longer exists.
    process.returncode = os.waitstatus_to_exitcode(status)
    if process.returncode != 0:
        raise subprocess.CalledProcessError(
            process.returncode, command, stderr=log_path.read_text(encoding="utf-8", errors="replace")
        )
    return RunOutcome(seconds=elapsed, peak_rss_bytes=usage.ru_maxrss * _RSS_UNIT_BYTES)


def time_runs(command: list[str], *, log_path: Path, repeats: int, warmup: int) -> list[RunOutcome]:
    """Run a command untimed ``warmup`` times, then return the outcome of ``repeats`` timed runs.

    Raises:
        subprocess.CalledProcessError: If any run fails.
    """
    for _ in range(warmup):
        run_once(command, log_path=log_path)
    return [run_once(command, log_path=log_path) for _ in range(repeats)]


def summarise(
    outcomes: list[RunOutcome], *, corpus: str, configuration: str, processes: int, variants: int
) -> EndToEndRow:
    """Return the row a configuration's timed runs add up to.

    Args:
        outcomes (list[RunOutcome]): The timed runs, at least one.
        corpus (str): The corpus kind annotated.
        configuration (str): The configuration's label.
        processes (int): Worker processes the runs were given.
        variants (int): Variants in the input.

    Returns:
        EndToEndRow: The median wall clock, what it works out to per variant, and the largest peak.
    """
    seconds = statistics.median(outcome.seconds for outcome in outcomes)
    return EndToEndRow(
        corpus=corpus,
        configuration=configuration,
        processes=processes,
        variants=variants,
        seconds=seconds,
        ms_per_1k=(seconds * 1000.0 * 1000.0 / variants) if variants else 0.0,
        variants_per_second=(variants / seconds) if seconds else 0.0,
        peak_rss_bytes=max(outcome.peak_rss_bytes for outcome in outcomes),
    )


def consequence_sets(output: Path) -> dict[tuple[str, str, str, str], frozenset[str]]:
    """Return each variant's consequence terms, keyed by the variant it was called on.

    Terms are compared as a set per variant rather than as file bytes: entry order can legitimately
    follow the worker count and the buffering, and neither of those is a wrong answer.

    Args:
        output (Path): An annotated VCF.

    Returns:
        dict[tuple[str, str, str, str], frozenset[str]]: Terms by
        ``(region, position, reference, alternates)``.
    """
    terms_by_variant: dict[tuple[str, str, str, str], frozenset[str]] = {}
    with output.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 8:
                continue
            key = (columns[0], columns[1], columns[3], columns[4])
            terms: set[str] = set()
            for entry in _csq_entries(columns[7]):
                terms.update(term for term in entry.split("&") if term)
            terms_by_variant[key] = frozenset(terms)
    return terms_by_variant


def _csq_entries(info: str) -> list[str]:
    """Return the consequence field of every CSQ entry in one INFO column."""
    entries: list[str] = []
    for field_text in info.split(";"):
        if not field_text.startswith(_CSQ_PREFIX):
            continue
        for entry in field_text[len(_CSQ_PREFIX) :].split(","):
            parts = entry.split("|")
            if len(parts) > 1:
                entries.append(parts[1])
    return entries


def output_mismatches(reference: Path, candidate: Path) -> list[str]:
    """Return how a candidate's output differs from the reference's, up to ``MISMATCH_LIMIT``.

    Args:
        reference (Path): The reference configuration's annotated VCF.
        candidate (Path): The configuration being checked.

    Returns:
        list[str]: One description per differing variant; empty when the two agree.
    """
    expected = consequence_sets(reference)
    actual = consequence_sets(candidate)
    mismatches: list[str] = []

    for key, terms in expected.items():
        if len(mismatches) >= MISMATCH_LIMIT:
            break
        found = actual.get(key)
        if found is None:
            mismatches.append(f"{_variant(key)}: missing from the output")
        elif found != terms:
            mismatches.append(f"{_variant(key)}: expected {sorted(terms)}, got {sorted(found)}")

    for key in actual:
        if len(mismatches) >= MISMATCH_LIMIT:
            break
        if key not in expected:
            mismatches.append(f"{_variant(key)}: not in the reference output")
    return mismatches


def _variant(key: tuple[str, str, str, str]) -> str:
    """Return a variant key as ``region:position ref/alt``."""
    region, position, reference, alternates = key
    return f"{region}:{position} {reference}/{alternates}"
