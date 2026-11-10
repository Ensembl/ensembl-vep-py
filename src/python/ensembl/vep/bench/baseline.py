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
"""Benchmark results on disk, and comparing a run against an earlier one.

A number without its context is not comparable, so a results file carries what produced it: the
corpus spec and its parameter hash, the git revision and whether the tree was dirty, the machine,
and the timing policy. A comparison across two different corpora is refused rather than reported.

Timings regress past a threshold. The declared cost per overlap is deterministic, so any change in
it is a change in behaviour -- reported as such, and not held to a threshold.
"""

from __future__ import annotations

import json
import os
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import attrs
from attrs import field, frozen

__all__ = [
    "COST_METRIC",
    "FORMAT_VERSION",
    "TIMING_METRICS",
    "BenchResults",
    "Comparison",
    "Regression",
    "RunInfo",
    "compare_results",
    "git_state",
    "machine_info",
    "utc_timestamp",
]

FORMAT_VERSION = 1
TIMING_METRICS = ("ms_per_1k", "batch_ms_per_1k", "ns_per_case")
"""The metrics a regression is judged on. Lower is better for all of them."""
COST_METRIC = "cost_per_context"


def machine_info() -> dict[str, object]:
    """Return what about this machine could explain a different timing."""
    return {
        "system": platform.system(),
        "release": platform.release(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "python": platform.python_version(),
        "cpus": os.cpu_count(),
    }


def git_state() -> dict[str, object]:
    """Return the git revision and whether the tree is dirty, or Nones outside a repository."""
    try:
        revision = subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True, timeout=10
        ).stdout.strip()
        status = subprocess.run(
            ["git", "status", "--porcelain"], capture_output=True, text=True, check=True, timeout=10
        ).stdout
    except (OSError, subprocess.SubprocessError):
        return {"revision": None, "dirty": None}
    return {"revision": revision, "dirty": bool(status.strip())}


def utc_timestamp() -> str:
    """Return the current UTC time, to the second."""
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


@frozen(kw_only=True, auto_attribs=True)
class RunInfo:
    """What produced a set of results.

    Attributes:
        verb (str): The bench verb that ran.
        corpus (dict[str, object]): The corpus spec.
        corpus_parameters (str): The spec's parameter hash.
        git (dict[str, object]): Revision and dirtiness.
        machine (dict[str, object]): See machine_info.
        timestamp (str): When the run finished.
        repeats (int): Timed samples per measurement.
        warmup (int): Untimed passes per measurement.
    """

    verb: str
    corpus: dict[str, object] = field(factory=dict)
    corpus_parameters: str = field(default="")
    git: dict[str, object] = field(factory=dict)
    machine: dict[str, object] = field(factory=dict)
    timestamp: str = field(default="")
    repeats: int = field(default=0)
    warmup: int = field(default=0)


@frozen(kw_only=True, auto_attribs=True)
class BenchResults:
    """A run's rows and what produced them.

    Attributes:
        info (RunInfo): The run's context.
        rows (tuple[dict[str, object], ...]): One row per measurement, each with ``corpus`` and
            ``configuration`` keys plus its metrics.
    """

    info: RunInfo
    rows: tuple[dict[str, object], ...] = field(default=())

    def to_json(self) -> str:
        """Return the results as indented JSON."""
        payload = {"format_version": FORMAT_VERSION, "info": attrs.asdict(self.info), "rows": list(self.rows)}
        return json.dumps(payload, indent=1, sort_keys=True) + "\n"

    @classmethod
    def from_json(cls, text: str) -> BenchResults:
        """Read results back.

        Raises:
            ValueError: If the file was written in a different format version.
        """
        payload = json.loads(text)
        if payload.get("format_version") != FORMAT_VERSION:
            raise ValueError(f"results format {payload.get('format_version')!r} is not {FORMAT_VERSION}")
        return cls(info=RunInfo(**payload["info"]), rows=tuple(payload["rows"]))

    def write(self, path: Path) -> None:
        """Write the results, creating the directory if needed."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.to_json(), encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> BenchResults:
        """Read results from a file."""
        return cls.from_json(path.read_text(encoding="utf-8"))


@frozen(kw_only=True, auto_attribs=True)
class Regression:
    """A timing that got slower than the threshold allows.

    Attributes:
        corpus (str): The corpus.
        configuration (str): The configuration.
        metric (str): The metric that regressed.
        baseline (float): Its value in the baseline.
        current (float): Its value now.
    """

    corpus: str
    configuration: str
    metric: str
    baseline: float
    current: float

    @property
    def change(self) -> float:
        """Return the fractional change, where 0.1 is ten percent slower."""
        return self.current / self.baseline - 1

    def line(self) -> str:
        """Return the regression as one readable line."""
        return (
            f"{self.corpus} / {self.configuration}: {self.metric} {self.baseline:.3f} -> {self.current:.3f} "
            f"({self.change:+.1%})"
        )


@frozen(kw_only=True, auto_attribs=True)
class Comparison:
    """A run compared against a baseline.

    Attributes:
        corpus_changed (bool): The two were run over different corpora, so nothing else is meaningful.
        regressions (tuple[Regression, ...]): Timings past the threshold.
        cost_changes (tuple[str, ...]): Configurations whose declared cost per overlap moved.
        unmatched (tuple[str, ...]): Rows of the current run the baseline has no counterpart for.
    """

    corpus_changed: bool = field(default=False)
    regressions: tuple[Regression, ...] = field(default=())
    cost_changes: tuple[str, ...] = field(default=())
    unmatched: tuple[str, ...] = field(default=())

    @property
    def ok(self) -> bool:
        """Return whether the run is comparable and no timing regressed."""
        return not (self.corpus_changed or self.regressions)


def _as_float(value: object) -> float | None:
    """Return a metric as a float, or None if it is not a number."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value)


def compare_results(baseline: BenchResults, current: BenchResults, *, threshold: float) -> Comparison:
    """Compare a run against a baseline, row by row.

    Args:
        baseline (BenchResults): The earlier run.
        current (BenchResults): This run.
        threshold (float): Fractional slowdown that counts as a regression.

    Returns:
        Comparison: What changed.
    """
    if baseline.info.corpus_parameters != current.info.corpus_parameters:
        return Comparison(corpus_changed=True)

    earlier = {(str(row["corpus"]), str(row["configuration"])): row for row in baseline.rows}
    regressions: list[Regression] = []
    cost_changes: list[str] = []
    unmatched: list[str] = []
    for row in current.rows:
        key = (str(row["corpus"]), str(row["configuration"]))
        previous = earlier.get(key)
        if previous is None:
            unmatched.append(f"{key[0]} / {key[1]}")
            continue
        for metric in TIMING_METRICS:
            before, after = _as_float(previous.get(metric)), _as_float(row.get(metric))
            if before is not None and after is not None and before > 0 and after > before * (1 + threshold):
                regressions.append(
                    Regression(
                        corpus=key[0], configuration=key[1], metric=metric, baseline=before, current=after
                    )
                )
        before_cost, after_cost = _as_float(previous.get(COST_METRIC)), _as_float(row.get(COST_METRIC))
        if before_cost is not None and after_cost is not None and before_cost != after_cost:
            cost_changes.append(f"{key[0]} / {key[1]}: cost/ctx {before_cost:.3f} -> {after_cost:.3f}")
    return Comparison(
        regressions=tuple(regressions), cost_changes=tuple(cost_changes), unmatched=tuple(unmatched)
    )
