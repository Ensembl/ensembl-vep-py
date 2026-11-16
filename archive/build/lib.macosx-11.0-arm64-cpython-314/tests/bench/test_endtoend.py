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
"""Tests for timing whole command-line runs."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from ensembl.vep.bench.endtoend import (
    MISMATCH_LIMIT,
    RunOutcome,
    child_environment,
    consequence_sets,
    output_mismatches,
    repository_path,
    run_command,
    run_once,
    summarise,
    time_runs,
)
from ensembl.vep.bench.suites import Configuration
from ensembl.vep.bench.workspace import Workspace

HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def annotated_vcf(path: Path, entries: dict[tuple[str, str, str, str], list[str]]) -> Path:
    """Write a minimal annotated VCF, one CSQ entry per term list given for a variant."""
    lines = [HEADER.rstrip("\n")]
    for (region, position, reference, alternate), terms in entries.items():
        csq = ",".join(f"{alternate}|{term}|MODIFIER" for term in terms)
        lines.append("\t".join((region, position, ".", reference, alternate, ".", ".", f"CSQ={csq}")))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


@pytest.fixture(name="workspace")
def fixture_workspace(tmp_path: Path) -> Workspace:
    """Return a workspace whose paths exist only as names; no corpus is written."""
    return Workspace(directory=tmp_path)


class TestRepositoryPath:
    """Pointing a run at this working copy."""

    def test_it_is_the_directory_holding_the_package(self) -> None:
        """The path is what ``ensembl.vep`` was imported from, so a run matches this process."""
        assert (repository_path() / "ensembl" / "vep" / "__init__.py").exists()

    def test_the_tree_goes_ahead_of_an_existing_python_path(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """An installed copy on PYTHONPATH cannot win over the tree being measured."""
        monkeypatch.setenv("PYTHONPATH", "/somewhere/else")

        value = child_environment()["PYTHONPATH"]

        assert value.startswith(str(repository_path()))
        assert value.endswith("/somewhere/else")

    def test_the_tree_is_set_when_nothing_was_there(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """With no PYTHONPATH set, the tree becomes the whole value."""
        monkeypatch.delenv("PYTHONPATH", raising=False)

        assert child_environment()["PYTHONPATH"] == str(repository_path())


class TestRunCommand:
    """The argv one configuration runs as."""

    def test_it_invokes_this_interpreter_and_the_cli_module(self, workspace: Workspace) -> None:
        """The module is run, not an installed entry point that could be another version."""
        command = run_command(
            Configuration(strategy="consequence.linear"),
            processes=1,
            workspace=workspace,
            output=workspace.directory / "out.vcf",
        )

        assert command[:5] == [sys.executable, "-m", "ensembl.vep.cli", "run", "default"]

    def test_caches_become_their_flags(self, workspace: Workspace) -> None:
        """A cache in the set is switched on, and one outside it explicitly off."""
        command = run_command(
            Configuration(strategy="consequence.columnar", caches=frozenset({"zonemap", "codon"})),
            processes=4,
            workspace=workspace,
            output=workspace.directory / "out.vcf",
        )

        assert "--zone_map" in command
        assert "--codon_cache" in command
        assert "--no-mapper_cache" in command
        assert command[command.index("--processes") + 1] == "4"

    def test_it_annotates_the_indexed_corpus(self, workspace: Workspace) -> None:
        """The run is given the compressed GFF3, which is the one with an index beside it."""
        command = run_command(
            Configuration(strategy="consequence.linear"),
            processes=1,
            workspace=workspace,
            output=workspace.directory / "out.vcf",
        )

        assert command[command.index("--gff") + 1] == str(workspace.indexed_gff3_path)


class TestRunOnce:
    """Running one command and reading what it cost."""

    def test_it_reports_wall_clock_and_peak_memory(self, tmp_path: Path) -> None:
        """Both numbers come back positive for a process that did something."""
        command = [sys.executable, "-c", "x = [0] * 1000000"]

        outcome = run_once(command, log_path=tmp_path / "run.log")

        assert outcome.seconds > 0
        assert outcome.peak_rss_bytes > 0

    def test_a_failing_run_raises_with_its_log(self, tmp_path: Path) -> None:
        """The error carries the run's stderr, so the runner can report why it failed."""
        command = [sys.executable, "-c", "import sys; print('boom', file=sys.stderr); sys.exit(3)"]

        with pytest.raises(subprocess.CalledProcessError) as raised:
            run_once(command, log_path=tmp_path / "run.log")

        assert raised.value.returncode == 3
        assert "boom" in raised.value.stderr

    def test_warmup_runs_are_not_timed(self, tmp_path: Path) -> None:
        """Only the repeats come back as outcomes."""
        command = [sys.executable, "-c", "pass"]

        outcomes = time_runs(command, log_path=tmp_path / "run.log", repeats=2, warmup=1)

        assert len(outcomes) == 2


class TestSummarise:
    """Turning repeats into a row."""

    def test_the_median_run_sets_the_rates(self) -> None:
        """Rates follow the median wall clock, and the peak is the largest seen."""
        outcomes = [
            RunOutcome(seconds=1.0, peak_rss_bytes=10),
            RunOutcome(seconds=2.0, peak_rss_bytes=30),
            RunOutcome(seconds=5.0, peak_rss_bytes=20),
        ]

        row = summarise(
            outcomes, corpus="default", configuration="consequence.linear", processes=1, variants=1000
        )

        assert row.seconds == 2.0
        assert row.ms_per_1k == pytest.approx(2000.0)
        assert row.variants_per_second == pytest.approx(500.0)
        assert row.peak_rss_bytes == 30

    def test_an_empty_input_rates_as_zero(self) -> None:
        """A corpus with no variants cannot divide, and reports zero rather than failing."""
        row = summarise(
            [RunOutcome(seconds=1.0, peak_rss_bytes=1)],
            corpus="default",
            configuration="consequence.linear",
            processes=1,
            variants=0,
        )

        assert row.ms_per_1k == 0.0
        assert row.variants_per_second == 0.0

    def test_the_row_serialises_every_measurement(self) -> None:
        """The results file carries what the table shows."""
        row = summarise(
            [RunOutcome(seconds=1.0, peak_rss_bytes=8)],
            corpus="heavy",
            configuration="consequence.columnar +zonemap",
            processes=2,
            variants=10,
        )

        payload = row.to_json()

        assert payload["corpus"] == "heavy"
        assert payload["configuration"] == "consequence.columnar +zonemap"
        assert payload["processes"] == 2
        assert payload["peak_rss_bytes"] == 8


class TestConsequenceSets:
    """Reading terms back out of an annotated VCF."""

    def test_terms_are_collected_per_variant(self, tmp_path: Path) -> None:
        """Every CSQ entry's terms land in one set for the variant."""
        path = annotated_vcf(
            tmp_path / "out.vcf", {("1", "100", "A", "G"): ["missense_variant", "splice_region_variant"]}
        )

        assert consequence_sets(path) == {
            ("1", "100", "A", "G"): frozenset({"missense_variant", "splice_region_variant"})
        }

    def test_ampersands_split_into_separate_terms(self, tmp_path: Path) -> None:
        """A compound consequence counts as the terms it names."""
        path = annotated_vcf(
            tmp_path / "out.vcf", {("1", "1", "A", "G"): ["splice_donor_variant&intron_variant"]}
        )

        assert consequence_sets(path)[("1", "1", "A", "G")] == frozenset(
            {"splice_donor_variant", "intron_variant"}
        )

    def test_other_info_fields_are_ignored(self, tmp_path: Path) -> None:
        """Only the CSQ field carries terms; anything else in INFO is not one."""
        path = tmp_path / "out.vcf"
        path.write_text(
            HEADER + "1\t100\t.\tA\tG\t.\t.\tAC=1;CSQ=G|missense_variant|MODIFIER;DP=9\n",
            encoding="utf-8",
        )

        assert consequence_sets(path)[("1", "100", "A", "G")] == frozenset({"missense_variant"})

    def test_a_row_without_enough_columns_is_skipped(self, tmp_path: Path) -> None:
        """A truncated line is ignored rather than raising."""
        path = tmp_path / "out.vcf"
        path.write_text(HEADER + "1\t100\t.\n", encoding="utf-8")

        assert not consequence_sets(path)


class TestOutputMismatches:
    """The end-to-end gate."""

    def test_identical_outputs_agree(self, tmp_path: Path) -> None:
        """Two runs that produced the same terms produce no mismatches."""
        entries = {("1", "100", "A", "G"): ["missense_variant"]}
        reference = annotated_vcf(tmp_path / "reference.vcf", entries)
        candidate = annotated_vcf(tmp_path / "candidate.vcf", entries)

        assert not output_mismatches(reference, candidate)

    def test_differing_terms_are_reported(self, tmp_path: Path) -> None:
        """The mismatch names the variant and both sides."""
        reference = annotated_vcf(tmp_path / "reference.vcf", {("1", "100", "A", "G"): ["missense_variant"]})
        candidate = annotated_vcf(
            tmp_path / "candidate.vcf", {("1", "100", "A", "G"): ["synonymous_variant"]}
        )

        mismatches = output_mismatches(reference, candidate)

        assert len(mismatches) == 1
        assert "1:100 A/G" in mismatches[0]
        assert "synonymous_variant" in mismatches[0]

    def test_a_missing_variant_is_reported(self, tmp_path: Path) -> None:
        """A variant the candidate never emitted is a mismatch, not a pass."""
        reference = annotated_vcf(tmp_path / "reference.vcf", {("1", "100", "A", "G"): ["missense_variant"]})
        candidate = annotated_vcf(tmp_path / "candidate.vcf", {})

        assert "missing from the output" in output_mismatches(reference, candidate)[0]

    def test_an_extra_variant_is_reported(self, tmp_path: Path) -> None:
        """A variant only the candidate emitted is a mismatch too."""
        reference = annotated_vcf(tmp_path / "reference.vcf", {})
        candidate = annotated_vcf(tmp_path / "candidate.vcf", {("1", "100", "A", "G"): ["missense_variant"]})

        assert "not in the reference output" in output_mismatches(reference, candidate)[0]

    def test_collection_stops_at_the_limit(self, tmp_path: Path) -> None:
        """One broken configuration cannot flood the log."""
        expected = {(str(position), "100", "A", "G"): ["missense_variant"] for position in range(50)}
        actual = {key: ["synonymous_variant"] for key in expected}
        reference = annotated_vcf(tmp_path / "reference.vcf", expected)
        candidate = annotated_vcf(tmp_path / "candidate.vcf", actual)

        assert len(output_mismatches(reference, candidate)) == MISMATCH_LIMIT
