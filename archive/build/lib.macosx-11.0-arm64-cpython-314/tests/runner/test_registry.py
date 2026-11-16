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
"""Unit tests for ensembl.vep.runner.registry."""

from __future__ import annotations

import pytest

from ensembl.vep.config import Config
from ensembl.vep.runner import registry as registry_module
from ensembl.vep.runner.base import ExitCode, Runner
from ensembl.vep.runner.registry import (
    available_runners,
    get_runner,
    register_runner,
    runner_group,
    runner_groups,
)


class _DummyRunner(Runner):
    """Minimal runner used to test registration and lookup."""

    name = "test.dummy"
    sections = ("logging",)
    description = "dummy runner"

    def __init__(self, config: Config) -> None:
        self.config = config

    def run(self) -> ExitCode:
        """Return a process exit code."""
        return ExitCode.OK


class TestAvailableRunners:
    """The registry always exposes the built-in runner."""

    def test_includes_the_built_in_runner(self) -> None:
        """Importing the package registers the default command."""
        assert "run.default" in available_runners()

    def test_the_run_command_uses_the_renamed_runner_class(self) -> None:
        """The built-in run command is provided by Runner."""
        assert get_runner("run.default").__name__ == "Runner"


class TestLookup:
    """Lookup behavior for registered and unregistered runners."""

    def test_unknown_name_raises_key_error_with_the_available_list(self) -> None:
        """Unknown names raise KeyError with the registry summary in the message."""
        with pytest.raises(KeyError, match="Unknown command 'nope'"):
            get_runner("nope")

    def test_registered_runner_can_be_resolved(self) -> None:
        """The registry can store and return a runner class."""
        register_runner("test.dummy")(_DummyRunner)
        try:
            assert get_runner("test.dummy") is _DummyRunner
        finally:
            registry_module._RUNNER_REGISTRY.pop("test.dummy", None)  # pylint: disable=protected-access


class TestGrouping:
    """Runner names split into groups the same way the CLI presents them."""

    def test_names_split_into_domain_and_verb(self) -> None:
        """A name without a dot is its own domain."""
        assert runner_group("bench.matrix") == "bench"
        assert runner_group("run") == "run"

    def test_groups_collect_the_verbs_of_a_domain(self) -> None:
        """A domain with one runner named after it is a subcommand in its own right."""
        assert runner_groups()["run"] == ("run.default", "run.haplo")

    def test_every_group_is_sorted(self) -> None:
        """Grouped runner names are returned in sorted order."""
        register_runner("test.beta")(_DummyRunner)
        register_runner("test.alpha")(_DummyRunner)
        try:
            assert runner_groups()["test"] == ("test.alpha", "test.beta")
        finally:
            registry_module._RUNNER_REGISTRY.pop("test.alpha", None)  # pylint: disable=protected-access
            registry_module._RUNNER_REGISTRY.pop("test.beta", None)  # pylint: disable=protected-access
