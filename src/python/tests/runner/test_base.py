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
"""Unit tests for ensembl.vep.runner.base."""

from __future__ import annotations

from typing import cast

from ensembl.vep.config import Config
from ensembl.vep.runner.base import ExitCode, Runner


class _StructuralRunner:
    """A concrete object that satisfies the runner protocol structurally."""

    name = "test.structural"
    sections = ("logging",)
    description = "structural test runner"

    def __init__(self, config: Config) -> None:
        self.config = config

    def run(self) -> ExitCode:
        """Return a success code."""
        return ExitCode.OK


class TestRunnerProtocol:
    """The runner base behaves like a runtime-checkable protocol."""

    def test_exit_code_values_are_stable(self) -> None:
        """The process exit values are part of the public contract."""
        assert ExitCode.OK == 0
        assert ExitCode.FAILURE == 1
        assert ExitCode.USAGE == 2

    def test_matching_objects_satisfy_the_protocol(self) -> None:
        """A structurally compatible object passes ``isinstance``."""
        runner = _StructuralRunner(config=cast(Config, object()))
        assert isinstance(runner, Runner)
        assert runner.name == "test.structural"
        assert runner.sections == ("logging",)
        assert runner.description == "structural test runner"
        assert runner.run() == ExitCode.OK
