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
"""Shared fixtures for the runner test suite."""

from __future__ import annotations

import logging
from collections.abc import Iterator

import pytest

from ensembl.vep.runner import registry


@pytest.fixture(autouse=True)
def restore_package_logging() -> Iterator[None]:
    """Restores the ``ensembl.vep`` logger after a test that runs the entry point.

    ``setup_logging`` configures that logger for the whole process, which is right for a command
    but leaks between tests: a level set here filters records before they reach the root logger,
    so a later test using ``caplog.at_level`` sees nothing and fails for a reason that has nothing
    to do with what it is testing.

    Yields:
        None: Control, while the test runs.
    """
    logger = logging.getLogger("ensembl.vep")
    level = logger.level
    handlers = list(logger.handlers)
    propagate = logger.propagate
    try:
        yield
    finally:
        logger.setLevel(level)
        logger.handlers = handlers
        logger.propagate = propagate


@pytest.fixture(autouse=True)
def restore_runner_registry() -> Iterator[None]:
    """Restore the runner registry after a test that registers one of its own.

    Registration is a module-level side effect, so a runner registered here would otherwise be
    seen by every later test that walks ``available_runners()``.

    Yields:
        None: Control, while the test runs.
    """
    registered = dict(registry._RUNNER_REGISTRY)  # pylint: disable=protected-access
    try:
        yield
    finally:
        registry._RUNNER_REGISTRY.clear()  # pylint: disable=protected-access
        registry._RUNNER_REGISTRY.update(registered)  # pylint: disable=protected-access
