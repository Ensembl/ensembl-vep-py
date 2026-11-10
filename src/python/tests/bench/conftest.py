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
"""Shared fixtures for the bench test suite."""

from __future__ import annotations

from collections.abc import Iterator

import pytest

from ensembl.vep.annotators.consequence.facts.base import REGISTRY


@pytest.fixture(autouse=True)
def restore_fact_registry() -> Iterator[None]:
    """Restore the fact registry after a test that instruments it.

    Counting fact computations swaps every fact's ``compute`` process-wide, relying on its own
    ``finally`` to put them back; this makes a failure partway through unable to leak.

    Yields:
        None: Control, while the test runs.
    """
    registered = dict(REGISTRY)
    try:
        yield
    finally:
        REGISTRY.clear()
        REGISTRY.update(registered)
