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
"""Run statistics: the term-level summary a run writes by default, and the opt-in fact profile.

Two planes, kept as two types on purpose. ``RunSummary`` costs about a millisecond a run and is
on by default; ``RunTrace`` (``ensembl.vep.trace.run``) is fact-level, measured at 1.3-1.7x on the
hottest path, and exists only when tracing is on. One merged type would stop the signature saying
what it costs to produce. Both satisfy ``Reportable``, so the renderers never branch on which they
were handed.
"""

from __future__ import annotations

from ensembl.vep.stats.report import Reportable, render_reportable
from ensembl.vep.stats.summary import RunSummary, RunSummaryCollector, RunSummaryDelta

__all__ = ["Reportable", "RunSummary", "RunSummaryCollector", "RunSummaryDelta", "render_reportable"]
