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
"""Report contracts and rendering."""

from __future__ import annotations

from typing import Protocol, runtime_checkable

from ensembl.vep.render.report import ReportFormat, ReportSection, render_report

__all__ = ["Reportable", "render_reportable"]


@runtime_checkable
class Reportable(Protocol):
    """Anything that can be rendered as a report.

    Members are properties rather than attributes: the implementers are frozen attrs classes, and
    a Protocol declaring a plain attribute demands a settable one.
    """

    @property
    def title(self) -> str:
        """Return the document title."""

    @property
    def subtitle(self) -> str:
        """Return one line shown under the title."""

    def sections(self) -> tuple[ReportSection, ...]:
        """Return the report sections in rendering order."""


def render_reportable(report: Reportable, report_format: str | ReportFormat) -> str:
    """Render any report in the requested format."""
    return render_report(report.title, report.subtitle, report.sections(), report_format)
