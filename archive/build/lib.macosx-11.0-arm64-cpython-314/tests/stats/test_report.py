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
"""Tests for the report contract and its two implementers."""

from __future__ import annotations

import json

from attrs import frozen

import ensembl.vep.annotators  # noqa: F401  pylint: disable=unused-import

# ^ keeps the annotators-before-caches import order.
from ensembl.vep.render.chart import Slice
from ensembl.vep.render.report import ChartKind, ReportSection
from ensembl.vep.stats import Reportable, RunSummary, render_reportable
from ensembl.vep.trace.run import RunTrace


@frozen(kw_only=True, auto_attribs=True)
class _TinyReport:
    """The smallest thing that satisfies Reportable."""

    @property
    def title(self) -> str:
        """Title."""
        return "Tiny"

    @property
    def subtitle(self) -> str:
        """Subtitle."""
        return "two sections"

    def sections(self) -> tuple[ReportSection, ...]:
        """Two sections: a table and a chart."""
        return (
            ReportSection(identifier="table", title="A table", headers=("Name", "Count"), rows=(("a", 1),)),
            ReportSection(
                identifier="chart",
                title="A chart",
                slices=(Slice(label="x", value=2), Slice(label="y", value=3)),
                chart=ChartKind.PIE,
            ),
        )


class TestReportable:
    """Both producers satisfy the one Protocol, structurally."""

    def test_a_run_summary_is_reportable(self) -> None:
        """The term-level plane."""
        assert isinstance(RunSummary(), Reportable)

    def test_a_run_trace_is_reportable(self) -> None:
        """The fact-level plane, without trace/ importing stats/."""
        assert isinstance(RunTrace.fold([]), Reportable)

    def test_an_empty_report_still_renders_in_every_format(self) -> None:
        """Empty collections, never None: the renderer draws nothing rather than branching."""
        for report in (RunSummary(), RunTrace.fold([])):
            for report_format in ("html", "text", "json"):
                assert render_reportable(report, report_format)


class TestRenderReportable:
    """The three formats, over one report."""

    def test_html_is_a_page_with_a_table_and_a_chart(self) -> None:
        """HTML carries the table and the pie."""
        rendered = render_reportable(_TinyReport(), "html")
        assert "<table" in rendered
        assert "<svg" in rendered
        assert "Tiny" in rendered

    def test_text_lists_slices_as_a_table(self) -> None:
        """Text has no charts, so the quantities become a two-column table."""
        rendered = render_reportable(_TinyReport(), "text")
        assert "## A chart" in rendered
        assert "x" in rendered and "2" in rendered

    def test_json_is_plain_data(self) -> None:
        """JSON round-trips the sections as data, with no rendering in it."""
        payload = json.loads(render_reportable(_TinyReport(), "json"))
        assert payload["title"] == "Tiny"
        assert payload["sections"][1]["slices"] == [{"label": "x", "value": 2}, {"label": "y", "value": 3}]
