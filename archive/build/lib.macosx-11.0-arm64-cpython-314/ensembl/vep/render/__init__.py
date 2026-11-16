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
"""Presentation. Domain-free: nothing here imports the rest of the package."""

from ensembl.vep.render.chart import Slice, render_bars, render_legend, render_pie, render_series
from ensembl.vep.render.page import Section, render_page
from ensembl.vep.render.table import TableFormat, render_table
from ensembl.vep.render.report import ChartKind, ReportFormat, ReportSection, render_report
from ensembl.vep.render.theme import Colour, Font, Metric, stylesheet

__all__ = [
    "ChartKind",
    "Colour",
    "Font",
    "Metric",
    "ReportFormat",
    "ReportSection",
    "Section",
    "Slice",
    "TableFormat",
    "render_bars",
    "render_legend",
    "render_page",
    "render_pie",
    "render_report",
    "render_series",
    "render_table",
    "stylesheet",
]
