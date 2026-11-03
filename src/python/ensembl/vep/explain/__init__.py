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
"""Overlap selection, annotation, comparison, and summaries for ``vep explain``.

Separate from ``runner/explain.py``, which is only the subcommand: reading the options, choosing
which rendering was asked for, and writing it out. Everything a second caller could want -- naming
and finding the overlaps of a real input, annotating one with tracing on, diffing two derivations,
summarising an overlap -- lives here.

The package sits above ``trace``: a trace is recorded during annotation, and these modules read
one afterwards.
"""

from __future__ import annotations

from ensembl.vep.explain.annotate import annotate_with_trace
from ensembl.vep.explain.compare import diff_terms, diff_traces
from ensembl.vep.explain.selection import (
    FeatureSelector,
    NamedOverlap,
    OverlapSelection,
    VariantSelector,
)
from ensembl.vep.explain.summary import format_summary, summary_html, summary_rows

__all__ = [
    "FeatureSelector",
    "NamedOverlap",
    "OverlapSelection",
    "VariantSelector",
    "annotate_with_trace",
    "diff_terms",
    "diff_traces",
    "format_summary",
    "summary_html",
    "summary_rows",
]
