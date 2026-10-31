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
"""Options for explaining a single overlap's derivation."""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

from attrs import define, field

from ensembl.vep.config.option_spec import option
from ensembl.vep.config.report import GRAPH_FORMATS
from ensembl.vep.config.sections.base import ConfigSection

__all__ = ["ExplainConfig"]

EXPLAIN_FORMATS = ("text",) + GRAPH_FORMATS
"""A derivation renders as text or as a diagram; it has no tabular form."""


@define(frozen=True, kw_only=True, auto_attribs=True)
class ExplainConfig(ConfigSection):
    """Which overlap to explain, and how to render the derivation.

    Attributes:
        name (str | None): Named overlap to select, for a synthetic corpus.
        index (int | None): Overlap position in the corpus, as an alternative to ``name``.
        variant (str | None): Variant to select from a real input, as ``CHR:POS:REF>ALT`` or
            ``CHR:POS:ID``.
        allele_num (int | None): Alternate allele a ``CHR:POS:ID`` variant selects, from 1; None for
            the first.
        transcript (str | None): Transcript stable id, to disambiguate when several overlap.
        compare (tuple[str, ...]): Two annotator names to diff against each other.
        list_overlaps (bool): List the selectable overlaps and exit.
        graph (bool): Draw the derivation as a diagram rather than printing it.
        output_format (str | None): The format to render in, or None for the default of the kind asked
            for. Spelled ``--output_format`` because ``--format`` is the input format.
        output_file (Path | None): Write the explanation here instead of to stdout.
    """

    section_name: ClassVar[str] = "explain"

    name: str | None = field(default=None, metadata=option("named overlap to explain"))
    index: int | None = field(default=None, metadata=option("overlap index to explain", metavar="N"))
    variant: str | None = field(
        default=None,
        metadata=option(
            "variant to explain, by alleles or by the input's own ID", metavar="CHR:POS:REF>ALT|CHR:POS:ID"
        ),
    )
    allele_num: int | None = field(
        default=None,
        metadata=option("alternate allele of a CHR:POS:ID variant to explain (default: 1)", metavar="N"),
    )
    transcript: str | None = field(default=None, metadata=option("transcript stable id", metavar="STABLE_ID"))
    compare: tuple[str, ...] = field(
        default=(), metadata=option("two annotators to diff", metavar="NAME,NAME")
    )
    list_overlaps: bool = field(
        default=False,
        metadata=option("list selectable overlaps and exit (see also: vep corpus describe)", flag="list"),
    )
    graph: bool = field(default=False, metadata=option("draw the derivation as a diagram"))
    output_format: str | None = field(
        default=None,
        metadata=option("report format (default: text, svg with --graph)", choices=EXPLAIN_FORMATS),
    )
    output_file: Path | None = field(default=None, metadata=option("write the explanation here"))
