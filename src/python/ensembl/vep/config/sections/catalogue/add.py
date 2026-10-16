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
"""Options for adding a catalogue term."""

from __future__ import annotations

from typing import ClassVar

from attrs import define, field

from ensembl.vep.common import ConsequenceImpact
from ensembl.vep.config.option_spec import option
from ensembl.vep.config.sections.base import ConfigSection

__all__ = ["CatalogueAddConfig"]


def _impact_names() -> tuple[str, ...]:
    """Return the impact buckets a term may be assigned."""
    return tuple(impact.value for impact in ConsequenceImpact)


@define(frozen=True, kw_only=True, auto_attribs=True)
class CatalogueAddConfig(ConfigSection):
    """The term being added.

    Attributes:
        term (str | None): The SO term to add.
        accession (str | None): Its SO accession.
        rank (int | None): 1-based rank to insert at, or None to append after the highest.
        tier (int): Tier grouping.
        impact (str | None): Impact bucket.
        so_parent (tuple[str, ...]): SO parent accessions, used to derive ``subsumes``.
        subsumes (tuple[str, ...]): Extra terms to subsume, on top of the derived ones.
        no_derive_subsumes (bool): Skip deriving from ``so_parent``.
        when (str | None): The ``when`` payload as JSON. Never inferred: a wrong formula is the one
            error conformance cannot catch, since every annotator would agree on it.
        allow_missing_when (bool): Write the entry as a placeholder rather than refusing.
        notes (str | None): Free-text notes.
    """

    section_name: ClassVar[str] = "catalogue_add"

    term: str | None = field(default=None, metadata=option("the SO term to add"))
    accession: str | None = field(default=None, metadata=option("SO accession", metavar="SO:XXXXXXX"))
    rank: int | None = field(default=None, metadata=option("rank to insert at", metavar="N"))
    tier: int = field(default=3, metadata=option("tier grouping"))
    impact: str | None = field(default=None, metadata=option("impact bucket", choices=_impact_names))
    so_parent: tuple[str, ...] = field(
        default=(), metadata=option("SO parent accessions, used to derive subsumes")
    )
    subsumes: tuple[str, ...] = field(default=(), metadata=option("extra terms to subsume"))
    no_derive_subsumes: bool = field(default=False, metadata=option("do not derive subsumes from so_parent"))
    when: str | None = field(default=None, metadata=option("the when payload as JSON", metavar="JSON"))
    allow_missing_when: bool = field(
        default=False, metadata=option("accept an entry with no when, as a placeholder")
    )
    notes: str | None = field(default=None, metadata=option("free-text notes"))
