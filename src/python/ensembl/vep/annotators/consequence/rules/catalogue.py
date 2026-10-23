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
"""Load the checked-in overlap consequence catalogue (``data/overlap_consequences.json``).

The catalogue is the shared source of term, accession, rank, tier, impact, formula, subsumption,
and ontology-parent metadata. ``when`` uses fact names and nested ``not``, ``all``, and ``any``
expressions. ``subsumes`` declares terms suppressed when both rules fire, while
``so_supported_parents`` records represented Sequence Ontology parents.

Fallback terms listed in ``EXCLUDED_TERMS`` intentionally have no ``when`` formula. All other
applicability conditions, including feature-kind checks, remain explicit in catalogue formulas.
"""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass, field
from functools import cache
from importlib.resources import files
from pathlib import Path
from typing import Any

from ensembl.vep.common import ConsequenceImpact

__all__ = [
    "ACCESSION_PATTERN",
    "CATALOGUE_RESOURCE",
    "ConsequenceDefinition",
    "DEFAULT_TERM",
    "ENTRY_KEY_ORDER",
    "EXCLUDED_TERMS",
    "SUPPORTED_SCHEMA_VERSIONS",
    "catalogue_digest",
    "catalogue_path",
    "dump_catalogue",
    "load_catalogue",
    "load_consequence_definitions",
]

CATALOGUE_RESOURCE = "data/overlap_consequences.json"
SUPPORTED_SCHEMA_VERSIONS = (5,)

ACCESSION_PATTERN = re.compile(r"^SO:\d{7}$")
"""Shape of a Sequence Ontology accession, e.g. ``SO:0001627``."""

ENTRY_KEY_ORDER = (
    "term",
    "accession",
    "rank",
    "tier",
    "impact",
    "so_supported_parents",
    "when",
    "subsumes",
    "notes",
)
"""Canonical per-entry key order, matching every existing schema-5 entry, so a freshly built entry
round-trips identically to a hand-edited one."""

DEFAULT_TERM = "intergenic_variant"
"""Fallback term applied by an annotator's caller when nothing else matched."""

EXCLUDED_TERMS = frozenset(
    {
        "intergenic_variant",
        "sequence_variant",
    }
)
"""Catalogue terms with no ``when`` documented."""


@dataclass(frozen=True)
class ConsequenceDefinition:
    """Represent one parsed entry from the overlap consequence catalogue.

    Attributes:
        term (str): The SO term, e.g. "missense_variant".
        accession (str): The SO accession, e.g. "SO:0001583".
        rank (int): Severity rank; 1 is most severe.
        tier (int): Tier grouping used by VEP's own evaluation order.
        impact (ConsequenceImpact | None): The IMPACT rating.
        when (object | None): A JSON boolean expression over fact atom names or None if this term
            is on EXCLUDED_TERMS.
        subsumes (tuple[str, ...]): Term names this consequence suppresses when when evaluated as true.
        so_supported_parents (tuple[str, ...]): SO ontology parent accessions that are also
            represented in this catalogue.
        notes (str | None): Free-text notes, if any.
    """

    term: str
    accession: str
    rank: int
    tier: int
    impact: ConsequenceImpact | None
    when: object | None
    subsumes: tuple[str, ...] = field(default_factory=tuple)
    so_supported_parents: tuple[str, ...] = field(default_factory=tuple)
    notes: str | None = None


def catalogue_path() -> Path:
    """Return the filesystem path of the packaged catalogue.

    Resolved as a package resource rather than relative to this file, so it survives being
    installed. Writing to it is refused unless the caller asks for it explicitly -- see
    ensembl.vep.annotators.consequence.rules.editing.
    """
    resource = files("ensembl.vep.annotators.consequence").joinpath(CATALOGUE_RESOURCE)
    return Path(str(resource))


def load_catalogue(path: Path | None = None) -> dict[str, Any]:
    """Return a catalogue's raw JSON payload, for editing rather than evaluation.

    ``load_consequence_definitions`` is the counterpart that parses the payload into rules; this
    keeps every key, including any the schema has that the parser ignores.

    Args:
        path (Path | None): The catalogue to read, or None for the packaged one.

    Returns:
        dict[str, Any]: The payload as parsed from JSON.

    Raises:
        ValueError: If the payload's ``schema_version`` is unsupported.
    """
    source = path or catalogue_path()
    payload = json.loads(source.read_text(encoding="utf-8"))
    schema_version = payload.get("schema_version")
    if schema_version not in SUPPORTED_SCHEMA_VERSIONS:
        raise ValueError(
            f"Unsupported overlap consequence catalogue schema version: {schema_version!r} "
            f"(supported: {SUPPORTED_SCHEMA_VERSIONS})"
        )
    return payload


def dump_catalogue(payload: dict[str, Any]) -> str:
    """Return a catalogue payload as JSON text, in the checked-in formatting."""
    return json.dumps(payload, indent=4) + "\n"


@cache
def load_consequence_definitions() -> tuple[ConsequenceDefinition, ...]:
    """Load and parse the package consequence catalogue.

    Returns:
        tuple[ConsequenceDefinition, ...]: Every catalogue entry, sorted by ``(rank, term)``.

    Raises:
        ValueError: If the catalogue's ``schema_version`` is unsupported, or a term appears more
            than once.
    """
    resource = files("ensembl.vep.annotators.consequence").joinpath(CATALOGUE_RESOURCE)
    payload = json.loads(resource.read_text(encoding="utf-8"))
    schema_version = payload.get("schema_version")
    if schema_version not in SUPPORTED_SCHEMA_VERSIONS:
        raise ValueError(
            f"Unsupported overlap consequence catalogue schema version: {schema_version!r} "
            f"(supported: {SUPPORTED_SCHEMA_VERSIONS})"
        )

    definitions: list[ConsequenceDefinition] = []
    seen_terms: set[str] = set()

    for entry in payload["consequences"]:
        term = entry["term"]
        if term in seen_terms:
            raise ValueError(f"Duplicate consequence catalogue term: {term}")
        seen_terms.add(term)

        definitions.append(
            ConsequenceDefinition(
                term=term,
                accession=entry["accession"],
                rank=int(entry.get("rank", 0)),
                tier=int(entry.get("tier", 0)),
                impact=ConsequenceImpact(entry["impact"]) if entry.get("impact") else None,
                when=entry.get("when"),
                subsumes=tuple(entry.get("subsumes", ())),
                so_supported_parents=tuple(entry.get("so_supported_parents", ())),
                notes=entry.get("notes"),
            )
        )

    return tuple(sorted(definitions, key=lambda definition: (definition.rank, definition.term)))


def catalogue_digest(path: Path | None = None) -> str:
    """Return a short digest of a catalogue file's bytes, the packaged catalogue by default.

    What a frozen expectation records beside its terms, so a later mismatch can be told apart:
    the catalogue changed under it, or the annotator did.
    """
    source = path if path is not None else catalogue_path()
    return hashlib.sha256(source.read_bytes()).hexdigest()[:16]
