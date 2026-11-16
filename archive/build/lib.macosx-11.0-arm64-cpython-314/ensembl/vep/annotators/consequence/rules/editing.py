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
"""Editing the consequence catalogue: adding and removing terms.

A ``when`` is never guessed. ``suggest_when_atoms`` offers candidates for a human to choose from,
and an entry written without one is marked as needing review rather than given a formula that looks
plausible -- a wrong formula is the one error the conformance suite cannot catch, since every
annotator would agree on the same wrong answer.
"""

from __future__ import annotations

import difflib
import json
from collections.abc import Sequence
from typing import Any

from ensembl.vep.annotators.consequence.facts.base import atom_names
from ensembl.vep.annotators.consequence.rules.catalogue import ACCESSION_PATTERN, ENTRY_KEY_ORDER
from ensembl.vep.annotators.consequence.rules.expressions import expr_from_json

__all__ = [
    "KEYWORD_TO_ATOM",
    "accession_index",
    "apply_add",
    "apply_remove",
    "build_entry",
    "derive_subsumes",
    "insert_at_rank",
    "next_rank",
    "remove_term",
    "resolve_when",
    "suggest_when_atoms",
]

# Ordered longest/most-specific keyword first, so e.g. "splice_donor_5th_base" is tried before the
# more general "splice_donor" for a term name that happens to contain both substrings. Built by
# reading which atom each existing catalogue term's own `when` centres on (see
# data/overlap_consequences.json) -- a best-effort seed list, not an exhaustive mapping.
KEYWORD_TO_ATOM: tuple[tuple[str, str], ...] = (
    ("splice_donor_5th_base", "overlap_donor_5th_base"),
    ("splice_donor_region", "overlap_donor_region"),
    ("splice_donor", "overlap_donor_site"),
    ("splice_acceptor", "overlap_acceptor_site"),
    ("polypyrimidine", "overlap_polypyrimidine_tract"),
    ("splice_region", "overlap_splice_region"),
    ("splice", "overlap_splice_region"),
    ("stop_gained", "alt_peptide_has_ter"),
    ("stop_retained", "cds_stop_retained"),
    ("stop_lost", "cds_stop_lost"),
    ("start_retained", "cds_start_retained"),
    ("start_lost", "cds_start_lost"),
    ("incomplete_terminal_codon", "partial_last_codon"),
    ("frameshift", "cds_length_change_inframe"),
    ("missense", "peptide_changed"),
    ("synonymous", "peptide_changed"),
    ("inframe_insertion", "is_insertion"),
    ("inframe_deletion", "is_deletion"),
    ("insertion", "is_insertion"),
    ("deletion", "is_deletion"),
    ("elongation", "is_structural_insertion"),
    ("truncation", "is_structural_deletion"),
    ("amplification", "is_structural_insertion"),
    ("5_prime_utr", "overlap_five_prime_utr"),
    ("three_prime_utr", "overlap_three_prime_utr"),
    ("3_prime_utr", "overlap_three_prime_utr"),
    ("five_prime_utr", "overlap_five_prime_utr"),
    ("intron", "overlap_intron"),
    ("exon", "overlap_exon"),
    ("upstream", "is_upstream"),
    ("downstream", "is_downstream"),
    ("nmd", "biotype_nmd"),
    ("non_coding", "biotype_non_coding"),
    ("coding_transcript", "biotype_protein_coding"),
    ("coding", "overlap_cds"),
    ("protein_altering", "peptide_changed"),
)
"""Curated keyword -> fact atom suggestions. Deliberately never applied automatically -- see
resolve_when and the module docstring's explanation of why a `when` is never guessed."""


def accession_index(entries: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """Build an accession -> entry lookup.

    Args:
        entries (list[dict[str, Any]]): The catalogue's ``consequences`` list.

    Returns:
        dict[str, dict[str, Any]]: Every entry, keyed by its SO accession.
    """
    return {entry["accession"]: entry for entry in entries}


def derive_subsumes(
    so_parent_accessions: Sequence[str], entries_by_accession: dict[str, dict[str, Any]]
) -> tuple[tuple[str, ...], list[str]]:
    """Walk ``so_supported_parents`` transitively to propose a ``subsumes`` set.

    Direction: if the new term A *is_a* an existing catalogue term B (B is one of A's SO parents),
    A is the more specific term, so A *subsumes* B -- suppressing B whenever A also fires. Walking
    stays inside the catalogue: an SO parent accession not represented here cannot be proposed
    (there is nothing to subsume), and a present-but-``when``-less ancestor (an EXCLUDED_TERMS-style
    term, see consequence/rules/catalogue.py) is skipped from the proposal itself -- it can never appear in a
    compiled FactRuleSet, so declaring a rule subsumes it would fail validation -- but is still
    noted, and its own further-up parents are still walked past it.

    Args:
        so_parent_accessions (list[str]): The new term's immediate SO parent accessions.
        entries_by_accession (dict[str, dict[str, Any]]): accession -> catalogue entry lookup.

    Returns:
        tuple[tuple[str, ...], list[str]]: The proposed, sorted ``subsumes`` term tuple, and a list
            of human-readable notes about the walk (skipped/unknown/excluded ancestors).
    """
    subsumed_terms: set[str] = set()
    notes: list[str] = []
    visited_accessions: set[str] = set()
    frontier = list(so_parent_accessions)

    while frontier:
        accession = frontier.pop()
        if accession in visited_accessions:
            continue
        visited_accessions.add(accession)

        entry = entries_by_accession.get(accession)
        if entry is None:
            notes.append(f"SO parent {accession} is not present in the catalogue -- skipped")
            continue

        if entry.get("when") is not None:
            subsumed_terms.add(entry["term"])
        else:
            notes.append(
                f"SO parent {entry['term']} ({accession}) has no 'when' in the catalogue -- not "
                "proposed as a subsumes target (it is never part of a compiled rule set), but its "
                "own SO parents are still walked"
            )
        frontier.extend(entry.get("so_supported_parents", ()))

    return tuple(sorted(subsumed_terms)), notes


def suggest_when_atoms(term: str) -> list[str]:
    """Suggest candidate fact atoms for a term name, by keyword and fuzzy matching.

    Never used to auto-populate a ``when`` -- purely a suggestion for a human author. See the
    module docstring for why a formula is never inferred automatically.

    Args:
        term (str): The SO term name being added, e.g. "stop_gained".

    Returns:
        list[str]: Candidate atom names, most-confident first (exact atom-name match, then
            curated keyword matches, then fuzzy name matches), deduplicated.
    """
    known_atoms = atom_names()
    lowered = term.lower()
    suggestions: list[str] = []

    if term in known_atoms:
        suggestions.append(term)
    for keyword, atom in KEYWORD_TO_ATOM:
        if keyword in lowered and atom in known_atoms and atom not in suggestions:
            suggestions.append(atom)
    for atom in difflib.get_close_matches(lowered, known_atoms, n=3, cutoff=0.5):
        if atom not in suggestions:
            suggestions.append(atom)

    return suggestions


def build_entry(
    *,
    term: str,
    accession: str,
    rank: int,
    tier: int,
    impact: str,
    so_supported_parents: Sequence[str],
    when: object | None,
    subsumes: tuple[str, ...],
    notes: str | None,
) -> dict[str, Any]:
    """Build a new catalogue entry dict in the canonical, round-trippable key order.

    Args:
        term (str): The SO term name.
        accession (str): The SO accession.
        rank (int): Severity rank.
        tier (int): Tier grouping.
        impact (str): One of ConsequenceImpact's values.
        so_supported_parents (list[str]): SO parent accessions also represented in the catalogue.
        when (object | None): A compiled-and-validated ``when`` payload, or None for a needs-review
            placeholder (see resolve_when).
        subsumes (tuple[str, ...]): Term names this entry suppresses when it also fires.
        notes (str | None): Free-text notes, if any.

    Returns:
        dict[str, Any]: The new entry, with keys present/ordered exactly as an existing schema-4
            catalogue entry would have them (``when``/``subsumes``/``notes`` omitted, not
            null/empty, when absent).
    """
    raw: dict[str, Any] = {
        "term": term,
        "accession": accession,
        "rank": rank,
        "tier": tier,
        "impact": impact,
        "so_supported_parents": list(so_supported_parents),
    }
    if when is not None:
        raw["when"] = when
    if subsumes:
        raw["subsumes"] = list(subsumes)
    if notes:
        raw["notes"] = notes
    return {key: raw[key] for key in ENTRY_KEY_ORDER if key in raw}


def next_rank(entries: list[dict[str, Any]]) -> int:
    """Return the rank one past the current highest rank.

    Args:
        entries (list[dict[str, Any]]): The catalogue's ``consequences`` list.

    Returns:
        int: ``max(rank for every entry) + 1``, or 1 if ``entries`` is empty.
    """
    return max((int(entry.get("rank", 0)) for entry in entries), default=0) + 1


def insert_at_rank(
    entries: list[dict[str, Any]], new_entry: dict[str, Any], requested_rank: int | None
) -> list[dict[str, Any]]:
    """Insert ``new_entry`` at ``requested_rank``, shifting later ranks down to make room.

    Args:
        entries (list[dict[str, Any]]): The catalogue's ``consequences`` list (not mutated).
        new_entry (dict[str, Any]): The entry to insert; its ``rank`` key is set by this function.
        requested_rank (int | None): The 1-based rank to insert at, or None to append after the
            current highest rank.

    Returns:
        list[dict[str, Any]]: A new list, sorted by rank, with ``new_entry`` inserted and every
            colliding/later rank shifted up by one.

    Raises:
        ValueError: If ``requested_rank`` is outside ``1..(current highest rank + 1)``.
    """
    if requested_rank is None:
        new_entry["rank"] = next_rank(entries)
        return sorted([*entries, new_entry], key=lambda entry: int(entry["rank"]))

    highest_rank = max((int(entry.get("rank", 0)) for entry in entries), default=0)
    if not 1 <= requested_rank <= highest_rank + 1:
        raise ValueError(f"rank {requested_rank} out of range (1..{highest_rank + 1})")

    shifted: list[dict[str, Any]] = []
    for entry in entries:
        if int(entry.get("rank", 0)) >= requested_rank:
            entry = dict(entry)
            entry["rank"] = int(entry["rank"]) + 1
        shifted.append(entry)
    new_entry["rank"] = requested_rank
    return sorted([*shifted, new_entry], key=lambda entry: int(entry["rank"]))


def remove_term(entries: list[dict[str, Any]], term: str) -> tuple[list[dict[str, Any]], list[str]]:
    """Remove ``term``, clean up every other entry's ``subsumes``, and close the rank gap.

    Args:
        entries (list[dict[str, Any]]): The catalogue's ``consequences`` list (not mutated).
        term (str): The term to remove. Caller has already verified it exists.

    Returns:
        tuple[list[dict[str, Any]], list[str]]: The remaining entries (rank-sorted, gap closed),
            and the sorted list of terms whose ``subsumes`` referenced ``term`` and were cleaned.
    """
    removed_entry = next(entry for entry in entries if entry["term"] == term)
    removed_rank = int(removed_entry.get("rank", 0))

    cleaned_from: list[str] = []
    remaining: list[dict[str, Any]] = []
    for entry in entries:
        if entry["term"] == term:
            continue
        entry = dict(entry)
        subsumes = entry.get("subsumes")
        if subsumes and term in subsumes:
            entry["subsumes"] = [subsumed for subsumed in subsumes if subsumed != term]
            if not entry["subsumes"]:
                del entry["subsumes"]
            cleaned_from.append(entry["term"])
        if int(entry.get("rank", 0)) > removed_rank:
            entry["rank"] = int(entry["rank"]) - 1
        remaining.append(entry)

    remaining.sort(key=lambda entry: int(entry["rank"]))
    return remaining, sorted(cleaned_from)


def resolve_when(
    term: str, when_argument: str | None, known_atoms: list[str]
) -> tuple[object | None, list[str], bool]:
    """Resolve ``--when`` into a validated payload, or a needs-review placeholder.

    Args:
        term (str): The term being added (used to generate suggestions when ``when_argument`` is
            not given).
        when_argument (str | None): The raw ``--when`` JSON text, or None if not given.
        known_atoms (list[str]): Every registered fact atom name (ensembl.vep.annotators.consequence.
            facts.base.atom_names()), used both to validate an explicit ``when`` and to generate
            suggestions for a missing one.

    Returns:
        tuple[object | None, list[str], bool]: ``(when_payload, messages, needs_review)``.
            ``when_payload`` is a compiled-valid JSON payload if ``when_argument`` was given, else
            None. ``needs_review`` is True whenever ``when_payload`` is None.

    Raises:
        ValueError: If ``when_argument`` is not valid JSON, does not compile to an expression, or
            references an atom name that is not registered.
    """
    if when_argument is None:
        suggestions = suggest_when_atoms(term)
        note = (
            f"no --when given; candidate atom(s) worth reviewing: {suggestions}"
            if suggestions
            else "no --when given; no confident candidate atom found by keyword/fuzzy matching"
        )
        return None, [note], True

    try:
        when_payload = json.loads(when_argument)
    except json.JSONDecodeError as decode_error:
        raise ValueError(f"when is not valid JSON: {decode_error}") from decode_error

    try:
        expression = expr_from_json(when_payload)
    except ValueError as compile_error:
        raise ValueError(f"when does not compile to an expression: {compile_error}") from compile_error

    unknown_atoms = expression.atoms() - set(known_atoms)
    if unknown_atoms:
        raise ValueError(f"when references unknown fact atom(s): {sorted(unknown_atoms)}")

    return when_payload, [f"using explicitly-authored --when: {json.dumps(when_payload)}"], False


def apply_add(
    payload: dict[str, Any],
    *,
    term: str,
    accession: str | None,
    impact: str | None,
    rank: int | None = None,
    tier: int = 3,
    so_parent: tuple[str, ...] = (),
    subsumes: tuple[str, ...] = (),
    no_derive_subsumes: bool = False,
    when: str | None = None,
    allow_missing_when: bool = False,
    notes: str | None = None,
) -> tuple[dict[str, Any], list[str]] | None:
    """Add a term to a catalogue payload, without writing anything.

    Args:
        payload (dict[str, Any]): The payload to add to, mutated in place.
        term (str): The SO term to add.
        accession (str | None): Its SO accession.
        impact (str | None): Its impact bucket.
        rank (int | None): 1-based rank to insert at, or None to append after the highest.
        tier (int): Tier grouping.
        so_parent (tuple[str, ...]): SO parent accessions, used to derive ``subsumes``.
        subsumes (tuple[str, ...]): Extra terms to subsume, on top of the derived ones.
        no_derive_subsumes (bool): Skip deriving from ``so_parent``.
        when (str | None): The ``when`` payload as JSON, or None to write a placeholder.
        allow_missing_when (bool): Accept the placeholder rather than refusing.
        notes (str | None): Free-text notes.

    Returns:
        tuple[dict[str, Any], list[str]] | None: ``(payload, messages)``, or None if the addition
            was rejected because no ``when`` was given and it was not explicitly allowed.

    Raises:
        ValueError: On a usage error -- missing accession or impact, duplicate term or accession,
            unknown subsumes target, or a malformed accession.
    """
    entries: list[dict[str, Any]] = payload["consequences"]

    if accession is None or impact is None:
        raise ValueError("adding a term requires an accession and an impact")
    if not ACCESSION_PATTERN.match(accession):
        raise ValueError(f"accession {accession!r} does not look like SO:NNNNNNN")
    if any(entry["term"] == term for entry in entries):
        raise ValueError(f"catalogue already has a term named {term!r}")
    if any(entry["accession"] == accession for entry in entries):
        raise ValueError(f"catalogue already has an entry with accession {accession!r}")

    known_terms = {entry["term"] for entry in entries}
    for target in subsumes:
        if target not in known_terms:
            raise ValueError(f"subsumes target {target!r} is not a catalogue term")

    derived: tuple[str, ...]
    derive_notes: list[str]
    if no_derive_subsumes:
        derived, derive_notes = (), []
    else:
        derived, derive_notes = derive_subsumes(so_parent, accession_index(entries))
    final_subsumes = tuple(dict.fromkeys((*derived, *subsumes)))

    when_payload, when_messages, needs_review = resolve_when(term, when, atom_names())

    messages = [
        f"deriving subsumes for {term!r} from {len(so_parent)} SO parent accession(s)",
        *derive_notes,
        f"derived subsumes: {list(derived)}",
        f"final subsumes (derived + explicit): {list(final_subsumes)}",
        *when_messages,
    ]

    if needs_review and not allow_missing_when:
        banner = "!" * 78
        messages.append(
            f"\n{banner}\n"
            "NEEDS REVIEW: no --when was given and none can be safely inferred.\n"
            f"  candidate fact atom(s) to consider: {suggest_when_atoms(term)}\n"
            "A boolean formula cannot be reliably inferred from a term name alone -- a wrong guess\n"
            "would pass conformance (every annotator agrees on the same wrong answer) and only show\n"
            "up as bad biology later. Author the `when` by hand and re-run with --when, or\n"
            "acknowledge the placeholder explicitly with --allow-missing-when.\n"
            f"{banner}"
        )
        for message in messages:
            print(message)
        return None

    new_entry = build_entry(
        term=term,
        accession=accession,
        rank=0,
        tier=tier,
        impact=impact,
        so_supported_parents=so_parent,
        when=when_payload,
        subsumes=final_subsumes,
        notes=notes,
    )
    payload["consequences"] = insert_at_rank(entries, new_entry, rank)
    return payload, messages


def apply_remove(payload: dict[str, Any], *, term: str) -> tuple[dict[str, Any], list[str]]:
    """Remove a term from a catalogue payload, without writing anything.

    Args:
        payload (dict[str, Any]): The payload to remove from, mutated in place.
        term (str): The SO term to remove.

    Returns:
        tuple[dict[str, Any], list[str]]: ``(payload, messages)``. Removal is never a needs-review
            case, so unlike ``apply_add`` this always succeeds when the term exists.

    Raises:
        ValueError: If no entry named ``term`` exists.
    """
    entries: list[dict[str, Any]] = payload["consequences"]
    if not any(entry["term"] == term for entry in entries):
        raise ValueError(f"no catalogue entry named {term!r}")

    remaining, cleaned_from = remove_term(entries, term)
    messages = [f"removed {term!r}"]
    if cleaned_from:
        messages.append(f"removed {term!r} from the 'subsumes' list of: {cleaned_from}")
    else:
        messages.append(f"{term!r} was not subsumed by any other entry -- no cleanup needed")

    payload["consequences"] = remaining
    return payload, messages
