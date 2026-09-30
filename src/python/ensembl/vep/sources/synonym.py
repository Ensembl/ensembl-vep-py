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
"""Sequence-region synonyms, and resolving a name against one source's own spellings."""

from __future__ import annotations

import logging
from collections.abc import Iterable, Mapping
from pathlib import Path

from attrs import define, field

__all__ = ["ChromosomeSynonyms", "SequenceRegionResolver"]

logger = logging.getLogger(__name__)


@define(kw_only=True, auto_attribs=True)
class ChromosomeSynonyms:
    """Sequence-region names that name the same region, held once for the whole run.

    Perl VEP keeps one synonym map per run and asks it for every source (``BaseVEP::
    chromosome_synonyms``); a name is not a synonym *of a source*, it is a synonym full stop, and
    which spelling a given source wants is decided by that source's own valid names. Keeping one
    map is what lets ``22 -> NC_000022.11`` serve the FASTA and the transcript file at once, in
    whichever direction each of them needs.

    Attributes:
        by_name (Mapping[str, frozenset[str]]): Every name, mapped to the other names for it.
    """

    by_name: Mapping[str, frozenset[str]] = field(factory=dict)

    @classmethod
    def from_file(cls, path: Path | str) -> ChromosomeSynonyms:
        """Read a synonyms file in Perl VEP's format.

        Each line lists whitespace-separated names for one region. The first name is the reference
        spelling and every later name is a synonym of it, recorded in both directions -- as Perl
        does, so two synonyms of one reference name are not made synonyms of each other.

        Args:
            path (Path | str): The file to read.

        Returns:
            ChromosomeSynonyms: The map the file describes.

        Raises:
            OSError: If the file cannot be read. Allowed to propagate: its message names the path,
                and a synonyms file the caller asked for and that is not there is not a run to
                continue silently.
        """

        pairs: dict[str, set[str]] = {}
        with open(path, encoding="utf-8") as synonyms_file:
            for line in synonyms_file:
                names = line.split()
                if len(names) < 2:
                    continue
                reference_name, synonyms = names[0], names[1:]
                for synonym in synonyms:
                    pairs.setdefault(reference_name, set()).add(synonym)
                    pairs.setdefault(synonym, set()).add(reference_name)
        logger.info("read %d sequence-region synonyms from %s", len(pairs), path)
        return cls(by_name={name: frozenset(synonyms) for name, synonyms in pairs.items()})

    def synonyms_of(self, name: str) -> frozenset[str]:
        """Return the other names for a region.

        Args:
            name (str): The name to look up.

        Returns:
            frozenset[str]: Its synonyms, empty when the map has none.
        """

        return self.by_name.get(name, frozenset())


@define(kw_only=True, auto_attribs=True)
class SequenceRegionResolver:
    """Resolve uploaded sequence-region names against one source.

    Perl's ``BaseVEP::get_source_chr_name`` in our shape: an exact source name wins, then a synonym
    the source declares, then the conventional case-insensitive ``chr`` prefix adjustment, and an
    unresolved name keeps its uploaded spelling. Resolution is memoised per instance, which is why
    each source gets its own resolver over the shared synonym map.

    Attributes:
        source_names (frozenset[str]): Names this source declares.
        synonyms (ChromosomeSynonyms): The run's synonym map, shared by every source.
    """

    source_names: frozenset[str]
    synonyms: ChromosomeSynonyms = field(factory=ChromosomeSynonyms)
    _resolved: dict[str, str] = field(factory=dict, init=False, repr=False)

    def __attrs_post_init__(self) -> None:
        """Reject source names that cannot be matched case-insensitively without ambiguity."""

        folded_names: dict[str, str] = {}
        for source_name in self.source_names:
            folded_name = source_name.casefold()
            previous_name = folded_names.get(folded_name)
            if previous_name is not None and previous_name != source_name:
                raise ValueError(f"source names {previous_name!r} and {source_name!r} differ only by case")
            folded_names[folded_name] = source_name

    def resolve(self, uploaded_name: str) -> str:
        """Return the source sequence-region name for an uploaded name."""

        cached = self._resolved.get(uploaded_name)
        if cached is not None:
            return cached

        resolved = self._resolve(uploaded_name)
        self._resolved[uploaded_name] = resolved
        return resolved

    def _resolve(self, uploaded_name: str) -> str:
        """Resolve one name without consulting the per-instance cache."""

        if uploaded_name in self.source_names:
            return uploaded_name
        synonym = self._declared_synonym(self.synonyms.synonyms_of(uploaded_name))
        if synonym is not None:
            return synonym
        chr_adjusted = self._chr_adjusted_name(uploaded_name)
        if chr_adjusted is not None:
            return chr_adjusted
        return uploaded_name

    def _declared_synonym(self, synonyms: Iterable[str]) -> str | None:
        """Return the synonym this source declares, lowest name first when several do.

        Perl walks an unordered hash here and takes whichever synonym comes out first. Sorting
        makes the same input give the same answer on every run, which output diffs depend on.
        """

        for synonym in sorted(synonyms):
            if synonym in self.source_names:
                return synonym
        return None

    def _chr_adjusted_name(self, uploaded_name: str) -> str | None:
        """Return a valid case-insensitive ``chr`` prefix adjustment."""

        if uploaded_name[:3].lower() == "chr":
            candidate = uploaded_name[3:]
            return self._source_name(candidate)
        return self._source_name(f"chr{uploaded_name}")

    def _source_name(self, candidate: str) -> str | None:
        """Return the valid source spelling of a candidate name."""

        if candidate in self.source_names:
            return candidate
        candidate_folded = candidate.casefold()
        for source_name in self.source_names:
            if source_name.casefold() == candidate_folded:
                return source_name
        return None
