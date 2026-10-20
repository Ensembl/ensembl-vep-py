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
"""Trace events and recording state for consequence annotations.

This module has no ``ensembl.vep`` imports because facts, contexts, and caches depend on its event
types. Higher-level reconstruction belongs in ``trace.reconstruct``.
"""

from __future__ import annotations

from collections.abc import Iterable
from enum import Enum

from attrs import define, field, frozen

__all__ = ["FactEvent", "FactOrigin", "PruneEvent", "RuleEvent", "TraceRecorder"]


class FactOrigin(Enum):
    """Identify how a fact value became known.

    Attributes:
        COMPUTED (FactOrigin): Computed by the registered fact function.
        PREFILLED (FactOrigin): Supplied by a zone map without computation.
        CACHED (FactOrigin): Supplied by the shared-codon cache.
        FORCED_FALSE (FactOrigin): Proved false through ``requires`` propagation.
    """

    COMPUTED = "computed"
    PREFILLED = "prefilled by zone map"
    CACHED = "codon cache hit"
    FORCED_FALSE = "forced False by requires"


@frozen(kw_only=True, auto_attribs=True)
class FactEvent:
    """Record a fact value and how it became known.

    Attributes:
        name (str): The fact's name.
        value (object): The fact's value.
        cost (int): Cost charged for a computed fact or avoided for a forced-false fact; zero for
            values supplied by a cache.
        origin (FactOrigin): How the value became known.
        reads (int): Number of ``get_fact`` requests for this fact.
        forced_by (tuple[str, ...]): False ``requires`` atoms that forced the value, when relevant.
    """

    name: str
    value: object
    cost: int
    origin: FactOrigin
    reads: int = field(default=0)
    forced_by: tuple[str, ...] = field(default=())


@frozen(kw_only=True, auto_attribs=True)
class RuleEvent:
    """Record a rule verdict for one overlap.

    Attributes:
        term (str): The SO term this rule decides.
        rank (int): The rule's severity rank.
        fired (bool): True if the rule's formula evaluated True, before subsumption pruning.
        settled_early (bool): Whether a guard or partial evaluation settled the rule.
        reason (str): Formula or partial-evaluation reason for the verdict.
    """

    term: str
    rank: int
    fired: bool
    settled_early: bool
    reason: str


@frozen(kw_only=True, auto_attribs=True)
class PruneEvent:
    """Record one subsumption decision.

    Attributes:
        dropped_term (str): Term removed from the emitted terms.
        subsumed_by (str): Firing term that subsumed ``dropped_term``.
    """

    dropped_term: str
    subsumed_by: str


@define(kw_only=True, auto_attribs=True)
class TraceRecorder:
    """Collect live state that cannot be reconstructed after annotation.

    Reconstruction derives values and rule outcomes from the completed fact context. The recorder
    retains only original read counts and cache origins, which replay cannot recover.

    Attributes:
        read_counts (dict[str, int]): Number of ``get_fact`` requests by fact name.
        served (dict[str, FactOrigin]): Cache origin by fact name.
    """

    read_counts: dict[str, int] = field(factory=dict)
    served: dict[str, FactOrigin] = field(factory=dict)

    def note_served(self, names: Iterable[str], origin: FactOrigin) -> None:
        """Record facts supplied by a cache without computation.

        Args:
            names (Iterable[str]): Newly memoised fact names supplied by the cache.
            origin (FactOrigin): Cache origin for the facts.
        """
        for name in names:
            self.served[name] = origin
