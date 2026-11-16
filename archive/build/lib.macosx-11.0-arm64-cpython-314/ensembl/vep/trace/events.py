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
"""The vocabulary a derivation trace is written in, and the one thing the live run must gather.

Deliberately dependency-free: ``annotators/consequence/facts/base.py`` and
``annotators/consequence/context/fact.py`` both import ``FactOrigin`` and ``TraceRecorder`` from
here, so anything this module imported from ``ensembl.vep`` would close a cycle back through
``annotators/__init__``.
"""

from __future__ import annotations

from collections.abc import Iterable
from enum import Enum

from attrs import define, field, frozen

__all__ = ["FactEvent", "FactOrigin", "PruneEvent", "RuleEvent", "TraceRecorder"]


class FactOrigin(Enum):
    """How a fact's value became known to a FactContext.

    Attributes:
        COMPUTED: ``get_fact`` called the fact's registered ``compute`` function and charged its
            declared cost.
        PREFILLED: Injected uncharged by a zone map (``ZoneMapCache.prefill``) via
            ``FactContext.prefill``, without ever calling the fact's ``compute`` function.
        CACHED: Injected uncharged by the shared-codon cache (``CodonCache.prefill``) on a hit.
        FORCED_FALSE: Proved False by ``requires`` propagation without ever being computed, so it
            never reached the memo table at all.
    """

    COMPUTED = "computed"
    PREFILLED = "prefilled by zone map"
    CACHED = "codon cache hit"
    FORCED_FALSE = "forced False"


@frozen(kw_only=True, auto_attribs=True)
class FactEvent:
    """One fact's value, how it became known, and how often it was read.

    Attributes:
        name (str): The fact's name (see ``consequence/facts/base.py``'s ``REGISTRY``).
        value (object): The fact's value.
        cost (int): The declared cost this fact charged, if COMPUTED; the cost it *would* have
            charged (cost avoided) if FORCED_FALSE; 0 for the uncharged origins.
        origin (FactOrigin): How the value became known.
        reads (int): How many times ``get_fact`` was asked for this fact. Zero for a fact that was
            forced False, or served by a cache nothing went on to read.
        forced_by (tuple[str, ...]): For FORCED_FALSE, the ``requires`` atoms found False. Empty
            for every other origin.
    """

    name: str
    value: object
    cost: int
    origin: FactOrigin
    reads: int = field(default=0)
    forced_by: tuple[str, ...] = field(default=())


@frozen(kw_only=True, auto_attribs=True)
class RuleEvent:
    """One rule's verdict for one overlap, and why it was reached.

    Attributes:
        term (str): The SO term this rule decides.
        rank (int): The rule's severity rank (1 is most severe).
        fired (bool): True if the rule's formula evaluated True, before subsumption pruning.
        settled_early (bool): True if the verdict was reached without evaluating the rule's full
            formula -- either a cheap-guard rejection or a partial (Kleene) evaluation settled by
            cheap and forced-False facts alone.
        reason (str): The settling atoms and their values when settled_early, the evaluated
            formula otherwise, or that the rule fell through to full evaluation.
    """

    term: str
    rank: int
    fired: bool
    settled_early: bool
    reason: str


@frozen(kw_only=True, auto_attribs=True)
class PruneEvent:
    """One subsumption decision: a firing term dropping another firing term.

    Attributes:
        dropped_term (str): The term removed from the emitted list.
        subsumed_by (str): The firing term whose ``subsumes`` declared ``dropped_term``.
    """

    dropped_term: str
    subsumed_by: str


@define(kw_only=True, auto_attribs=True)
class TraceRecorder:
    """The only thing the annotation path writes to when tracing is on.

    Everything else in a trace is reconstructed afterwards from the memo table and the rule set
    (see ``trace/reconstruct.py``). These two are what reconstruction cannot recover: a replay sees
    only memo hits, so it cannot observe the original read pattern, and it cannot tell a value a
    cache served from one the fact layer computed, because the cache has since been written to.

    Both are cheap by construction. ``read_counts`` is one dict bump inside the single
    ``is not None`` branch ``get_fact`` already pays; ``served`` is written at most twice per
    overlap, by ``FactContext.prefill``.

    Attributes:
        read_counts (dict[str, int]): Fact name -> how many times ``get_fact`` was asked for it.
        served (dict[str, FactOrigin]): Fact name -> the uncharged origin that supplied it, for
            the facts a cache put into the memo rather than the fact layer computing them.
    """

    read_counts: dict[str, int] = field(factory=dict)
    served: dict[str, FactOrigin] = field(factory=dict)

    def note_served(self, names: Iterable[str], origin: FactOrigin) -> None:
        """Record that a cache supplied these facts uncharged.

        Args:
            names (Iterable[str]): The facts the cache put into the memo. The caller passes only
                the names that were *not* already memoised, so a cache miss that computed the
                value the normal way stays COMPUTED.
            origin (FactOrigin): Which cache supplied them.
        """
        for name in names:
            self.served[name] = origin
