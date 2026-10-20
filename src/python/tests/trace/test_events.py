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
"""Tests for the trace value types and the live-run recorder."""

from __future__ import annotations

from ensembl.vep.trace.events import FactEvent, FactOrigin, PruneEvent, RuleEvent, TraceRecorder


class TestTraceRecorder:
    """What the annotation path gathers, and what it deliberately does not."""

    def test_starts_empty(self) -> None:
        """A fresh recorder has counted nothing and been served nothing."""
        recorder = TraceRecorder()
        assert recorder.read_counts == {}
        assert recorder.served == {}

    def test_note_served_attributes_every_name_to_one_origin(self) -> None:
        """note_served records the cache each name came from."""
        recorder = TraceRecorder()
        recorder.note_served(["overlap_exon", "overlap_feature"], FactOrigin.PREFILLED)
        recorder.note_served(["_codon_alleles"], FactOrigin.CACHED)
        assert recorder.served == {
            "overlap_exon": FactOrigin.PREFILLED,
            "overlap_feature": FactOrigin.PREFILLED,
            "_codon_alleles": FactOrigin.CACHED,
        }

    def test_a_later_origin_wins_for_the_same_name(self) -> None:
        """The caller filters out already-memoised names, so a second call is a real re-serve."""
        recorder = TraceRecorder()
        recorder.note_served(["overlap_exon"], FactOrigin.PREFILLED)
        recorder.note_served(["overlap_exon"], FactOrigin.CACHED)
        assert recorder.served["overlap_exon"] is FactOrigin.CACHED


class TestEventDefaults:
    """The value types' defaults, which are what an unfilled field reads as."""

    def test_fact_event_defaults_to_unread_and_unforced(self) -> None:
        """reads and forced_by default to nothing rather than None."""
        event = FactEvent(name="overlap_cds", value=True, cost=20, origin=FactOrigin.COMPUTED)
        assert event.reads == 0
        assert event.forced_by == ()

    def test_events_are_frozen(self) -> None:
        """The event types are value objects, so a consumer cannot edit a trace it was handed."""
        for event in (
            FactEvent(name="a", value=True, cost=0, origin=FactOrigin.COMPUTED),
            RuleEvent(term="t", rank=1, fired=True, settled_early=False, reason=""),
            PruneEvent(dropped_term="a", subsumed_by="b"),
        ):
            assert getattr(type(event), "__setattr__", None) is not object.__setattr__

    def test_memo_hit_is_not_an_origin(self) -> None:
        """Read counts superseded it: a fact is recorded once, with how often it was read."""
        assert not hasattr(FactOrigin, "MEMO_HIT")
