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
"""The cross-payload feature pool contract, and the transcript pool that implements it.

``BinnedFeatureCache`` (the worker's per-payload transcript cache) is cleared before every new
payload, so a worker that revisits the same gene across many payloads -- a dense or targeted input,
or simply a long run of variants on one chromosome -- rebuilds the same ``Transcript`` from GXF
records every time. Building one is not free: grouping records, constructing ``Exon``/``Translation``
objects, and resolving the shared ``Gene``.

``Transcript`` objects are never mutated after ``GXFTranscriptSource`` builds them, so reusing the
same object across payloads is safe -- the only growth over a transcript's pooled lifetime is
whatever ``cached_property`` values get resolved on it (introns, spliced sequence, CDS, translated
sequence, and so on), which is why the pool is bounded by bytes rather than kept forever.

``TranscriptPool`` is keyed by ``(region, stable_id, start, end)`` -- the same fields
``BinnedFeatureCache`` already uses as its own dedup key -- and mirrors ``ZoneMapCache``'s
byte-budgeted LRU shape: build outside the pool, insert and evict under a lock.
"""

from __future__ import annotations

from collections import OrderedDict
from threading import RLock
from typing import Protocol, TypeVar, runtime_checkable

from attrs import define, field

from ensembl.vep.models import Feature
from ensembl.vep.models.genic import Transcript

__all__ = [
    "DEFAULT_TRANSCRIPT_POOL_BYTES",
    "FeatureKey",
    "FeaturePool",
    "TranscriptKey",
    "TranscriptPool",
]

FeatureT = TypeVar("FeatureT", bound=Feature)

FeatureKey = tuple[str, str, int, int]
"""``(region, stable_id, start, end)`` -- a built feature's pool identity.

Matches ``BinnedFeatureCache.add``'s own dedup key, so a feature is addressed the same way whether
the caller means "the cache for this payload" or "the pool across payloads".
"""


@runtime_checkable
class FeaturePool(Protocol[FeatureT]):
    """Define an interchangeable cache of built features, reused across payloads.

    A source builds features from its own records; a pool lets the next payload reuse those
    objects instead of rebuilding them. Implementations decide what to retain and when to evict,
    so a pool is a hint, never a source of truth: ``get`` may miss at any time and ``put`` may
    discard what it was given.

    Attributes:
        max_bytes (int): Byte budget for pooled features. ``<= 0`` disables the pool.
        hit_count (int): Number of ``get`` calls served from the pool so far.
        miss_count (int): Number of ``get`` calls that found nothing so far.
    """

    max_bytes: int
    hit_count: int
    miss_count: int

    def get(self, key: FeatureKey) -> FeatureT | None:
        """Return the pooled feature for a key, or None on a miss."""

        raise NotImplementedError

    def put(self, key: FeatureKey, feature: FeatureT) -> None:
        """Offer a freshly built feature to the pool."""

        raise NotImplementedError

    def clear(self) -> None:
        """Clear pooled features."""

        raise NotImplementedError


DEFAULT_TRANSCRIPT_POOL_BYTES = 128 * 1024 * 1024
"""Default byte budget for ``TranscriptPool``, shared by the config default and every call site."""

TranscriptKey = FeatureKey
"""``(region, stable_id, start, end)`` -- the ``FeatureKey`` spelling a transcript pool uses."""

_BYTES_PER_EXON = 600
"""Estimated per-exon overhead: the ``Exon`` object, its ``Slice`` reference and attrs bookkeeping.

Calibrated loosely, not measured with ``tracemalloc`` the way ``ZoneMap.estimated_bytes`` was --
see the module docstring on why an exact figure is not the point. It is one term of two so that a
transcript with many small exons and one with few huge exons are both charged something realistic.
"""

_BYTES_PER_GENOMIC_BASE = 0.5
"""Estimated per-base overhead for the transcript's own genomic span.

Stands in for the ``cached_property`` values a pooled transcript accumulates once fully read
(spliced sequence, CDS, translated sequence, introns, mapper state): each is at most on the order
of one byte per base, and several (translated sequence at a third the length; CDS restricted to
the coding portion; non-coding transcripts with no CDS/translation at all) are shorter, so half a
byte per base approximates the mix without enumerating every property.

Calibrated with ``tracemalloc`` over a real dense region (GRCh38 chr1:1-20,000,000, 8,021
transcripts, 70,786 exons): building alone grew the heap by 35.0 MB; reading ``spliced_seq``,
``translated_seq``, ``introns`` and ``cds`` on every transcript afterwards grew it to 102.2 MB.
At this value the estimate for that region is 155.1 MB, 1.52x the fully-read figure. A first
attempt at 3 B/base put the estimate at 719 MB for the unread population alone (roughly 20x too
high), capping the pool at a few hundred transcripts regardless of ``max_bytes``; plain 1 B/base
(the next value tried) still overshot fully-read at 267.8 MB, 2.62x.
"""


def _estimate_transcript_bytes(transcript: Transcript) -> int:
    """Estimate a built transcript's heap footprint, for LRU byte accounting.

    ``600 * len(spliced_exons) + 0.5 * (end - start)``: a per-exon term for the ``Exon`` objects
    themselves, plus a per-base term standing in for cached derived state (introns, spliced
    sequence, CDS, translated sequence, ...) once read. This is deliberately approximate -- an
    exact accounting would have to walk every ``cached_property`` a transcript might have resolved
    by the time it is evicted, which depends on what the workflow happened to read from it. The
    budget is a soft cap for the same reason: bound pool memory to roughly the right order of
    magnitude, not to an exact byte count. It also does not count memory a transcript keeps alive
    indirectly -- an evicted transcript can still be reachable (and so not actually freed) through
    a pooled sibling's ``gene.transcripts``.

    Args:
        transcript (Transcript): The transcript to size.

    Returns:
        int: The estimated size in bytes, clamped to be non-negative.
    """
    span = max(0, transcript.end - transcript.start)
    return _BYTES_PER_EXON * len(transcript.spliced_exons) + int(_BYTES_PER_GENOMIC_BASE * span)


@define(kw_only=True, auto_attribs=True, eq=False)
class TranscriptPool(FeaturePool[Transcript]):
    """A byte-budgeted, least-recently-used cache of built ``Transcript`` models.

    Reuse across payloads within one worker -- a payload's own ``BinnedFeatureCache`` is still
    cleared every time (``BufferedWorker._begin_payload``); this pool is what survives that clear.
    Construct one per worker and hand it to ``GXFTranscriptSource.transcript_pool``.

    Only one sequence region is held at a time: inserting a transcript from a different region than
    the one currently pooled drops every entry from the old region first. A sorted VCF input visits
    one region for a long run of payloads before moving to the next, so this keeps the common case
    (many payloads over the same region) cheap without holding memory for regions no longer being
    annotated. It is a size optimisation, not a correctness requirement -- a transcript is always
    addressed by its full key, so nothing breaks if the same region is revisited later; it is just
    a pool miss.

    ``eq=False`` keeps identity semantics: this is a shared, mutable instance, not a value object,
    so two pools with equal counters must not compare equal.

    Attributes:
        max_bytes (int): Byte budget for pooled transcripts. ``<= 0`` disables the pool: ``get``
            always misses and ``put`` never stores anything. Default is
            ``DEFAULT_TRANSCRIPT_POOL_BYTES``.
        hit_count (int): Number of ``get`` calls served from the pool so far.
        miss_count (int): Number of ``get`` calls that found nothing so far.
        eviction_count (int): Number of entries evicted for being over the byte budget so far.
        region_eviction_count (int): Number of entries dropped so far because the pooled region
            changed.
    """

    max_bytes: int = DEFAULT_TRANSCRIPT_POOL_BYTES
    hit_count: int = field(default=0, init=False)
    miss_count: int = field(default=0, init=False)
    eviction_count: int = field(default=0, init=False)
    region_eviction_count: int = field(default=0, init=False)
    _region: str | None = field(default=None, init=False, repr=False)
    _transcripts: OrderedDict[TranscriptKey, Transcript] = field(factory=OrderedDict, init=False, repr=False)
    _sizes: dict[TranscriptKey, int] = field(factory=dict, init=False, repr=False)
    _total_bytes: int = field(default=0, init=False, repr=False)
    _lock: RLock = field(factory=RLock, init=False, repr=False)

    def get(self, key: TranscriptKey) -> Transcript | None:
        """Return the pooled transcript for ``key``, moving it to most-recently-used.

        Args:
            key (TranscriptKey): The transcript's pool identity.

        Returns:
            Transcript | None: The cached transcript, or None on a miss (including when the pool
            is disabled).
        """
        if self.max_bytes <= 0:
            return None

        with self._lock:
            transcript = self._transcripts.get(key)
            if transcript is None:
                self.miss_count += 1
                return None
            self._transcripts.move_to_end(key)
            self.hit_count += 1
            return transcript

    def put(self, key: TranscriptKey, transcript: Transcript) -> None:
        """Insert a freshly built transcript, evicting by region change and then by byte budget.

        A no-op when the pool is disabled (``max_bytes <= 0``).

        Args:
            key (TranscriptKey): The transcript's pool identity. ``key[0]`` is the region.
            transcript (Transcript): The built transcript to store.
        """
        if self.max_bytes <= 0:
            return

        region = key[0]
        size = _estimate_transcript_bytes(transcript)

        with self._lock:
            if self._region is not None and self._region != region:
                self._drop_region_locked(self._region)
            self._region = region

            if key in self._sizes:
                self._total_bytes -= self._sizes[key]
            self._transcripts[key] = transcript
            self._transcripts.move_to_end(key)
            self._sizes[key] = size
            self._total_bytes += size

            while self._total_bytes > self.max_bytes and len(self._transcripts) > 1:
                evicted_key, _ = self._transcripts.popitem(last=False)
                self._total_bytes -= self._sizes.pop(evicted_key, 0)
                self.eviction_count += 1

    def clear(self) -> None:
        """Clear the pool contents."""
        with self._lock:
            self._transcripts.clear()
            self._sizes.clear()
            self._total_bytes = 0
            self._region = None

    def __len__(self) -> int:
        """Return the number of pooled transcripts."""
        return len(self._transcripts)

    def _drop_region_locked(self, region: str) -> None:
        """Remove every entry for ``region``. Caller must hold ``_lock``."""
        stale_keys = [key for key in self._transcripts if key[0] == region]
        for key in stale_keys:
            del self._transcripts[key]
            self._total_bytes -= self._sizes.pop(key, 0)
        self.region_eviction_count += len(stale_keys)
