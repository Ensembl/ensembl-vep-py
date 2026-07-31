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
"""Abstract base for sequence sources."""

from __future__ import annotations

from abc import ABC, abstractmethod
from threading import RLock

from attrs import define, field

from ensembl.vep.caches import SequenceCache
from ensembl.vep.common import Strand, reverse_complement


@define(kw_only=True, auto_attribs=True)
class SequenceSource(ABC):
    """Unified interface for sources that provide sequence for slices."""

    cache: SequenceCache | None = field(default=None)
    _seq_region_lengths: dict[str, int] = field(factory=dict, init=False, repr=False)
    _lock: RLock = field(factory=RLock, init=False, repr=False)

    def get_seq_region_length(self, seq_region_name: str) -> int | None:
        """Return the length of a sequence region, caching when available.

        Args:
            seq_region_name (str): The sequence region name.

        Returns:
            int | None: The sequence region length, or None if unavailable.
        """

        with self._lock:
            cached = self._seq_region_lengths.get(seq_region_name)
            if cached is not None:
                return cached

            length = self._get_seq_region_length(seq_region_name)
            if length is None:
                return None
        
            self._seq_region_lengths[seq_region_name] = length

        return length

    def get_sequence_by_start_end_strand(
        self,
        seq_region_name: str,
        start: int,
        end: int,
        *,
        strand: Strand = Strand.FORWARD,
    ) -> str | None:
        """Return sequence from a sequence region for the specified coordinates and strand.

        Args:
            seq_region_name (str): The sequence region name.
            start (int): Start coordinate for the sequence span.
            end (int): End coordinate for the sequence span.
            strand (Strand): Requested strand orientation.

        Returns:
            str | None: The requested sequence, or None if unavailable.
        """

        if start < 1 or end < start:
            raise ValueError("sequence spans must be 1-based and start must be <= end")

        requested_length = end - start + 1

        if self.cache is None:
            sequence = self._load_sequence_span(seq_region_name, start, end)
        else:
            sequence = self.cache.get_sequence(
                seq_region_name,
                start,
                end,
                self._load_sequence_span,
            )

        if sequence is None:
            sequence = "N" * requested_length
        else:
            sequence = self._pad_to_requested_length(sequence, start, end)

        if strand == Strand.REVERSE:
            return reverse_complement(sequence)

        return sequence

    @abstractmethod
    def close(self) -> None:
        """Release any resources owned by the source."""

        raise NotImplementedError

    def __enter__(self) -> SequenceSource:
        """Return the source instance for use as a context manager.

        Returns:
            SequenceSource: The source instance.
        """

        return self

    def __exit__(self, *_exc_info: object) -> None:
        """Close the source when the context exits."""

        self.close()

    def _load_sequence_span(self, seq_region_name: str, start: int, end: int) -> str | None:
        """Load a span directly from the concrete source without using cache.

        Args:
            seq_region_name (str): The sequence region name.
            start (int): Start coordinate for the sequence span.
            end (int): End coordinate for the sequence span.

        Returns:
            str | None: The requested sequence, or None if unavailable.
        """

        sequence = self._fetch_sequence_chunk(seq_region_name, start, end)
        return self._pad_to_requested_length(sequence, start, end)

    def _pad_to_requested_length(self, sequence: str, start: int, end: int) -> str:
        """Pad a span with Ns when the source cannot provide the full range.

        Args:
            sequence (str): The sequence to pad.
            start (int): Start coordinate for the sequence span.
            end (int): End coordinate for the sequence span.

        Returns:
            str: The padded sequence.
        """

        requested_length = end - start + 1
        if len(sequence) < requested_length:
            sequence += "N" * (requested_length - len(sequence))

        return sequence

    def clear_cache(self) -> None:
        """Clear derived sequence caches."""

        with self._lock:
            self._seq_region_lengths.clear()

        if self.cache is not None:
            self.cache.clear()

    @abstractmethod
    def _fetch_sequence_chunk(self, seq_region_name: str, start: int, end: int) -> str:
        """Return a forward-oriented chunk of sequence.

        Args:
            seq_region_name (str): The sequence region name.
            start (int): Start coordinate for the sequence span.
            end (int): End coordinate for the sequence span.

        Returns:
            str: The requested sequence chunk, padded or truncated as needed.
        """

        raise NotImplementedError

    @abstractmethod
    def _get_seq_region_length(self, seq_region_name: str) -> int | None:
        """Return the sequence region length if available.

        Args:
            seq_region_name (str): The sequence region name.

        Returns:
            int | None: The sequence region length, or None if unknown.
        """

        raise NotImplementedError
