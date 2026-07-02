"""Abstract base for sequence sources."""

from __future__ import annotations

from abc import ABC, abstractmethod
from threading import RLock

from attrs import define, field

from ensembl.vep.caches import SequenceChunkCache
from ensembl.vep.common import Strand, reverse_complement


@define(kw_only=True, auto_attribs=True)
class SequenceSource(ABC):
    """Unified interface for sources that provide sequence for slices."""

    cache: SequenceChunkCache | None = field(default=None)
    _seq_region_lengths: dict[str, int] = field(factory=dict, init=False, repr=False)
    _lock: RLock = field(factory=RLock, init=False, repr=False)

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

    def get_seq_region_length(self, seq_region_name: str) -> int | None:
        """Return the length of a sequence region, caching the answer when available.

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

        with self._lock:
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
        """Return sequence for a linear span, on the requested strand.

        Args:
            seq_region_name (str): The sequence region name.
            start (int): 1-based inclusive start coordinate.
            end (int): 1-based inclusive end coordinate.
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
            sequence = self._load_sequence_span_cached(seq_region_name, start, end)

        if sequence is None:
            sequence = "N" * requested_length

        if strand == Strand.REVERSE:
            return reverse_complement(sequence)

        return sequence

    def _load_sequence_span_cached(self, seq_region_name: str, start: int, end: int) -> str | None:
        """Load a span using the configured chunk cache."""

        assert self.cache is not None
        chunk_size = self.cache.chunk_size
        chunk_start = self.cache.chunk_start(start)
        sequence = ""

        for current_start in range(chunk_start, end + 1, chunk_size):
            current_end = min(end, current_start + chunk_size - 1)
            chunk = self._get_chunk_cached(seq_region_name, current_start)
            local_start = max(start, current_start) - current_start
            local_end = current_end - current_start + 1
            sequence += chunk[local_start:local_end]

        return self._pad_to_requested_length(sequence, start, end)

    def _get_chunk_cached(self, seq_region_name: str, chunk_start: int) -> str:
        """Return a cached chunk or load it from the concrete source."""

        assert self.cache is not None
        cache_key = (seq_region_name, chunk_start)

        def loader() -> str:
            chunk_end = chunk_start + self.cache.chunk_size - 1
            return self._fetch_sequence_chunk(seq_region_name, chunk_start, chunk_end)

        return self.cache.get_or_load(cache_key, loader)

    def _load_sequence_span(self, seq_region_name: str, start: int, end: int) -> str | None:
        """Load a span directly from the concrete source without caching."""

        sequence = ""
        current = start
        while current <= end:
            current_end = min(end, current + self.chunk_size - 1)
            chunk = self._fetch_sequence_chunk(seq_region_name, current, current_end)
            sequence += chunk
            current = current_end + 1

        return self._pad_to_requested_length(sequence, start, end)

    def _pad_to_requested_length(self, sequence: str, start: int, end: int) -> str:
        """Pad a span with Ns when the source cannot provide the full range."""

        requested_length = end - start + 1
        if len(sequence) < requested_length:
            sequence += "N" * (requested_length - len(sequence))

        return sequence

    @property
    def chunk_size(self) -> int:
        """Return the active chunk size."""

        if self.cache is not None:
            return self.cache.chunk_size

        return 1 << 12

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
            start (int): 1-based inclusive start coordinate.
            end (int): 1-based inclusive end coordinate.

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
