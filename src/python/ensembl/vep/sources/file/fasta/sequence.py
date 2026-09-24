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
"""FASTA-backed sequence source."""

from __future__ import annotations

import gzip
import logging
import os
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import TextIO

from attrs import define, field

from ensembl.vep.parsers.fasta import FastaParser
from ensembl.vep.sources.sequence_source import SequenceSource

logger = logging.getLogger(__name__)


@define(kw_only=True, auto_attribs=True)
class FastaSequenceSource(SequenceSource):
    """Sequence source backed by a FASTA file.

    Attributes:
        file_path (str | Path | None): FASTA file path.
    """

    file_path: str | Path | None = field(default=None)
    adaptor: FastaParser = field(factory=FastaParser)
    _metadata_lengths: dict[str, int] | None = field(default=None, init=False, repr=False)
    _warned_sequential: bool = field(default=False, init=False, repr=False)
    _adaptor_pid: int | None = field(default=None, init=False, repr=False)

    def close(self) -> None:
        """Release FASTA metadata and derived sequence caches."""

        self._close_adaptor()
        self._metadata_lengths = None
        self.clear_cache()

    def _close_adaptor(self) -> None:
        """Close the parser and forget that it is open for span fetches."""

        self.adaptor.close()
        self._adaptor_pid = None

    def _fetch_sequence_chunk(self, seq_region_name: str, start: int, end: int) -> str:
        """Returns a FASTA sequence chunk.

        Args:
            seq_region_name (str): The sequence region name.
            start (int): 1-based inclusive start coordinate.
            end (int): 1-based inclusive end coordinate.

        Returns:
            str: The requested chunk or Ns if unavailable.
        """

        if self.file_path is None:
            return "N" * (end - start + 1)

        # Features fetch many short spans, so the reader stays open between them rather than
        # re-reading the index for every exon. A handle inherited through fork shares its file
        # offset with the parent, so a forked worker opens its own.
        if self._adaptor_pid != os.getpid():
            self.adaptor.open(self.file_path, mode="auto")
            self._adaptor_pid = os.getpid()

        if self.adaptor.reader_mode == "query":
            if not self.adaptor.seek(seq_region_name=seq_region_name, start=start, end=end):
                return "N" * (end - start + 1)

            record = self.adaptor.next()
            if record is None:
                return "N" * (end - start + 1)

            return record.sequence

        self._close_adaptor()
        self._warn_about_sequential_access()
        return self._scan_sequence_chunk(seq_region_name, start, end)

    def _get_seq_region_length(self, seq_region_name: str) -> int | None:
        """Returns the FASTA sequence length for a region.

        Args:
            seq_region_name (str): The sequence region name.

        Query-mode metadata comes from the Rust reader's ``.fai`` index. When auto mode falls back
        to sequential reading, a bounded scan derives the same authoritative lengths without
        retaining reference sequences.

        Returns:
            int | None: The sequence length or None if unavailable.
        """

        self._load_metadata_lengths()

        assert self._metadata_lengths is not None
        length = self._metadata_lengths.get(seq_region_name)
        return int(length) if length is not None else None

    @property
    def sequence_region_names(self) -> frozenset[str]:
        """Return names advertised by the FASTA metadata index."""

        self._load_metadata_lengths()
        assert self._metadata_lengths is not None
        return frozenset(self._metadata_lengths)

    def _load_metadata_lengths(self) -> None:
        """Load FASTA metadata lengths once."""

        if self._metadata_lengths is not None:
            return
        if self.file_path is None:
            metadata = self.adaptor.metadata
        else:
            with self.adaptor.open(self.file_path, mode="auto"):
                metadata = self.adaptor.metadata
            self._adaptor_pid = None
        self._metadata_lengths = {
            str(name): int(length) for name, length in metadata.get("sequence_lengths", {}).items()
        }
        if not self._metadata_lengths and self.file_path is not None:
            self._warn_about_sequential_access()
            self._metadata_lengths = self._scan_sequence_lengths()

    # TBD: instead of falling back to _RustFastaAdaptor sequental mode which does not have
    # seek method and reading sequence chunk efficiently (it reads the whole sequence) we are
    # using this Python fallback. This is not an important branch at this moment because we
    # expect a indexed fasta but we should implement a more efficient sequential mode in Rust.
    def _scan_sequence_lengths(self) -> dict[str, int]:
        """Return sequence lengths from a bounded sequential FASTA scan."""

        lengths: dict[str, int] = {}
        sequence_region_name: str | None = None
        with self._open_sequential_text() as handle:
            for line in handle:
                if line.startswith(">"):
                    sequence_region_name = self._header_name(line)
                    lengths.setdefault(sequence_region_name, 0)
                    continue
                if sequence_region_name is not None:
                    lengths[sequence_region_name] += len("".join(line.split()))
        return lengths

    def _scan_sequence_chunk(self, seq_region_name: str, start: int, end: int) -> str:
        """Return one span without retaining its complete unindexed reference contig."""

        chunk: list[str] = []
        position = 0
        in_requested_region = False
        with self._open_sequential_text() as handle:
            for line in handle:
                if line.startswith(">"):
                    if in_requested_region:
                        break
                    in_requested_region = self._header_name(line) == seq_region_name
                    position = 0
                    continue
                if not in_requested_region:
                    continue

                sequence_line = "".join(line.split()).upper()
                line_start = position + 1
                line_end = position + len(sequence_line)
                overlap_start = max(start, line_start)
                overlap_end = min(end, line_end)
                if overlap_start <= overlap_end:
                    local_start = overlap_start - line_start
                    local_end = overlap_end - line_start + 1
                    chunk.append(sequence_line[local_start:local_end])
                position = line_end
                if position >= end:
                    break

        sequence = "".join(chunk)
        return sequence + "N" * (end - start + 1 - len(sequence))

    @contextmanager
    def _open_sequential_text(self) -> Iterator[TextIO]:
        """Open plain or gzip-compatible FASTA text for a sequential fallback."""

        assert self.file_path is not None
        path = Path(self.file_path)
        with path.open("rb") as probe:
            compressed = probe.read(2) == b"\x1f\x8b"
        if compressed:
            with gzip.open(path, mode="rt", encoding="utf-8") as handle:
                yield handle
            return
        with path.open("r", encoding="utf-8") as handle:
            yield handle

    @staticmethod
    def _header_name(line: str) -> str:
        """Return the sequence-region token from a FASTA header."""

        return line[1:].strip().split(maxsplit=1)[0]

    def _warn_about_sequential_access(self) -> None:
        """Warn once that an index is required for scalable random access."""

        if self._warned_sequential:
            return
        logger.warning(
            "%s has no usable FASTA index; sequential fallback is bounded but may be slow",
            self.file_path,
        )
        self._warned_sequential = True
