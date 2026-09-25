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

import logging
import os
from pathlib import Path

from attrs import define, field
from attrs.setters import frozen

from ensembl.vep.parsers.fasta import FastaParser
from ensembl.vep.sources.sequence_source import SequenceSource

logger = logging.getLogger(__name__)


@define(kw_only=True, auto_attribs=True)
class FastaSequenceSource(SequenceSource):
    """Sequence source backed by a FASTA file.

    Attributes:
        file_path (str | Path | None): FASTA file path.
    """

    file_path: str | Path | None = field(default=None, on_setattr=frozen)
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

        if self.adaptor.reader_mode == "sequential":
            self._warn_about_sequential_access()

        if not self.adaptor.seek(seq_region_name=seq_region_name, start=start, end=end):
            return "N" * (end - start + 1)

        record = self.adaptor.next()
        if record is None:
            return "N" * (end - start + 1)

        return record.sequence

    def _get_seq_region_length(self, seq_region_name: str) -> int | None:
        """Returns the FASTA sequence length for a region.

        Args:
            seq_region_name (str): The sequence region name.

        Query-mode metadata comes from the Rust reader's ``.fai`` index. Without an index, the
        Rust reader derives the same lengths from one scan that retains no sequence.

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
                if self.adaptor.reader_mode == "sequential":
                    self._warn_about_sequential_access()
            self._adaptor_pid = None
        self._metadata_lengths = {
            str(name): int(length) for name, length in metadata.get("sequence_lengths", {}).items()
        }

    def _warn_about_sequential_access(self) -> None:
        """Warn once that an index is required for scalable random access."""

        if self._warned_sequential:
            return
        logger.warning(
            "%s has no usable FASTA index; sequential fallback is bounded but may be slow",
            self.file_path,
        )
        self._warned_sequential = True
