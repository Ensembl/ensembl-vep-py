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

from pathlib import Path

from attrs import define, field

from ensembl.vep.parsers.fasta import FastaParser
from ensembl.vep.sources.sequence_source import SequenceSource


@define(kw_only=True, auto_attribs=True)
class FastaSequenceSource(SequenceSource):
    """Sequence source backed by a FASTA file.

    Attributes:
        file_path (str | Path | None): FASTA file path.
        mode (str): Reader mode used by the parser.
    """

    file_path: str | Path | None = field(default=None)
    mode: str = field(default="auto")
    adaptor: FastaParser = field(init=False, factory=FastaParser)
    _sequences: dict[str, str] | None = field(default=None, init=False, repr=False)

    def close(self) -> None:
        """Release loaded FASTA sequence data."""

        self.adaptor.close()
        self.clear_cache()
    
    def _fetch_sequence_chunk(self, seq_region_name: str, start: int, end: int) -> str:
        """Return a FASTA sequence chunk.

        Args:
            seq_region_name (str): The sequence region name.
            start (int): 1-based inclusive start coordinate.
            end (int): 1-based inclusive end coordinate.

        Returns:
            str: The requested chunk or Ns if unavailable.
        """

        with self.adaptor.open(self.file_path, mode="auto"):
            if self.adaptor.reader_mode == "query":
                if not self.adaptor.seek(seq_region_name=seq_region_name, start=start, end=end):
                    return "N" * (end - start + 1)

                record = self.adaptor.next()
                if record is None:
                    return "N" * (end - start + 1)

                return record.sequence

        sequence = self._load_sequences().get(seq_region_name)
        if sequence is None:
            return "N" * (end - start + 1)

        raw_start = start - 1
        raw_end = end
        chunk = sequence[raw_start:raw_end]
        if len(chunk) < (end - start + 1):
            chunk += "N" * ((end - start + 1) - len(chunk))

        return chunk

    def _get_seq_region_length(self, seq_region_name: str) -> int | None:
        """Return the FASTA sequence length for a region.

        Args:
            seq_region_name (str): The sequence region name.

        Returns:
            int | None: The sequence length or None if unavailable.
        """

        lengths = self.adaptor.metadata.get("sequence_lengths", {})
        length = lengths.get(seq_region_name)
        return int(length) if length is not None else None

    def _load_sequences(self) -> dict[str, str]:
        """Load all FASTA records into an in-memory lookup map."""

        if self._sequences is None:
            with self.adaptor.open(self.file_path, mode="sequential"):
                self._sequences = {record.seqid: record.sequence for record in self.adaptor.records(limit=None)}

        return self._sequences
