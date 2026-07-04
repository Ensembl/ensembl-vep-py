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
    _parser: FastaParser = field(init=False, repr=False)
    _sequences: dict[str, str] | None = field(default=None, init=False, repr=False)

    def __attrs_post_init__(self) -> None:
        """Initialize the underlying FASTA parser."""

        self._parser = FastaParser(
            file_path=str(self.file_path) if self.file_path is not None else None,
            mode=self.mode,
        )

    def open(self, file_path: str | Path, mode: str = "auto") -> FastaSequenceSource:
        """Attach a FASTA file to this source.

        Args:
            file_path (str | Path): FASTA file path.
            mode (str): Reader mode used by the parser.

        Returns:
            FastaSequenceSource: The source instance.
        """

        self.file_path = file_path
        self.mode = mode
        self._parser.open(str(file_path), mode)
        self.clear_cache()
        return self

    def close(self) -> None:
        """Release loaded FASTA sequence data."""

        self._parser.close()
        self.clear_cache()

    def clear_cache(self) -> None:
        """Clear cached FASTA sequence data."""

        super().clear_cache()
        self._sequences = None

    def _fetch_sequence_chunk(self, seq_region_name: str, start: int, end: int) -> str:
        """Return a FASTA sequence chunk.

        Args:
            seq_region_name (str): The sequence region name.
            start (int): 1-based inclusive start coordinate.
            end (int): 1-based inclusive end coordinate.

        Returns:
            str: The requested chunk or Ns if unavailable.
        """

        if self._parser.reader_mode == "query":
            if not self._parser.seek(seq_region_name=seq_region_name, start=start, end=end):
                return "N" * (end - start + 1)

            record = self._parser.next()
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

        if self._parser.reader_mode == "query":
            if not self._parser.seek(seq_region_name=seq_region_name):
                return None

            record = self._parser.next()
            if record is None:
                return None

            return len(record.sequence)

        sequence = self._load_sequences().get(seq_region_name)
        return len(sequence) if sequence is not None else None

    def _load_sequences(self) -> dict[str, str]:
        """Load all FASTA records into an in-memory lookup map."""

        if self._sequences is None:
            self._sequences = {record.seqid: record.sequence for record in self._parser.records(limit=None)}

        return self._sequences
