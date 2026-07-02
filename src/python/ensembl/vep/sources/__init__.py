"""Sequence source interfaces."""

from .sequence_source import SequenceSource
from .file.fasta.sequence_source import FastaSequenceSource

__all__ = ["SequenceSource", "FastaSequenceSource"]
