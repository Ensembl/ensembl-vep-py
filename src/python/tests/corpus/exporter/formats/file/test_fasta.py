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
"""Tests for FASTA corpus export."""

from __future__ import annotations

from pathlib import Path

import pytest

from ensembl.vep.common import CoordSystemLevel, Strand
from ensembl.vep.corpus import Corpus
from ensembl.vep.corpus.exporter.formats.file.fasta import FASTA_NAME, fasta_lines
from ensembl.vep.corpus.exporter.manifest import CorpusManifest
from ensembl.vep.models.sequence import CoordSystem, Slice
from ensembl.vep.parsers.fasta import FastaParser
from ensembl.vep.sources.file.fasta.sequence import FastaSequenceSource


class TestFasta:
    """Chromosome sequences exported as FASTA."""

    def test_writes_each_sequence_region_in_full(
        self, source: Corpus, saved: tuple[Path, CorpusManifest]
    ) -> None:
        """Every corpus sequence region is represented by its complete sequence."""
        directory, _manifest = saved
        records = {record.seqid: record.sequence for record in FastaParser.parse_file(directory / FASTA_NAME)}
        assert records == {region.seq_region_name: region.seq for region in source.sequence_regions}

    def test_saved_fasta_is_a_sequence_source(
        self, source: Corpus, saved: tuple[Path, CorpusManifest]
    ) -> None:
        """A normal sequence source reads back the first exported region unchanged."""
        directory, _manifest = saved
        original = source.sequence_regions[0]
        restored = Slice(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name=original.seq_region_name,
            start=1,
            end=original.end,
            strand=Strand.FORWARD,
            sequence_source=FastaSequenceSource(file_path=directory / FASTA_NAME),
        )
        assert restored.seq == original.seq

    def test_exporting_empty_sequence_is_error(self, source: Corpus) -> None:
        """An empty sequence does not get output."""
        original = source.sequence_regions[0]
        with original.unsafe() as unsafe_slice:
            setattr(unsafe_slice, "_seq", "")
        with pytest.raises(ValueError, match="No sequence for"):
            list(fasta_lines([source]))
