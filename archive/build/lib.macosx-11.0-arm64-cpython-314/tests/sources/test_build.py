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
"""Tests for building a run's sources from its configuration."""

from __future__ import annotations

from pathlib import Path

import pytest

from ensembl.vep.config.sections.common.source import SourceConfig
from ensembl.vep.parsers.vcf import VCFParser
from ensembl.vep.sources.build import (
    build_annotation_sources,
    build_chromosome_synonyms,
    build_input_source,
    build_sequence_source,
    transcript_source_of,
)
from ensembl.vep.sources.file.fasta.sequence import FastaSequenceSource
from ensembl.vep.sources.file.gff3.transcript import GFF3TranscriptSource
from ensembl.vep.sources.file.gtf.transcript import GTFTranscriptSource
from ensembl.vep.sources.file.vcf.variation import VCFVariationSource
from ensembl.vep.sources.synonym import ChromosomeSynonyms


def _config(**overrides: object) -> SourceConfig:
    """Return a configuration naming a GFF3 and a FASTA, with the given fields replaced."""

    fields: dict[str, object] = {
        "input_file": Path("variants.vcf"),
        "gff": Path("transcripts.gff3"),
        "fasta": Path("reference.fa"),
    }
    fields.update(overrides)
    return SourceConfig(**fields)  # type: ignore[arg-type]


class TestBuildAnnotationSources:
    """Choosing and wiring what a run annotates against."""

    def test_builds_a_gff3_source_sharing_the_reference(self) -> None:
        """The transcript source reads sequence from the same source the input does."""

        sequence_source = build_sequence_source(_config())
        sources = build_annotation_sources(
            _config(), sequence_source=sequence_source, synonyms=ChromosomeSynonyms()
        )

        assert isinstance(sources[0], GFF3TranscriptSource)
        assert sources[0].sequence_source is sequence_source

    def test_builds_a_gtf_source_when_one_is_configured(self) -> None:
        """Select the GTF transcript source when requested."""

        config = _config(gff=None, gtf=Path("transcripts.gtf"))
        sources = build_annotation_sources(
            config,
            sequence_source=build_sequence_source(config),
            synonyms=ChromosomeSynonyms(),
        )

        assert isinstance(sources[0], GTFTranscriptSource)

    def test_wires_name_resolution_in_both_directions(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Records read from the file get reference names, and queries go out in the file's names."""
        fasta = tmp_path / "reference.fa"
        fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
        gff = tmp_path / "transcripts.gff3"
        gff.write_text(
            "##gff-version 3\n"
            "##sequence-region 1 1 4\n"
            "1\tensembl\tmRNA\t1\t4\t.\t+\t.\tID=transcript:1\n",
            encoding="utf-8",
        )
        monkeypatch.setattr(FastaSequenceSource, "_get_seq_region_length", lambda _source, _name: 4)
        config = _config(gff=gff, fasta=fasta)
        synonyms = ChromosomeSynonyms(
            by_name={"chr1": frozenset({"1"}), "1": frozenset({"chr1"})},
        )

        source = transcript_source_of(
            build_annotation_sources(config, sequence_source=build_sequence_source(config), synonyms=synonyms)
        )

        assert source.reference_region_resolver is not None
        assert source.reference_region_resolver.resolve("1") == "chr1"
        assert source.source_region_resolver is not None
        assert source.source_region_resolver.resolve("chr1") == "1"

    @pytest.mark.parametrize(
        ("overrides", "message"),
        [
            pytest.param({"gtf": Path("a.gtf")}, "give --gff or --gtf", id="both_formats"),
            pytest.param({"gff": None}, "needs transcripts", id="no_format"),
        ],
    )
    def test_rejects_ambiguous_or_absent_transcripts(
        self, overrides: dict[str, object], message: str
    ) -> None:
        """Reject a transcript configuration that names two sources, or none."""

        config = _config(**overrides)

        with pytest.raises(ValueError, match=message):
            build_annotation_sources(
                config,
                sequence_source=build_sequence_source(config),
                synonyms=ChromosomeSynonyms(),
            )

    @pytest.mark.parametrize(
        ("overrides", "message"),
        [
            pytest.param({"cache": True, "dir_cache": Path("~/.vep")}, "VEP cache", id="cache"),
            pytest.param({"database": True}, "Ensembl database", id="database"),
            pytest.param({"custom": ("clinvar.vcf.gz,ClinVar,vcf",)}, "custom annotation", id="custom"),
        ],
    )
    def test_names_the_source_it_cannot_build_yet(self, overrides: dict[str, object], message: str) -> None:
        """A configured source with no reader says which one, rather than being ignored.

        Perl VEP annotates from a cache, a database and custom files as well as from GFF3/GTF. The
        options exist here so a run asking for one is told plainly, instead of silently annotating
        against something else.
        """

        config = _config(**overrides)

        with pytest.raises(NotImplementedError, match=message):
            build_annotation_sources(
                config,
                sequence_source=build_sequence_source(config),
                synonyms=ChromosomeSynonyms(),
            )


class TestBuildInputSource:
    """The source holding the variants being annotated."""

    def test_builds_a_vcf_source_sharing_the_reference(self) -> None:
        """The input source reads sequence from the shared reference source."""

        config = _config()
        sequence_source = build_sequence_source(config)

        source = build_input_source(config, sequence_source=sequence_source, synonyms=ChromosomeSynonyms())

        assert isinstance(source, VCFVariationSource)
        assert source.sequence_source is sequence_source

    def test_rejects_a_configuration_with_no_variants(self) -> None:
        """There is nothing to annotate without an input file."""

        config = _config(input_file=None)

        with pytest.raises(ValueError, match="needs variants"):
            build_input_source(
                config,
                sequence_source=build_sequence_source(config),
                synonyms=ChromosomeSynonyms(),
            )

    def test_resolves_slices_without_changing_vcf_records(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Name resolution happens at the source boundary; the record keeps its own spelling."""

        monkeypatch.setattr(
            FastaSequenceSource,
            "sequence_region_names",
            property(lambda _source: frozenset({"chr1"})),
        )
        monkeypatch.setattr(FastaSequenceSource, "_get_seq_region_length", lambda _source, _name: 100)
        config = _config()
        source = build_input_source(
            config,
            sequence_source=build_sequence_source(config),
            synonyms=ChromosomeSynonyms(),
        )
        record = VCFParser().parse_text(
            "##fileformat=VCFv4.3\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "1\t10\t.\tA\tT\t.\tPASS\t.\n"
        )[0]

        variation = source.records_to_features([record])[0]

        assert record.chrom == "1"
        assert variation.slice is not None
        assert variation.slice.seq_region_name == "chr1"


class TestBuildSequenceSource:
    """The reference every other source shares."""

    def test_rejects_a_configuration_with_no_reference(self) -> None:
        """Consequences need reference sequence, so a FASTA is not optional here."""

        with pytest.raises(ValueError, match="needs a reference"):
            build_sequence_source(_config(fasta=None))


class TestBuildChromosomeSynonyms:
    """Reading the run's synonym map from configuration."""

    def test_reads_the_configured_file(self, tmp_path: Path) -> None:
        """``--synonyms`` names the file the map comes from."""

        path = tmp_path / "synonyms.txt"
        path.write_text("chr1\tNC_000001.11\n", encoding="utf-8")

        synonyms = build_chromosome_synonyms(_config(synonyms=path))

        assert synonyms.synonyms_of("NC_000001.11") == frozenset({"chr1"})

    def test_has_no_synonyms_without_the_option(self) -> None:
        """A run that names no synonyms file resolves names on exact match and ``chr`` alone."""

        assert not build_chromosome_synonyms(_config()).by_name


class TestTranscriptSourceOf:
    """Picking the transcript source out of a run's annotation sources."""

    def test_returns_the_transcript_source(self) -> None:
        """Consequence calling asks for transcripts by what they are, not where they came from."""

        config = _config()
        sources = build_annotation_sources(
            config,
            sequence_source=build_sequence_source(config),
            synonyms=ChromosomeSynonyms(),
        )

        assert transcript_source_of(sources) is sources[0]

    def test_reports_a_set_with_no_transcripts(self) -> None:
        """Say so plainly rather than failing later with no transcript to annotate against."""

        with pytest.raises(ValueError, match="no transcript source"):
            transcript_source_of([])
