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

from ensembl.vep.common.constants import INPUT_FORMAT_EXTENSIONS
from ensembl.vep.config.sections.common.source import SourceConfig, _input_formats
from ensembl.vep.sources.build import (
    build_chromosome_synonyms,
    build_input_source,
    build_sequence_source,
    build_transcript_source,
    detect_input_format,
)
from ensembl.vep.sources.file.fasta.sequence import FastaSequenceSource
from ensembl.vep.sources.file.gff3.transcript import GFF3TranscriptSource
from ensembl.vep.sources.file.gtf.transcript import GTFTranscriptSource
from ensembl.vep.sources.file.vcf.variation import PicklableVCFVariationSource
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


class TestBuildTranscriptSource:
    """Choosing and wiring the transcripts a run annotates against."""

    def test_builds_a_gff3_source_sharing_the_reference(self) -> None:
        """The transcript source reads sequence from the same source the input does."""

        sequence_source = build_sequence_source(_config())
        source = build_transcript_source(
            _config(), sequence_source=sequence_source, synonyms=ChromosomeSynonyms()
        )

        assert isinstance(source, GFF3TranscriptSource)
        assert source.sequence_source is sequence_source

    def test_builds_a_gtf_source_when_one_is_configured(self) -> None:
        """Select the GTF transcript source when requested."""

        config = _config(gff=None, gtf=Path("transcripts.gtf"))
        source = build_transcript_source(
            config,
            sequence_source=build_sequence_source(config),
            synonyms=ChromosomeSynonyms(),
        )

        assert isinstance(source, GTFTranscriptSource)

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

        source = build_transcript_source(
            config, sequence_source=build_sequence_source(config), synonyms=synonyms
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
            build_transcript_source(
                config,
                sequence_source=build_sequence_source(config),
                synonyms=ChromosomeSynonyms(),
            )

    @pytest.mark.parametrize(
        ("overrides", "message"),
        [
            pytest.param({"cache": True, "dir_cache": Path("~/.vep")}, "VEP cache", id="cache"),
            pytest.param({"database": True}, "Ensembl database", id="database"),
        ],
    )
    def test_names_the_source_it_cannot_build_yet(self, overrides: dict[str, object], message: str) -> None:
        """A configured source with no reader says which one, rather than being ignored.

        Perl VEP annotates from a cache and a database as well as from GFF3/GTF. The options exist
        here so a run asking for one is told plainly, instead of silently annotating against
        something else. ``--custom`` is refused by its own track builder, which is what a run asks
        for custom annotation through.
        """

        config = _config(**overrides)

        with pytest.raises(NotImplementedError, match=message):
            build_transcript_source(
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

        assert isinstance(source, PicklableVCFVariationSource)
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

    def test_resolves_slices_without_changing_the_record(self, monkeypatch: pytest.MonkeyPatch) -> None:
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
        record = "1\t10\t.\tA\tT\t.\tPASS\t."

        assert isinstance(source, PicklableVCFVariationSource)
        variation = source.records_to_features([record])[0]

        assert variation.uploaded_locus is not None
        assert variation.uploaded_locus.seq_region_name == "1"
        assert variation.slice is not None
        assert variation.slice.seq_region_name == "chr1"


class TestRegionResolution:
    """What names a source resolves regions into."""

    def test_a_source_built_without_a_reference_resolves_no_names(self) -> None:
        """A run given no FASTA leaves every region name as the input spells it."""

        source = build_input_source(_config(fasta=None), sequence_source=None, synonyms=ChromosomeSynonyms())

        assert source.sequence_region_resolver is not None
        assert source.sequence_region_resolver.source_names == frozenset()


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


class TestDetectInputFormat:
    """Which reader an input file is read with."""

    @pytest.mark.parametrize(
        "name",
        [
            pytest.param("variants.vcf", id="plain"),
            pytest.param("variants.VCF", id="upper_case"),
            pytest.param("variants.bcf", id="bcf"),
            pytest.param("variants.vcf.gz", id="gzip"),
            pytest.param("variants.vcf.bgz", id="bgzip"),
        ],
    )
    def test_a_vcf_is_recognised_however_it_is_named(self, name: str) -> None:
        """Read the extension that names the content, not the compression."""

        assert detect_input_format(Path(name)) == "vcf"

    def test_an_unrecognised_name_is_read_as_vcf(self) -> None:
        """Keep reading an input of any name as VCF, which is what VEP has always done."""

        assert detect_input_format(Path("variant_effect_output.txt")) == "vcf"

    def test_a_configured_format_wins_over_the_name(self) -> None:
        """Only the caller knows what an unqualified file holds."""

        assert detect_input_format(Path("variants.txt"), "VCF") == "vcf"

    def test_a_format_with_no_reader_is_refused(self) -> None:
        """Say what can be read rather than reading the file as something else."""

        with pytest.raises(ValueError, match="no reader for input format 'bed'"):
            detect_input_format(Path("variants.bed"), "bed")

    def test_the_format_names_come_from_one_table(self) -> None:
        """``--format`` offers what the build conditional reads, which is one closed set."""

        assert _input_formats() == tuple(sorted(INPUT_FORMAT_EXTENSIONS))

    def test_the_configured_format_selects_the_source(self) -> None:
        """``--format`` reaches the source that is built, not only the name reported."""

        config = _config(format="vcf")

        source = build_input_source(config, sequence_source=None, synonyms=ChromosomeSynonyms())

        assert isinstance(source, PicklableVCFVariationSource)

    def test_an_unreadable_format_stops_the_build(self) -> None:
        """A run refuses an input format it has no reader for."""

        config = _config(format="bed")

        with pytest.raises(ValueError, match="no reader for input format 'bed'"):
            build_input_source(config, sequence_source=None, synonyms=ChromosomeSynonyms())
