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
"""Tests for file-backed corpus loading."""

from __future__ import annotations

from pathlib import Path

import pytest

from ensembl.vep.corpus import LoadedCorpus, build_corpus
from ensembl.vep.corpus.describe import variant_class
from ensembl.vep.corpus.exporter import get_exporter
from ensembl.vep.corpus.exporter.formats.file.fasta import FASTA_NAME
from ensembl.vep.corpus.exporter.formats.file.gff3 import GFF3_NAME
from ensembl.vep.corpus.exporter.formats.file.vcf import VCF_NAME
from ensembl.vep.corpus.spec import CorpusSpec
from ensembl.vep.sources.file.fasta.sequence import FastaSequenceSource


def _file_spec(directory: Path) -> CorpusSpec:
    """Return the specification for a corpus file bundle."""
    return CorpusSpec(
        kind="loaded.files",
        parameters={
            "fasta": str(directory / FASTA_NAME),
            "gff3": str(directory / GFF3_NAME),
            "vcf": str(directory / VCF_NAME),
        },
    )


@pytest.fixture(name="saved_default_corpus")
def fixture_saved_default_corpus(tmp_path: Path) -> CorpusSpec:
    """Write a small corpus and return the specification that reloads it."""
    synthetic_corpus = build_corpus(CorpusSpec(kind="synthetic.default", parameters={"fuzz": 0}))
    get_exporter("file").export(synthetic_corpus, tmp_path)
    return _file_spec(tmp_path)


class TestLoadedCorpus:
    """A file-backed corpus loaded through the registry."""

    def test_rejects_non_tuple_membership(self, saved_default_corpus: CorpusSpec) -> None:
        """Require loaders to materialise membership before construction."""
        with build_corpus(saved_default_corpus) as loaded_corpus:
            with pytest.raises(TypeError, match="tuple"):
                LoadedCorpus(
                    spec=loaded_corpus.spec,
                    sequence_regions=list(loaded_corpus.sequence_regions),  # type: ignore[arg-type]
                    named_variants=loaded_corpus.named_variants,
                    transcripts=loaded_corpus.transcripts,
                    named_overlaps=loaded_corpus.named_overlaps,
                )

    def test_reloads_variants_transcripts_and_overlaps(
        self,
        saved_default_corpus: CorpusSpec,
    ) -> None:
        """Retain every exported corpus relationship on reload."""
        loaded_corpus = build_corpus(saved_default_corpus)

        assert isinstance(loaded_corpus, LoadedCorpus)
        assert loaded_corpus.named_variants
        assert loaded_corpus.transcripts
        assert loaded_corpus.named_overlaps
        assert any(
            name not in {overlap_name for overlap_name, _overlap in loaded_corpus.named_overlaps}
            for name, _allele in loaded_corpus.named_variants
        )

    def test_duplicate_vcf_ids_receive_globally_unique_names(
        self,
        saved_default_corpus: CorpusSpec,
    ) -> None:
        """Generated suffixes cannot collide with IDs that were already present in the VCF."""
        vcf_path = Path(str(saved_default_corpus.parameter("vcf")))
        lines = vcf_path.read_text(encoding="utf-8").splitlines()
        record_indexes = [index for index, line in enumerate(lines) if line and not line.startswith("#")]
        for index, stable_id in zip(record_indexes[:3], ("x:2", "x", "x"), strict=True):
            columns = lines[index].split("\t")
            columns[2] = stable_id
            lines[index] = "\t".join(columns)
        vcf_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

        with build_corpus(saved_default_corpus) as loaded_corpus:
            names = [name for name, _allele in loaded_corpus.named_variants]
        assert len(names) == len(set(names))
        assert {"x", "x:2", "x:3"} <= set(names)

    def test_shares_each_sequence_region_instance(
        self,
        saved_default_corpus: CorpusSpec,
    ) -> None:
        """Attach loaded transcripts and variants to the same region object."""
        loaded_corpus = build_corpus(saved_default_corpus)

        assert len(loaded_corpus.sequence_regions) == 1
        sequence_region = loaded_corpus.sequence_regions[0]
        assert all(transcript.slice is sequence_region for transcript in loaded_corpus.transcripts)
        assert all(allele.slice is sequence_region for _name, allele in loaded_corpus.named_variants)

    def test_resolves_uploaded_transcript_and_vcf_region_names(
        self,
        saved_default_corpus: CorpusSpec,
    ) -> None:
        """Reconcile uploaded ``chr`` prefixes with the FASTA region identity."""
        gff3_path = Path(str(saved_default_corpus.parameter("gff3")))
        gff3_lines = gff3_path.read_text(encoding="utf-8").splitlines()
        sequence_line = next(line for line in gff3_lines if line.startswith("##sequence-region "))
        reference_name = sequence_line.split()[1]
        uploaded_name = f"chr{reference_name}"
        gff3_path.write_text(
            "\n".join(
                (
                    line.replace(
                        f"##sequence-region {reference_name} ", f"##sequence-region {uploaded_name} "
                    )
                    if line.startswith("##sequence-region ")
                    else line.replace(f"{reference_name}\t", f"{uploaded_name}\t", 1)
                )
                for line in gff3_lines
            )
            + "\n",
            encoding="utf-8",
        )

        vcf_path = Path(str(saved_default_corpus.parameter("vcf")))
        vcf_path.write_text(
            vcf_path.read_text(encoding="utf-8")
            .replace(
                f"ID={reference_name},",
                f"ID={uploaded_name},",
            )
            .replace(f"{reference_name}\t", f"{uploaded_name}\t"),
            encoding="utf-8",
        )

        with build_corpus(saved_default_corpus) as loaded_corpus:
            assert [region.seq_region_name for region in loaded_corpus.sequence_regions] == [reference_name]
            assert loaded_corpus.named_overlaps
            assert all(
                transcript.slice is loaded_corpus.sequence_regions[0]
                for transcript in loaded_corpus.transcripts
            )
            assert all(
                allele.slice is loaded_corpus.sequence_regions[0]
                for _name, allele in loaded_corpus.named_variants
            )

    @pytest.mark.xfail(
        strict=True,
        reason="unresolved FASTA contigs currently receive transcript-variant overlaps",
    )
    def test_does_not_annotate_variants_on_unresolved_fasta_contig(
        self,
        saved_default_corpus: CorpusSpec,
    ) -> None:
        """Keep variants on an unknown reference separate from uploaded transcripts."""
        gff3_path = Path(str(saved_default_corpus.parameter("gff3")))
        gff3_lines = gff3_path.read_text(encoding="utf-8").splitlines()
        sequence_line = next(line for line in gff3_lines if line.startswith("##sequence-region "))
        reference_name = sequence_line.split()[1]
        missing_reference = "missing_reference"
        gff3_path.write_text(
            "\n".join(
                (
                    line.replace(
                        f"##sequence-region {reference_name} ",
                        f"##sequence-region {missing_reference} ",
                    )
                    if line.startswith("##sequence-region ")
                    else line.replace(f"{reference_name}\t", f"{missing_reference}\t", 1)
                )
                for line in gff3_lines
            )
            + "\n",
            encoding="utf-8",
        )

        vcf_path = Path(str(saved_default_corpus.parameter("vcf")))
        vcf_path.write_text(
            vcf_path.read_text(encoding="utf-8")
            .replace(f"ID={reference_name},", f"ID={missing_reference},")
            .replace(f"{reference_name}\t", f"{missing_reference}\t"),
            encoding="utf-8",
        )

        with build_corpus(saved_default_corpus) as loaded_corpus:
            assert loaded_corpus.transcripts
            assert loaded_corpus.named_variants
            assert not loaded_corpus.named_overlaps

    def test_vcf_empty_allele_spelling_preserves_variant_classes(
        self,
        saved_default_corpus: CorpusSpec,
    ) -> None:
        """The VCF reader's dash sentinel remains a deletion/insertion, not a substitution."""
        with build_corpus(saved_default_corpus) as loaded_corpus:
            by_name = dict(loaded_corpus.named_variants)
            assert variant_class(by_name["ENST0001_start_codon_del1"]) == "deletion"
            assert variant_class(by_name["ENST0001_start_codon_ins3"]) == "insertion"

    def test_gff3_transcript_ids_are_stable_across_round_trips(
        self,
        saved_default_corpus: CorpusSpec,
        tmp_path: Path,
    ) -> None:
        """Do not add a ``transcript:`` prefix each time a corpus is re-exported."""
        gff3_path = saved_default_corpus.parameter("gff3")
        assert isinstance(gff3_path, str)
        original_directory = Path(gff3_path).parent
        original_gff3 = (original_directory / GFF3_NAME).read_text(encoding="utf-8")
        first_corpus = build_corpus(saved_default_corpus)
        first_directory = tmp_path / "first-round-trip"
        get_exporter("file").export(first_corpus, first_directory)

        second_corpus = build_corpus(_file_spec(first_directory))
        second_directory = tmp_path / "second-round-trip"
        get_exporter("file").export(second_corpus, second_directory)

        first_gff3 = (first_directory / GFF3_NAME).read_text(encoding="utf-8")
        second_gff3 = (second_directory / GFF3_NAME).read_text(encoding="utf-8")
        assert [transcript.stable_id for transcript in first_corpus.transcripts] == [
            transcript.stable_id for transcript in second_corpus.transcripts
        ]
        assert original_gff3 == first_gff3 == second_gff3
        assert "ID=CDS:" in first_gff3
        assert "ID=CDS:;" not in first_gff3

    def test_gff3_translation_ids_are_stable_across_round_trips(
        self,
        saved_default_corpus: CorpusSpec,
        tmp_path: Path,
    ) -> None:
        """Check that CDS records are exported and reloaded properly"""
        first_corpus = build_corpus(saved_default_corpus)
        first_directory = tmp_path / "first-round-trip"
        get_exporter("file").export(first_corpus, first_directory)
        first_corpus_coding_translations = [
            transcript.translation
            for transcript in first_corpus.transcripts
            if transcript.biotype != "lincRNA"
        ]
        for translation in first_corpus_coding_translations:
            assert translation is not None

        second_corpus = build_corpus(_file_spec(first_directory))
        second_corpus_coding_translations = [
            transcript.translation
            for transcript in second_corpus.transcripts
            if transcript.biotype != "lincRNA"
        ]
        for translation in second_corpus_coding_translations:
            assert translation is not None

        assert [translation.stable_id for translation in first_corpus.transcripts] == [
            translation.stable_id for translation in second_corpus.transcripts
        ]

    def test_close_is_idempotent_and_context_managed(
        self,
        saved_default_corpus: CorpusSpec,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Close the owned source once across explicit and contextual cleanup."""
        loaded_corpus = build_corpus(saved_default_corpus)
        assert isinstance(loaded_corpus, LoadedCorpus)
        assert loaded_corpus.sequence_source is not None
        close_calls = 0
        sequence_source = loaded_corpus.sequence_source
        original_close = FastaSequenceSource.close

        def count_close(source: FastaSequenceSource) -> None:
            """Count resource cleanup while retaining the normal source behaviour."""
            nonlocal close_calls
            if source is sequence_source:
                close_calls += 1
            original_close(source)

        monkeypatch.setattr(FastaSequenceSource, "close", count_close)
        loaded_corpus.close()
        loaded_corpus.close()
        with loaded_corpus as context_corpus:
            assert context_corpus is loaded_corpus

        assert close_calls == 1

    @pytest.mark.parametrize(
        "parameter",
        [
            pytest.param("fasta", id="fasta"),
            pytest.param("gff3", id="gff3"),
            pytest.param("vcf", id="vcf"),
        ],
    )
    def test_rejects_missing_bundle_file(
        self,
        saved_default_corpus: CorpusSpec,
        parameter: str,
    ) -> None:
        """Name the missing input file instead of constructing a partial corpus."""
        parameters = saved_default_corpus.parameters
        parameters[parameter] = f"{parameters[parameter]}.missing"

        with pytest.raises(FileNotFoundError, match="missing"):
            build_corpus(CorpusSpec(kind="loaded.files", parameters=parameters))

    def test_rejects_ambiguous_transcript_format(self, saved_default_corpus: CorpusSpec) -> None:
        """Require exactly one supported transcript annotation format."""
        parameters = saved_default_corpus.parameters
        parameters["gtf"] = parameters["gff3"]

        with pytest.raises(ValueError, match="exactly one"):
            build_corpus(CorpusSpec(kind="loaded.files", parameters=parameters))
