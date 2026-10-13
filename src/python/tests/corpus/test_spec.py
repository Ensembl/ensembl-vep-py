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
"""Tests for corpus specifications and synthetic compatibility options."""

from __future__ import annotations

from pathlib import Path

import pytest

from ensembl.vep.corpus.base import Corpus
from ensembl.vep.corpus.registry import build_corpus
from ensembl.vep.corpus.spec import CorpusSpec
from ensembl.vep.corpus.synthetic.default import default_corpus
from ensembl.vep.corpus.synthetic.heavy import heavy_corpus


def _alternates(corpus: Corpus) -> list[tuple[str, str | None]]:
    """Return each overlap's name and alternate allele, in order."""
    return [(name, overlap.overlapping_feature.alternate_allele) for name, overlap in corpus.named_overlaps]


def _positions(corpus: Corpus) -> list[tuple[str, int]]:
    """Return each overlap's name and start, in order."""
    return [(name, overlap.overlapping_feature.start) for name, overlap in corpus.named_overlaps]


class TestSpec:
    """The value itself."""

    def test_defaults_are_the_generator_defaults(self) -> None:
        """A default spec builds the corpus the fixtures have always used."""
        spec = CorpusSpec()
        assert (spec.genome_seed, spec.boundary_seed, spec.fuzz_seed) == (7, 11, 23)
        assert (spec.heavy_genome_seed, spec.sweep_seed) == (31, 41)

    def test_from_options_reads_attributes_by_name(self) -> None:
        """Anything shaped like the corpus section will do; a missing attribute keeps its default."""

        class _Options:
            kind = "heavy"
            step = 600
            sweep_seed = 99

        spec = CorpusSpec.from_options(_Options())
        assert (spec.kind, spec.step, spec.sweep_seed) == ("heavy", 600, 99)
        assert spec.fuzz == 4000

    def test_kinds_expands_legacy_synthetic_names(self) -> None:
        """The compatibility aliases still expose their original members."""
        assert CorpusSpec(kind="both").kinds == ("default", "heavy")
        assert CorpusSpec(kind="heavy").kinds == ("heavy",)

    def test_registry_owns_corpus_names(self) -> None:
        """A spec accepts any non-empty factory name, including undotted names."""
        assert CorpusSpec(kind="custom").kind == "custom"
        with pytest.raises(ValueError, match="must not be empty"):
            CorpusSpec(kind="")

    def test_generic_parameters_are_factory_owned(self) -> None:
        """Synthetic parameter names have no global meaning for another factory."""
        spec = CorpusSpec(kind="custom", parameters={"fuzz": -1, "step": "wide"})
        assert spec.parameters == {"fuzz": -1, "step": "wide"}

    @pytest.mark.parametrize("field,value", [("fuzz", -1), ("step", 0), ("step", -1)])
    def test_sizes_are_validated_at_construction(self, field: str, value: int) -> None:
        """Invalid sizes fail at the component boundary, not inside ``range`` later."""
        with pytest.raises(ValueError, match=field):
            if field == "fuzz":
                CorpusSpec(fuzz=value)
            else:
                CorpusSpec(step=value)

    def test_seeds_are_reported_only_for_built_corpora(self) -> None:
        """A heavy-only spec does not claim boundary or fuzz seeds."""
        assert set(CorpusSpec(kind="heavy").seeds_by_kind()) == {"heavy"}
        assert CorpusSpec(kind="both").seeds_by_kind()["default"] == {"genome": 7, "boundary": 11, "fuzz": 23}

    def test_the_parameter_hash_is_stable_and_seed_sensitive(self) -> None:
        """Same parameters, same hash; any seed change, a different one."""
        assert CorpusSpec().parameter_hash() == CorpusSpec().parameter_hash()
        assert CorpusSpec().parameter_hash() != CorpusSpec(fuzz_seed=24).parameter_hash()

    def test_the_legacy_default_parameter_hash_is_preserved(self) -> None:
        """Existing benchmark baselines remain comparable after the spec generalisation."""
        assert CorpusSpec().parameter_hash() == "243a78d0956960d1"

    def test_non_synthetic_parameters_are_not_polluted_by_generator_defaults(self) -> None:
        """A generic spec carries only the parameters owned by its registered factory."""
        spec = CorpusSpec(kind="loaded.files", parameters={"fasta": "corpus.fa"})
        assert spec.parameters == {"fasta": "corpus.fa"}

    def test_generic_options_decode_factory_parameters(self) -> None:
        """Registered factories receive repeatable command parameters with scalar types."""

        class _Options:
            kind = "loaded.files"
            parameters = ("fasta=reference.fa", "window=250", "strict=true")
            corpus_dir = None
            fuzz = 99

        assert CorpusSpec.from_options(_Options()).parameters == {
            "fasta": "reference.fa",
            "strict": True,
            "window": 250,
        }

    def test_loaded_bundle_directory_expands_to_standard_file_names(self) -> None:
        """The existing corpus directory option can select a coordinated file bundle."""

        class _Options:
            kind = "loaded.files"
            parameters: tuple[str, ...] = ()
            corpus_dir = Path("bundle")

        assert CorpusSpec.from_options(_Options()).parameters == {
            "fasta": "bundle/corpus.fa",
            "gff3": "bundle/corpus.gff3",
            "vcf": "bundle/corpus.vcf",
        }

    def test_a_corpus_preserves_its_spec(self) -> None:
        """The registry passes its specification through to the built corpus."""
        spec = CorpusSpec(kind="default", fuzz=5, boundary_seed=12)
        assert build_corpus(spec).spec == spec


class TestSeedsReachTheirStreams:
    """Each seed moves its own stream, and only that one."""

    def test_the_default_seeds_change_nothing(self) -> None:
        """Forwarding the default seeds builds what the unseeded generator builds."""
        seeded = build_corpus(CorpusSpec(kind="default", fuzz=20))
        assert _positions(seeded) == [
            (name, overlap.overlapping_feature.start)
            for name, overlap in default_corpus(fuzz_count=20).named_overlaps
        ]
        heavy = build_corpus(CorpusSpec(kind="heavy", step=600))
        assert _alternates(heavy) == [
            (name, overlap.overlapping_feature.alternate_allele)
            for name, overlap in heavy_corpus(step=600).named_overlaps
        ]

    def test_the_genome_seed_changes_the_chromosome(self) -> None:
        """Only the sequence underneath the variants moves."""
        before = build_corpus(CorpusSpec(kind="default", fuzz=0)).sequence_regions[0].seq
        after = build_corpus(CorpusSpec(kind="default", fuzz=0, genome_seed=8)).sequence_regions[0].seq
        assert before != after

    def test_the_boundary_seed_changes_the_boundary_alternates(self) -> None:
        """Boundary SNVs draw their alternate base from the boundary stream."""
        before = build_corpus(CorpusSpec(kind="default", fuzz=0))
        after = build_corpus(CorpusSpec(kind="default", fuzz=0, boundary_seed=12))
        assert [name for name, _ in _alternates(before)] == [name for name, _ in _alternates(after)]
        assert _alternates(before) != _alternates(after)
        assert before.sequence_regions[0].seq == after.sequence_regions[0].seq

    def test_the_fuzz_seed_changes_the_fuzz_positions_only(self) -> None:
        """The boundary variants stay put; the random ones move."""
        before = _positions(build_corpus(CorpusSpec(kind="default", fuzz=30)))
        after = _positions(build_corpus(CorpusSpec(kind="default", fuzz=30, fuzz_seed=24)))
        boundary_before = [entry for entry in before if not entry[0].startswith("fuzz")]
        boundary_after = [entry for entry in after if not entry[0].startswith("fuzz")]
        assert boundary_before == boundary_after
        assert [entry for entry in before if entry[0].startswith("fuzz")] != [
            entry for entry in after if entry[0].startswith("fuzz")
        ]

    def test_the_sweep_seed_changes_the_heavy_locus_alternates(self) -> None:
        """The sweep's SNVs draw from the sweep stream, which used to ignore its seed."""
        before = build_corpus(CorpusSpec(kind="heavy", step=600))
        after = build_corpus(CorpusSpec(kind="heavy", step=600, sweep_seed=42))
        assert _alternates(before) != _alternates(after)
