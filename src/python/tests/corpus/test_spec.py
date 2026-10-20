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
"""Tests for corpus specifications and synthetic parameters."""

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

    def test_default_uses_the_canonical_combined_name(self) -> None:
        """The default specification expands to both canonical synthetic corpora."""
        specification = CorpusSpec()
        assert specification.kind == "synthetic.all"
        assert specification.kinds == ("synthetic.default", "synthetic.heavy")

    def test_parameters_are_returned_as_a_copy(self) -> None:
        """Callers cannot mutate a specification through the public mapping."""
        specification = CorpusSpec(kind="custom", parameters={"seed": 7})
        parameters = specification.parameters
        parameters["seed"] = 8
        assert specification.parameters == {"seed": 7}

    def test_defaults_are_the_generator_defaults(self) -> None:
        """A default spec builds the corpus the fixtures have always used."""
        spec = CorpusSpec()
        assert (spec.genome_seed, spec.boundary_seed, spec.fuzz_seed) == (7, 11, 23)
        assert (spec.heavy_genome_seed, spec.sweep_seed) == (31, 41)

    def test_from_options_reads_attributes_by_name(self) -> None:
        """Anything shaped like the corpus section will do; a missing attribute keeps its default."""

        class _Options:
            kind = "synthetic.heavy"
            step = 600
            sweep_seed = 99

        spec = CorpusSpec.from_options(_Options())
        assert (spec.kind, spec.step, spec.sweep_seed) == ("synthetic.heavy", 600, 99)
        assert spec.fuzz == 4000

    def test_kinds_expands_combined_synthetic_corpus(self) -> None:
        """The combined corpus exposes both concrete registered members."""
        assert CorpusSpec(kind="synthetic.all").kinds == (
            "synthetic.default",
            "synthetic.heavy",
        )
        assert CorpusSpec(kind="synthetic.heavy").kinds == ("synthetic.heavy",)

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
        assert set(CorpusSpec(kind="synthetic.heavy").seeds_by_kind()) == {"synthetic.heavy"}
        assert CorpusSpec(kind="synthetic.all").seeds_by_kind()["synthetic.default"] == {
            "genome": 7,
            "boundary": 11,
            "fuzz": 23,
        }

    def test_the_parameter_hash_is_stable_and_seed_sensitive(self) -> None:
        """Same parameters, same hash; any seed change, a different one."""
        assert CorpusSpec().parameter_hash() == CorpusSpec().parameter_hash()
        assert CorpusSpec().parameter_hash() != CorpusSpec(fuzz_seed=24).parameter_hash()

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

    def test_generic_mapping_options_preserve_scalar_parameters(self) -> None:
        """Configuration mappings retain each supported scalar type."""
        options = type(
            "Options",
            (),
            {
                "kind": "custom",
                "parameters": {"text": "value", "integer": 7, "float": 1.5, "flag": True, "empty": None},
            },
        )()

        assert CorpusSpec.from_options(options).parameters == {
            "text": "value",
            "integer": 7,
            "float": 1.5,
            "flag": True,
            "empty": None,
        }

    @pytest.mark.parametrize(
        ("configured", "exception", "message"),
        [
            pytest.param({1: "value"}, TypeError, "names must be strings", id="non-string-name"),
            pytest.param({"value": []}, TypeError, "must be a scalar", id="non-scalar-value"),
            pytest.param(("missing-separator",), ValueError, "KEY=VALUE", id="malformed-assignment"),
            pytest.param(("=value",), ValueError, "needs a name", id="empty-name"),
            pytest.param(("labels=[1, 2]",), ValueError, "must be scalars", id="non-scalar-json"),
            pytest.param(7, TypeError, "mapping or a sequence", id="unsupported-parameters"),
        ],
    )
    def test_generic_options_reject_invalid_parameters(
        self,
        configured: object,
        exception: type[Exception],
        message: str,
    ) -> None:
        """Reject malformed, non-scalar, and unsupported generic parameter configuration."""
        options = type("Options", (), {"kind": "custom", "parameters": configured})()

        with pytest.raises(exception, match=message):
            CorpusSpec.from_options(options)

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

    def test_golden_bundle_directory_becomes_a_loader_parameter(self) -> None:
        """The golden loader receives the configured bundle directory without synthetic defaults."""

        class _Options:
            kind = "golden.file"
            corpus_dir = Path("golden-bundle")

        assert CorpusSpec.from_options(_Options()).parameters == {"directory": "golden-bundle"}

    def test_from_dict_supports_legacy_synthetic_parameters(self) -> None:
        """Older flat representations retain their explicit synthetic overrides."""
        specification = CorpusSpec.from_dict({"kind": "synthetic.default", "fuzz": 12})
        assert specification.fuzz == 12

    @pytest.mark.parametrize(
        ("value", "message"),
        [
            pytest.param({"kind": 1}, "kind must be a string", id="non-string-kind"),
            pytest.param({"parameters": ()}, "parameters must be a mapping", id="non-mapping-parameters"),
        ],
    )
    def test_from_dict_rejects_invalid_structure(self, value: dict[str, object], message: str) -> None:
        """Deserialised specifications validate their kind and parameter container."""
        with pytest.raises(TypeError, match=message):
            CorpusSpec.from_dict(value)

    def test_with_parameters_returns_an_updated_copy(self) -> None:
        """Replacing parameters leaves the original frozen specification unchanged."""
        original = CorpusSpec(kind="custom", parameters={"window": 100, "strict": False})
        updated = original.with_parameters(window=250, strict=True)

        assert original.parameters == {"window": 100, "strict": False}
        assert updated.parameters == {"window": 250, "strict": True}

    def test_seed_properties_reject_boolean_values(self) -> None:
        """Seed accessors do not treat booleans as integer seeds."""
        with pytest.raises(TypeError, match="genome_seed.*integer"):
            _ = CorpusSpec(genome_seed=True).genome_seed

    @pytest.mark.parametrize(
        "parameters",
        [
            pytest.param({1: "value"}, id="non-string-name"),
            pytest.param({"value": []}, id="non-scalar-value"),
        ],
    )
    def test_constructor_rejects_invalid_factory_parameters(self, parameters: dict[object, object]) -> None:
        """Factories can only receive string names and scalar values."""
        with pytest.raises(TypeError):
            CorpusSpec(kind="custom", parameters=parameters)  # type: ignore[arg-type]

    def test_a_corpus_preserves_its_spec(self) -> None:
        """The registry passes its specification through to the built corpus."""
        spec = CorpusSpec(kind="synthetic.default", fuzz=5, boundary_seed=12)
        assert build_corpus(spec).spec == spec


class TestSeedsReachTheirStreams:
    """Each seed moves its own stream, and only that one."""

    def test_the_default_seeds_change_nothing(self) -> None:
        """Forwarding the default seeds builds what the unseeded generator builds."""
        seeded = build_corpus(CorpusSpec(kind="synthetic.default", fuzz=20))
        assert _positions(seeded) == [
            (name, overlap.overlapping_feature.start)
            for name, overlap in default_corpus(fuzz_count=20).named_overlaps
        ]
        heavy = build_corpus(CorpusSpec(kind="synthetic.heavy", step=600))
        assert _alternates(heavy) == [
            (name, overlap.overlapping_feature.alternate_allele)
            for name, overlap in heavy_corpus(step=600).named_overlaps
        ]

    def test_the_genome_seed_changes_the_chromosome(self) -> None:
        """Only the sequence underneath the variants moves."""
        before = build_corpus(CorpusSpec(kind="synthetic.default", fuzz=0)).sequence_regions[0].seq
        after = (
            build_corpus(CorpusSpec(kind="synthetic.default", fuzz=0, genome_seed=8)).sequence_regions[0].seq
        )
        assert before != after

    def test_the_boundary_seed_changes_the_boundary_alternates(self) -> None:
        """Boundary SNVs draw their alternate base from the boundary stream."""
        before = build_corpus(CorpusSpec(kind="synthetic.default", fuzz=0))
        after = build_corpus(CorpusSpec(kind="synthetic.default", fuzz=0, boundary_seed=12))
        assert [name for name, _ in _alternates(before)] == [name for name, _ in _alternates(after)]
        assert _alternates(before) != _alternates(after)
        assert before.sequence_regions[0].seq == after.sequence_regions[0].seq

    def test_the_fuzz_seed_changes_the_fuzz_positions_only(self) -> None:
        """The boundary variants stay put; the random ones move."""
        before = _positions(build_corpus(CorpusSpec(kind="synthetic.default", fuzz=30)))
        after = _positions(build_corpus(CorpusSpec(kind="synthetic.default", fuzz=30, fuzz_seed=24)))
        boundary_before = [entry for entry in before if not entry[0].startswith("fuzz")]
        boundary_after = [entry for entry in after if not entry[0].startswith("fuzz")]
        assert boundary_before == boundary_after
        assert [entry for entry in before if entry[0].startswith("fuzz")] != [
            entry for entry in after if entry[0].startswith("fuzz")
        ]

    def test_the_sweep_seed_changes_the_heavy_locus_alternates(self) -> None:
        """The sweep's SNVs draw from the sweep stream, which used to ignore its seed."""
        before = build_corpus(CorpusSpec(kind="synthetic.heavy", step=600))
        after = build_corpus(CorpusSpec(kind="synthetic.heavy", step=600, sweep_seed=42))
        assert _alternates(before) != _alternates(after)
