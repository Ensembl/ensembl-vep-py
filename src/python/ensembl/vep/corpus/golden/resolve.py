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
"""Resolve a loaded golden corpus or explicitly generate its fallback."""

from __future__ import annotations

from contextlib import AbstractContextManager
from importlib.resources import as_file, files
from pathlib import Path

from ensembl.vep.corpus.base import Corpus
from ensembl.vep.corpus.golden.base import GoldenCorpus
from ensembl.vep.corpus.golden.terms import GoldenTerms, overlap_key
from ensembl.vep.corpus.registry import build_corpus, register_corpus
from ensembl.vep.corpus.spec import CorpusSpec

__all__ = ["load_golden_corpus"]

_PACKAGE = "ensembl.vep.corpus.golden"
_DEFAULT_DATASET = "default"


def _loaded_spec(directory: Path) -> CorpusSpec:
    """Return the file-loader specification for a golden directory."""
    return CorpusSpec(
        kind="loaded.files",
        parameters={
            "fasta": str(directory / "corpus.fa"),
            "gff3": str(directory / "corpus.gff3"),
            "vcf": str(directory / "corpus.vcf"),
        },
    )


def _identity_spec(specification: CorpusSpec) -> CorpusSpec:
    """Return dataset identity without directory or generation policy parameters."""
    return CorpusSpec(
        kind="golden.file",
        parameters={"dataset": specification.parameter("dataset", _DEFAULT_DATASET)},
    )


def _validate_expected_terms(corpus: Corpus, expected: GoldenTerms) -> None:
    """Require one expected-term entry for every loaded corpus overlap."""
    overlap_keys = {overlap_key(name, overlap) for name, overlap in corpus.named_overlaps}
    term_keys = set(expected.terms)
    missing_terms = sorted(overlap_keys - term_keys)
    extra_terms = sorted(term_keys - overlap_keys)
    if missing_terms or extra_terms:
        details = []
        if missing_terms:
            details.append(f"missing terms for {missing_terms!r}")
        if extra_terms:
            details.append(f"terms for unknown overlaps {extra_terms!r}")
        raise ValueError("golden terms must match loaded corpus overlaps: " + "; ".join(details))


def load_golden_corpus(
    specification: CorpusSpec,
    *,
    directory: Path | None = None,
) -> GoldenCorpus:
    """Load pinned corpus files and terms.

    This corpus-layer operation deliberately has no annotator dependency. Recording terms and
    generating a missing fallback belong to the annotator integration that consumes this model.

    Args:
        specification (CorpusSpec): Identity of the requested golden dataset.
        directory (Path | None): Explicit bundle directory, or ``None`` for packaged data.
    Returns:
        GoldenCorpus: Loaded corpus with pinned expected terms.

    Raises:
        FileNotFoundError: If required files are absent.
        TypeError: If the dataset parameter has an invalid type.
        ValueError: If loaded terms identify a different corpus specification.
    """
    identity_specification = _identity_spec(specification)
    resource_context: AbstractContextManager[object] | None = None
    if directory is None:
        dataset = identity_specification.parameter("dataset", _DEFAULT_DATASET)
        if not isinstance(dataset, str):
            raise TypeError("golden corpus dataset must be a string")
        resource_context = as_file(files(_PACKAGE).joinpath("data", dataset))
        directory = Path(resource_context.__enter__())  # pylint: disable=unnecessary-dunder-call

    if not directory.is_dir():
        if resource_context is not None:
            resource_context.__exit__(None, None, None)
        raise FileNotFoundError(directory)

    terms_path = directory / "terms.json"
    required = (directory / "corpus.fa", directory / "corpus.gff3", directory / "corpus.vcf", terms_path)
    missing = [path for path in required if not path.is_file()]
    if missing:
        if resource_context is not None:
            resource_context.__exit__(None, None, None)
        raise FileNotFoundError(missing[0])

    corpus: Corpus | None = None
    try:
        corpus = build_corpus(_loaded_spec(directory))
        expected = GoldenTerms.load(terms_path)
        if expected.spec != identity_specification:
            raise ValueError(
                f"golden bundle identifies {expected.spec.as_dict()!r}, "
                f"not requested {identity_specification.as_dict()!r}"
            )
        _validate_expected_terms(corpus, expected)
    except Exception:
        if corpus is not None:
            corpus.close()
        if resource_context is not None:
            resource_context.__exit__(None, None, None)
        raise
    return GoldenCorpus(
        spec=expected.spec,
        corpus=corpus,
        expected_terms=expected,
        resource_context=resource_context,
    )


@register_corpus("golden.file")
def _registered_golden_corpus(specification: CorpusSpec) -> Corpus:
    """Build the registered packaged golden corpus."""
    directory_value = specification.parameter("directory")
    directory = Path(directory_value) if isinstance(directory_value, str) else None
    return load_golden_corpus(
        _identity_spec(specification),
        directory=directory,
    )
