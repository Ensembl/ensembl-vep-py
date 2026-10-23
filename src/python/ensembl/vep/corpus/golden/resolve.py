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
from ensembl.vep.corpus.golden.terms import GoldenTerms, record_golden
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


def _build_golden_from_synthetic(
    specification: CorpusSpec, synthetic_specification: CorpusSpec
) -> GoldenCorpus:
    """Build an in-memory golden corpus from a synthetic corpus specification."""
    # Import lazily so registering the golden corpus does not couple the corpus package import to
    # annotator discovery. Generation is an explicit fallback, not part of ordinary loading.
    from ensembl.vep.annotators.factory import build_annotator  # pylint: disable=import-outside-toplevel

    corpus = build_corpus(synthetic_specification)
    try:
        annotator = build_annotator(annotator="consequence.linear")
        recorded = record_golden(corpus, annotator, "consequence.linear")
    except Exception:
        corpus.close()
        raise
    expected = GoldenTerms(
        spec=specification,
        annotator=recorded.annotator,
        catalogue=recorded.catalogue,
        terms=recorded.terms,
    )
    return GoldenCorpus(spec=specification, corpus=corpus, expected_terms=expected)


def load_golden_corpus(
    specification: CorpusSpec,
    *,
    directory: Path | None = None,
    fallback: CorpusSpec | None = None,
    allow_generate: bool = False,
) -> GoldenCorpus:
    """Load pinned corpus files and terms, or explicitly generate a missing fallback.

    Existing but malformed data always raises. Generation is limited to absent data so a corrupt
    golden reference cannot silently validate against newly generated expectations.

    Args:
        specification (CorpusSpec): Identity of the requested golden dataset.
        directory (Path | None): Explicit bundle directory, or ``None`` for packaged data.
        fallback (CorpusSpec | None): Corpus to annotate when missing data may be generated.
        allow_generate (bool): Whether an absent dataset may use ``fallback``.

    Returns:
        GoldenCorpus: Loaded or generated corpus with pinned expected terms.

    Raises:
        FileNotFoundError: If required corpus or term files are absent without an allowed fallback.
        TypeError: If the dataset parameter has an invalid type.
        ValueError: If loaded terms identify a different corpus specification.
    """
    resource_context: AbstractContextManager[object] | None = None
    if directory is None:
        dataset = specification.parameter("dataset", _DEFAULT_DATASET)
        if not isinstance(dataset, str):
            raise TypeError("golden corpus dataset must be a string")
        resource_context = as_file(files(_PACKAGE).joinpath("data", dataset))
        try:
            # The context has to remain open for the lifetime of the returned GoldenCorpus.
            directory = Path(resource_context.__enter__())  # pylint: disable=unnecessary-dunder-call
        except FileNotFoundError:
            if not allow_generate or fallback is None:
                raise
            return _build_golden_from_synthetic(specification, fallback)

    if not directory.is_dir():
        if resource_context is not None:
            resource_context.__exit__(None, None, None)
        if allow_generate and fallback is not None:
            return _build_golden_from_synthetic(specification, fallback)
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
        # Load the corpus files through the regular file-backed factory
        corpus = build_corpus(_loaded_spec(directory))
        expected = GoldenTerms.load(terms_path)
        if expected.spec != specification:
            raise ValueError(
                f"golden bundle identifies {expected.spec.as_dict()!r}, "
                f"not requested {specification.as_dict()!r}"
            )
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
    allow_generate = specification.parameter("allow_generate", False) is True
    fallback_kind = specification.parameter("fallback", "synthetic.default")
    fallback = None
    if isinstance(fallback_kind, str):
        fallback = CorpusSpec(kind=fallback_kind, parameters=specification.parameters)
    return load_golden_corpus(
        _identity_spec(specification),
        directory=directory,
        fallback=fallback,
        allow_generate=allow_generate,
    )
