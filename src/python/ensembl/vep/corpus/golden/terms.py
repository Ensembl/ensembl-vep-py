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
"""Frozen expected terms emitted for a corpus and pinned to a file.

Conformance proves every annotator agrees with ``consequence.linear``. It cannot catch all of them
agreeing on the same wrong answer; a golden file can. It freezes the terms one annotator produced
for one corpus, together with the corpus spec and the catalogue that produced them, so a later run
is checked against a fixed expectation rather than only against itself.

Named ``golden`` on purpose: in this repository "oracle" already means the external Perl
``VariationEffect.pm``, whose whole value is that it is *not* us. This is the opposite thing.

A mismatch is reported three distinct ways, because each calls for a different response. The
catalogue changed under the file -- regenerate it, deliberately. The file was recorded for another
corpus -- it is the wrong file. Or the terms moved with neither changed -- a regression.
"""

from __future__ import annotations

import json
from pathlib import Path

from attrs import field, frozen

from ensembl.vep.annotators.base import Annotator
from ensembl.vep.annotators.consequence.rules.catalogue import catalogue_digest
from ensembl.vep.corpus.base import Corpus
from ensembl.vep.corpus.spec import CorpusSpec
from ensembl.vep.models.overlap import OverlapConsequence, TranscriptVariationAllele

__all__ = ["FORMAT_VERSION", "GoldenReport", "GoldenTerms", "overlap_key", "record_golden"]

FORMAT_VERSION = 1
"""Bumped whenever the file layout changes, so an old file is rejected rather than misread."""

_REPORT_LIMIT = 20


def overlap_key(name: str, overlap: TranscriptVariationAllele) -> str:
    """Return the key an overlap is filed under: its feature, then its variant's name.

    A name identifies the variant, not the overlap -- one variant pairs with every feature in its
    flank window -- so the feature has to be part of the key.
    """
    return f"{overlap.feature.stable_id}::{name}"


@frozen(kw_only=True, auto_attribs=True)
class GoldenTerms:
    """The terms one annotator emitted for every overlap of one corpus.

    Attributes:
        spec (CorpusSpec): The corpus the terms were recorded over.
        annotator (str): The annotator that emitted them.
        catalogue (str): Digest of the catalogue they were emitted against.
        terms (dict[str, tuple[str, ...]]): ``feature::name`` -> the emitted terms, in order.
    """

    spec: CorpusSpec
    annotator: str
    catalogue: str
    terms: dict[str, tuple[str, ...]] = field(factory=dict)

    def to_json(self) -> str:
        """Return the file text: a header, then one overlap per line.

        One line per overlap is deliberate. Regenerating after a catalogue change then diffs as
        exactly the overlaps whose terms moved, which is the review a regeneration needs.
        """
        header = {
            "format_version": FORMAT_VERSION,
            "spec": self.spec.as_dict(),
            "annotator": self.annotator,
            "catalogue": self.catalogue,
        }
        lines = ["{"]
        for name, value in header.items():
            lines.append(f" {json.dumps(name)}: {json.dumps(value, sort_keys=True)},")
        lines.append(' "terms": {')
        entries = sorted(self.terms.items())
        for index, (key, value) in enumerate(entries):
            separator = "," if index < len(entries) - 1 else ""
            lines.append(f"  {json.dumps(key)}: {json.dumps(list(value))}{separator}")
        lines.append(" }")
        lines.append("}")
        return "\n".join(lines) + "\n"

    @classmethod
    def from_json(cls, text: str) -> GoldenTerms:
        """Read terms back from file text.

        Raises:
            ValueError: If the file was written in a different format version.
        """
        payload = json.loads(text)
        version = payload.get("format_version")
        if version != FORMAT_VERSION:
            raise ValueError(f"golden file format {version!r} is not the supported {FORMAT_VERSION}")
        terms = {key: tuple(value) for key, value in payload["terms"].items()}
        return cls(
            spec=CorpusSpec.from_dict(payload["spec"]),
            annotator=payload["annotator"],
            catalogue=payload["catalogue"],
            terms=terms,
        )

    def write(self, path: Path) -> None:
        """Write the terms to a file, creating its directory if needed."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.to_json(), encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> GoldenTerms:
        """Read terms from a file."""
        return cls.from_json(path.read_text(encoding="utf-8"))


def record_golden(
    corpus: Corpus, annotator: Annotator[OverlapConsequence], annotator_name: str
) -> GoldenTerms:
    """Annotate every overlap of a corpus and freeze the terms.

    Args:
        corpus (Corpus): The corpus to record over.
        annotator (Annotator[OverlapConsequence]): The annotator whose terms are frozen.
        annotator_name (str): Its registry name, recorded in the file.

    Returns:
        GoldenTerms: The frozen terms.

    Raises:
        ValueError: If two overlaps would be filed under the same key.
    """
    named_overlaps = list(corpus.named_overlaps)
    results = annotator.annotate_feature_overlaps([overlap for _name, overlap in named_overlaps])
    terms: dict[str, tuple[str, ...]] = {}
    for (name, overlap), consequences in zip(named_overlaps, results, strict=True):
        key = overlap_key(name, overlap)
        if key in terms:
            raise ValueError(f"two overlaps share the key {key!r}; a golden file needs them unique")
        terms[key] = tuple(consequence.term for consequence in consequences)
    return GoldenTerms(spec=corpus.spec, annotator=annotator_name, catalogue=catalogue_digest(), terms=terms)


@frozen(kw_only=True, auto_attribs=True)
class GoldenReport:
    """How recorded terms differ from golden ones.

    Attributes:
        catalogue_changed (bool): The catalogue is not the one the golden file was recorded against.
        spec_changed (bool): The corpus is not the one the golden file was recorded over.
        annotator_changed (bool): A different annotator produced the recording. Informational --
            checking another strategy against one annotator's golden file is a legitimate use.
        missing (tuple[str, ...]): Keys the golden file has and the recording does not.
        unexpected (tuple[str, ...]): Keys the recording has and the golden file does not.
        mismatches (tuple[tuple[str, tuple[str, ...], tuple[str, ...]], ...]): Key, golden terms,
            recorded terms, for every overlap whose terms moved.
    """

    catalogue_changed: bool = field(default=False)
    spec_changed: bool = field(default=False)
    annotator_changed: bool = field(default=False)
    missing: tuple[str, ...] = field(default=())
    unexpected: tuple[str, ...] = field(default=())
    mismatches: tuple[tuple[str, tuple[str, ...], tuple[str, ...]], ...] = field(default=())

    @classmethod
    def compare(cls, golden: GoldenTerms, recorded: GoldenTerms) -> GoldenReport:
        """Compare a recording against golden terms."""
        golden_keys = set(golden.terms)
        recorded_keys = set(recorded.terms)
        mismatches = []
        for key in sorted(golden_keys & recorded_keys):
            if golden.terms[key] != recorded.terms[key]:
                mismatches.append((key, golden.terms[key], recorded.terms[key]))
        return cls(
            catalogue_changed=golden.catalogue != recorded.catalogue,
            spec_changed=golden.spec != recorded.spec,
            annotator_changed=golden.annotator != recorded.annotator,
            missing=tuple(sorted(golden_keys - recorded_keys)),
            unexpected=tuple(sorted(recorded_keys - golden_keys)),
            mismatches=tuple(mismatches),
        )

    @property
    def ok(self) -> bool:
        """Whether the recording matches, against the same catalogue and corpus."""
        return not (
            self.catalogue_changed or self.spec_changed or self.missing or self.unexpected or self.mismatches
        )

    def lines(self, limit: int = _REPORT_LIMIT) -> list[str]:
        """Return the report as readable lines, the reason first, then at most ``limit`` examples."""
        lines = []
        if self.catalogue_changed:
            lines.append(
                "the catalogue changed since the golden file was recorded; if that change is intended, "
                "re-record with --update_golden and review the diff"
            )
        if self.spec_changed:
            lines.append("the golden file was recorded over a different corpus spec")
        if self.annotator_changed:
            lines.append("note: checking a different annotator than the one that recorded the file")
        if self.missing:
            lines.append(
                f"{len(self.missing)} overlaps missing from the recording, e.g. {list(self.missing[:3])}"
            )
        if self.unexpected:
            lines.append(
                f"{len(self.unexpected)} overlaps not in the golden file, e.g. {list(self.unexpected[:3])}"
            )
        if self.mismatches:
            lines.append(f"{len(self.mismatches)} overlaps emit different terms:")
            for key, golden_terms, recorded_terms in self.mismatches[:limit]:
                lines.append(f"  {key}: golden {list(golden_terms)} != recorded {list(recorded_terms)}")
            if len(self.mismatches) > limit:
                lines.append(f"  ... and {len(self.mismatches) - limit} more")
        if not lines:
            lines.append("recorded terms match the golden file")
        return lines
