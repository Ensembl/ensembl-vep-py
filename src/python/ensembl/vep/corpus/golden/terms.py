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
"""Pinned consequence terms for a golden corpus.

The serialised terms are part of the corpus layer. Recording new terms is an annotator concern and
is intentionally implemented by the consumer-side integration rather than this package.
"""

from __future__ import annotations

import json
from pathlib import Path
from attrs import field, frozen

from ensembl.vep.corpus.spec import CorpusSpec
from ensembl.vep.models.overlap import TranscriptVariationAllele

__all__ = ["FORMAT_VERSION", "GoldenReport", "GoldenTerms", "overlap_key"]

FORMAT_VERSION = 1
"""Bumped whenever the file layout changes, so an old file is rejected rather than misread."""

_REPORT_LIMIT = 20


def overlap_key(name: str, overlap: TranscriptVariationAllele) -> str:
    """Return the stable ``feature::variant`` key for an overlap.

    Args:
        name (str): Corpus-level variant name.
        overlap (TranscriptVariationAllele): Transcript overlap to identify.

    Returns:
        str: Key containing both feature and variant identities.
    """
    return f"{overlap.feature.stable_id}::{name}"


@frozen(kw_only=True, auto_attribs=True)
class GoldenTerms:
    """The terms one annotator emitted for every overlap of one corpus.

    Attributes:
        spec (CorpusSpec): The corpus the terms were recorded over.
        annotator (str): The annotator that emitted them.
        catalogue (str): Digest of the catalogue they were emitted against.
        terms (dict[str, tuple[str, ...]]): ``feature::name`` to emitted terms, in order.
    """

    spec: CorpusSpec
    annotator: str
    catalogue: str
    terms: dict[str, tuple[str, ...]] = field(factory=dict)

    def to_json(self) -> str:
        """Return deterministic file text with one overlap per line."""
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
        """Read terms from JSON text.

        Args:
            text (str): Serialised golden terms.

        Returns:
            GoldenTerms: Deserialised terms.

        Raises:
            ValueError: If the file uses an unsupported format version.
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
        """Write the terms, creating the destination directory when needed."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.to_json(), encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> GoldenTerms:
        """Read terms from a file."""
        return cls.from_json(path.read_text(encoding="utf-8"))


@frozen(kw_only=True, auto_attribs=True)
class GoldenReport:
    """Describe how recorded terms differ from golden terms.

    Attributes:
        catalogue_changed (bool): Whether catalogue provenance changed.
        spec_changed (bool): Whether the recording used another corpus.
        annotator_changed (bool): Whether another annotator made the recording.
        missing (tuple[str, ...]): Keys absent from the recording.
        unexpected (tuple[str, ...]): Keys absent from the golden terms.
        mismatches (tuple[tuple[str, tuple[str, ...], tuple[str, ...]], ...]): Key, expected
            terms, and recorded terms for changed overlaps.
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
        """Return whether terms match against the same catalogue and corpus."""
        return not (
            self.catalogue_changed or self.spec_changed or self.missing or self.unexpected or self.mismatches
        )

    def lines(self, limit: int = _REPORT_LIMIT) -> list[str]:
        """Return a readable report with at most ``limit`` mismatch examples."""
        lines = []
        if self.catalogue_changed:
            lines.append(
                "the catalogue changed since the golden file was recorded; if that change is intended, "
                "re-record and review the diff"
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
