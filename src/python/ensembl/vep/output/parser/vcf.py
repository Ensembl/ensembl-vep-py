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
"""Annotated VCF consequence parsing.

VEP writes them as ``CSQ`` and SnpEff as ``ANN``: one INFO field whose layout the header's
``Description`` declares, one entry per (allele, feature) separated by commas, fields by pipes,
several terms by ``&``. Both are read the same way; only the field names differ.

Two mechanics that must stay right, both confirmed on this codebase:

* ``seed_attrs`` takes the *public* name and prepends the underscore itself:
  ``seed_attrs(overlap_consequences=...)``, never ``_overlap_consequences``.
* ``transcript_variation_alleles`` carries ``on_setattr=invalidates("_overlap_consequences")``, so
  assigning it after seeding silently wipes the seed. This reader never assigns it, which is also
  why it needs no ``--gff``/``--fasta``: the terms live on the allele, and
  ``Variation.most_severe_consequence`` and ``.class_so_term`` answer with zero transcripts.
"""

from __future__ import annotations

import re
from collections.abc import Iterator, Sequence
from pathlib import Path

from attrs import define, field

from ensembl.vep.annotators.consequence.rules.rulesets import FactRuleSet, build_fact_ruleset
from ensembl.vep.common import CoordSystemLevel, Locus, Strand
from ensembl.vep.common.constants import EMPTY_ALLELE
from ensembl.vep.models.overlap import OverlapConsequence
from ensembl.vep.models.sequence import CoordSystem, Slice
from ensembl.vep.models.variation import Variation
from ensembl.vep.output.parser.base import OutputParser
from ensembl.vep.parsers.vcf import VCFParser, VCFRecord

__all__ = ["ANNOTATION_FIELDS", "UNKNOWN_RANK", "VCFOutputParser", "csq_alleles", "parser_for"]

ANNOTATION_FIELDS = {"CSQ": ("Allele", "Consequence"), "ANN": ("Allele", "Annotation")}
"""INFO key -> the names of its allele and consequence fields, as each tool spells them."""

UNKNOWN_RANK = 999
"""Rank given to a term the catalogue does not model. Rank 1 is the most severe and the model's
default is 0, so an unranked term must sit *below* every real one rather than above."""

_FORMAT_PATTERN = re.compile(r"Format:\s*(.+)$")


def annotation_layouts(metadata: dict) -> dict[str, list[str]]:
    """Return each annotation INFO key present in a header with the field names it declares.

    Args:
        metadata (dict): ``VCFParser.metadata``.

    Returns:
        dict[str, list[str]]: INFO key -> field names, for every key in ``ANNOTATION_FIELDS`` the
            header declares.
    """
    layouts: dict[str, list[str]] = {}
    for key in ANNOTATION_FIELDS:
        declaration = metadata.get("info", {}).get(key)
        if declaration is None:
            continue
        description = declaration.get("description", "")
        match = _FORMAT_PATTERN.search(description)
        # VEP says "Format: A|B|C"; SnpEff quotes the list: 'A | B | C'.
        field_text = match.group(1) if match else description.split("'")[1] if "'" in description else ""
        layouts[key] = [name.strip() for name in field_text.strip("'\" ").split("|")]
    return layouts


def csq_alleles(reference: str, alternates: Sequence[str]) -> list[str]:
    """Return how VEP spells each alternate allele in the CSQ ``Allele`` field.

    An SNV is its base. When every allele shares the reference's first base -- an indel with its
    VCF padding -- that base is dropped and an allele left empty becomes ``-``.
    """
    shares_padding = bool(reference) and all(
        alternate and alternate[0] == reference[0] for alternate in alternates
    )
    if shares_padding and (len(reference) > 1 or any(len(alternate) > 1 for alternate in alternates)):
        return [alternate[1:] or EMPTY_ALLELE for alternate in alternates]
    return list(alternates)


@define(kw_only=True, auto_attribs=True)
class VCFOutputParser(OutputParser):
    """Parser for annotated VCF files.

    Attributes:
        ruleset (FactRuleSet): Where a known term's rank, accession and impact come from. A term
            the catalogue does not model is kept, with ``UNKNOWN_RANK``.
    """

    ruleset: FactRuleSet = field(factory=build_fact_ruleset)

    def parse(self, path: str | Path) -> Iterator[Variation]:
        """Yield one variation per record.

        Args:
            path (str | Path): The annotated VCF.

        Yields:
            Variation: Each record, its alleles seeded with the terms the file assigned them.
        """
        parser = VCFParser()
        with parser.open(path):
            layouts = annotation_layouts(parser.metadata)
            while (record := parser.next()) is not None:
                yield self._variation(record, layouts)

    def _variation(self, record: VCFRecord, layouts: dict[str, list[str]]) -> Variation:
        """Build one variation from a record and seed its alleles' consequences."""
        reference = record.reference
        alternates = list(record.alt)
        end = record.pos + len(reference) - 1
        chromosome = Slice(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name=record.chrom,
            start=1,
            end=end,
            strand=Strand.FORWARD,
        )
        variation = Variation(
            start=record.pos,
            end=end,
            strand=Strand.FORWARD,
            slice=chromosome,
            stable_id=record.id[0] if record.id else None,
            allele_string="/".join([reference, *alternates]),
            # The record read back is the uploaded one: nothing resolved its name on the way here,
            # so the file's own spelling is both the slice's and the uploaded locus's.
            uploaded_locus=Locus(
                seq_region_name=record.chrom,
                start=record.pos,
                end=end,
                allele_string="/".join([reference, *alternates]),
            ),
        )

        terms_by_index = self._terms_by_allele_index(record, layouts, reference, alternates)
        for index, allele in enumerate(variation.variation_alleles):
            consequences = [self._consequence(term) for term in terms_by_index.get(index, [])]
            allele.seed_attrs(overlap_consequences=consequences)
        return variation

    def _terms_by_allele_index(
        self, record: VCFRecord, layouts: dict[str, list[str]], reference: str, alternates: list[str]
    ) -> dict[int, list[str]]:
        """Return the terms each alternate allele carries, by its position in ``ALT``."""
        terms: dict[int, list[str]] = {}
        for key, names in layouts.items():
            raw = record.info.get(key)
            if not raw:
                continue
            allele_field, term_field = ANNOTATION_FIELDS[key]
            if allele_field not in names or term_field not in names:
                continue
            allele_index, term_index = names.index(allele_field), names.index(term_field)
            spellings = csq_alleles(reference, alternates) if key == "CSQ" else alternates
            for entry in str(raw).split(","):
                fields = entry.split("|")
                if len(fields) <= max(allele_index, term_index):
                    continue
                if fields[allele_index] not in spellings:
                    continue
                position = spellings.index(fields[allele_index])
                terms.setdefault(position, []).extend(term for term in fields[term_index].split("&") if term)
        return terms

    def _consequence(self, term: str) -> OverlapConsequence:
        """Return the catalogue's consequence for a term, or an unranked stand-in."""
        rule = self.ruleset.by_term.get(term)
        if rule is not None:
            return rule.to_overlap_consequence()
        return OverlapConsequence(term=term, accession="", rank=UNKNOWN_RANK, impact=None)


def parser_for(path: str | Path, input_format: str | None = None) -> OutputParser:
    """Return the reader for a file, by its declared format or its extension.

    Raises:
        ValueError: If the format is not one this package reads.
    """
    # Local imports: the seams import this module's base, and this is the one place all three meet.
    # pylint: disable=import-outside-toplevel
    from ensembl.vep.output.parser.json import JSONOutputParser
    from ensembl.vep.output.parser.tab import TabOutputParser

    name = str(path).lower()
    if name.endswith(".gz"):
        name = name[:-3]
    chosen = input_format or name.rsplit(".", 1)[-1]
    readers: dict[str, type[OutputParser]] = {
        "vcf": VCFOutputParser,
        "tab": TabOutputParser,
        "txt": TabOutputParser,
        "json": JSONOutputParser,
    }
    if chosen not in readers:
        raise ValueError(f"no reader for {chosen!r} output; expected one of {sorted(readers)}")
    return readers[chosen]()
