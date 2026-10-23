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
"""Tests for the columnar consequence annotator."""

from __future__ import annotations

import pytest

from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.annotators.consequence.rules.expressions import Expr
from ensembl.vep.annotators.consequence.rules.rulesets import build_fact_ruleset
from ensembl.vep.annotators.consequence.strategies.columnar import (
    BitColumnBuffer,
    ColumnarFactsAnnotator,
    RulePlan,
    _AndPlanNode,
    _AtomPlanNode,
    _NotPlanNode,
    _OrPlanNode,
)
from ensembl.vep.caches.codon_cache import CodonEditLRUCache
from ensembl.vep.caches.zone_map import ZoneMapCache
from ensembl.vep.common import Strand
from ensembl.vep.models.genic import Exon, Transcript
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.models.variation import VariationAllele
from ensembl.vep.models.sequence import CoordSystemLevel, CoordSystem, Slice


class _UnsupportedExpression(Expr):
    """Expression node unsupported by the columnar plan compiler."""


class _NotATranscriptOverlap:
    """Object that is not a transcript variation allele."""


class TestRulePlan:
    """Tests for compiling columnar rule plans."""

    def test_unsupported_expression_is_rejected(self) -> None:
        """Test an unknown expression node raises a type error."""
        expression = _UnsupportedExpression()

        with pytest.raises(TypeError, match="unsupported expression node type"):
            RulePlan.build(expression)


class TestColumnarFactsAnnotator:
    """Tests for columnar batch handling and validation."""

    def test_empty_helpers_return_empty_results(self) -> None:
        """Test empty evaluation and tracing chunks are no-ops."""
        annotator = ColumnarFactsAnnotator()

        assert annotator._evaluate_chunk([]) == []  # pylint: disable=protected-access
        assert not annotator._trace_chunk([])  # pylint: disable=protected-access

    def test_non_transcript_overlap_is_rejected(self) -> None:
        """Test a batch containing a non-transcript overlap raises a type error."""
        annotator = ColumnarFactsAnnotator()

        with pytest.raises(TypeError, match="only annotates TranscriptVariationAllele"):
            annotator.annotate_feature_overlaps([_NotATranscriptOverlap()])  # type: ignore[list-item]


class TestOptionalCaches:
    """Tests for the optional strategy caches."""

    def test_zone_map_is_used_without_changing_terms(
        self, boundary_overlaps: dict[str, TranscriptVariationAllele]
    ) -> None:
        """Test the enabled zone map serves facts without changing the result."""
        overlap = boundary_overlaps["ENST0001::ENST0001_cds3_del1"]
        expected = ColumnarFactsAnnotator().annotate_feature_overlaps([overlap])
        zone_map = ZoneMapCache()
        annotator = ColumnarFactsAnnotator(zone_map=zone_map)

        actual = annotator.annotate_feature_overlaps([overlap])

        assert actual == expected
        assert annotator.annotator_context.zone_map_cache is zone_map
        assert zone_map.maps_built == 1
        assert zone_map.facts_served > 0

    def test_codon_cache_is_reused_without_changing_terms(
        self, boundary_overlaps: dict[str, TranscriptVariationAllele]
    ) -> None:
        """Test the enabled codon cache serves a repeated coding edit."""
        overlap = boundary_overlaps["ENST0001::ENST0001_cds3_del1"]
        codon_cache = CodonEditLRUCache()
        annotator = ColumnarFactsAnnotator(codon_cache=codon_cache)

        consequences = annotator.annotate_feature_overlaps([overlap, overlap])

        assert consequences[1] == consequences[0]
        assert annotator.annotator_context.codon_cache is codon_cache
        assert codon_cache.miss_count == 1
        assert codon_cache.hit_count == 1


def _find_and_node_with_atoms(plan_node: object, names: frozenset[str]) -> _AndPlanNode | None:
    """Depth-first search for an ``_AndPlanNode`` whose direct atom children are exactly ``names``."""
    if isinstance(plan_node, _AndPlanNode):
        direct_atom_names = frozenset(
            child.fact_name for child in plan_node.children if isinstance(child, _AtomPlanNode)
        )
        if direct_atom_names == names:
            return plan_node
    children: tuple[object, ...]
    if isinstance(plan_node, (_AndPlanNode, _OrPlanNode)):
        children = plan_node.children
    elif isinstance(plan_node, _NotPlanNode):
        children = (plan_node.child,)
    else:
        children = ()
    for child in children:
        found = _find_and_node_with_atoms(child, names)
        if found is not None:
            return found
    return None


class DummySlice(Slice):
    """Slice with a sequence supplied through the unsafe test hook."""

    def __init__(self, sequence: str) -> None:
        """Initialize the slice with a cached sequence."""
        super().__init__(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name="chrPin",
            start=1,
            end=len(sequence),
            strand=Strand.FORWARD,
        )
        with self.unsafe() as unsafe_slice:
            setattr(unsafe_slice, "_seq", sequence)


def _build_partial_codon_stop_gain_overlap() -> TranscriptVariationAllele:
    """A coding SNV inside an N-filled exon: its narrow codon window translates with an ambiguous
    ("X") residue, so ``peptide_has_X`` is True for it -- unlike an ordinary boundary-corpus SNV.
    """
    exon_one_seq, exon_three_seq = "ATGGGA", "CCCTAA"
    slice_ = DummySlice(exon_one_seq + "N" * 9 + exon_three_seq)
    transcript = Transcript(start=1, end=21, strand=Strand.FORWARD, stable_id="ENST_XSTOPGAIN", slice=slice_)
    transcript.spliced_exons = [
        Exon(start=1, end=6, strand=Strand.FORWARD, stable_id="ENST_XSTOPGAIN_e0", slice=slice_),
        Exon(start=10, end=12, strand=Strand.FORWARD, stable_id="ENST_XSTOPGAIN_e1"),  # no slice -> N-filled
        Exon(start=16, end=21, strand=Strand.FORWARD, stable_id="ENST_XSTOPGAIN_e2", slice=slice_),
    ]
    transcript.seed_attrs(
        first_coding_exon_index=0, last_coding_exon_index=2, exon_coding_start=1, exon_coding_end=6
    )
    variation_allele = VariationAllele(
        start=10,
        end=10,
        strand=Strand.FORWARD,
        slice=slice_,
        reference_allele="N",
        alternate_allele="A",
    )
    return TranscriptVariationAllele(feature=transcript, overlapping_feature=variation_allele)


class TestStopGainedPlanOrder:
    """Stage 2's ordering claim: ``peptide_has_X`` sorts before ``cds_stop_gain`` inside
    ``stop_gained``'s compiled plan, so the expensive full-CDS translation
    (``_full_cds_peptide_alleles``) is reached only for overlaps ``peptide_has_X`` already
    settled True.
    """

    def test_peptide_has_x_sorts_before_cds_stop_gain(self) -> None:
        """The compiled plan orders peptide_has_X ahead of cds_stop_gain within their shared And."""
        rule = build_fact_ruleset().by_term["stop_gained"]
        plan = RulePlan.build(rule.when)

        and_node = _find_and_node_with_atoms(
            plan.root_plan_node, frozenset({"peptide_has_X", "cds_stop_gain"})
        )

        assert and_node is not None
        names_in_order = [child.fact_name for child in and_node.children if isinstance(child, _AtomPlanNode)]
        assert names_in_order.index("peptide_has_X") < names_in_order.index("cds_stop_gain")

    def test_full_cds_peptide_alleles_only_computed_where_peptide_has_x_is_true(
        self, boundary_overlaps: dict[str, TranscriptVariationAllele]
    ) -> None:
        """Evaluating stop_gained's compiled plan alone over a mixed batch never computes
        ``_full_cds_peptide_alleles`` for an overlap where ``peptide_has_X`` is False."""
        rule = build_fact_ruleset().by_term["stop_gained"]
        plan = RulePlan.build(rule.when)

        x_overlap = _build_partial_codon_stop_gain_overlap()
        overlaps = [x_overlap, *boundary_overlaps.values()]
        fact_contexts = [FactContext(overlap=overlap) for overlap in overlaps]
        bit_columns = BitColumnBuffer(fact_contexts)
        live_mask = (1 << len(fact_contexts)) - 1

        plan.evaluate(live_mask, bit_columns)

        assert fact_contexts[0].memo.get("peptide_has_X") is True
        computed_full_cds_count = 0
        for fact_context in fact_contexts:
            if "_full_cds_peptide_alleles" in fact_context.memo:
                computed_full_cds_count += 1
                assert fact_context.memo.get("peptide_has_X") is True, fact_context.overlap
        assert computed_full_cds_count >= 1
