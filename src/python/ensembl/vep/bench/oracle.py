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
"""The Perl oracle: the same cases through our annotators and through ``VariationEffect.pm``.

The oracle is the external reference -- the whole point is that it is *not* us -- so a
disagreement is reported, never gated on. The Perl script tests only the coding, intronic, upstream
and downstream categories, so the comparison asks whether our terms *cover* the Perl ones, not
whether they are equal.

Coverage is counted twice, because the direct count undercounts agreement. The Perl script reports
raw predicates, while our annotators emit terms after subsumption: an intronic variant in the
polypyrimidine tract emits ``splice_polypyrimidine_tract_variant``, which suppresses the
``intron_variant`` the Perl side reports. Counting the terms our emitted ones suppress, through the
catalogue's ``subsumes`` chains, separates that from a real disagreement.

The case list here and ``build_cases`` in ``bench/variation_effect.pl`` must stay aligned by name
and order; ``CASE_NAMES`` is what a test pins against the Perl file. The transcript is built with
the corpus's own ``build_transcript``, so its translation is derived exactly as the suite's is --
the hand-rolled builder this replaces had fallen out of step with the model and no longer ran.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import time
from collections.abc import Mapping, Sequence
from pathlib import Path

from attrs import frozen

from ensembl.vep.annotators.base import Annotator
from ensembl.vep.common import CoordSystemLevel, Strand
from ensembl.vep.corpus.synthetic.components.transcript import build_transcript
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.overlap import OverlapConsequence, TranscriptVariationAllele
from ensembl.vep.models.sequence import CoordSystem, Slice
from ensembl.vep.models.variation import VariationAllele

__all__ = [
    "CASE_NAMES",
    "OracleCase",
    "covered_cases",
    "covered_through_subsumption",
    "oracle_cases",
    "oracle_transcript",
    "perl_seconds",
    "perl_terms",
    "perl_unavailable_reason",
    "python_seconds",
    "python_terms",
    "subsumption_closure",
]

CASE_NAMES = (
    "coding_exact",
    "coding_near_boundary",
    "intron_primary",
    "intron_secondary",
    "upstream_near",
    "upstream_far",
    "downstream_near",
    "downstream_far",
    "no_hit",
)
"""The cases, in the order the Perl script builds and reports them."""

_CHROMOSOME_LENGTH = 10000
_PERL_TIMEOUT_SECONDS = 600


@frozen(kw_only=True, auto_attribs=True)
class OracleCase:
    """A named overlap, paired with the Perl case of the same name.

    Attributes:
        name (str): The case name, shared with the Perl script.
        overlap (TranscriptVariationAllele): The overlap to annotate.
    """

    name: str
    overlap: TranscriptVariationAllele


def oracle_transcript() -> Transcript:
    """Return a small forward-strand coding transcript: two exons around one intron, all of both coding.

    Exon 1 is 1000-1029 and exon 2 1050-1088, so the CDS is 69 bases -- 23 codons, a multiple of
    three, so no partial-codon edge leaks into a case that is not about one. The background is a
    plain ``A`` repeat: the comparison needs only the category a variant lands in to be
    unambiguous, not a realistic open reading frame.
    """
    chromosome = Slice(
        coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
        seq_region_name="1",
        start=1,
        end=_CHROMOSOME_LENGTH,
        strand=Strand.FORWARD,
    )
    with chromosome.unsafe() as unsafe_slice:
        setattr(unsafe_slice, "_seq", "A" * _CHROMOSOME_LENGTH)
    return build_transcript(
        stable_id="ENST00000000001",
        synthetic_slice=chromosome,
        strand=Strand.FORWARD,
        exons_genomic_ascending=[(1000, 1029), (1050, 1088)],
        cds_genomic_bounds=(1000, 1088),
        gene=None,
        biotype="protein_coding",
    )


def _snv(transcript: Transcript, position: int) -> TranscriptVariationAllele:
    """Return an A>T SNV at a genomic position, against the transcript."""
    allele = VariationAllele(
        start=position,
        end=position,
        strand=Strand.FORWARD,
        slice=transcript.slice,
        reference_allele="A",
        alternate_allele="T",
    )
    return TranscriptVariationAllele(feature=transcript, overlapping_feature=allele)


def oracle_cases() -> tuple[OracleCase, ...]:
    """Return the cases, named and ordered as the Perl script's."""
    transcript = oracle_transcript()
    positions = (
        1010,  # mid exon 1, well inside the CDS
        1028,  # two bases from the donor site
        1039,  # deep intronic, the intron's midpoint
        1032,  # near the donor site, offset +3
        transcript.start - 1,
        transcript.start - 4000,
        transcript.end + 1,
        transcript.end + 4000,
        transcript.end + 6000,  # beyond the 5 kb flank window
    )
    return tuple(
        OracleCase(name=name, overlap=_snv(transcript, position))
        for name, position in zip(CASE_NAMES, positions, strict=True)
    )


def python_terms(
    annotator: Annotator[OverlapConsequence], cases: Sequence[OracleCase]
) -> tuple[tuple[str, ...], ...]:
    """Return one annotator's terms for every case, in case order."""
    return tuple(
        tuple(consequence.term for consequence in annotator.annotate_feature_overlap(case.overlap))
        for case in cases
    )


def perl_unavailable_reason(perl_modules: Path | None, script: Path) -> str | None:
    """Return why the Perl oracle cannot run here, or None if it can."""
    if shutil.which("perl") is None:
        return "perl is not on PATH"
    if not script.is_file():
        return f"oracle script not found: {script}"
    if perl_modules is None:
        return "no ensembl-variation modules directory given (--perl_modules)"
    if not perl_modules.is_dir():
        return f"ensembl-variation modules not found at {perl_modules}"
    return None


def _perl_command(perl_modules: Path, script: Path, *arguments: str) -> list[str]:
    """Return the command line that runs the oracle script."""
    return ["perl", "-I", str(perl_modules), str(script), *arguments]


def perl_terms(perl_modules: Path, script: Path, case_count: int) -> tuple[tuple[str, ...], ...]:
    """Run the oracle for its terms, one tuple per case in its own order.

    Raises:
        ValueError: If the oracle reports a different number of cases, meaning the two case lists
            have drifted apart.
        subprocess.CalledProcessError: If the script fails.
    """
    result = subprocess.run(
        _perl_command(perl_modules, script, "--mode", "json"),
        check=True,
        capture_output=True,
        text=True,
        timeout=_PERL_TIMEOUT_SECONDS,
    )
    parsed = json.loads(result.stdout)
    if len(parsed) != case_count:
        raise ValueError(
            f"the Perl oracle reported {len(parsed)} cases and this side has {case_count}: "
            "build_cases in the Perl script and CASE_NAMES here have drifted apart"
        )
    return tuple(tuple(terms) for terms in parsed)


def perl_seconds(perl_modules: Path, script: Path, *, number: int, repeats: int) -> list[float]:
    """Return the elapsed seconds of each timed oracle run of ``number`` workload executions.

    Each run is a whole process, so the figure includes Perl's start-up; the per-case number it
    yields is an upper bound, and says so wherever it is reported.
    """
    durations = []
    for _ in range(max(1, repeats)):
        started = time.perf_counter()
        subprocess.run(
            _perl_command(perl_modules, script, "--mode", "bench", "--number", str(number)),
            check=True,
            capture_output=True,
            timeout=_PERL_TIMEOUT_SECONDS,
        )
        durations.append(time.perf_counter() - started)
    return durations


def python_seconds(
    annotator: Annotator[OverlapConsequence], cases: Sequence[OracleCase], *, number: int, repeats: int
) -> list[float]:
    """Return the elapsed seconds of each timed run of ``number`` workload executions."""
    durations = []
    for _ in range(max(1, repeats)):
        started = time.perf_counter()
        for _ in range(number):
            for case in cases:
                annotator.annotate_feature_overlap(case.overlap)
        durations.append(time.perf_counter() - started)
    return durations


def covered_cases(perl: Sequence[Sequence[str]], ours: Sequence[Sequence[str]]) -> int:
    """Return how many cases our terms cover every Perl-reported term of."""
    return sum(1 for perl_case, our_case in zip(perl, ours, strict=True) if set(perl_case) <= set(our_case))


def subsumption_closure(subsumes: Mapping[str, Sequence[str]]) -> dict[str, frozenset[str]]:
    """Return every term each term suppresses, directly or through a chain of ``subsumes``."""
    closure: dict[str, frozenset[str]] = {}
    for term in subsumes:
        reached: set[str] = set()
        pending = list(subsumes.get(term, ()))
        while pending:
            current = pending.pop()
            if current in reached:
                continue
            reached.add(current)
            pending.extend(subsumes.get(current, ()))
        closure[term] = frozenset(reached)
    return closure


def covered_through_subsumption(
    perl: Sequence[Sequence[str]], ours: Sequence[Sequence[str]], subsumes: Mapping[str, Sequence[str]]
) -> int:
    """Return how many cases our terms cover, counting the terms they suppress as covered.

    Args:
        perl (Sequence[Sequence[str]]): The Perl-reported terms, per case.
        ours (Sequence[Sequence[str]]): Our emitted terms, per case.
        subsumes (Mapping[str, Sequence[str]]): Term -> the terms it directly suppresses.

    Returns:
        int: Cases whose every Perl-reported term is emitted by us or suppressed by a term we emit.
    """
    closure = subsumption_closure(subsumes)
    covered = 0
    for perl_case, our_case in zip(perl, ours, strict=True):
        explained = set(our_case)
        for term in our_case:
            explained |= closure.get(term, frozenset())
        if set(perl_case) <= explained:
            covered += 1
    return covered
