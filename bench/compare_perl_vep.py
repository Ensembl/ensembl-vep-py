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
"""Compare ensembl-vep-py with Perl VEP: speed over input sizes, and consequence agreement.

Both tools annotate the first N records of the same VCF against the same GFF3 and FASTA, for each
N in ``--sizes``. Wall time is plotted against N. The outputs of the largest size are then compared
per (variant, allele, transcript). Two kinds of difference are kept apart, because they have
different causes: a transcript only one tool annotated (the two GFF readers loaded different
transcripts) and a transcript both annotated with different terms (the consequence logic differs).
The report gives agreement for each, the most common term-set patterns, a per-term table over the
shared transcripts, and the first ``--show`` disagreements with the ``vep explain`` command that
derives ours. Every disagreement is also written to a TSV.

Run with matplotlib available, e.g.::

    uv run --with matplotlib python compare_perl_vep.py \
        --input ~/.vep/input/clinvar_20260905_chr22.vcf.gz \
        --gff ~/.vep/gff/genes.gff3.bgz \
        --fasta ~/.vep/fasta/Homo_sapiens.GRCh38.dna.toplevel.fa \
        --ensembl_root ~/ensembl-repos --sizes 1000,10000,all --fork 1,4
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import statistics
import subprocess
import sys
import time
from collections import Counter
from collections.abc import Iterator
from pathlib import Path

from attrs import define, field

PERL_LIBRARIES = (
    "ensembl/modules",
    "ensembl-variation/modules",
    "ensembl-io/modules",
    "ensembl-funcgen/modules",
    "bioperl-live/lib",
    "ensembl-vep",
    "ensembl-vep/modules",
)
"""Checkouts under ``--ensembl_root`` that Perl VEP needs on PERL5LIB."""

TOOLS = ("perl", "python")

AnnotationKey = tuple[str, int, str, str, str, str]
"""(chrom, pos, id, alt column, CSQ allele, transcript id) -- one compared annotation."""


@define(kw_only=True, auto_attribs=True)
class Timing:
    """Wall times of one tool on one input size.

    Attributes:
        tool (str): ``perl`` or ``python``.
        variants (int): Records in the input.
        fork (int): Worker processes given to the tool.
        seconds (list[float]): Wall time of each repeat.
    """

    tool: str
    variants: int
    fork: int
    seconds: list[float] = field(factory=list)

    @property
    def median(self) -> float:
        """Return the median wall time over the repeats."""
        return statistics.median(self.seconds)

    @property
    def per_second(self) -> float:
        """Return variants annotated per second from the median."""
        return self.variants / self.median if self.median else 0.0


@define(kw_only=True, auto_attribs=True)
class Record:
    """What one VCF record's CSQ says, by allele and transcript.

    Attributes:
        chrom (str): CHROM column.
        pos (int): POS column.
        identifier (str): ID column.
        ref (str): REF column.
        alts (list[str]): ALT column, split.
        terms (dict[tuple[str, str], frozenset[str]]): (CSQ allele, transcript) -> consequence terms.
    """

    chrom: str
    pos: int
    identifier: str
    ref: str
    alts: list[str]
    terms: dict[tuple[str, str], frozenset[str]] = field(factory=dict)


@define(kw_only=True, auto_attribs=True)
class Agreement:
    """What comparing the two outputs found.

    Attributes:
        counts (Counter[str]): Overall tallies, keyed by what they count.
        per_term (Counter[tuple[str, str]]): (term, side) -> count over transcripts both tools
            annotated; side is ``both``, ``perl_only`` or ``python_only``.
        patterns (Counter[tuple[str, str]]): (Perl terms, our terms) -> count, for shared
            transcripts whose term sets differ.
        one_sided (Counter[tuple[str, str]]): (kind, transcript) -> count, for transcripts only one
            tool annotated.
        disagreements (list[Disagreement]): Every disagreement, in input order.
    """

    counts: Counter[str] = field(factory=Counter)
    per_term: Counter[tuple[str, str]] = field(factory=Counter)
    patterns: Counter[tuple[str, str]] = field(factory=Counter)
    one_sided: Counter[tuple[str, str]] = field(factory=Counter)
    disagreements: list[Disagreement] = field(factory=list)


@define(kw_only=True, auto_attribs=True)
class Disagreement:
    """One annotation the two tools did not agree on.

    Attributes:
        key (AnnotationKey): Which variant, allele and transcript.
        perl (frozenset[str] | None): Perl's terms, None when Perl did not annotate the pair.
        python (frozenset[str] | None): Our terms, None when we did not annotate the pair.
        explain (str): The ``vep explain`` command that derives our side.
    """

    key: AnnotationKey
    perl: frozenset[str] | None
    python: frozenset[str] | None
    explain: str

    @property
    def kind(self) -> str:
        """Return ``perl_only``, ``python_only`` or ``terms_differ``."""
        if self.python is None:
            return "perl_only"
        if self.perl is None:
            return "python_only"
        return "terms_differ"


def _open_text(path: Path) -> Iterator[str]:
    """Yield the lines of a plain or gzip/bgzip-compressed text file."""
    opener = gzip.open if path.suffix in {".gz", ".bgz"} else open
    with opener(path, "rt", encoding="utf-8") as handle:
        yield from handle


def write_subset(source: Path, destination: Path, size: int | None) -> int:
    """Write the header and the first ``size`` records of a VCF, or all of them for None.

    Returns:
        int: The number of records written.
    """
    written = 0
    with destination.open("w", encoding="utf-8") as out:
        for line in _open_text(source):
            if line.startswith("#"):
                out.write(line)
                continue
            if size is not None and written >= size:
                break
            out.write(line)
            written += 1
    return written


def perl_command(arguments: argparse.Namespace, input_file: Path, output_file: Path, fork: int) -> list[str]:
    """Return the Perl VEP command line for one run."""
    command = [
        "perl",
        str(arguments.ensembl_root / "ensembl-vep" / "vep"),
        "--input_file", str(input_file),
        "--gff", str(arguments.gff),
        "--fasta", str(arguments.fasta),
        "--vcf",
        "--output_file", str(output_file),
        "--force_overwrite",
        "--no_stats",
    ]  # fmt: skip
    if fork > 1:
        command.extend(["--fork", str(fork)])
    return command + arguments.perl_args.split()


def python_command(
    arguments: argparse.Namespace, input_file: Path, output_file: Path, fork: int
) -> list[str]:
    """Return the ensembl-vep-py command line for one run."""
    command = [
        sys.executable, "-m", "ensembl.vep.cli", "run", "default",
        "--input_file", str(input_file),
        "--gff", str(arguments.gff),
        "--fasta", str(arguments.fasta),
        "--output_file", str(output_file),
        "--output_format", "vcf",
        "--force_overwrite",
        "--no_stats",
        "--no_progress",
        "--fork", str(fork),
    ]  # fmt: skip
    print(sys.executable)
    return command + arguments.python_args.split()


def perl_environment(ensembl_root: Path) -> dict[str, str]:
    """Return the environment with the Ensembl checkouts prepended to PERL5LIB."""
    environment = dict(os.environ)
    libraries = [str(ensembl_root / library) for library in PERL_LIBRARIES]
    if environment.get("PERL5LIB"):
        libraries.append(environment["PERL5LIB"])
    environment["PERL5LIB"] = os.pathsep.join(libraries)
    return environment


def run_timed(command: list[str], environment: dict[str, str], log_file: Path) -> float:
    """Run a command, append its output to a log, and return its wall time in seconds.

    Raises:
        SystemExit: If the command fails; the log names why.
    """
    started = time.perf_counter()
    with log_file.open("a", encoding="utf-8") as log:
        log.write(f"\n$ {' '.join(command)}\n")
        log.flush()
        completed = subprocess.run(
            command, env=environment, stdout=log, stderr=subprocess.STDOUT, check=False
        )
    elapsed = time.perf_counter() - started
    if completed.returncode != 0:
        raise SystemExit(f"command failed ({completed.returncode}), see {log_file}:\n  {' '.join(command)}")
    return elapsed


def csq_layout(line: str) -> list[str]:
    """Return the CSQ field names an ``##INFO=<ID=CSQ`` header line declares."""
    description = line.split("Format:", 1)[1]
    return description.strip().rstrip('">').strip().split("|")


def feature_id(value: str) -> str:
    """Return a transcript id without a GFF3 type prefix: ours keeps ``transcript:``, Perl drops it."""
    return value.split(":", 1)[-1]


def read_csq(path: Path) -> dict[tuple[str, int, str, str], Record]:
    """Read every record's CSQ into term sets, keyed by (chrom, pos, id, alt column).

    Only Transcript rows are kept, so intergenic and regulatory rows do not count as disagreements
    about transcripts.
    """
    layout: list[str] = []
    records: dict[tuple[str, int, str, str], Record] = {}
    for line in _open_text(path):
        if line.startswith("##INFO=<ID=CSQ"):
            layout = csq_layout(line)
            continue
        if line.startswith("#"):
            continue
        columns = line.rstrip("\n").split("\t")
        chrom, pos, identifier, ref, alt, info = (
            columns[0],
            int(columns[1]),
            columns[2],
            columns[3],
            columns[4],
            columns[7],
        )
        record = Record(chrom=chrom, pos=pos, identifier=identifier, ref=ref, alts=alt.split(","))
        records[(chrom, pos, identifier, alt)] = record

        csq = next((entry[4:] for entry in info.split(";") if entry.startswith("CSQ=")), None)
        if csq is None or not layout:
            continue
        for annotation in csq.split(","):
            values = dict(zip(layout, annotation.split("|")))
            if values.get("Feature_type") != "Transcript" or not values.get("Feature"):
                continue
            terms = frozenset(values.get("Consequence", "").split("&")) - {""}
            record.terms[(values.get("Allele", ""), feature_id(values["Feature"]))] = terms
    return records


def csq_alleles(ref: str, alts: list[str]) -> list[str]:
    """Return each ALT as VEP spells it in CSQ: shared first base stripped, ``-`` when empty."""
    alleles = [ref, *alts]
    bases = set("ACGTN")
    if all(allele and set(allele.upper()) <= bases for allele in alleles):
        if len({allele[0] for allele in alleles}) == 1:
            alleles = [allele[1:] or "-" for allele in alleles]
    return alleles[1:]


def explain_command(arguments: argparse.Namespace, record: Record, allele: str, transcript: str) -> str:
    """Return the ``vep explain`` command that derives our terms for one annotation."""
    alleles = csq_alleles(record.ref, record.alts)
    allele_num = alleles.index(allele) + 1 if allele in alleles else 1
    if record.identifier not in {"", "."}:
        variant = f"{record.chrom}:{record.pos}:{record.identifier.split(';')[0]}"
        allele_option = f" --allele_num {allele_num}" if len(record.alts) > 1 else ""
    else:
        trimmed = alleles != record.alts
        start = record.pos + 1 if trimmed else record.pos
        reference = (record.ref[1:] or "-") if trimmed else record.ref
        variant = f"{record.chrom}:{start}:{reference}>{allele}"
        allele_option = ""
    return (
        f"vep explain --input_file {arguments.input} --gff {arguments.gff} --fasta {arguments.fasta} "
        f"--variant '{variant}'{allele_option} --transcript {transcript}"
    )


def compare(arguments: argparse.Namespace, perl_output: Path, python_output: Path) -> Agreement:
    """Compare the two outputs annotation by annotation."""
    perl_records = read_csq(perl_output)
    python_records = read_csq(python_output)
    agreement = Agreement()
    counts = agreement.counts

    counts["records_in_one_output_only"] = len(perl_records.keys() ^ python_records.keys())
    counts["records_compared"] = len(perl_records.keys() & python_records.keys())

    for record_key, perl_record in perl_records.items():
        python_record = python_records.get(record_key)
        if python_record is None:
            continue
        for pair in sorted(perl_record.terms.keys() | python_record.terms.keys()):
            perl_terms = perl_record.terms.get(pair)
            python_terms = python_record.terms.get(pair)
            counts["annotations"] += 1
            if perl_terms is not None and python_terms is not None:
                counts["shared"] += 1
                for term in perl_terms & python_terms:
                    agreement.per_term[(term, "both")] += 1
                for term in perl_terms - python_terms:
                    agreement.per_term[(term, "perl_only")] += 1
                for term in python_terms - perl_terms:
                    agreement.per_term[(term, "python_only")] += 1
            if perl_terms == python_terms:
                counts["agree"] += 1
                continue

            allele, transcript = pair
            ref = python_record.ref
            disagreement = Disagreement(
                key=(*record_key, ref, transcript),
                perl=perl_terms,
                python=python_terms,
                explain=explain_command(arguments, python_record, allele, transcript),
            )
            counts[disagreement.kind] += 1
            if disagreement.kind == "terms_differ":
                agreement.patterns[(_terms_text(perl_terms), _terms_text(python_terms))] += 1
            else:
                agreement.one_sided[(disagreement.kind, transcript)] += 1
            agreement.disagreements.append(disagreement)
    return agreement


def _terms_text(terms: frozenset[str] | None) -> str:
    """Return a term set as VEP writes it, or ``(not annotated)``."""
    return "&".join(sorted(terms)) if terms is not None else "(not annotated)"


def _percent(part: int, whole: int) -> str:
    """Return ``part`` as a percentage of ``whole``, or ``-`` when there is nothing to divide."""
    return f"{100 * part / whole:.3f}%" if whole else "-"


def print_agreement(agreement: Agreement, *, show: int, term: str | None, kind: str | None, top: int) -> None:
    """Print agreement by kind, the common patterns, the per-term table, and the first disagreements.

    Args:
        agreement (Agreement): What the comparison found.
        show (int): Disagreements to print.
        term (str | None): Print only disagreements whose term sets differ in this term.
        kind (str | None): Print only disagreements of this kind.
        top (int): Rows to print of each ranked table.
    """
    counts = agreement.counts
    annotations = counts["annotations"]
    shared = counts["shared"]
    agree_shared = shared - counts["terms_differ"]
    print("\n== consequence agreement (per variant, allele and transcript) ==")
    print(f"records compared                    {counts['records_compared']:>10,}")
    print(f"records in one output only          {counts['records_in_one_output_only']:>10,}")
    print(f"annotations, either tool            {annotations:>10,}")
    identical = counts["agree"]
    print(f"  identical                         {identical:>10,}  {_percent(identical, annotations)}")
    print(f"  transcript annotated by Perl only {counts['perl_only']:>10,}")
    print(f"  transcript annotated by us only   {counts['python_only']:>10,}")
    print(f"annotations, both tools             {shared:>10,}")
    print(f"  identical terms                   {agree_shared:>10,}  {_percent(agree_shared, shared)}")
    print(f"  different terms                   {counts['terms_differ']:>10,}")

    print(f"\n== top {top} term-set patterns on shared transcripts ==")
    print(f"{'count':>7}  perl  ->  python")
    for (perl_text, python_text), count in agreement.patterns.most_common(top):
        print(f"{count:>7,}  {perl_text}  ->  {python_text}")

    print(f"\n== top {top} transcripts annotated by one tool only ==")
    for (side, transcript), count in agreement.one_sided.most_common(top):
        print(f"{count:>7,}  {transcript:<20} {side}")

    terms = sorted({name for name, _side in agreement.per_term})
    print("\n== per term, over transcripts both tools annotated ==")
    print(f"{'term':<40}{'both':>10}{'perl only':>11}{'py only':>10}{'agree':>9}")
    for name in terms:
        both = agreement.per_term[(name, "both")]
        perl_only = agreement.per_term[(name, "perl_only")]
        python_only = agreement.per_term[(name, "python_only")]
        total = both + perl_only + python_only
        flag = "" if perl_only == python_only == 0 else "  <--"
        print(f"{name:<40}{both:>10,}{perl_only:>11,}{python_only:>10,}{100 * both / total:>8.2f}%{flag}")

    selected = []
    for item in agreement.disagreements:
        if kind is not None and item.kind != kind:
            continue
        if term is not None and term not in (item.perl or frozenset()) ^ (item.python or frozenset()):
            continue
        selected.append(item)
    filters = ", ".join(f"{name}={value}" for name, value in (("kind", kind), ("term", term)) if value)
    heading = f" ({filters})" if filters else ""
    print(f"\n== first {min(show, len(selected))} of {len(selected):,} disagreements{heading} ==")
    for item in selected[:show]:
        chrom, pos, identifier, alt, allele, transcript = item.key
        print(f"\n{chrom}:{pos} {identifier} ALT={alt} allele={allele} {transcript} [{item.kind}]")
        print(f"  perl:   {_terms_text(item.perl)}")
        print(f"  python: {_terms_text(item.python)}")
        print(f"  {item.explain}")


def write_disagreements(path: Path, disagreements: list[Disagreement]) -> None:
    """Write every disagreement as a TSV, for filtering beyond what is printed."""
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ["chrom", "pos", "id", "alt", "ref", "transcript", "kind", "perl", "python", "explain"]
        )
        for item in disagreements:
            writer.writerow(
                [*item.key, item.kind, _terms_text(item.perl), _terms_text(item.python), item.explain]
            )


def print_timings(timings: list[Timing]) -> None:
    """Print the median wall time and throughput of every run, with Perl's time over ours."""
    print("\n== speed (median wall time) ==")
    columns = ("perl s", "python s", "speed-up")
    print(
        f"{'variants':>10}{'fork':>6}"
        + "".join(f"{name:>10}" for name in columns)
        + f"{'perl v/s':>11}{'py v/s':>10}"
    )
    by_run = {(timing.tool, timing.variants, timing.fork): timing for timing in timings}
    for variants, fork in sorted({(timing.variants, timing.fork) for timing in timings}):
        perl = by_run.get(("perl", variants, fork))
        python = by_run.get(("python", variants, fork))
        perl_seconds = f"{perl.median:.2f}" if perl else "-"
        python_seconds = f"{python.median:.2f}" if python else "-"
        speedup = f"{perl.median / python.median:.2f}x" if perl and python and python.median else "-"
        perl_rate = f"{perl.per_second:,.0f}" if perl else "-"
        python_rate = f"{python.per_second:,.0f}" if python else "-"
        timing_columns = f"{perl_seconds:>10}{python_seconds:>10}{speedup:>10}"
        print(f"{variants:>10,}{fork:>6}{timing_columns}{perl_rate:>11}{python_rate:>10}")


def plot_timings(timings: list[Timing], path: Path) -> None:
    """Plot wall time and throughput against input size, one line per tool and fork count."""
    # Optional: only the plot needs it, so the timings and the comparison run without it.
    import matplotlib  # pylint: disable=import-outside-toplevel,import-error

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # pylint: disable=import-outside-toplevel,import-error

    figure, (time_axis, rate_axis) = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    colours = {"perl": "#2a6fdb", "python": "#e8743b"}
    styles = ["-", "--", ":", "-."]
    forks = sorted({timing.fork for timing in timings})
    for tool in TOOLS:
        for fork_index, fork in enumerate(forks):
            runs = sorted(
                (timing for timing in timings if timing.tool == tool and timing.fork == fork),
                key=lambda timing: timing.variants,
            )
            if not runs:
                continue
            sizes = [timing.variants for timing in runs]
            label = f"{'Perl VEP' if tool == 'perl' else 'ensembl-vep-py'} (fork {fork})"
            style = {
                "color": colours[tool],
                "linestyle": styles[fork_index % len(styles)],
                "marker": "o",
                "linewidth": 2,
            }
            time_axis.plot(sizes, [timing.median for timing in runs], label=label, **style)
            low = [min(timing.seconds) for timing in runs]
            high = [max(timing.seconds) for timing in runs]
            time_axis.fill_between(sizes, low, high, color=colours[tool], alpha=0.12)
            rate_axis.plot(sizes, [timing.per_second for timing in runs], label=label, **style)

    # pylint: disable-next=import-outside-toplevel,import-error
    from matplotlib.ticker import FuncFormatter, NullFormatter

    plain = FuncFormatter(lambda value, _position: f"{value:,.10g}")
    for axis in (time_axis, rate_axis):
        axis.set_xscale("log")
        axis.xaxis.set_major_formatter(plain)
        axis.xaxis.set_minor_formatter(NullFormatter())
        axis.set_xlabel("variants in input")
        axis.grid(True, which="both", alpha=0.3)
    time_axis.set_yscale("log")
    time_axis.yaxis.set_major_formatter(plain)
    time_axis.yaxis.set_minor_formatter(plain)
    time_axis.set_ylabel("wall time (s), median; band = min to max")
    time_axis.set_title("Run time")
    rate_axis.set_ylabel("variants per second")
    rate_axis.set_title("Throughput")
    rate_axis.legend(frameon=False)
    figure.savefig(path, dpi=150)
    plt.close(figure)


def parse_sizes(text: str) -> list[int | None]:
    """Parse ``1000,10000,all`` into record counts, None meaning the whole input."""
    sizes: list[int | None] = []
    for piece in text.split(","):
        piece = piece.strip().lower()
        sizes.append(None if piece == "all" else int(piece))
    return sizes


def parse_arguments(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse the command line."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("--input", type=Path, required=True, help="VCF to annotate (plain, gz or bgz)")
    parser.add_argument("--gff", type=Path, required=True, help="GFF3 transcripts, given to both tools")
    parser.add_argument("--fasta", type=Path, required=True, help="reference FASTA, given to both tools")
    parser.add_argument(
        "--ensembl_root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
        help="directory holding ensembl-vep, ensembl, ensembl-variation, ensembl-io and bioperl-live",
    )
    parser.add_argument(
        "--sizes", default="1000,10000,all", help="record counts to time; 'all' is the whole input"
    )
    parser.add_argument("--fork", default="1", help="worker counts to time each size with, e.g. 1,4")
    parser.add_argument("--repeats", type=int, default=1, help="runs per tool, size and fork count")
    parser.add_argument(
        "--outdir", type=Path, default=Path("perl_comparison"), help="where outputs and the plot go"
    )
    parser.add_argument("--perl_args", default="", help="extra arguments for Perl VEP, as one string")
    parser.add_argument("--python_args", default="", help="extra arguments for ensembl-vep-py, as one string")
    parser.add_argument("--tools", default="perl,python", help="which tools to run; compare needs both")
    parser.add_argument(
        "--compare_only", action="store_true", help="skip timing; compare existing outputs in --outdir"
    )
    parser.add_argument("--show", type=int, default=20, help="disagreements to print")
    parser.add_argument("--term", default=None, help="print only disagreements involving this SO term")
    parser.add_argument(
        "--kind",
        choices=("terms_differ", "perl_only", "python_only"),
        default=None,
        help="print only disagreements of this kind",
    )
    parser.add_argument("--top", type=int, default=15, help="rows of each ranked table")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Time both tools, plot the timings, and compare the largest outputs."""
    arguments = parse_arguments(argv)
    arguments.input = arguments.input.expanduser().resolve()
    arguments.gff = arguments.gff.expanduser().resolve()
    arguments.fasta = arguments.fasta.expanduser().resolve()
    arguments.ensembl_root = arguments.ensembl_root.expanduser().resolve()
    outdir: Path = arguments.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    log_file = outdir / "runs.log"

    sizes = parse_sizes(arguments.sizes)
    forks = [int(piece) for piece in arguments.fork.split(",")]
    tools = [tool for tool in arguments.tools.split(",") if tool in TOOLS]
    builders = {"perl": perl_command, "python": python_command}
    environments = {"perl": perl_environment(arguments.ensembl_root), "python": dict(os.environ)}
    largest = sizes[-1]
    largest_label = "all" if largest is None else str(largest)

    timings: list[Timing] = []
    if not arguments.compare_only:
        for size in sizes:
            label = "all" if size is None else str(size)
            subset = outdir / f"input_{label}.vcf"
            variants = write_subset(arguments.input, subset, size)
            for fork in forks:
                for tool in tools:
                    timing = Timing(tool=tool, variants=variants, fork=fork)
                    output = outdir / f"{tool}_{label}_fork{fork}.vcf"
                    for repeat in range(arguments.repeats):
                        command = builders[tool](arguments, subset, output, fork)
                        timing.seconds.append(run_timed(command, environments[tool], log_file))
                        print(
                            f"{tool:<7} {variants:>9,} variants  fork {fork}  run {repeat + 1}: "
                            f"{timing.seconds[-1]:.2f}s",
                            flush=True,
                        )
                    timings.append(timing)

        print_timings(timings)
        (outdir / "timings.json").write_text(
            json.dumps(
                [
                    {"tool": t.tool, "variants": t.variants, "fork": t.fork, "seconds": t.seconds}
                    for t in timings
                ],
                indent=2,
            ),
            encoding="utf-8",
        )
        try:
            plot_timings(timings, outdir / "speed.png")
            print(f"\nplot: {outdir / 'speed.png'}")
        except ImportError:
            print("\nmatplotlib is not installed; run with 'uv run --with matplotlib' for the plot")

    perl_output = outdir / f"perl_{largest_label}_fork{forks[0]}.vcf"
    python_output = outdir / f"python_{largest_label}_fork{forks[0]}.vcf"
    if not (perl_output.exists() and python_output.exists()):
        print(f"\nno pair of outputs to compare ({perl_output.name}, {python_output.name})")
        return 0

    agreement = compare(arguments, perl_output, python_output)
    print(f"\ncomparing {perl_output.name} with {python_output.name}")
    print_agreement(
        agreement, show=arguments.show, term=arguments.term, kind=arguments.kind, top=arguments.top
    )
    write_disagreements(outdir / "disagreements.tsv", agreement.disagreements)
    summary = {
        "counts": dict(agreement.counts),
        "per_term": {f"{term}\t{side}": count for (term, side), count in sorted(agreement.per_term.items())},
        "patterns": [
            {"perl": perl_text, "python": python_text, "count": count}
            for (perl_text, python_text), count in agreement.patterns.most_common()
        ],
    }
    (outdir / "agreement.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"\nall disagreements: {outdir / 'disagreements.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
