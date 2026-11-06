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
"""Compare Python, Rust, and augmented-list interval indexes."""

from __future__ import annotations

import argparse
import random
import timeit
from collections.abc import Sequence
from datetime import datetime
from operator import attrgetter
from typing import Any

import numpy as np

import intervaltree  # pylint: disable=import-error
import intervaltree_rs  # pylint: disable=import-error
from ailist import AIList

from ensembl.vep.common import CoordSystemLevel, Strand
from ensembl.vep.models.sequence import CoordSystem, Slice
from ensembl.vep.models.variation.variation import Variation


def genomic_slice(seq_region_name: str = "1", seq_region_end: int = 250_000_000) -> Slice:
    """Create a Slice instance for the benchmark variants."""
    return Slice(
        coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
        seq_region_name=seq_region_name,
        start=1,
        end=seq_region_end,
        strand=Strand.FORWARD,
    )


def create_variation_objects(
    length: int, seq_region_name: str = "1", var_region_end: int = 500_000
) -> list[Variation]:
    """Create variation objects for exercising an interval index."""
    random.seed(82736401937474632)
    var_ends = random.sample(range(1, var_region_end), k=length)
    var_lengths = random.choices([1, 2, 100], cum_weights=[80, 95, 100], k=length)
    return [
        Variation(
            start=max(1, end - length + 1),
            end=end,
            strand=Strand.FORWARD,
            slice=genomic_slice(seq_region_name),
            allele_string="A/T",
        )
        for end, length in zip(var_ends, var_lengths)
    ]


_TEST_VARIANTS = [
    Variation(start=9, end=10, strand=Strand.FORWARD),
    Variation(start=10, end=10, strand=Strand.FORWARD),
    Variation(start=10, end=11, strand=Strand.FORWARD),
    Variation(start=10, end=12, strand=Strand.FORWARD),
    Variation(start=11, end=11, strand=Strand.FORWARD),
    Variation(start=11, end=12, strand=Strand.FORWARD),
    Variation(start=12, end=14, strand=Strand.FORWARD),
    Variation(start=13, end=13, strand=Strand.FORWARD),
    Variation(start=14, end=15, strand=Strand.FORWARD),
    Variation(start=14, end=16, strand=Strand.FORWARD),
    Variation(start=15, end=15, strand=Strand.FORWARD),
    Variation(start=15, end=16, strand=Strand.FORWARD),
    Variation(start=16, end=16, strand=Strand.FORWARD),
    Variation(start=16, end=17, strand=Strand.FORWARD),
]

_TEST_RANGES = [(11, 15)]


def assert_overlapping_variants(var_list: Sequence[Variation], query_range: tuple[int, int]) -> bool:
    """Return whether the query found the expected variation objects."""
    var_list = sorted(var_list, key=attrgetter("start", "end"))
    if query_range == (11, 15):
        try:
            assert var_list == _TEST_VARIANTS[2:12]
        except AssertionError:
            return False
    return True


def check_overlap_correctness() -> None:
    """Report whether each interval implementation finds the expected overlaps."""
    py_tree = intervaltree.IntervalTree()
    for var in _TEST_VARIANTS:
        py_tree.add(intervaltree.Interval(var.start, var.end + 1, var))
    for interval in _TEST_RANGES:
        var_list = [overlap[2] for overlap in py_tree.overlap(interval[0], interval[1] + 1)]
        if assert_overlapping_variants(var_list, interval):
            print(f"PyTree overlaps: {interval} = CORRECT")
        else:
            print(f"PyTree overlaps: {interval} = INCORRECT!")
            for var in sorted(var_list, key=attrgetter("start", "end")):
                print(f"{var.start=}, {var.end=}")

    rust_tree = intervaltree_rs.IntervalTree()
    for var in _TEST_VARIANTS:
        rust_tree.insert((var.start, var.end + 1, var))
    for interval in _TEST_RANGES:
        var_list = [overlap[2] for overlap in rust_tree.search(interval[0] + 1, interval[1], inclusive=True)]
        if assert_overlapping_variants(var_list, interval):
            print(f"RustTree overlaps: {interval} = CORRECT")
        else:
            print(f"RustTree overlaps: {interval} = INCORRECT!")
            for var in sorted(var_list, key=attrgetter("start", "end")):
                print(f"{var.start=}, {var.end=}")

    ailist = AIList()
    for i, var in enumerate(_TEST_VARIANTS):
        ailist.add(var.start, var.end + 1, i)
    for interval in _TEST_RANGES:
        var_list = [_TEST_VARIANTS[idx] for idx in ailist.intersect_ids(interval[0], interval[1] + 1)]
        if assert_overlapping_variants(var_list, interval):
            print(f"AIList overlaps: {interval} = CORRECT")
        else:
            print(f"AIList overlaps: {interval} = INCORRECT!")
            for var in sorted(var_list, key=attrgetter("start", "end")):
                print(f"{var.start=}, {var.end=}")


class PyIntTree:
    """Benchmark fixture for the Python ``intervaltree`` package.

    Attributes:
        num_vars (int): Number of generated variations.
        variants (list[Variation]): Generated variations to index.
    """

    def __init__(self, num_vars: int = 100) -> None:
        """Create a fixture containing ``num_vars`` variations."""
        self.num_vars = num_vars
        self.variants = create_variation_objects(num_vars)
        self.tree: Any = None
        self.overlaps: Any = None
        self.query_start = 0
        self.query_end = 0

    def intervaltree_create(self) -> None:
        """Build an interval tree by adding one variation at a time."""
        tree = intervaltree.IntervalTree()
        for var in self.variants:
            tree.add(intervaltree.Interval(var.start, var.end + 1, var))
        self.tree = tree

    def intervaltree_create_from_tuples(self) -> None:
        """Build an interval tree from coordinate tuples."""
        self.tree = intervaltree.IntervalTree.from_tuples(
            [(var.start, var.end + 1, var) for var in self.variants]
        )

    def intervaltree_overlap(self, start: int = 99_100, end: int = 200_000) -> None:
        """Query the interval tree for overlapping variations."""
        self.query_start = start
        self.query_end = end
        self.overlaps = self.tree.overlap(start, end + 1)
        # assert len(overlaps) == 1005, "Incorrect number of overlaps"

    def check_overlap_correctness(self) -> None:
        """Report variations that do not overlap the recorded query."""
        for _left, _right, var in sorted(self.overlaps):
            overlaps = (
                (var.start >= self.query_start and var.start <= self.query_end)
                or (var.end >= self.query_start and var.end <= self.query_end)
                or (var.start < self.query_start and var.start > self.query_end)
            )
            if not overlaps:
                print(f"Overlap error: {self.query_start}-{self.query_end}, {var}")


class RustIntTree:
    """Benchmark fixture for the Rust-backed ``intervaltree_rs`` package.

    Attributes:
        num_vars (int): Number of generated variations.
        variants (list[Variation]): Generated variations to index.
    """

    def __init__(self, num_vars: int = 100) -> None:
        """Create a fixture containing ``num_vars`` variations."""
        self.num_vars = num_vars
        self.variants = create_variation_objects(num_vars)
        self.tree: Any = None
        self.overlaps: Any = None
        self.query_start = 0
        self.query_end = 0

    def intervaltree_create(self) -> None:
        """Build an interval tree by adding one variation at a time."""
        tree = intervaltree_rs.IntervalTree()
        for var in self.variants:
            tree.insert((var.start, var.end + 1, var))
        self.tree = tree

    def intervaltree_create_from_tuples(self) -> None:
        """Build an interval tree from coordinate tuples."""
        self.tree = intervaltree_rs.IntervalTree.from_tuples(
            [(var.start, var.end + 1, var) for var in self.variants]
        )

    def intervaltree_overlap(self, start: int = 99_100, end: int = 200_000) -> None:
        """Query the interval tree for overlapping variations."""
        self.query_start = start
        self.query_end = end
        self.overlaps = self.tree.search(start + 1, end, inclusive=True)
        # assert len(overlaps) == 1005, "Incorrect number of overlaps"

    def check_overlap_correctness(self) -> None:
        """Report variations that do not overlap the recorded query."""
        for _left, _right, var in sorted(self.overlaps):
            overlaps = (
                (var.start >= self.query_start and var.start <= self.query_end)
                or (var.end >= self.query_start and var.end <= self.query_end)
                or (var.start < self.query_start and var.start > self.query_end)
            )
            if not overlaps:
                print(f"Overlap error: {self.query_start}-{self.query_end}, {var}")


class AIListTree:
    """Benchmark fixture for an augmented interval list.

    Attributes:
        num_vars (int): Number of generated variations.
        variants (list[Variation]): Generated variations to index.
    """

    def __init__(self, num_vars: int = 100) -> None:
        """Create a fixture containing ``num_vars`` variations."""
        self.num_vars = num_vars
        self.variants = create_variation_objects(num_vars)
        self.tree: Any = None
        self.overlaps: Any = None
        self.query_start = 0
        self.query_end = 0

    def intervaltree_create(self) -> None:
        """Build an augmented interval list one variation at a time."""
        tree = AIList()
        for i, var in enumerate(self.variants):
            tree.add(var.start, var.end + 1, i)
        self.tree = tree

    def intervaltree_create_from_tuples(self) -> None:
        """Build an augmented interval list from coordinate arrays."""
        starts = np.asarray([var.start for var in self.variants])
        ends = np.asarray([var.end for var in self.variants])
        index = np.arange(len(self.variants))
        tree = AIList()
        self.tree = tree.from_array(starts, ends, index)

    def intervaltree_overlap(self, start: int = 99_100, end: int = 200_000) -> None:
        """Query the augmented interval list for overlapping variations."""
        self.query_start = start
        self.query_end = end
        self.overlaps = [self.variants[idx] for idx in self.tree.intersect_ids(start, end + 1)]
        # assert len(overlaps) == 1005, "Incorrect number of overlaps"

    def nearest(self, start: int = 99_100, end: int = 200_000) -> None:
        """Query the nearest indexed variation."""
        self.tree.closest(start=start, end=end + 1, k=1)

    def check_overlap_correctness(self) -> None:
        """Report variations that do not overlap the recorded query."""
        for idx in self.overlaps:
            var = self.variants[idx]
            overlaps = (
                (var.start >= self.query_start and var.start <= self.query_end)
                or (var.end >= self.query_start and var.end <= self.query_end)
                or (var.start < self.query_start and var.start > self.query_end)
            )
            if not overlaps:
                print(f"Overlap error: {self.query_start}-{self.query_end}, {var}")


def bench(
    cls: type[PyIntTree] | type[RustIntTree] | type[AIListTree],
    num_vars: int,
    n: int = 10_000_000,
) -> tuple[object, str, str, str, str, str, str]:
    """Benchmark construction and queries for one interval-index class."""
    obj = cls(num_vars)
    times = timeit.repeat("obj.intervaltree_create()", number=n, repeat=5, globals={"obj": obj})
    min_setup_s = min(times) / n * 1e6

    times = timeit.repeat(
        "obj.intervaltree_create_from_tuples()",
        number=n,
        repeat=5,
        globals={"obj": obj},
    )
    min_setup_tuple_s = min(times) / n * 1e6

    obj.intervaltree_create()
    times = timeit.repeat(
        "obj.intervaltree_overlap()",
        number=n,
        repeat=5,
        globals={"obj": obj},
    )
    min_query_s = min(times) / n * 1e6

    min_nearest_s = None
    if cls.__name__ == "AIListTree":
        times = timeit.repeat(
            "obj.nearest()",
            number=n,
            repeat=5,
            globals={"obj": obj},
        )
        min_nearest_s = min(times) / n * 1e6

    obj.intervaltree_overlap()
    print(f"{cls.__name__:22}:", f"{'tree setup':22}", f"{min_setup_s:.2f} µs/call")
    print(f"{cls.__name__:22}:", f"{'tree setup (tuples)':22}", f"{min_setup_tuple_s:.2f} µs/call")
    print(f"{cls.__name__:22}:", f"{'overlap calc':22}", f"{min_query_s:.2f} µs/call")
    if min_nearest_s:
        print(f"{cls.__name__:22}:", f"{'nearest calc':22}", f"{min_nearest_s:.2f} µs/call")
    return (
        obj,
        cls.__name__,
        str(num_vars),
        str(n),
        f"{min_setup_s}",
        f"{min_setup_tuple_s}",
        f"{min_query_s}",
    )


def write_out_memory_stats(snapshot: Any, class_name: str, args: argparse.Namespace) -> None:
    """Write line-level memory statistics for one implementation."""
    stats = snapshot.statistics("lineno")
    with open(f"{args.memory_file}-{class_name}.txt", mode="w") as mem_fh:
        print(f"Memory use: {class_name}", file=mem_fh)
        for stat in stats:
            print(f"{stat}", file=mem_fh)


def main(args: argparse.Namespace) -> None:
    """Run the configured interval-index benchmarks."""
    print("Benchmarking IntervalTree creation and overlap calc:")
    now = datetime.now()
    print(f"Starting, {now.strftime('%d/%m/%Y %H:%M:%S')}")

    with open(args.output_file, mode="w") as out_fh:
        print(
            "\t".join(
                [
                    "Implementation",
                    "NumVariants",
                    "NumTimes",
                    "Min[Setup] (µs)",
                    "Min[Setup - Tuples] (µs)",
                    "Min[Overlap] (µs)",
                ]
            ),
            file=out_fh,
        )
        for num_vars in args.num_vars:
            print(f"Number of variants in tree = {num_vars}\nNumber of iterations = {args.times}")
            now = datetime.now()
            print(f"Time = {now.strftime('%d/%m/%Y %H:%M:%S')}")
            _py_tree, *results = bench(PyIntTree, num_vars, n=args.times)
            print("\t".join(results), file=out_fh)

            _rust_tree, *results = bench(RustIntTree, num_vars, n=args.times)
            print("\t".join(results), file=out_fh)

            _ai_list, *results = bench(AIListTree, num_vars, n=args.times)
            print("\t".join(results), file=out_fh)

    now = datetime.now()
    print(f"Finished, {now.strftime('%d/%m/%Y %H:%M:%S')}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script to benchmark different implementation of interval trees",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--output_file",
        metavar="FILE",
        type=str,
        default="bench/interval_tree_bench.txt",
        help="File name for the output file",
    )
    parser.add_argument(
        "--memory_file",
        metavar="MEM_FILE",
        type=str,
        default="bench/interval_tree_bench_mem",
        help="Base name for the memory output file",
    )
    parser.add_argument(
        "--num_vars", type=int, nargs="*", default=[10, 100], help="Set the number of variants to test"
    )
    parser.add_argument(
        "--times", type=int, default=1000, help="Set the number of times to run each expression"
    )
    parser.add_argument("--correctness_check", action="store_true", help="Just do correctness checks")
    parser.add_argument("--debug", action="count", default=0, help="Prints debugging information")
    params = parser.parse_args()
    if params.debug > 0:
        print(params)
    if params.correctness_check:
        check_overlap_correctness()
    else:
        main(params)
