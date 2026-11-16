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
"""Stream VCF rows into sequence-region-safe primitive batches."""

from __future__ import annotations

import gzip
from collections.abc import Iterator
from pathlib import Path

from ensembl.vep.buffer.records import VCFBufferPayload

__all__ = ["count_vcf_records", "input_header_lines", "iter_vcf_buffers", "iter_vcf_columns"]


def count_vcf_records(path: Path) -> int:
    """Return the number of variant records without constructing parser objects."""
    with path.open("rb") as probe:
        compressed = probe.read(2) == b"\x1f\x8b"

    if compressed:
        with gzip.open(path, "rb") as handle:
            return sum(1 for line in handle if line.strip() and not line.startswith(b"#"))

    with path.open("rb") as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith(b"#"))


def _text_lines(path: Path) -> Iterator[str]:
    """Yield decoded input lines from plain or gzip-compressed VCF."""
    with path.open("rb") as probe:
        compressed = probe.read(2) == b"\x1f\x8b"

    if compressed:
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            yield from handle
        return

    with path.open("r", encoding="utf-8", newline="") as handle:
        yield from handle


def input_header_lines(path: Path) -> list[str]:
    """Return an input VCF's leading header lines without line endings."""
    headers: list[str] = []
    for line in _text_lines(path):
        if not line.startswith("#"):
            break
        headers.append(line.rstrip("\r\n"))
    return headers


def iter_vcf_columns(path: Path) -> Iterator[tuple[str, ...]]:
    """Yield every input VCF column without normalising its text."""
    for line in _text_lines(path):
        if not line.strip() or line.startswith("#"):
            continue
        columns = tuple(line.rstrip("\r\n").split("\t"))
        if len(columns) < 8:
            raise ValueError("a VCF record must contain at least eight columns")
        yield columns


def iter_vcf_buffers(
    path: Path,
    *,
    batch_size: int,
) -> Iterator[VCFBufferPayload]:
    """Yield ordered VCF batches without mixing sequence regions.

    Args:
        path (Path): Input VCF path.
        batch_size (int): Maximum record count in a same-region batch.

    Yields:
        Iterator[VCFBufferPayload]: Primitive record batches in input order.
    """
    if batch_size < 1:
        raise ValueError("batch_size must be at least one")

    index = 0
    region: str | None = None
    columns: list[tuple[str, ...]] = []
    for record_columns in iter_vcf_columns(path):
        record_region = record_columns[0]
        if region is not None and (record_region != region or len(columns) == batch_size):
            yield VCFBufferPayload(index=index, sequence_region=region, columns=tuple(columns))
            index += 1
            columns = []
        region = record_region
        columns.append(record_columns)
    if region is not None:
        yield VCFBufferPayload(index=index, sequence_region=region, columns=tuple(columns))
