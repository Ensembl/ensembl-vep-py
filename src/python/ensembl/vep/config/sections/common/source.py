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
"""Where variants and features are read from."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import ClassVar

from attrs import define, field

from ensembl.vep.common.constants import DEFAULT_DISTANCE
from ensembl.vep.config.option_spec import option
from ensembl.vep.config.sections.base import ConfigSection

__all__ = ["SourceConfig"]


def _as_distances(value: int | str | Iterable[int]) -> tuple[int, ...]:
    """Return a distance given as one number, a comma-separated string, or a pair, as a tuple."""

    if isinstance(value, str):
        return tuple(int(piece.strip()) for piece in value.split(",") if piece.strip())
    if isinstance(value, int):
        return (value,)
    return tuple(int(piece) for piece in value)


@define(frozen=True, kw_only=True, auto_attribs=True)
class SourceConfig(ConfigSection):
    """Input variants and the annotation sources to overlap them with.

    Attributes:
        input_file (Path | None): Variant input, or None to use a synthetic corpus instead.
        format (str | None): Input format, or None to infer from the file extension. Perl VEP
            spells it ``--format``; ``--input_format`` is accepted as a synonym.
        gff (Path | None): GFF3 file of transcripts.
        gtf (Path | None): GTF file of transcripts.
        fasta (Path | None): Reference FASTA, required alongside ``gff``/``gtf`` for sequence.
        synonyms (Path | None): File of sequence-region synonyms, one region per line: a reference
            name followed by its synonyms, whitespace-separated. Perl VEP's ``--synonyms`` format.
        cache (Path | None): VEP cache directory to annotate from. Accepted so the option exists
            where Perl VEP has it; no reader is implemented yet.
        database (bool): Annotate from an Ensembl database. Accepted, not implemented yet.
        custom (tuple[str, ...]): Custom annotation files, in Perl VEP's ``--custom`` spelling.
            Accepted, not implemented yet.
        distance (tuple[int, ...]): How far outside a feature a variant may fall and still be
            paired with it, as Perl VEP's ``--distance``: one value for both directions, or
            ``upstream,downstream``. ``--flank`` is accepted as a synonym. Read through
            ``upstream_distance`` and ``downstream_distance`` rather than directly.
    """

    section_name: ClassVar[str] = "source"

    input_file: Path | None = field(default=None, metadata=option("variant input file", short="i"))
    format: str | None = field(
        default=None,
        metadata=option(
            "input format, inferred from the file extension when omitted",
            synonyms=("input_format",),
            choices=("vcf",),
        ),
    )
    gff: Path | None = field(default=None, metadata=option("GFF3 file of transcripts"))
    gtf: Path | None = field(default=None, metadata=option("GTF file of transcripts"))
    fasta: Path | None = field(default=None, metadata=option("reference FASTA"))
    synonyms: Path | None = field(
        default=None, metadata=option("file of sequence-region synonyms, one region per line")
    )
    cache: Path | None = field(
        default=None, metadata=option("VEP cache directory to annotate from (not implemented yet)")
    )
    database: bool = field(
        default=False, metadata=option("annotate from an Ensembl database (not implemented yet)")
    )
    custom: tuple[str, ...] = field(
        default=(), metadata=option("custom annotation files (not implemented yet)", metavar="SPEC")
    )
    distance: tuple[int, ...] = field(
        default=(DEFAULT_DISTANCE,),
        converter=_as_distances,
        metadata=option(
            "distance outside a feature still counted as overlapping, as one value or up,down",
            synonyms=("flank",),
            metavar="BP[,DOWNSTREAM]",
        ),
    )

    @property
    def upstream_distance(self) -> int:
        """Bases before a feature's 5' end that still pair with it."""

        return self.distance[0]

    @property
    def downstream_distance(self) -> int:
        """Bases past a feature's 3' end that still pair with it, the upstream value when unset."""

        return self.distance[-1]

    def __attrs_post_init__(self) -> None:
        """Reject a distance that names neither one nor two non-negative values."""

        if not 1 <= len(self.distance) <= 2:
            raise ValueError("distance takes one value or upstream,downstream")
        if any(value < 0 for value in self.distance):
            raise ValueError("distance must be zero or greater")

    def has_features(self) -> bool:
        """Returns whether a transcript source has been configured.

        Returns:
            bool: True if either ``gff`` or ``gtf`` is set.
        """
        return self.gff is not None or self.gtf is not None
