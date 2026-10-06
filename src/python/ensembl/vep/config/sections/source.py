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

from pathlib import Path
from typing import ClassVar

from attrs import define, field

from ensembl.vep.config.option_spec import option
from ensembl.vep.config.sections.base import ConfigSection

__all__ = ["SourceConfig"]


@define(frozen=True, kw_only=True, auto_attribs=True)
class SourceConfig(ConfigSection):
    """Input variants and the annotation sources to overlap them with.

    Attributes:
        input_file (Path | None): Variant input, or None to use a synthetic corpus instead.
        input_format (str | None): Input format, or None to infer from the file extension.
        gff (Path | None): GFF3 file of transcripts.
        gtf (Path | None): GTF file of transcripts.
        fasta (Path | None): Reference FASTA, required alongside ``gff``/``gtf`` for sequence.
        flank (int): How far outside a feature a variant may fall and still be paired with it.
    """

    section_name: ClassVar[str] = "source"

    input_file: Path | None = field(default=None, metadata=option("variant input file"))
    input_format: str | None = field(
        default=None,
        metadata=option("input format, inferred from the file extension when omitted", choices=("vcf",)),
    )
    gff: Path | None = field(default=None, metadata=option("GFF3 file of transcripts"))
    gtf: Path | None = field(default=None, metadata=option("GTF file of transcripts"))
    fasta: Path | None = field(default=None, metadata=option("reference FASTA"))
    flank: int = field(
        default=5000, metadata=option("distance outside a feature still counted as overlapping")
    )

    def __attrs_post_init__(self) -> None:
        """Reject a negative feature flank."""

        if self.flank < 0:
            raise ValueError("flank must be zero or greater")

    def has_features(self) -> bool:
        """Returns whether a transcript source has been configured.

        Returns:
            bool: True if either ``gff`` or ``gtf`` is set.
        """
        return self.gff is not None or self.gtf is not None
