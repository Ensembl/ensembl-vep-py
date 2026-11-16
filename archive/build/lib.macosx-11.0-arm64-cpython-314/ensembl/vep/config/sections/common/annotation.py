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
"""Annotation strategy and cache options.

The default is the columnar staged strategy with every cache off. The zone map, codon and
mapper caches are all opt-in, because their benefit depends on the workload -- the zone map
pays for itself once a transcript is touched by more than one variant, which dense or targeted
panels do and sparse genome-wide input does not. ``readthrough`` changes the modelled answer;
the remaining options only affect the cost of reaching it.
"""

from __future__ import annotations

from typing import ClassVar

from attrs import define, field

from ensembl.vep.config.option_spec import option
from ensembl.vep.config.sections.base import ConfigSection

__all__ = ["AnnotationConfig"]

READTHROUGH_LEVELS = ("none", "utr", "beyond_transcript")

# Mirrors ensembl.vep.caches.zone_map.DEFAULT_ZONE_MAP_BYTES (the source of truth every other
# call site imports). Not imported here so this config layer stays free of the domain-layer
# dependencies caches.zone_map pulls in transitively (annotators.consequence.*, models.*) -- the
# same reason _annotator_names() below resolves its own dependency lazily rather than at import
# time.
_DEFAULT_ZONE_MAP_BYTES = 256 * 1024 * 1024


def _annotator_names() -> tuple[str, ...]:
    """Returns the registered annotator names.

    Resolved lazily so this module never imports ensembl.vep.annotators, which would make the
    configuration layer depend on the domain it configures.

    Returns:
        tuple[str, ...]: Every registered annotator name, sorted.
    """
    # pylint: disable=import-outside-toplevel
    from ensembl.vep.annotators.registry import available_annotators

    return available_annotators()


@define(frozen=True, kw_only=True, auto_attribs=True)
class AnnotationConfig(ConfigSection):
    """Which annotator to run and which caches to give it.

    Attributes:
        annotator (str): Registered annotator name. Defaults to the columnar staged strategy.
        batch_size (int): Overlaps per chunk, for the batch-capable strategies.
        zone_map (bool): Enable the per-transcript positional-fact zone map. Disabled by default;
            worth enabling for inputs with many variants per transcript.
        zone_map_bytes (int): Byte budget for the zone map cache.
        codon_cache (bool): Enable the shared-codon edit cache.
        codon_cache_entries (int): Maximum codon cache entries.
        mapper_cache (bool): Enable the genomic-to-CDS mapper index cache.
        mapper_cache_bytes (int): Byte budget for the mapper index cache.
        readthrough (str): How far translation may run past a destroyed stop codon.
        early_subsumption_skip (bool): Skip rules already suppressed by a hit, for the staged
            strategy. Measured a small net loss; off by default.
        fact_cost_override (bool): Use the packaged measured fact costs instead of declared tiers.
    """

    section_name: ClassVar[str] = "annotation"

    annotator: str = field(
        default="consequence.columnar-staged",
        metadata=option("annotator to run", choices=_annotator_names, metavar="NAME"),
    )
    batch_size: int = field(default=256, metadata=option("overlaps per batch", metavar="N"))
    zone_map: bool = field(
        default=False,
        metadata=option("enable the positional zone map for inputs with many variants per transcript"),
    )
    zone_map_bytes: int = field(
        default=_DEFAULT_ZONE_MAP_BYTES, metadata=option("zone map byte budget", metavar="BYTES")
    )
    codon_cache: bool = field(default=False, metadata=option("enable the shared-codon cache"))
    codon_cache_entries: int = field(default=65536, metadata=option("codon cache size", metavar="N"))
    mapper_cache: bool = field(default=False, metadata=option("enable the mapper index cache"))
    mapper_cache_bytes: int = field(
        default=64 * 1024 * 1024, metadata=option("mapper cache byte budget", metavar="BYTES")
    )
    readthrough: str = field(
        default="none",
        metadata=option("how far to translate past a destroyed stop codon", choices=READTHROUGH_LEVELS),
    )
    early_subsumption_skip: bool = field(
        default=False, metadata=option("skip rules already suppressed by a hit (staged only)")
    )
    fact_cost_override: bool = field(
        default=False,
        metadata=option("use measured fact costs", env="ENSEMBL_VEP_FACT_COST_OVERRIDE"),
    )

    def __attrs_post_init__(self) -> None:
        """Reject batch and enabled-cache limits that cannot be used."""

        if self.batch_size < 1:
            raise ValueError("batch_size must be at least one")
        if self.zone_map and self.zone_map_bytes < 1:
            raise ValueError("zone_map_bytes must be at least one when zone_map is enabled")
        if self.codon_cache and self.codon_cache_entries < 1:
            raise ValueError("codon_cache_entries must be at least one when codon_cache is enabled")
        if self.mapper_cache and self.mapper_cache_bytes < 1:
            raise ValueError("mapper_cache_bytes must be at least one when mapper_cache is enabled")
