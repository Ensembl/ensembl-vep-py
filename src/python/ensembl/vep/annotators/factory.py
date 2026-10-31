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
"""Build annotators and run-scoped context from plain values.

Registered strategies accept different keyword arguments, so construction derives supported
options from each factory signature instead of maintaining a duplicate table.
"""

from __future__ import annotations

import inspect

from ensembl.vep.annotators.base import Annotator
from ensembl.vep.annotators.consequence.context.annotator import AnnotatorContext
from ensembl.vep.annotators.consequence.facts import apply_packaged_fact_cost_overrides
from ensembl.vep.annotators.registry import annotator_factory, get_annotator
from ensembl.vep.caches.codon_cache import CodonEditLRUCache
from ensembl.vep.caches.mapper_cache import MapperIndexLRUCache
from ensembl.vep.caches.zone_map import DEFAULT_ZONE_MAP_BYTES, ZoneMapCache
from ensembl.vep.common import DEFAULT_DISTANCE, Readthrough
from ensembl.vep.models.overlap import OverlapConsequence

__all__ = ["build_annotator", "build_annotator_context"]


def build_annotator_context(
    *,
    zone_map: bool = False,
    zone_map_bytes: int = DEFAULT_ZONE_MAP_BYTES,
    codon_cache: bool = False,
    codon_cache_entries: int = 65536,
    mapper_cache: bool = False,
    mapper_cache_bytes: int = 64 * 1024 * 1024,
    readthrough: str | Readthrough = Readthrough.NONE,
    fact_cost_override: bool = False,
    upstream_distance: int = DEFAULT_DISTANCE,
    downstream_distance: int = DEFAULT_DISTANCE,
) -> AnnotatorContext:
    """Return the run-scoped context described by a set of cache choices.

    Args:
        zone_map (bool): Enable the per-transcript positional-fact zone map.
        zone_map_bytes (int): Byte budget for the zone map cache.
        codon_cache (bool): Enable the shared-codon edit cache.
        codon_cache_entries (int): Maximum codon cache entries.
        mapper_cache (bool): Enable the genomic-to-CDS mapper index cache.
        mapper_cache_bytes (int): Byte budget for the mapper index cache.
        readthrough (str | Readthrough): How far translation may run past a destroyed stop codon,
            as a level or its lowercased name.
        fact_cost_override (bool): Order fact evaluation by the packaged measured costs instead of the
            declared tiers. The table is process-wide, so this changes every annotator built after it.
        upstream_distance (int): Bases before a feature's 5' end that still count as upstream.
        downstream_distance (int): Bases past a feature's 3' end that still count as downstream.

    Returns:
        AnnotatorContext: The context, with only the enabled caches attached.
    """
    if isinstance(readthrough, str):
        readthrough = Readthrough[readthrough.upper()]
    if fact_cost_override:
        apply_packaged_fact_cost_overrides()
    return AnnotatorContext(
        zone_map_cache=(
            ZoneMapCache(
                max_bytes=zone_map_bytes,
                upstream_distance=upstream_distance,
                downstream_distance=downstream_distance,
            )
            if zone_map
            else None
        ),
        codon_cache=CodonEditLRUCache(max_entries=codon_cache_entries) if codon_cache else None,
        mapper_index_cache=MapperIndexLRUCache(max_bytes=mapper_cache_bytes) if mapper_cache else None,
        readthrough=readthrough,
        upstream_distance=upstream_distance,
        downstream_distance=downstream_distance,
    )


def build_annotator(
    *,
    annotator: str = "consequence.linear",
    batch_size: int = 256,
    early_subsumption_skip: bool = False,
    zone_map: bool = False,
    zone_map_bytes: int = DEFAULT_ZONE_MAP_BYTES,
    codon_cache: bool = False,
    codon_cache_entries: int = 65536,
    mapper_cache: bool = False,
    mapper_cache_bytes: int = 64 * 1024 * 1024,
    readthrough: str | Readthrough = Readthrough.NONE,
    fact_cost_override: bool = False,
    upstream_distance: int = DEFAULT_DISTANCE,
    downstream_distance: int = DEFAULT_DISTANCE,
    annotator_context: AnnotatorContext | None = None,
) -> Annotator[OverlapConsequence]:
    """Return a configured annotator, with its caches attached.

    Args:
        annotator (str): Registered annotator name.
        batch_size (int): Overlaps per chunk, for the strategies that batch.
        early_subsumption_skip (bool): Skip rules already suppressed by a hit, where supported.
        zone_map (bool): Enable the positional-fact zone map.
        zone_map_bytes (int): Byte budget for the zone map cache.
        codon_cache (bool): Enable the shared-codon edit cache.
        codon_cache_entries (int): Maximum codon cache entries.
        mapper_cache (bool): Enable the mapper index cache.
        mapper_cache_bytes (int): Byte budget for the mapper index cache.
        readthrough (str | Readthrough): How far translation may run past a destroyed stop codon.
        fact_cost_override (bool): Use the packaged measured fact costs, when building the context.
        upstream_distance (int): Bases before a feature's 5' end that still count as upstream.
        downstream_distance (int): Bases past a feature's 3' end that still count as downstream.
        annotator_context (AnnotatorContext | None): A context to share across annotators, or None
            to build one from the cache arguments.

    Returns:
        Annotator[OverlapConsequence]: The annotator instance.

    Raises:
        KeyError: If no annotator is registered under ``annotator``.
    """
    context = annotator_context or build_annotator_context(
        zone_map=zone_map,
        zone_map_bytes=zone_map_bytes,
        codon_cache=codon_cache,
        codon_cache_entries=codon_cache_entries,
        mapper_cache=mapper_cache,
        mapper_cache_bytes=mapper_cache_bytes,
        readthrough=readthrough,
        fact_cost_override=fact_cost_override,
        upstream_distance=upstream_distance,
        downstream_distance=downstream_distance,
    )
    candidate_kwargs = {
        "annotator_context": context,
        "batch_size": batch_size,
        "early_subsumption_skip": early_subsumption_skip,
    }
    accepted = inspect.signature(annotator_factory(annotator)).parameters
    return get_annotator(
        annotator,
        **{name: value for name, value in candidate_kwargs.items() if name in accepted},
    )
