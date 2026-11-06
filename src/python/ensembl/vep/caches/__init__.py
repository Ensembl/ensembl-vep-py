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
"""Cache implementations."""

from __future__ import annotations

from .binned_cache import BinnedCache, BinnedFeatureCache
from .codon_cache import (
    CodonCache,
    CodonEdit,
    CodonEditKey,
    CodonEditLRUCache,
    CodonEditLoader,
)
from .feature_pool import (
    DEFAULT_TRANSCRIPT_POOL_BYTES,
    FeatureKey,
    FeaturePool,
    TranscriptKey,
    TranscriptPool,
)
from .interval_tree import (
    FeatureIntervalTreeCache,
    Interval,
    IntervalTree,
    IntervalTreeCache,
    IntervalTreeKey,
)
from .mapper_cache import MapperIndexCache, MapperIndexLRUCache, MapperIndexLoader
from .sequence_cache import SequenceBinnedLRUCache, SequenceCache
from .slice_cache import SliceCache, SliceKey, WeakSliceCache
from .zone_map import SERVED_FACT_NAMES, ZoneMap, ZoneMapCache

__all__ = [
    "BinnedCache",
    "BinnedFeatureCache",
    "CodonCache",
    "CodonEdit",
    "CodonEditKey",
    "CodonEditLRUCache",
    "CodonEditLoader",
    "DEFAULT_TRANSCRIPT_POOL_BYTES",
    "FeatureIntervalTreeCache",
    "FeatureKey",
    "FeaturePool",
    "Interval",
    "IntervalTree",
    "IntervalTreeCache",
    "IntervalTreeKey",
    "MapperIndexCache",
    "MapperIndexLRUCache",
    "MapperIndexLoader",
    "SequenceCache",
    "SequenceBinnedLRUCache",
    "SERVED_FACT_NAMES",
    "TranscriptKey",
    "TranscriptPool",
    "ZoneMap",
    "ZoneMapCache",
    "SliceCache",
    "SliceKey",
    "WeakSliceCache",
]
