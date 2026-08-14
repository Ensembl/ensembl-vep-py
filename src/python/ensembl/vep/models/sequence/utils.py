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
"""Sequence utils"""

from ensembl.vep.common.enums import CoordSystemLevel, Strand, ClassSOTerm
from ensembl.vep.models.sequence import CoordSystem, Slice
from ensembl.vep.sources.sequence_source import SequenceSource

_CHROMOSOME_COORD_SYSTEM = CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME, default=True)
_PLACEHOLDER_SEQ_REGION_END = 2**31 - 1


def get_or_build_slice(
    seq_region_name: str, slice_cache: dict[str, Slice], sequence_source: SequenceSource | None = None
) -> Slice:
    """Return the shared Slice for a sequence region, creating it when needed.

    Args:
        seq_region_name: The sequence-region name to look up.
        slice_cache: Previously built Slice cache keyed by sequence-region name.

    Returns:
        Slice: The shared Slice object for the sequence region.
    """

    if seq_region_name not in slice_cache:
        seq_region_length = None
        if sequence_source is not None:
            seq_region_length = sequence_source.get_seq_region_length(seq_region_name)

        slice_cache[seq_region_name] = Slice(
            coord_system=_CHROMOSOME_COORD_SYSTEM,
            seq_region_name=seq_region_name,
            start=1,
            end=seq_region_length if seq_region_length is not None else _PLACEHOLDER_SEQ_REGION_END,
            strand=Strand.FORWARD,
            seq_region_length=seq_region_length,
            sequence_source=sequence_source,
        )
    return slice_cache[seq_region_name]
