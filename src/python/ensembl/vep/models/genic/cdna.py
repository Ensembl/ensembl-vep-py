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
"""cDNA class."""

from __future__ import annotations

from typing import TYPE_CHECKING
from attrs import define, field

from ensembl.vep.models import Feature

if TYPE_CHECKING:
    from ensembl.vep.models.genic import Transcript, Translation


@define(kw_only=True, auto_attribs=True)
class cDNA(Feature):
    """
    Represents a CDS.

    Attributes:
        transcript: The transcript to which this cDNA belongs.
    """

    transcript: Transcript | None = field(default=None, init=False)
