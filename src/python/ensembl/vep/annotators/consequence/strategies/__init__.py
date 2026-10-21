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
"""Interchangeable ways to evaluate the consequence rule set.

``ConsequenceAnnotator`` is the shared implementation each strategy subclasses. The modules
beside it are the strategies themselves. Importing this package registers every ``consequence.*``
annotator as an import side effect.

Registry names represent evaluation strategies rather than strategy-cache combinations. Callers
configure zone-map and codon caches independently through the annotator context.

* ``linear`` -- rank-ordered scan with a cheap-guard digest cache. The conformance reference.
* ``staged`` -- cheap facts first, then ``requires`` propagation; expensive tier only on demand.
* ``columnar`` -- bit-parallel over a batch; cuts dispatch, not fact cost.
* ``columnar-staged`` -- both of the above: bit-parallel with forced-false staging masks.
* ``columnar-kernels`` -- columnar staging with selected columns filled by native batch kernels.

What was ``regular.columnar-staged-zonemap`` is now
``ColumnarStagedAnnotator(zone_map=ZoneMapCache(), codon_cache=CodonEditLRUCache())``.
"""

from __future__ import annotations

from .base import *
from .linear import *
from .staged import *
from .columnar import *
from .columnar_staged import *
from .columnar_kernels import *
