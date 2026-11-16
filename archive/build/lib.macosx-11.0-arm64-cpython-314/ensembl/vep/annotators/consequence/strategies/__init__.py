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

`base.ConsequenceAnnotator` is the shared implementation each strategy subclasses; the modules
beside it are the strategies themselves. Importing this package registers every ``consequence.*``
annotator as an import side effect.

**Strategies only, no cache cross-product.** There used to also be ``regular.zonemap`` and
``regular.columnar-staged-zonemap`` (registered under the pre-rename ``regular.`` namespace),
which were not different ways of evaluating rules at all -- they were ``linear`` and
``columnar-staged`` with a cache switched on. That made the strategy list a cross-product: every
new optional cache would have doubled it, and the names had to spell out the combination. The zone
map and the shared-codon cache are now caller-selected utilities living in
``ensembl/vep/caches/``, passed to any strategy via ``zone_map=``/``codon_cache=`` (or the
annotator context directly), so this list holds only genuinely distinct evaluation strategies:

* ``linear`` -- rank-ordered scan with a cheap-guard digest cache. The conformance reference.
* ``staged`` -- cheap facts first, then ``requires`` propagation; expensive tier only on demand.
* ``columnar`` -- bit-parallel over a batch; cuts dispatch, not fact cost.
* ``columnar-staged`` -- both of the above: bit-parallel with forced-false staging masks.
* ``columnar-kernels`` -- ``columnar-staged`` with the cheap fact columns filled by native batch
  kernels; only flat integer arrays cross FFI.

What was ``regular.columnar-staged-zonemap`` is now
``ColumnarStagedAnnotator(zone_map=ZoneMapCache(), codon_cache=CodonEditLRUCache())``.
"""

from .base import *
from .linear import *
from .staged import *
from .columnar import *
from .columnar_staged import *
from .columnar_kernels import *
