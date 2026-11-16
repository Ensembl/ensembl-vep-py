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
"""The cost-tagged facts layer.

Importing this package is what populates REGISTRY (facts/base.py): each fact module below
registers its facts as an import side effect via the fact decorator (also in facts/base.py).

FactContext is **not** re-exported here. It lives in the sibling ``consequence/context/`` package
alongside AnnotatorContext, because the two contexts are a pair distinguished by lifetime and
splitting them across packages would hide that. Import it from
``ensembl.vep.annotators.consequence.context``.
"""

from .base import *

# Imported for their fact-registration side effects; also re-exported for direct use.
# helper first: kind/position/sequence all import from it.
from .helper import *
from .kind import *
from .position import *
from .sequence import *
