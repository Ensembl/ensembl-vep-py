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
"""The consequence annotator: what's true about a variant, which SO term follows, and how to work it out.

Three layers, in dependency order:

* `facts` -- cost-tagged boolean facts over one (variant, feature) overlap, plus the `context`
  they memoise into.
* `rules` -- SO terms as declarative boolean formulas over fact atoms, compiled from the
  checked-in consequence catalogue (`rules.catalogue`); see that submodule's `__init__.py` for why
  the loader lives there rather than as a third sibling package.
* `strategies` -- interchangeable ways of evaluating those rules against those facts. Every
  strategy must produce *identical* output; `consequence.linear` is the conformance reference.

Importing this package registers every fact and every ``consequence.*`` annotator as an import
side effect. `strategies` is imported last because it depends on the two layers above it.
"""

from . import facts
from . import rules
from . import strategies
