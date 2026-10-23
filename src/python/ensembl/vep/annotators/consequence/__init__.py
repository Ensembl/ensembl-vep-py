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
"""Consequence annotation facts, rules, and strategies.

The package has three dependency-ordered layers:

* ``facts`` defines cost-tagged properties of a variant-feature overlap.
* ``rules`` compiles Sequence Ontology terms into Boolean formulas over facts.
* ``strategies`` provides interchangeable evaluators that produce identical consequences.

Importing this package registers every fact and every ``consequence.*`` annotator as an import
side effect. ``strategies`` is imported last because it depends on facts and rules.
"""

from __future__ import annotations

from . import facts
from . import rules
from . import strategies
