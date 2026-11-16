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
"""The boolean-rules layer: SO terms as formulas over fact atoms.

``catalogue`` (the checked-in consequence catalogue loader) lives here rather than as a sibling
of ``facts``/``rules``: nothing in ``facts/`` reads it -- the class-membership facts it used to
feed (``feature_kind_of``, before the ``exclude_kind`` denylist was replaced with ordinary
formula conjuncts) are gone -- and every remaining consumer (``rulesets.build_fact_ruleset``,
``validation``'s catalogue-drift checks) is already inside ``rules/``.
"""

from .catalogue import *
from .expressions import *
from .rulesets import *
from .validation import *
