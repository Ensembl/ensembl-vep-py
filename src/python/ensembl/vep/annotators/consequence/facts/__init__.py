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
"""Cost-tagged facts used by consequence rules.

Importing this package populates ``REGISTRY`` through ``register_fact`` decorators.

``FactContext`` remains in ``consequence.context`` beside ``AnnotatorContext`` because the two
types differ by lifetime rather than by responsibility.
"""

from __future__ import annotations

from .base import *

# Imported for their fact-registration side effects; also re-exported for direct use.
# helper first: kind/position/sequence all import from it.
from .helper import *
from .kind import *
from .position import *
from .sequence import *
