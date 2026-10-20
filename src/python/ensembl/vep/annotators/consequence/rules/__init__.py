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
"""Boolean consequence rules over fact atoms.

The checked-in catalogue belongs here because only rule construction, validation, and editing
consume it; the fact layer remains independent of catalogue metadata.
"""

from __future__ import annotations

from .catalogue import *
from .expressions import *
from .rulesets import *
from .validation import *
