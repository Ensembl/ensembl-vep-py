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
"""Golden corpora and pinned annotation terms."""

from .base import GoldenCorpus
from .resolve import load_golden_corpus
from .terms import FORMAT_VERSION, GoldenReport, GoldenTerms, overlap_key

__all__ = [
    "FORMAT_VERSION",
    "GoldenCorpus",
    "GoldenReport",
    "GoldenTerms",
    "load_golden_corpus",
    "overlap_key",
]
