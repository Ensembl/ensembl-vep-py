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
"""Benchmarking: the conformance gate, timing, the strategy matrix, sweeps, results and the Perl oracle.

Measurement only. The corpus a benchmark runs over is ``ensembl.vep.corpus``'s, not something this
package owns; the runners in ``ensembl.vep.runner.bench`` are the command line on top.
"""

from __future__ import annotations
