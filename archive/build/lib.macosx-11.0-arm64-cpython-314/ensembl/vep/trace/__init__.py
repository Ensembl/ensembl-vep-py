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
"""Explainable derivations: what a consequence annotator did for one overlap, and why.

Promoted out of ``annotators/consequence/`` because the fact layer, the rule layer and the caches
all touch it, and a cross-cutting concern owned by one of them reads as belonging to that one. It
sits below ``annotators`` in the dependency order, which is also what lets ``caches/`` name a fact
origin without reaching up into the annotators.

Six modules, in dependency order:

* ``events`` -- the value types, and ``TraceRecorder``, the only thing the live run writes to.
* ``annotation`` -- the finished, immutable ``AnnotationTrace`` and how it is rendered.
* ``reconstruct`` -- rebuilds one from a context the run has finished with.
* ``run`` -- folds many into what a whole run cost.
* ``explain`` -- diffs two of them, or one against terms somebody else assigned.
* ``summary`` -- describes the overlap a trace is about: variant, transcript, codons, whole CDS.

Only the first two are re-exported here. The others import the fact, rule or output layers, which
import this package, so importing them from here would close the cycle -- reach them by module
(``from ensembl.vep.trace import reconstruct``) instead.
"""

from ensembl.vep.trace.annotation import (
    AnnotationTrace,
    describe_overlap,
    format_explanation,
    maybe_start_recorder,
    tracing_enabled,
)
from ensembl.vep.trace.events import FactEvent, FactOrigin, PruneEvent, RuleEvent, TraceRecorder

__all__ = [
    "AnnotationTrace",
    "FactEvent",
    "FactOrigin",
    "PruneEvent",
    "RuleEvent",
    "TraceRecorder",
    "describe_overlap",
    "format_explanation",
    "maybe_start_recorder",
    "tracing_enabled",
]
