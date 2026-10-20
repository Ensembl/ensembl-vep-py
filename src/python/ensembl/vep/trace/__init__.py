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
"""Derivation traces for consequence annotations.

A trace records the facts, rule verdicts, subsumption decisions, and emitted terms for one
overlap. The package sits below ``annotators`` because fact evaluation and caches both depend on
its event types.

Only annotation and event types are re-exported here. Reconstruction and aggregation depend on
higher layers that already import this package, so import those features from their modules to
avoid dependency cycles. Comparing and summarising a finished trace is not tracing at all: that
lives in ``ensembl.vep.explain``.
"""

from __future__ import annotations

from ensembl.vep.trace.annotation import AnnotationTrace, annotation_tracing, describe_overlap
from ensembl.vep.trace.events import FactEvent, FactOrigin, PruneEvent, RuleEvent, TraceRecorder

__all__ = [
    "AnnotationTrace",
    "FactEvent",
    "FactOrigin",
    "PruneEvent",
    "RuleEvent",
    "TraceRecorder",
    "annotation_tracing",
    "describe_overlap",
]
