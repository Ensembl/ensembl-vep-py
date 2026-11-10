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
"""Traced annotation of one overlap."""

from __future__ import annotations

from typing import Any

from attrs import asdict

from ensembl.vep.annotators.factory import build_annotator
from ensembl.vep.config import Config
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.trace import AnnotationTrace, annotation_tracing

__all__ = ["NO_TRACE", "annotate_with_trace"]

NO_TRACE = (
    "the native rust.* annotators never record a trace. That is not a gap: every annotator is "
    "proven byte-identical to consequence.linear by the conformance suite, so the same overlap "
    "run through consequence.linear has the same derivation to show."
)


def annotate_with_trace(
    config: Config, annotator_name: str, overlap: TranscriptVariationAllele
) -> tuple[list[str], AnnotationTrace | None]:
    """Run one annotator over one overlap with tracing on.

    Args:
        config (Config): The run's configuration, for the annotator and cache options.
        annotator_name (str): A registered annotator name.
        overlap (TranscriptVariationAllele): The overlap to annotate.

    Returns:
        tuple[list[str], AnnotationTrace | None]: The emitted terms, and the derivation -- None
            for a native rust.* annotator, which never traces.
    """
    options: dict[str, Any] = asdict(config.annotation)
    options["annotator"] = annotator_name
    annotator = build_annotator(
        **options,
        upstream_distance=config.source.upstream_distance,
        downstream_distance=config.source.downstream_distance,
    )
    with annotation_tracing(quiet=config.logging.effective_level().upper() != "DEBUG"):
        consequences = annotator.annotate_feature_overlap(overlap)
    return [consequence.term for consequence in consequences], getattr(annotator, "last_trace", None)
