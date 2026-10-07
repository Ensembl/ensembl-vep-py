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
"""Tests for the section declarations."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from ensembl.vep.config import AnnotationConfig, Config, ExecutionConfig, LoggingConfig, SourceConfig
from ensembl.vep.config.option_spec import section_options
from ensembl.vep.config.parser.common import resolve_option_key
from ensembl.vep.config.sections import BenchEndToEndConfig, BenchMatrixConfig
from ensembl.vep.config.sections.base import ConfigSection


class TestSectionBehaviour:
    """Convenience methods on sections."""

    def test_quiet_wins_over_verbose(self) -> None:
        """Quiet wins over verbose."""

        assert LoggingConfig(verbose=True, quiet=True).effective_level() == "error"

    def test_verbose_lowers_the_level(self) -> None:
        """Verbose lowers the log level."""

        assert LoggingConfig(verbose=True).effective_level() == "debug"

    def test_an_explicit_level_is_used_when_neither_shorthand_is_set(self) -> None:
        """An explicit level is used as-is."""

        assert LoggingConfig(log_level="warning").effective_level() == "warning"

    @pytest.mark.parametrize(
        "kwargs,expected",
        [({}, False), ({"gff": Path("a.gff3")}, True), ({"gtf": Path("a.gtf")}, True)],
    )
    def test_a_transcript_source_is_detected(self, kwargs: dict, expected: bool) -> None:
        """Either format counts as features."""

        assert SourceConfig(**kwargs).has_features() is expected

    def test_fork_is_a_spelling_of_processes_not_a_field(self) -> None:
        """The Perl-compatible spelling resolves to the one process-count field."""
        spelling_owner = resolve_option_key("fork", (ExecutionConfig,))

        assert spelling_owner == ("execution", "processes")
        assert not hasattr(ExecutionConfig(), "fork")

    def test_transcript_pool_bytes_zero_is_allowed(self) -> None:
        """A zero byte budget disables pooling rather than being rejected as invalid."""
        assert ExecutionConfig(transcript_pool_bytes=0).transcript_pool_bytes == 0

    def test_transcript_pool_bytes_default_matches_the_caches_module(self) -> None:
        """The config layer keeps its own literal (to avoid a domain import), so it must agree
        with ``caches.feature_pool.DEFAULT_TRANSCRIPT_POOL_BYTES``, the source of truth every
        other call site imports."""
        feature_pool = pytest.importorskip("ensembl.vep.caches.feature_pool")

        assert ExecutionConfig().transcript_pool_bytes == feature_pool.DEFAULT_TRANSCRIPT_POOL_BYTES

    def test_annotation_defaults_select_columnar_staging_with_every_cache_off(self) -> None:
        """Select columnar staging by default; every cache, including the zone map, is opt-in."""

        annotation = AnnotationConfig()

        assert annotation.annotator == "consequence.columnar-staged"
        assert not annotation.zone_map
        assert not annotation.codon_cache
        assert not annotation.mapper_cache

    @pytest.mark.parametrize(
        "config_type,kwargs,message",
        [
            (AnnotationConfig, {"batch_size": 0}, "batch_size"),
            (
                AnnotationConfig,
                {"zone_map": True, "zone_map_bytes": 0},
                "zone_map_bytes",
            ),
            (
                AnnotationConfig,
                {"codon_cache": True, "codon_cache_entries": 0},
                "codon_cache_entries",
            ),
            (
                AnnotationConfig,
                {"mapper_cache": True, "mapper_cache_bytes": 0},
                "mapper_cache_bytes",
            ),
            (SourceConfig, {"distance": -1}, "distance"),
            (ExecutionConfig, {"processes": 0}, "processes"),
            (ExecutionConfig, {"threads": 0}, "threads"),
            (ExecutionConfig, {"buffer_size": 0}, "buffer_size"),
            (ExecutionConfig, {"annotation_window": 0}, "annotation_window"),
            (ExecutionConfig, {"max_in_flight": 0}, "max_in_flight"),
            (ExecutionConfig, {"transcript_pool_bytes": -1}, "transcript_pool_bytes"),
        ],
    )
    def test_invalid_worker_settings_are_rejected_before_execution(
        self,
        config_type: type[AnnotationConfig] | type[SourceConfig] | type[ExecutionConfig],
        kwargs: dict[str, object],
        message: str,
    ) -> None:
        """Invalid worker inputs fail while the parent constructs configuration."""

        with pytest.raises(ValueError, match=message):
            config_type(**kwargs)


class TestBenchEndToEndSection:
    """The end-to-end benchmark's own options."""

    def test_a_single_process_run_is_the_default(self) -> None:
        """One worker, so the default run measures the pipeline rather than the scheduler."""
        assert BenchEndToEndConfig().processes == (1,)

    def test_the_default_cache_sets_match_the_matrix_verb(self) -> None:
        """A configuration label has to mean the same run in both benchmark tables."""
        assert BenchEndToEndConfig().caches == BenchMatrixConfig().caches

    @pytest.mark.parametrize(
        ("kwargs", "message"),
        [
            pytest.param({"processes": ()}, "at least one worker count", id="no_processes"),
            pytest.param({"processes": (1, 0)}, "process count", id="zero_process_count"),
            pytest.param({"batch_size": 0}, "batch_size", id="zero_batch_size"),
        ],
    )
    def test_unrunnable_counts_are_rejected(self, kwargs: dict[str, Any], message: str) -> None:
        """A configuration that could not be run fails while it is built."""
        with pytest.raises(ValueError, match=message):
            BenchEndToEndConfig(**kwargs)


class TestEverySection:
    """Properties every section must satisfy."""

    @pytest.mark.parametrize(
        "section_type", list(Config.section_types.values()), ids=lambda t: t.section_name
    )
    def test_a_section_is_constructible_with_no_arguments(self, section_type: type[ConfigSection]) -> None:
        """Every section can be built with no arguments."""
        assert section_type().section_name == section_type.section_name

    @pytest.mark.parametrize(
        "section_type", list(Config.section_types.values()), ids=lambda t: t.section_name
    )
    def test_every_option_declares_help_text(self, section_type: type[ConfigSection]) -> None:
        """Every option declares help text."""
        for field_name, (_field, spec) in section_options(section_type).items():
            assert spec.help, f"{section_type.section_name}.{field_name}"

    @pytest.mark.parametrize(
        "section_type", list(Config.section_types.values()), ids=lambda t: t.section_name
    )
    def test_every_section_names_itself(self, section_type: type[ConfigSection]) -> None:
        """Every section names itself."""
        assert section_type.section_name

    def test_section_names_are_unique(self) -> None:
        """Section names stay unique."""
        names = [section_type.section_name for section_type in Config.section_types.values()]
        assert len(names) == len(set(names))
