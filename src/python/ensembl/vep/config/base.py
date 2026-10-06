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
"""Configuration roots and layer loading."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, ClassVar

import attrs
from attrs import define, field, fields

from ensembl.vep.config.option_spec import section_options
from ensembl.vep.config.parser.common import resolve_option
from ensembl.vep.config.parser import build_parser, parse_args, read_config_file, read_environment
from ensembl.vep.config.sections import (
    AnnotationConfig,
    BenchConfig,
    BenchEndToEndConfig,
    BenchMatrixConfig,
    BenchOracleConfig,
    BenchSweepConfig,
    CatalogueAddConfig,
    CatalogueConfig,
    CatalogueCostsConfig,
    CatalogueRemoveConfig,
    CatalogueShowConfig,
    CatalogueValidateConfig,
    CatalogueWriteConfig,
    ConfigSection,
    CorpusConfig,
    CorpusDescribeConfig,
    CorpusGoldenConfig,
    CorpusSaveConfig,
    ExecutionConfig,
    ExplainConfig,
    LoggingConfig,
    RunConfig,
    SourceConfig,
    StatsConfig,
)
from ensembl.vep.config.sources import ConfigLayer, ConfigSource, coerce_layer, merge_layers

__all__ = ["Config", "load"]


@define(frozen=True, kw_only=True, auto_attribs=True)
class Config:
    """The assembled configuration and its value sources."""

    annotation: AnnotationConfig = field(factory=AnnotationConfig)
    bench: BenchConfig = field(factory=BenchConfig)
    bench_endtoend: BenchEndToEndConfig = field(factory=BenchEndToEndConfig)
    bench_matrix: BenchMatrixConfig = field(factory=BenchMatrixConfig)
    bench_oracle: BenchOracleConfig = field(factory=BenchOracleConfig)
    bench_sweep: BenchSweepConfig = field(factory=BenchSweepConfig)
    catalogue: CatalogueConfig = field(factory=CatalogueConfig)
    catalogue_add: CatalogueAddConfig = field(factory=CatalogueAddConfig)
    catalogue_costs: CatalogueCostsConfig = field(factory=CatalogueCostsConfig)
    catalogue_remove: CatalogueRemoveConfig = field(factory=CatalogueRemoveConfig)
    catalogue_show: CatalogueShowConfig = field(factory=CatalogueShowConfig)
    catalogue_validate: CatalogueValidateConfig = field(factory=CatalogueValidateConfig)
    catalogue_write: CatalogueWriteConfig = field(factory=CatalogueWriteConfig)
    corpus: CorpusConfig = field(factory=CorpusConfig)
    corpus_describe: CorpusDescribeConfig = field(factory=CorpusDescribeConfig)
    corpus_golden: CorpusGoldenConfig = field(factory=CorpusGoldenConfig)
    corpus_save: CorpusSaveConfig = field(factory=CorpusSaveConfig)
    explain: ExplainConfig = field(factory=ExplainConfig)
    execution: ExecutionConfig = field(factory=ExecutionConfig)
    logging: LoggingConfig = field(factory=LoggingConfig)
    run: RunConfig = field(factory=RunConfig)
    source: SourceConfig = field(factory=SourceConfig)
    stats: StatsConfig = field(factory=StatsConfig)
    command: str = field(default="")
    sources: dict[str, dict[str, ConfigSource]] = field(factory=dict)

    section_types: ClassVar[dict[str, type[ConfigSection]]] = {
        "annotation": AnnotationConfig,
        "bench": BenchConfig,
        "bench_endtoend": BenchEndToEndConfig,
        "bench_matrix": BenchMatrixConfig,
        "bench_oracle": BenchOracleConfig,
        "bench_sweep": BenchSweepConfig,
        "catalogue": CatalogueConfig,
        "catalogue_add": CatalogueAddConfig,
        "catalogue_costs": CatalogueCostsConfig,
        "catalogue_remove": CatalogueRemoveConfig,
        "catalogue_show": CatalogueShowConfig,
        "catalogue_validate": CatalogueValidateConfig,
        "catalogue_write": CatalogueWriteConfig,
        "corpus": CorpusConfig,
        "corpus_describe": CorpusDescribeConfig,
        "corpus_golden": CorpusGoldenConfig,
        "corpus_save": CorpusSaveConfig,
        "explain": ExplainConfig,
        "execution": ExecutionConfig,
        "logging": LoggingConfig,
        "run": RunConfig,
        "source": SourceConfig,
        "stats": StatsConfig,
    }
    """Map section attribute names to their classes."""

    @classmethod
    def section_types_for(cls, names: tuple[str, ...]) -> tuple[type[ConfigSection], ...]:
        """Return the section classes for the named attributes."""
        section_types: list[type[ConfigSection]] = []
        for name in names:
            section_types.append(cls.section_types[name])
        return tuple(section_types)

    @classmethod
    def build_section(
        cls,
        section_type: type[ConfigSection],
        values: dict[str, object],
    ) -> ConfigSection:
        """Build one section from its merged values."""
        known: dict[str, object] = {}
        for name, value in values.items():
            if name in section_options(section_type):
                known[name] = value
        return section_type(**known)  # type: ignore[call-arg]

    @classmethod
    def sources_for(
        cls, layers: tuple[tuple[ConfigSource, ConfigLayer], ...]
    ) -> dict[str, dict[str, ConfigSource]]:
        """Return the winning source for each option."""
        sources: dict[str, dict[str, ConfigSource]] = {}
        for source, layer in layers:
            for section, values in layer.items():
                for field_name in values:
                    sources.setdefault(section, {})[field_name] = source
        return sources

    @classmethod
    def for_command(cls, command: str, **options: object) -> Config:
        """Return the configuration for one command, named as the command line names its options.

        The command's own sections decide where each option belongs, so a caller gives flat option
        names and never has to know which section declares which. Values are coerced as the command
        line coerces them, so a path may be given as a string.

        Args:
            command (str): A registered command, such as ``run.default``.
            **options (object): Option values, named by field or by flag.

        Returns:
            Config: The command's sections built from these options, every other section left at
                its defaults.

        Raises:
            ValueError: If the command is not registered, or an option is not one it reads.
        """

        section_names = cls._command_section_names(command)
        section_types = cls.section_types_for(section_names)

        given: ConfigLayer = {}
        for key, value in options.items():
            resolved = resolve_option(key, value, section_types)
            if resolved is None:
                raise ValueError(f"{key!r} is not an option of {command!r}")
            section_name, field_name, resolved_value = resolved
            given.setdefault(section_name, {})[field_name] = resolved_value

        coerced = coerce_layer(given, section_types)
        built: dict[str, ConfigSection] = {}
        for name in section_names:
            section_type = cls.section_types[name]
            built[name] = cls.build_section(section_type, coerced.get(section_type.section_name, {}))
        return cls(**built, command=command)  # type: ignore[arg-type]

    @staticmethod
    def _command_section_names(command: str) -> tuple[str, ...]:
        """Return the sections a registered command reads.

        The registry is imported here rather than at module scope because importing it registers
        the runners, and every runner module imports this one.

        Raises:
            ValueError: If no command is registered under the name.
        """

        # pylint: disable=import-outside-toplevel,import-error,no-name-in-module
        from ensembl.vep.runner import get_runner

        try:
            return get_runner(command).sections
        except KeyError as unknown:
            raise ValueError(str(unknown).strip("\"'")) from unknown

    def sections(self) -> tuple[ConfigSection, ...]:
        """Return the configured sections in declaration order."""
        sections: list[ConfigSection] = []
        for attribute in fields(type(self)):
            value = getattr(self, attribute.name)
            if isinstance(value, ConfigSection):
                sections.append(value)
        return tuple(sections)

    def section(self, name: str) -> ConfigSection:
        """Return one section by its declared name."""
        for section in self.sections():
            if section.section_name == name:
                return section
        raise KeyError(f"Unknown configuration section '{name}'")

    def source_of(self, section: str, field_name: str) -> ConfigSource:
        """Return the layer that supplied one option."""
        return self.sources.get(section, {}).get(field_name, ConfigSource.DEFAULT)

    def as_mapping(self) -> dict[str, dict[str, object]]:
        """Return the declared options as ``section -> field -> value``."""
        mapping: dict[str, dict[str, object]] = {}
        for section in self.sections():
            values: dict[str, object] = {}
            for name in section_options(type(section)):
                values[name] = getattr(section, name)
            mapping[section.section_name] = values
        return mapping

    def evolve(self, **changes: Any) -> Config:
        """Return a copy with selected sections replaced."""
        return attrs.evolve(self, **changes)


def load(
    argv: list[str] | None = None,
    *,
    command: str = "",
    section_names: tuple[str, ...] = (),
    parser: argparse.ArgumentParser | None = None,
    config_file: Path | None = None,
    overrides: ConfigLayer | None = None,
    environ: dict[str, str] | None = None,
) -> Config:
    """Load configuration from the available layers."""
    section_names = section_names or tuple(Config.section_types)
    section_types = Config.section_types_for(section_names)

    layers: list[tuple[ConfigSource, ConfigLayer]] = []
    if config_file is not None:
        layers.append((ConfigSource.FILE, read_config_file(config_file, section_types)))
    layers.append((ConfigSource.ENVIRONMENT, read_environment(section_types, environ)))
    if argv is not None:
        command_line_parser = parser or build_parser(f"vep {command}".strip(), "", section_types)
        layers.append((ConfigSource.COMMAND_LINE, parse_args(command_line_parser, argv, section_types)))
    if overrides:
        layers.append((ConfigSource.COMMAND_LINE, overrides))

    merged = merge_layers(*(layer for _, layer in layers))
    built: dict[str, ConfigSection] = {}
    for name in section_names:
        section_type = Config.section_types[name]
        built[name] = Config.build_section(section_type, merged.get(section_type.section_name, {}))
    sources = Config.sources_for(tuple(layers))
    return Config(**built, command=command, sources=sources)  # type: ignore[arg-type]
