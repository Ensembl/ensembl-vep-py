# Developer Guide

Conventions and architectural notes for contributing to `ensembl-vep-py`. Per-method/per-class detail lives in docstrings and the auto-generated [code reference](../reference/index.md) — these pages cover the conventions and reasoning that aren't obvious from reading a single file in isolation.

- [Style Conventions](style_conventions.md) — typing, logging, and docstring conventions.
- [Writing Classes](writing_classes.md) — class naming, attribute categories (core/derived/externally-set), `attrs` configuration, caching and cache invalidation, mixins, and containers.
- [Position in a Sequence Region](coordinates_and_slices.md) — coordinate systems, slices, and coordinate conventions (including the insertion-coordinate special case).
- [Getting away with circular dependency](circular_dependencies.md) — import patterns that avoid circular imports across the package structure.
- [Parsers and Feature Sources](parsers_and_feature_sources.md) — the `Parser` interface and the `FeatureSource`/`FileSource`/`DBSource` hierarchy used to turn raw input files (and, eventually, databases) into domain `Feature` objects.
