---
myst:
  html_meta:
    "description lang=en": "Guidance for developers contributing to the Ensembl VEP project."
---

# Developer Guide

Conventions and architectural notes for contributing to `ensembl-vep-py`. Per-method/per-class detail lives in docstrings and the auto-generated [API reference](../auto/reference/api/ensembl.vep.rst) — these pages cover the conventions and reasoning that aren't obvious from reading a single file in isolation. The design of the main components is explained in [Concepts](../concepts/index.md).

```{toctree}
:maxdepth: 3

setup
style_conventions
writing_classes
circular_dependencies
benchmarking
documentation
coverage_report
writing_plugins
adding_annotation_sources
```
