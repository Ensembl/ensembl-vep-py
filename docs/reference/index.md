---
myst:
  html_meta:
    "description lang=en": |
      Reference documentation for the ensembl-vep command line, data formats and Python API.
---

# Reference

Look up the exact commands, options, formats and terms Ensembl VEP supports. The CLI reference, the consequence terms and the API reference are generated from the source code by `make clidoc`, `make termdoc` and `make apidoc`; `make docs` regenerates all three before building the site.

```{toctree}
:maxdepth: 2

CLI reference <../auto/reference/cli/index>
data_formats
Consequence terms <../auto/reference/consequence_terms>
API reference <../auto/reference/api/ensembl.vep>
```
