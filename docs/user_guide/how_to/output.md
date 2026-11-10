# Writing output and run summaries

This guide shows how to choose where annotations are written, in which format, and how to create a summary of a run.

## Output file

`--output_file` (or `-o`) names the annotation output. The default is `variant_effect_output.txt` in the current directory.

Ensembl VEP does not replace an existing file. Add `--force_overwrite` (or `--force`) to replace it:

```bash
vep run --config vep.toml --output_file annotated.vcf --force_overwrite
```

The output and the run summary are written to temporary files first and moved into place together when the run succeeds, so a failed run leaves no partial output behind.

To write the annotations to standard output instead, for example to pipe them into another tool, give `STDOUT` or `-` as the output file:

```bash
vep run --config vep.toml --output_file STDOUT | grep -v '^##' | head
```

## Output format

`--output_format` selects the format:

- `vcf` (the default) writes the input VCF records with their annotations added in the `CSQ` INFO field. `--vcf` is a shorthand for `--output_format vcf`.
- `json` writes one JSON object per variant, one object per line. `--json` is a shorthand for `--output_format json`.

See [Data formats](../../reference/data_formats.md) for the content of each format.

## Run summary

`vep run` writes a summary of the run beside the output, named after the output file with a `_summary` suffix: `annotated.vcf` gives `annotated.vcf_summary.html`. The summary counts the variants, their classes and their consequences.

- `--stats_file` names the summary file instead.
- `--stats_format` selects `html` (the default), `text` or `json`. `--stats_html` and `--stats_text` are shorthands.
- `--no_stats` skips the summary.

No summary is derived when the output goes to standard output; give `--stats_file` to write one.

## Summarising an existing output file

`vep stats` builds the same summary from a file Ensembl VEP has already annotated, without running the annotation again:

```bash
vep stats --input_file annotated.vcf --output_format html --output_file annotated_summary.html
```

`vep stats` reads the consequence annotations from the `CSQ` (or `ANN`) INFO field described in the VCF header, so it does not need the transcript or reference files. See the [`vep stats` reference](../../auto/reference/cli/stats.md) for its options.
