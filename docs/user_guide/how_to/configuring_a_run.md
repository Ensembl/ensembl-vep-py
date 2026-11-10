# Configuring a run

Every option of a `vep` command can be set in three places: on the command line, in an environment variable, or in a configuration file. This guide shows each one and how they combine.

## Command-line flags

Flags are the most direct way to set an option. Most options take a value:

```bash
vep run --input_file variants.vcf --gtf genes.gtf.bgz --fasta genome.fa.bgz
```

Switches take no value. Each switch also has a `--no_` form that turns it off again, which is useful to override a value set in a configuration file or an environment variable:

```bash
vep run --config vep.toml --no_stats
```

Some options accept a second spelling. For example, `--format` can also be given as `--input_format`, and `--distance` as `--flank`. The [CLI reference](../../auto/reference/cli/index.md) lists every flag of every command with its other spellings and its default.

## Environment variables

Every option can also be set by an environment variable named `ENSEMBL_VEP_<SECTION>_<OPTION>`, in upper case. The section is the option table the option is listed under in the CLI reference. For example, `--output_format` belongs to the `run` options:

```bash
export ENSEMBL_VEP_RUN_OUTPUT_FORMAT=json
vep run --input_file variants.vcf --gtf genes.gtf.bgz --fasta genome.fa.bgz
```

Switches accept `1`, `true`, `yes` or `on`, and `0`, `false`, `no` or `off`. The CLI reference shows the variable name of each option.

## Configuration files

`--config FILE` reads options from a file. The file can be written in [TOML](https://toml.io/), with one table per option section:

```toml
[source]
input_file = "examples/homo_sapiens/homo_sapiens_GRCh38.vcf"
gtf = "examples/homo_sapiens/genes.gtf.bgz"
fasta = "examples/homo_sapiens/unmasked.fa.bgz"

[run]
output_format = "json"
no_stats = true
```

```bash
vep run --config vep.toml --output_file annotated.json
```

`--config` can appear anywhere on the command line, and one file can hold the sections of several commands: a command reads the tables it uses and ignores the others.

A file without any `[table]` header is read as a Perl VEP configuration file instead, with one option and its value per line:

```text
input_file examples/homo_sapiens/homo_sapiens_GRCh38.vcf
gtf examples/homo_sapiens/genes.gtf.bgz
fasta examples/homo_sapiens/unmasked.fa.bgz
```

Options this implementation does not recognise are skipped in a Perl VEP configuration file. Inside a TOML table, an unrecognised option is skipped with a warning.

## Which value wins

When an option is set in more than one place, the command line wins over the environment, and the environment wins over the configuration file. An option that is not set anywhere takes its default.

For example, with the `vep.toml` above, this run writes VCF rather than JSON, because the environment variable overrides the file:

```bash
ENSEMBL_VEP_RUN_OUTPUT_FORMAT=vcf vep run --config vep.toml --output_file annotated.vcf
```
