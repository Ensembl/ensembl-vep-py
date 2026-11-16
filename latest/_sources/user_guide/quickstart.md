# Quickstart guide

This guide annotates a small set of human variants against the GRCh38 reference assembly and Ensembl gene annotation. The example data is included in the repository, so the run does not need to download any additional files.

Follow the [source installation instructions](install.md#source), then run the commands below from the root of the `ensembl-vep-py` repository.

## Annotate the variants

Run Ensembl VEP with the example VCF, GTF and FASTA files:

```bash
vep run \
  --input_file examples/homo_sapiens/homo_sapiens_GRCh38.vcf \
  --gtf examples/homo_sapiens/genes.gtf.bgz \
  --fasta examples/homo_sapiens/unmasked.fa.bgz \
  --output_file quickstart_output.vcf \
  --stats_file quickstart_summary.html
```

The command annotates the example variants and creates two files in the current directory:

- `quickstart_output.vcf` contains the input variants and their consequence annotations in the `CSQ` INFO field;
- `quickstart_summary.html` contains a summary of the run and its results.

Open `quickstart_summary.html` in a web browser, or inspect the first annotated records on the command line:

```bash
grep -v '^##' quickstart_output.vcf | head
```

Ensembl VEP does not replace an existing output file by default. Add `--force_overwrite` to the command when you want to repeat the example and replace both output files.

See the [`vep run` reference](../reference/cli/run.md) for the full set of input, output and execution options.
