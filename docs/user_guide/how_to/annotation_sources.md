# Choosing annotation sources

Ensembl VEP annotates variants against a transcript file and a reference sequence. This guide shows which files a run needs, how to prepare them, and how to match sequence-region names between them.

## Transcripts and reference sequence

A run needs both of the following:

- a transcript file, given with either `--gff` (GFF3) or `--gtf` (GTF), but not both;
- a reference sequence, given with `--fasta`.

```bash
vep run \
  --input_file variants.vcf \
  --gtf genes.gtf.bgz \
  --fasta genome.fa.bgz
```

Use the transcript file and the reference sequence of the same assembly as your variants.

## Indexing the files

All three files can be plain text or compressed with `bgzip`. They can be read without an index, but a run is much faster when they are indexed, because Ensembl VEP can then read only the regions around your variants instead of scanning each file:

- **GFF3 and GTF:** sort the file, compress it with `bgzip` and index it with `tabix`. A transcript file read without an index must keep all the records of a gene together.
- **FASTA:** index the file with `samtools faidx`. A plain FASTA file needs a `.fai` index; a `bgzip` compressed FASTA file needs both a `.fai` and a `.gzi` index.

The example data in `examples/homo_sapiens` is compressed and indexed this way.

## Matching sequence-region names

Different files can name the same chromosome differently, for example `22`, `chr22` or `NC_000022.11`. Ensembl VEP matches the names in your variants and transcripts to the names in the reference sequence in this order:

1. the exact name;
2. a synonym from the synonyms file, when one is given;
3. the same name with the `chr` prefix added or removed, ignoring case.

When the names in your files differ in some other way, give a synonyms file with `--synonyms`. Each line holds two whitespace-separated names for the same region: the name used by the reference sequence first, then its synonym.

```text
22	NC_000022.11
X	NC_000023.11
```

The output always reports a variant on the sequence region named in your input file.

## Upstream and downstream distance

A variant near a transcript but outside it is reported as `upstream_gene_variant` or `downstream_gene_variant` when it lies within the distance set by `--distance`. Give one value for both directions, or two comma-separated values for upstream and downstream:

```bash
vep run --config vep.toml --distance 5000,2000
```

See the {ref}`source options <cli-section-source>` for the default.

## Sources that are not available yet

The `--cache`, `--dir_cache` and `--database` options are accepted, but annotating from a VEP cache or an Ensembl database is not implemented yet. A run given either option stops with an error saying so.
