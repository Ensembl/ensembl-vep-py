(cli-stats)=
# vep stats

Commands to explain and summarise annotation.

Summarise an annotated file, or profile the facts over a corpus.

```text
usage: vep stats [-h] [--log_level {debug,info,warning,error,critical}] [--log_file LOG_FILE] [-v]
                 [-q] [-i INPUT_FILE] [--format {vcf}] [--gff GFF] [--gtf GTF] [--fasta FASTA]
                 [--synonyms SYNONYMS] [--cache] [--dir_cache DIR_CACHE] [--database]
                 [--custom SPEC] [--distance BP[,DOWNSTREAM]]
                 [--kind {golden.file,loaded.files,synthetic.all,synthetic.default,synthetic.heavy}]
                 [--corpus_parameter KEY=VALUE] [--fuzz N] [--step N] [--genome_seed GENOME_SEED]
                 [--boundary_seed BOUNDARY_SEED] [--fuzz_seed FUZZ_SEED]
                 [--heavy_genome_seed HEAVY_GENOME_SEED] [--sweep_seed SWEEP_SEED]
                 [--corpus_dir CORPUS_DIR] [--annotator NAME] [--batch_size N] [--zone_map]
                 [--zone_map_bytes BYTES] [--codon_cache] [--codon_cache_entries N]
                 [--mapper_cache] [--mapper_cache_bytes BYTES]
                 [--readthrough {none,utr,beyond_transcript}] [--early_subsumption_skip]
                 [--fact_cost_override] [--output_file OUTPUT_FILE]
                 [--output_format {html,text,json}] [--profile] [--config FILE]
```

Common options: {ref}`Logging options <cli-section-logging>`, {ref}`Source options <cli-section-source>`, {ref}`Corpus options <cli-section-corpus>`, {ref}`Annotation options <cli-section-annotation>`.

Global options: {ref}`--config FILE <cli-global>`.

(cli-section-stats)=
## Stats options

What ``vep stats`` reads and what it writes.

In a TOML config file these options go under the `[stats]` table.

| Option | Default | Environment variable | Description |
| ------ | ------- | -------------------- | ----------- |
| `--output_file OUTPUT_FILE` | – | `ENSEMBL_VEP_STATS_OUTPUT_FILE` | Report output file. |
| `--output_format {html,text,json}` | `html` | `ENSEMBL_VEP_STATS_OUTPUT_FORMAT` | Report format. |
| `--profile, --no_profile` | `false` | `ENSEMBL_VEP_STATS_PROFILE` | Collect fact-level profiling (slower). |
