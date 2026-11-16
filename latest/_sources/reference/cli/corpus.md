(cli-corpus)=
# vep corpus

Commands to build and check the synthetic corpus.

(cli-corpus-describe)=
## vep corpus describe

Describe a registered corpus and list its selectable overlaps.

```text
usage: vep corpus describe [-h] [--log_level {debug,info,warning,error,critical}]
                           [--log_file LOG_FILE] [-v] [-q]
                           [--kind {golden.file,loaded.files,synthetic.all,synthetic.default,synthetic.heavy}]
                           [--corpus_parameter KEY=VALUE] [--fuzz N] [--step N]
                           [--genome_seed GENOME_SEED] [--boundary_seed BOUNDARY_SEED]
                           [--fuzz_seed FUZZ_SEED] [--heavy_genome_seed HEAVY_GENOME_SEED]
                           [--sweep_seed SWEEP_SEED] [--corpus_dir CORPUS_DIR] [--annotator NAME]
                           [--batch_size N] [--zone_map] [--zone_map_bytes BYTES] [--codon_cache]
                           [--codon_cache_entries N] [--mapper_cache] [--mapper_cache_bytes BYTES]
                           [--readthrough {none,utr,beyond_transcript}] [--early_subsumption_skip]
                           [--fact_cost_override] [--select SELECT]
                           [--output_format {text,markdown,confluence,html,svg,dot,mermaid}]
                           [--output_file OUTPUT_FILE] [--annotate] [--config FILE]
```

Common options: {ref}`Logging options <cli-section-logging>`, {ref}`Corpus options <cli-section-corpus>`, {ref}`Annotation options <cli-section-annotation>`.

Global options: {ref}`--config FILE <cli-global>`.

(cli-section-corpus-describe)=
### Corpus describe options

What to show about a corpus, and how.

In a TOML config file these options go under the `[corpus_describe]` table.

| Option | Default | Environment variable | Description |
| ------ | ------- | -------------------- | ----------- |
| `--select SELECT` | – | `ENSEMBL_VEP_CORPUS_DESCRIBE_SELECT` | List only this variant's overlaps. |
| `--output_format {text,markdown,confluence,html,svg,dot,mermaid}` | – | `ENSEMBL_VEP_CORPUS_DESCRIBE_OUTPUT_FORMAT` | Report format (default: markdown). |
| `--output_file OUTPUT_FILE` | – | `ENSEMBL_VEP_CORPUS_DESCRIBE_OUTPUT_FILE` | Write output here. |
| `--annotate, --no_annotate` | `false` | `ENSEMBL_VEP_CORPUS_DESCRIBE_ANNOTATE` | Also annotate the corpus and count its terms. |

(cli-corpus-golden)=
## vep corpus golden

Record or check the terms an annotator emits for a corpus.

```text
usage: vep corpus golden [-h] [--log_level {debug,info,warning,error,critical}]
                         [--log_file LOG_FILE] [-v] [-q]
                         [--kind {golden.file,loaded.files,synthetic.all,synthetic.default,synthetic.heavy}]
                         [--corpus_parameter KEY=VALUE] [--fuzz N] [--step N]
                         [--genome_seed GENOME_SEED] [--boundary_seed BOUNDARY_SEED]
                         [--fuzz_seed FUZZ_SEED] [--heavy_genome_seed HEAVY_GENOME_SEED]
                         [--sweep_seed SWEEP_SEED] [--corpus_dir CORPUS_DIR] [--annotator NAME]
                         [--batch_size N] [--zone_map] [--zone_map_bytes BYTES] [--codon_cache]
                         [--codon_cache_entries N] [--mapper_cache] [--mapper_cache_bytes BYTES]
                         [--readthrough {none,utr,beyond_transcript}] [--early_subsumption_skip]
                         [--fact_cost_override] [--golden_file GOLDEN_FILE] [--update_golden]
                         [--config FILE]
```

Common options: {ref}`Logging options <cli-section-logging>`, {ref}`Corpus options <cli-section-corpus>`, {ref}`Annotation options <cli-section-annotation>`.

Global options: {ref}`--config FILE <cli-global>`.

(cli-section-corpus-golden)=
### Corpus golden options

Where golden terms live, and whether to record them or check against them.

In a TOML config file these options go under the `[corpus_golden]` table.

| Option | Default | Environment variable | Description |
| ------ | ------- | -------------------- | ----------- |
| `--golden_file GOLDEN_FILE` | – | `ENSEMBL_VEP_CORPUS_GOLDEN_GOLDEN_FILE` | Golden terms file to check or record. |
| `--update_golden, --no_update_golden` | `false` | `ENSEMBL_VEP_CORPUS_GOLDEN_UPDATE_GOLDEN` | Record the terms and overwrite the golden file. |

(cli-corpus-save)=
## vep corpus save

Write a registered corpus as FASTA, GFF3 and VCF.

```text
usage: vep corpus save [-h] [--log_level {debug,info,warning,error,critical}]
                       [--log_file LOG_FILE] [-v] [-q]
                       [--kind {golden.file,loaded.files,synthetic.all,synthetic.default,synthetic.heavy}]
                       [--corpus_parameter KEY=VALUE] [--fuzz N] [--step N]
                       [--genome_seed GENOME_SEED] [--boundary_seed BOUNDARY_SEED]
                       [--fuzz_seed FUZZ_SEED] [--heavy_genome_seed HEAVY_GENOME_SEED]
                       [--sweep_seed SWEEP_SEED] [--corpus_dir CORPUS_DIR] [--exporter {file}]
                       [--force] [--config FILE]
```

Common options: {ref}`Logging options <cli-section-logging>`, {ref}`Corpus options <cli-section-corpus>`.

Global options: {ref}`--config FILE <cli-global>`.

(cli-section-corpus-save)=
### Corpus save options

How to export a corpus. The destination is ``corpus.corpus_dir``.

In a TOML config file these options go under the `[corpus_save]` table.

| Option | Default | Environment variable | Description |
| ------ | ------- | -------------------- | ----------- |
| `--exporter {file}` | `file` | `ENSEMBL_VEP_CORPUS_SAVE_EXPORTER` | Registered corpus exporter. |
| `--force, --no_force` | `false` | `ENSEMBL_VEP_CORPUS_SAVE_FORCE` | Overwrite an existing saved corpus. |
