(cli-run)=
# vep run

Commands to annotate variants.

Without a verb, `vep run` runs `vep run default`.

Without any command, `vep` runs `vep run default`.

(cli-run-default)=
## vep run default

Annotate variants against a set of annotation sources.

```text
usage: vep run default [-h] [--log_level {debug,info,warning,error,critical}]
                       [--log_file LOG_FILE] [-v] [-q] [-i INPUT_FILE] [--format {vcf}]
                       [--gff GFF] [--gtf GTF] [--fasta FASTA] [--synonyms SYNONYMS] [--cache]
                       [--dir_cache DIR_CACHE] [--database] [--custom SPEC]
                       [--distance BP[,DOWNSTREAM]] [--annotator NAME] [--batch_size N]
                       [--zone_map] [--zone_map_bytes BYTES] [--codon_cache]
                       [--codon_cache_entries N] [--mapper_cache] [--mapper_cache_bytes BYTES]
                       [--readthrough {none,utr,beyond_transcript}] [--early_subsumption_skip]
                       [--fact_cost_override] [--processes PROCESSES]
                       [--start_method {auto,spawn,fork,forkserver}] [--threads THREADS]
                       [--buffer_size BUFFER_SIZE] [--annotation_window ANNOTATION_WINDOW]
                       [--max_in_flight MAX_IN_FLIGHT] [--transcript_pool_bytes BYTES]
                       [-o OUTPUT_FILE] [--output_format OUTPUT_FORMAT] [--vcf] [--json]
                       [--force_overwrite] [--no_progress] [--no_stats] [--stats_file STATS_FILE]
                       [--stats_format {html,text,json}] [--stats_html] [--stats_text]
                       [--config FILE]
```

Common options: {ref}`Logging options <cli-section-logging>`, {ref}`Source options <cli-section-source>`, {ref}`Annotation options <cli-section-annotation>`.

Shared options: {ref}`Execution options <cli-section-execution>`, {ref}`Run options <cli-section-run>`.

Global options: {ref}`--config FILE <cli-global>`.

(cli-run-haplo)=
## vep run haplo

Annotate phased haplotypes.

```text
usage: vep run haplo [-h] [--log_level {debug,info,warning,error,critical}] [--log_file LOG_FILE]
                     [-v] [-q] [-i INPUT_FILE] [--format {vcf}] [--gff GFF] [--gtf GTF]
                     [--fasta FASTA] [--synonyms SYNONYMS] [--cache] [--dir_cache DIR_CACHE]
                     [--database] [--custom SPEC] [--distance BP[,DOWNSTREAM]] [--annotator NAME]
                     [--batch_size N] [--zone_map] [--zone_map_bytes BYTES] [--codon_cache]
                     [--codon_cache_entries N] [--mapper_cache] [--mapper_cache_bytes BYTES]
                     [--readthrough {none,utr,beyond_transcript}] [--early_subsumption_skip]
                     [--fact_cost_override] [--processes PROCESSES]
                     [--start_method {auto,spawn,fork,forkserver}] [--threads THREADS]
                     [--buffer_size BUFFER_SIZE] [--annotation_window ANNOTATION_WINDOW]
                     [--max_in_flight MAX_IN_FLIGHT] [--transcript_pool_bytes BYTES]
                     [-o OUTPUT_FILE] [--output_format OUTPUT_FORMAT] [--vcf] [--json]
                     [--force_overwrite] [--no_progress] [--no_stats] [--stats_file STATS_FILE]
                     [--stats_format {html,text,json}] [--stats_html] [--stats_text]
                     [--config FILE]
```

Common options: {ref}`Logging options <cli-section-logging>`, {ref}`Source options <cli-section-source>`, {ref}`Annotation options <cli-section-annotation>`.

Shared options: {ref}`Execution options <cli-section-execution>`, {ref}`Run options <cli-section-run>`.

Global options: {ref}`--config FILE <cli-global>`.

## Options shared by `vep run` commands

(cli-section-execution)=
### Execution options

Configure worker processes and optional independent-track threads.

Used by {ref}`vep run default <cli-run-default>`, {ref}`vep run haplo <cli-run-haplo>`.

In a TOML config file these options go under the `[execution]` table.

| Option | Default | Environment variable | Description |
| ------ | ------- | -------------------- | ----------- |
| `--processes PROCESSES, --fork PROCESSES` | `1` | `ENSEMBL_VEP_EXECUTION_PROCESSES` | Worker process count. |
| `--start_method {auto,spawn,fork,forkserver}` | `auto` | `ENSEMBL_VEP_EXECUTION_START_METHOD` | Multiprocessing start method. |
| `--threads THREADS` | `1` | `ENSEMBL_VEP_EXECUTION_THREADS` | Independent annotation-track thread count. |
| `--buffer_size BUFFER_SIZE` | `5000` | `ENSEMBL_VEP_EXECUTION_BUFFER_SIZE` | Records per input buffer. |
| `--annotation_window ANNOTATION_WINDOW` | `1024` | `ENSEMBL_VEP_EXECUTION_ANNOTATION_WINDOW` | Variations annotated together within an input buffer. |
| `--max_in_flight MAX_IN_FLIGHT` | – | `ENSEMBL_VEP_EXECUTION_MAX_IN_FLIGHT` | Maximum submitted input buffers. |
| `--transcript_pool_bytes BYTES` | `134217728` | `ENSEMBL_VEP_EXECUTION_TRANSCRIPT_POOL_BYTES` | Transcript pool byte budget (0 disables). |

(cli-section-run)=
### Run options

Where annotation output and the run summary go.

Used by {ref}`vep run default <cli-run-default>`, {ref}`vep run haplo <cli-run-haplo>`.

In a TOML config file these options go under the `[run]` table.

| Option | Default | Environment variable | Description |
| ------ | ------- | -------------------- | ----------- |
| `-o OUTPUT_FILE, --output_file OUTPUT_FILE` | `variant_effect_output.txt` | `ENSEMBL_VEP_RUN_OUTPUT_FILE` | Annotation output file. |
| `--output_format OUTPUT_FORMAT` | `vcf` | `ENSEMBL_VEP_RUN_OUTPUT_FORMAT` | Annotation output format (Perl VEP defaults to its own 'vep' format, which has no writer here yet). `--vcf` is shorthand for `vcf`. `--json` is shorthand for `json`. |
| `--force_overwrite, --force, --no_force_overwrite, --no_force` | `false` | `ENSEMBL_VEP_RUN_FORCE_OVERWRITE` | Overwrite an existing output file. |
| `--stats_file STATS_FILE, --sf STATS_FILE` | – | `ENSEMBL_VEP_RUN_STATS_FILE` | Run summary output file. |
| `--stats_format {html,text,json}` | `html` | `ENSEMBL_VEP_RUN_STATS_FORMAT` | Run summary format. `--stats_html` is shorthand for `html`. `--stats_text` is shorthand for `text`. |
