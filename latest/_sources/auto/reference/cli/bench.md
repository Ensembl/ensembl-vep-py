(cli-bench)=
# vep bench

Commands to measure annotation speed.

(cli-bench-endtoend)=
## vep bench endtoend

Time whole command-line runs over a corpus on disk.

```text
usage: vep bench endtoend [-h] [--log_level {debug,info,warning,error,critical}]
                          [--log_file LOG_FILE] [-v] [-q]
                          [--kind {golden.file,loaded.files,synthetic.all,synthetic.default,synthetic.heavy}]
                          [--corpus_parameter KEY=VALUE] [--fuzz N] [--step N]
                          [--genome_seed GENOME_SEED] [--boundary_seed BOUNDARY_SEED]
                          [--fuzz_seed FUZZ_SEED] [--heavy_genome_seed HEAVY_GENOME_SEED]
                          [--sweep_seed SWEEP_SEED] [--corpus_dir CORPUS_DIR] [--repeats N]
                          [--warmup N] [--no_gate] [--results_file RESULTS_FILE]
                          [--baseline_file BASELINE_FILE]
                          [--regression_threshold REGRESSION_THRESHOLD]
                          [--output_format {text,markdown,confluence,html,svg,dot,mermaid}]
                          [--strategies NAME,...] [--caches SET,...] [--processes N,...]
                          [--batch_size N] [--work_dir WORK_DIR] [--clean] [--config FILE]
```

Common options: {ref}`Logging options <cli-section-logging>`, {ref}`Corpus options <cli-section-corpus>`.

Shared options: {ref}`Bench options <cli-section-bench>`.

Global options: {ref}`--config FILE <cli-global>`.

(cli-section-bench-endtoend)=
### Bench endtoend options

Which configurations to run end to end, where to run them, and what to leave behind.

In a TOML config file these options go under the `[bench_endtoend]` table.

| Option | Default | Environment variable | Description |
| ------ | ------- | -------------------- | ----------- |
| `--strategies NAME,...` | – | `ENSEMBL_VEP_BENCH_ENDTOEND_STRATEGIES` | Strategies to time (default: every consequence strategy). |
| `--caches SET,...` | `none,zonemap,zonemap+codon` | `ENSEMBL_VEP_BENCH_ENDTOEND_CACHES` | Cache sets to cross each strategy with. |
| `--processes N,...` | `1` | `ENSEMBL_VEP_BENCH_ENDTOEND_PROCESSES` | Worker process counts to time each configuration with. |
| `--batch_size N` | `256` | `ENSEMBL_VEP_BENCH_ENDTOEND_BATCH_SIZE` | Overlaps per batch. |
| `--work_dir WORK_DIR` | – | `ENSEMBL_VEP_BENCH_ENDTOEND_WORK_DIR` | Directory for the corpus, its indexes and the run outputs. |
| `--clean, --no_clean` | `true` | `ENSEMBL_VEP_BENCH_ENDTOEND_CLEAN` | Remove the files this run created. |

(cli-bench-matrix)=
## vep bench matrix

Time every strategy and cache configuration over a corpus.

```text
usage: vep bench matrix [-h] [--log_level {debug,info,warning,error,critical}]
                        [--log_file LOG_FILE] [-v] [-q]
                        [--kind {golden.file,loaded.files,synthetic.all,synthetic.default,synthetic.heavy}]
                        [--corpus_parameter KEY=VALUE] [--fuzz N] [--step N]
                        [--genome_seed GENOME_SEED] [--boundary_seed BOUNDARY_SEED]
                        [--fuzz_seed FUZZ_SEED] [--heavy_genome_seed HEAVY_GENOME_SEED]
                        [--sweep_seed SWEEP_SEED] [--corpus_dir CORPUS_DIR] [--repeats N]
                        [--warmup N] [--no_gate] [--results_file RESULTS_FILE]
                        [--baseline_file BASELINE_FILE]
                        [--regression_threshold REGRESSION_THRESHOLD]
                        [--output_format {text,markdown,confluence,html,svg,dot,mermaid}]
                        [--strategies NAME,...] [--caches SET,...] [--batch_size N]
                        [--config FILE]
```

Common options: {ref}`Logging options <cli-section-logging>`, {ref}`Corpus options <cli-section-corpus>`.

Shared options: {ref}`Bench options <cli-section-bench>`.

Global options: {ref}`--config FILE <cli-global>`.

(cli-section-bench-matrix)=
### Bench matrix options

Which strategies, crossed with which cache sets.

In a TOML config file these options go under the `[bench_matrix]` table.

| Option | Default | Environment variable | Description |
| ------ | ------- | -------------------- | ----------- |
| `--strategies NAME,...` | – | `ENSEMBL_VEP_BENCH_MATRIX_STRATEGIES` | Strategies to time (default: every consequence strategy). |
| `--caches SET,...` | `none,zonemap,zonemap+codon` | `ENSEMBL_VEP_BENCH_MATRIX_CACHES` | Cache sets to cross each strategy with. |
| `--batch_size N` | `256` | `ENSEMBL_VEP_BENCH_MATRIX_BATCH_SIZE` | Chunk size for the batch timing. |

(cli-bench-oracle)=
## vep bench oracle

Compare our strategies with the Perl VariationEffect.pm oracle.

```text
usage: vep bench oracle [-h] [--log_level {debug,info,warning,error,critical}]
                        [--log_file LOG_FILE] [-v] [-q] [--repeats N] [--warmup N] [--no_gate]
                        [--results_file RESULTS_FILE] [--baseline_file BASELINE_FILE]
                        [--regression_threshold REGRESSION_THRESHOLD]
                        [--output_format {text,markdown,confluence,html,svg,dot,mermaid}]
                        [--perl_modules PERL_MODULES] [--perl_oracle_script PERL_ORACLE_SCRIPT]
                        [--strategies NAME,...] [--oracle_number N] [--config FILE]
```

Common options: {ref}`Logging options <cli-section-logging>`.

Shared options: {ref}`Bench options <cli-section-bench>`.

Global options: {ref}`--config FILE <cli-global>`.

(cli-section-bench-oracle)=
### Bench oracle options

Where the Perl oracle lives, and what to compare against it.

In a TOML config file these options go under the `[bench_oracle]` table.

| Option | Default | Environment variable | Description |
| ------ | ------- | -------------------- | ----------- |
| `--perl_modules PERL_MODULES` | – | `ENSEMBL_VEP_BENCH_ORACLE_PERL_MODULES` | Ensembl-variation modules directory. |
| `--perl_oracle_script PERL_ORACLE_SCRIPT` | `bench/variation_effect.pl` | `ENSEMBL_VEP_BENCH_ORACLE_PERL_ORACLE_SCRIPT` | The Perl oracle script. |
| `--strategies NAME,...` | – | `ENSEMBL_VEP_BENCH_ORACLE_STRATEGIES` | Strategies to compare (default: linear and the kernels). |
| `--oracle_number N` | `200` | `ENSEMBL_VEP_BENCH_ORACLE_ORACLE_NUMBER` | Workload executions per timed repeat. |

(cli-bench-sweep)=
## vep bench sweep

Sweep one numeric annotation setting and time each value.

```text
usage: vep bench sweep [-h] [--log_level {debug,info,warning,error,critical}]
                       [--log_file LOG_FILE] [-v] [-q]
                       [--kind {golden.file,loaded.files,synthetic.all,synthetic.default,synthetic.heavy}]
                       [--corpus_parameter KEY=VALUE] [--fuzz N] [--step N]
                       [--genome_seed GENOME_SEED] [--boundary_seed BOUNDARY_SEED]
                       [--fuzz_seed FUZZ_SEED] [--heavy_genome_seed HEAVY_GENOME_SEED]
                       [--sweep_seed SWEEP_SEED] [--corpus_dir CORPUS_DIR] [--repeats N]
                       [--warmup N] [--no_gate] [--results_file RESULTS_FILE]
                       [--baseline_file BASELINE_FILE]
                       [--regression_threshold REGRESSION_THRESHOLD]
                       [--output_format {text,markdown,confluence,html,svg,dot,mermaid}]
                       [--strategies NAME,...]
                       [--sweep_knob {batch_size,codon_cache_entries,mapper_cache_bytes,zone_map_bytes}]
                       [--sweep_values N,N,...] [--config FILE]
```

Common options: {ref}`Logging options <cli-section-logging>`, {ref}`Corpus options <cli-section-corpus>`.

Shared options: {ref}`Bench options <cli-section-bench>`.

Global options: {ref}`--config FILE <cli-global>`.

(cli-section-bench-sweep)=
### Bench sweep options

Which setting to sweep, over which values, for which strategies.

In a TOML config file these options go under the `[bench_sweep]` table.

| Option | Default | Environment variable | Description |
| ------ | ------- | -------------------- | ----------- |
| `--strategies NAME,...` | – | `ENSEMBL_VEP_BENCH_SWEEP_STRATEGIES` | Strategies to sweep (default: the batch strategies). |
| `--sweep_knob {batch_size,codon_cache_entries,mapper_cache_bytes,zone_map_bytes}` | `batch_size` | `ENSEMBL_VEP_BENCH_SWEEP_SWEEP_KNOB` | Setting to vary. |
| `--sweep_values N,N,...` | `32,64,128,256,512,1024` | `ENSEMBL_VEP_BENCH_SWEEP_SWEEP_VALUES` | Values to try. |

## Options shared by `vep bench` commands

(cli-section-bench)=
### Bench options

How every benchmark is measured and reported. Read by every ``bench`` verb.

Used by {ref}`vep bench endtoend <cli-bench-endtoend>`, {ref}`vep bench matrix <cli-bench-matrix>`, {ref}`vep bench oracle <cli-bench-oracle>`, {ref}`vep bench sweep <cli-bench-sweep>`.

In a TOML config file these options go under the `[bench]` table.

| Option | Default | Environment variable | Description |
| ------ | ------- | -------------------- | ----------- |
| `--repeats N` | `3` | `ENSEMBL_VEP_BENCH_REPEATS` | Timed samples per configuration. |
| `--warmup N` | `1` | `ENSEMBL_VEP_BENCH_WARMUP` | Untimed passes before timing. |
| `--results_file RESULTS_FILE` | – | `ENSEMBL_VEP_BENCH_RESULTS_FILE` | Write results JSON here. |
| `--baseline_file BASELINE_FILE` | – | `ENSEMBL_VEP_BENCH_BASELINE_FILE` | Compare against this results JSON. |
| `--regression_threshold REGRESSION_THRESHOLD` | `0.05` | `ENSEMBL_VEP_BENCH_REGRESSION_THRESHOLD` | Fractional slowdown counted as a regression. |
| `--output_format {text,markdown,confluence,html,svg,dot,mermaid}` | – | `ENSEMBL_VEP_BENCH_OUTPUT_FORMAT` | Report format (default: text). |
