(cli-catalogue)=
# vep catalogue

Commands to maintain the consequence catalogue.

(cli-catalogue-add)=
## vep catalogue add

Add a term to the consequence catalogue.

```text
usage: vep catalogue add [-h] [--log_level {debug,info,warning,error,critical}]
                         [--log_file LOG_FILE] [-v] [-q] [--catalogue CATALOGUE]
                         [--catalogue_overwrite] [--dry_run] [--term TERM]
                         [--accession SO:XXXXXXX] [--rank N] [--tier TIER]
                         [--impact {HIGH,MODERATE,LOW,MODIFIER}] [--so_parent SO_PARENT]
                         [--subsumes SUBSUMES] [--no_derive_subsumes] [--when JSON]
                         [--allow_missing_when] [--notes NOTES] [--config FILE]
```

Common options: {ref}`Logging options <cli-section-logging>`.

Shared options: {ref}`Catalogue options <cli-section-catalogue>`, {ref}`Catalogue write options <cli-section-catalogue-write>`.

Global options: {ref}`--config FILE <cli-global>`.

(cli-section-catalogue-add)=
### Catalogue add options

The term being added.

In a TOML config file these options go under the `[catalogue_add]` table.

| Option | Default | Environment variable | Description |
| ------ | ------- | -------------------- | ----------- |
| `--term TERM` | – | `ENSEMBL_VEP_CATALOGUE_ADD_TERM` | The SO term to add. |
| `--accession SO:XXXXXXX` | – | `ENSEMBL_VEP_CATALOGUE_ADD_ACCESSION` | SO accession. |
| `--rank N` | – | `ENSEMBL_VEP_CATALOGUE_ADD_RANK` | Rank to insert at. |
| `--tier TIER` | `3` | `ENSEMBL_VEP_CATALOGUE_ADD_TIER` | Tier grouping. |
| `--impact {HIGH,MODERATE,LOW,MODIFIER}` | – | `ENSEMBL_VEP_CATALOGUE_ADD_IMPACT` | Impact bucket. |
| `--so_parent SO_PARENT` | – | `ENSEMBL_VEP_CATALOGUE_ADD_SO_PARENT` | SO parent accessions, used to derive subsumes. |
| `--subsumes SUBSUMES` | – | `ENSEMBL_VEP_CATALOGUE_ADD_SUBSUMES` | Extra terms to subsume. |
| `--when JSON` | – | `ENSEMBL_VEP_CATALOGUE_ADD_WHEN` | The when payload as JSON. |
| `--allow_missing_when, --no_allow_missing_when` | `false` | `ENSEMBL_VEP_CATALOGUE_ADD_ALLOW_MISSING_WHEN` | Accept an entry with no when, as a placeholder. |
| `--notes NOTES` | – | `ENSEMBL_VEP_CATALOGUE_ADD_NOTES` | Free-text notes. |

(cli-catalogue-costs)=
## vep catalogue costs

Measure each fact's real cost against its declared tier.

```text
usage: vep catalogue costs [-h] [--log_level {debug,info,warning,error,critical}]
                           [--log_file LOG_FILE] [-v] [-q] [--catalogue CATALOGUE]
                           [--output_format {text,markdown,confluence,html,svg,dot,mermaid}]
                           [--output_file OUTPUT_FILE]
                           [--kind {golden.file,loaded.files,synthetic.all,synthetic.default,synthetic.heavy}]
                           [--corpus_parameter KEY=VALUE] [--fuzz N] [--step N]
                           [--genome_seed GENOME_SEED] [--boundary_seed BOUNDARY_SEED]
                           [--fuzz_seed FUZZ_SEED] [--heavy_genome_seed HEAVY_GENOME_SEED]
                           [--sweep_seed SWEEP_SEED] [--corpus_dir CORPUS_DIR] [--config FILE]
```

Common options: {ref}`Logging options <cli-section-logging>`, {ref}`Corpus options <cli-section-corpus>`.

Shared options: {ref}`Catalogue options <cli-section-catalogue>`.

Global options: {ref}`--config FILE <cli-global>`.

(cli-section-catalogue-costs)=
### Catalogue costs options

How the measured costs are reported and persisted.

In a TOML config file these options go under the `[catalogue_costs]` table.

| Option | Default | Environment variable | Description |
| ------ | ------- | -------------------- | ----------- |
| `--output_format {text,markdown,confluence,html,svg,dot,mermaid}` | – | `ENSEMBL_VEP_CATALOGUE_COSTS_OUTPUT_FORMAT` | Report format (default: text). |
| `--output_file OUTPUT_FILE` | – | `ENSEMBL_VEP_CATALOGUE_COSTS_OUTPUT_FILE` | Write the measured cost table here. |

(cli-catalogue-remove)=
## vep catalogue remove

Remove a term from the consequence catalogue.

```text
usage: vep catalogue remove [-h] [--log_level {debug,info,warning,error,critical}]
                            [--log_file LOG_FILE] [-v] [-q] [--catalogue CATALOGUE]
                            [--catalogue_overwrite] [--dry_run] [--term TERM] [--config FILE]
```

Common options: {ref}`Logging options <cli-section-logging>`.

Shared options: {ref}`Catalogue options <cli-section-catalogue>`, {ref}`Catalogue write options <cli-section-catalogue-write>`.

Global options: {ref}`--config FILE <cli-global>`.

(cli-section-catalogue-remove)=
### Catalogue remove options

The term being removed.

In a TOML config file these options go under the `[catalogue_remove]` table.

| Option | Default | Environment variable | Description |
| ------ | ------- | -------------------- | ----------- |
| `--term TERM` | – | `ENSEMBL_VEP_CATALOGUE_REMOVE_TERM` | The SO term to remove. |

(cli-catalogue-show)=
## vep catalogue show

List catalogue terms or registered facts.

```text
usage: vep catalogue show [-h] [--log_level {debug,info,warning,error,critical}]
                          [--log_file LOG_FILE] [-v] [-q] [--catalogue CATALOGUE]
                          [--subject {terms,facts}] [--select SELECT] [--graph]
                          [--output_format {text,markdown,confluence,html,svg,dot,mermaid}]
                          [--sort {name,impact,cost_tier}] [--output_file OUTPUT_FILE]
                          [--config FILE]
```

Common options: {ref}`Logging options <cli-section-logging>`.

Shared options: {ref}`Catalogue options <cli-section-catalogue>`.

Global options: {ref}`--config FILE <cli-global>`.

(cli-section-catalogue-show)=
### Catalogue show options

What to list, and how to render it.

In a TOML config file these options go under the `[catalogue_show]` table.

| Option | Default | Environment variable | Description |
| ------ | ------- | -------------------- | ----------- |
| `--subject {terms,facts}` | `terms` | `ENSEMBL_VEP_CATALOGUE_SHOW_SUBJECT` | What to list. |
| `--select SELECT` | – | `ENSEMBL_VEP_CATALOGUE_SHOW_SELECT` | Show only this term or fact. |
| `--graph, --no_graph` | `false` | `ENSEMBL_VEP_CATALOGUE_SHOW_GRAPH` | Draw as a diagram. |
| `--output_format {text,markdown,confluence,html,svg,dot,mermaid}` | – | `ENSEMBL_VEP_CATALOGUE_SHOW_OUTPUT_FORMAT` | Report format (default: markdown, svg with --graph). |
| `--sort {name,impact,cost_tier}` | – | `ENSEMBL_VEP_CATALOGUE_SHOW_SORT` | Row order. |
| `--output_file OUTPUT_FILE` | – | `ENSEMBL_VEP_CATALOGUE_SHOW_OUTPUT_FILE` | Write output here. |

(cli-catalogue-validate)=
## vep catalogue validate

Validate the consequence rule set.

```text
usage: vep catalogue validate [-h] [--log_level {debug,info,warning,error,critical}]
                              [--log_file LOG_FILE] [-v] [-q] [--catalogue CATALOGUE]
                              [--solver {none,z3}] [--verify_requires]
                              [--kind {golden.file,loaded.files,synthetic.all,synthetic.default,synthetic.heavy}]
                              [--corpus_parameter KEY=VALUE] [--fuzz N] [--step N]
                              [--genome_seed GENOME_SEED] [--boundary_seed BOUNDARY_SEED]
                              [--fuzz_seed FUZZ_SEED] [--heavy_genome_seed HEAVY_GENOME_SEED]
                              [--sweep_seed SWEEP_SEED] [--corpus_dir CORPUS_DIR] [--config FILE]
```

Common options: {ref}`Logging options <cli-section-logging>`, {ref}`Corpus options <cli-section-corpus>`.

Shared options: {ref}`Catalogue options <cli-section-catalogue>`.

Global options: {ref}`--config FILE <cli-global>`.

(cli-section-catalogue-validate)=
### Catalogue validate options

Validation tiers to run.

In a TOML config file these options go under the `[catalogue_validate]` table.

| Option | Default | Environment variable | Description |
| ------ | ------- | -------------------- | ----------- |
| `--solver {none,z3}` | `none` | `ENSEMBL_VEP_CATALOGUE_VALIDATE_SOLVER` | Whole-ruleset SAT backend. |
| `--verify_requires, --no_verify_requires` | `false` | `ENSEMBL_VEP_CATALOGUE_VALIDATE_VERIFY_REQUIRES` | Check every declared requires against a corpus. |

## Options shared by `vep catalogue` commands

(cli-section-catalogue)=
### Catalogue options

Which catalogue to act on. Every verb reads this one.

Used by {ref}`vep catalogue add <cli-catalogue-add>`, {ref}`vep catalogue costs <cli-catalogue-costs>`, {ref}`vep catalogue remove <cli-catalogue-remove>`, {ref}`vep catalogue show <cli-catalogue-show>`, {ref}`vep catalogue validate <cli-catalogue-validate>`.

In a TOML config file these options go under the `[catalogue]` table.

| Option | Default | Environment variable | Description |
| ------ | ------- | -------------------- | ----------- |
| `--catalogue CATALOGUE` | – | `ENSEMBL_VEP_CATALOGUE_CATALOGUE_FILE` | Catalogue JSON to read or write. |

(cli-section-catalogue-write)=
### Catalogue write options

Read only by the verbs that change the catalogue, so the listing verbs never offer it.

Used by {ref}`vep catalogue add <cli-catalogue-add>`, {ref}`vep catalogue remove <cli-catalogue-remove>`.

In a TOML config file these options go under the `[catalogue_write]` table.

| Option | Default | Environment variable | Description |
| ------ | ------- | -------------------- | ----------- |
| `--catalogue_overwrite, --no_catalogue_overwrite` | `false` | `ENSEMBL_VEP_CATALOGUE_WRITE_CATALOGUE_OVERWRITE` | Write the packaged catalogue in place. |
| `--dry_run, --no_dry_run` | `false` | `ENSEMBL_VEP_CATALOGUE_WRITE_DRY_RUN` | Report the change without writing it. |
