(cli-common)=
# Common options

Option sections read by commands from more than one group. Each command's page links to
the ones it accepts.

(cli-section-logging)=
## Logging options

How much the run says, and where it says it.

Used by {ref}`vep bench endtoend <cli-bench-endtoend>`, {ref}`vep bench matrix <cli-bench-matrix>`, {ref}`vep bench oracle <cli-bench-oracle>`, {ref}`vep bench sweep <cli-bench-sweep>`, {ref}`vep catalogue add <cli-catalogue-add>`, {ref}`vep catalogue costs <cli-catalogue-costs>`, {ref}`vep catalogue remove <cli-catalogue-remove>`, {ref}`vep catalogue show <cli-catalogue-show>`, {ref}`vep catalogue validate <cli-catalogue-validate>`, {ref}`vep corpus describe <cli-corpus-describe>`, {ref}`vep corpus golden <cli-corpus-golden>`, {ref}`vep corpus save <cli-corpus-save>`, {ref}`vep explain <cli-explain>`, {ref}`vep run default <cli-run-default>`, {ref}`vep run haplo <cli-run-haplo>`, {ref}`vep stats <cli-stats>`.

In a TOML config file these options go under the `[logging]` table.

| Option | Default | Environment variable | Description |
| ------ | ------- | -------------------- | ----------- |
| `--log_level {debug,info,warning,error,critical}, --level {debug,info,warning,error,critical}` | `info` | `ENSEMBL_VEP_LOGGING_LOG_LEVEL` | Logging level. |
| `--log_file LOG_FILE, --warning_file LOG_FILE` | – | `ENSEMBL_VEP_LOGGING_LOG_FILE` | Write log records to this file. |
| `-v, --verbose, --no_verbose` | `false` | `ENSEMBL_VEP_LOGGING_VERBOSE` | Log at debug level. |
| `-q, --quiet, --no_quiet` | `false` | `ENSEMBL_VEP_LOGGING_QUIET` | Log errors only. |

(cli-section-corpus)=
## Corpus options

Which registered corpus to build, and how.

Used by {ref}`vep bench endtoend <cli-bench-endtoend>`, {ref}`vep bench matrix <cli-bench-matrix>`, {ref}`vep bench sweep <cli-bench-sweep>`, {ref}`vep catalogue costs <cli-catalogue-costs>`, {ref}`vep catalogue validate <cli-catalogue-validate>`, {ref}`vep corpus describe <cli-corpus-describe>`, {ref}`vep corpus golden <cli-corpus-golden>`, {ref}`vep corpus save <cli-corpus-save>`, {ref}`vep explain <cli-explain>`, {ref}`vep stats <cli-stats>`.

In a TOML config file these options go under the `[corpus]` table.

| Option | Default | Environment variable | Description |
| ------ | ------- | -------------------- | ----------- |
| `--kind {golden.file,loaded.files,synthetic.all,synthetic.default,synthetic.heavy}` | `synthetic.all` | `ENSEMBL_VEP_CORPUS_KIND` | Registered corpus to use. |
| `--corpus_parameter KEY=VALUE` | – | `ENSEMBL_VEP_CORPUS_PARAMETERS` | Factory parameter as KEY=VALUE (repeatable). |
| `--fuzz N` | `4000` | `ENSEMBL_VEP_CORPUS_FUZZ` | Fuzz variant count. |
| `--step N` | `13` | `ENSEMBL_VEP_CORPUS_STEP` | Heavy-locus variant sweep step. |
| `--genome_seed GENOME_SEED` | `7` | `ENSEMBL_VEP_CORPUS_GENOME_SEED` | Seed for the synthetic genome. |
| `--boundary_seed BOUNDARY_SEED` | `11` | `ENSEMBL_VEP_CORPUS_BOUNDARY_SEED` | Seed for the boundary variants. |
| `--fuzz_seed FUZZ_SEED` | `23` | `ENSEMBL_VEP_CORPUS_FUZZ_SEED` | Seed for the fuzz variants. |
| `--heavy_genome_seed HEAVY_GENOME_SEED` | `31` | `ENSEMBL_VEP_CORPUS_HEAVY_GENOME_SEED` | Seed for the heavy-locus genome. |
| `--sweep_seed SWEEP_SEED` | `41` | `ENSEMBL_VEP_CORPUS_SWEEP_SEED` | Seed for the heavy-locus variant sweep. |
| `--corpus_dir CORPUS_DIR` | – | `ENSEMBL_VEP_CORPUS_CORPUS_DIR` | Directory to load or save a corpus. |

(cli-section-annotation)=
## Annotation options

Which annotator to run and which caches to give it.

Used by {ref}`vep corpus describe <cli-corpus-describe>`, {ref}`vep corpus golden <cli-corpus-golden>`, {ref}`vep explain <cli-explain>`, {ref}`vep run default <cli-run-default>`, {ref}`vep run haplo <cli-run-haplo>`, {ref}`vep stats <cli-stats>`.

In a TOML config file these options go under the `[annotation]` table.

| Option | Default | Environment variable | Description |
| ------ | ------- | -------------------- | ----------- |
| `--annotator NAME` | `consequence.columnar-staged` | `ENSEMBL_VEP_ANNOTATION_ANNOTATOR` | Annotator to run. |
| `--batch_size N` | `256` | `ENSEMBL_VEP_ANNOTATION_BATCH_SIZE` | Overlaps per batch. |
| `--zone_map, --no_zone_map` | `false` | `ENSEMBL_VEP_ANNOTATION_ZONE_MAP` | Enable the positional zone map for inputs with many variants per transcript. |
| `--zone_map_bytes BYTES` | `268435456` | `ENSEMBL_VEP_ANNOTATION_ZONE_MAP_BYTES` | Zone map byte budget. |
| `--codon_cache, --no_codon_cache` | `false` | `ENSEMBL_VEP_ANNOTATION_CODON_CACHE` | Enable the shared-codon cache. |
| `--codon_cache_entries N` | `65536` | `ENSEMBL_VEP_ANNOTATION_CODON_CACHE_ENTRIES` | Codon cache size. |
| `--mapper_cache, --no_mapper_cache` | `false` | `ENSEMBL_VEP_ANNOTATION_MAPPER_CACHE` | Enable the mapper index cache. |
| `--mapper_cache_bytes BYTES` | `67108864` | `ENSEMBL_VEP_ANNOTATION_MAPPER_CACHE_BYTES` | Mapper cache byte budget. |
| `--readthrough {none,utr,beyond_transcript}` | `none` | `ENSEMBL_VEP_ANNOTATION_READTHROUGH` | How far to translate past a destroyed stop codon. |
| `--early_subsumption_skip, --no_early_subsumption_skip` | `false` | `ENSEMBL_VEP_ANNOTATION_EARLY_SUBSUMPTION_SKIP` | Skip rules already suppressed by a hit (staged only). |
| `--fact_cost_override, --no_fact_cost_override` | `false` | `ENSEMBL_VEP_FACT_COST_OVERRIDE` | Use measured fact costs. |

(cli-section-source)=
## Source options

Input variant and annotation source configuration.

Used by {ref}`vep explain <cli-explain>`, {ref}`vep run default <cli-run-default>`, {ref}`vep run haplo <cli-run-haplo>`, {ref}`vep stats <cli-stats>`.

In a TOML config file these options go under the `[source]` table.

| Option | Default | Environment variable | Description |
| ------ | ------- | -------------------- | ----------- |
| `-i INPUT_FILE, --input_file INPUT_FILE` | – | `ENSEMBL_VEP_SOURCE_INPUT_FILE` | Variant input file. |
| `--format {vcf}, --input_format {vcf}` | – | `ENSEMBL_VEP_SOURCE_FORMAT` | Input format, inferred from the file extension when omitted. |
| `--gff GFF` | – | `ENSEMBL_VEP_SOURCE_GFF` | GFF3 file of transcripts. |
| `--gtf GTF` | – | `ENSEMBL_VEP_SOURCE_GTF` | GTF file of transcripts. |
| `--fasta FASTA` | – | `ENSEMBL_VEP_SOURCE_FASTA` | Reference FASTA. |
| `--synonyms SYNONYMS` | – | `ENSEMBL_VEP_SOURCE_SYNONYMS` | File of sequence-region synonyms, one region per line. |
| `--cache, --no_cache` | `false` | `ENSEMBL_VEP_SOURCE_CACHE` | Annotate from a VEP cache (not implemented yet). |
| `--dir_cache DIR_CACHE` | – | `ENSEMBL_VEP_SOURCE_DIR_CACHE` | Directory of the VEP cache. |
| `--database, --no_database` | `false` | `ENSEMBL_VEP_SOURCE_DATABASE` | Annotate from an Ensembl database (not implemented yet). |
| `--custom SPEC` | – | `ENSEMBL_VEP_SOURCE_CUSTOM` | Custom annotation spec, repeatable (not implemented yet). |
| `--distance BP[,DOWNSTREAM], --flank BP[,DOWNSTREAM]` | `5000` | `ENSEMBL_VEP_SOURCE_DISTANCE` | Distance outside a feature still counted as overlapping, as one value or up,down. |
