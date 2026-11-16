(cli-explain)=
# vep explain

Commands to explain and summarise annotation.

Explain one overlap's derivation, fact by fact and rule by rule.

```text
usage: vep explain [-h] [--log_level {debug,info,warning,error,critical}] [--log_file LOG_FILE]
                   [-v] [-q] [-i INPUT_FILE] [--format {vcf}] [--gff GFF] [--gtf GTF]
                   [--fasta FASTA] [--synonyms SYNONYMS] [--cache] [--dir_cache DIR_CACHE]
                   [--database] [--custom SPEC] [--distance BP[,DOWNSTREAM]]
                   [--kind {golden.file,loaded.files,synthetic.all,synthetic.default,synthetic.heavy}]
                   [--corpus_parameter KEY=VALUE] [--fuzz N] [--step N]
                   [--genome_seed GENOME_SEED] [--boundary_seed BOUNDARY_SEED]
                   [--fuzz_seed FUZZ_SEED] [--heavy_genome_seed HEAVY_GENOME_SEED]
                   [--sweep_seed SWEEP_SEED] [--corpus_dir CORPUS_DIR] [--annotator NAME]
                   [--batch_size N] [--zone_map] [--zone_map_bytes BYTES] [--codon_cache]
                   [--codon_cache_entries N] [--mapper_cache] [--mapper_cache_bytes BYTES]
                   [--readthrough {none,utr,beyond_transcript}] [--early_subsumption_skip]
                   [--fact_cost_override] [--name NAME] [--index N]
                   [--variant CHR:POS:REF>ALT|CHR:POS:ID] [--allele_num N]
                   [--transcript STABLE_ID] [--compare NAME,NAME] [--list] [--graph]
                   [--output_format {text,svg,dot,mermaid}] [--output_file OUTPUT_FILE]
                   [--config FILE]
```

Common options: {ref}`Logging options <cli-section-logging>`, {ref}`Source options <cli-section-source>`, {ref}`Corpus options <cli-section-corpus>`, {ref}`Annotation options <cli-section-annotation>`.

Global options: {ref}`--config FILE <cli-global>`.

(cli-section-explain)=
## Explain options

Which overlap to explain, and how to render the derivation.

In a TOML config file these options go under the `[explain]` table.

| Option | Default | Environment variable | Description |
| ------ | ------- | -------------------- | ----------- |
| `--name NAME` | – | `ENSEMBL_VEP_EXPLAIN_NAME` | Named overlap to explain. |
| `--index N` | – | `ENSEMBL_VEP_EXPLAIN_INDEX` | Overlap index to explain. |
| `--variant CHR:POS:REF>ALT\|CHR:POS:ID` | – | `ENSEMBL_VEP_EXPLAIN_VARIANT` | Variant to explain, by alleles or by the input's own ID. |
| `--allele_num N` | – | `ENSEMBL_VEP_EXPLAIN_ALLELE_NUM` | Alternate allele of a CHR:POS:ID variant to explain (default: 1). |
| `--transcript STABLE_ID` | – | `ENSEMBL_VEP_EXPLAIN_TRANSCRIPT` | Transcript stable id. |
| `--compare NAME,NAME` | – | `ENSEMBL_VEP_EXPLAIN_COMPARE` | Two annotators to diff. |
| `--list, --no_list` | `false` | `ENSEMBL_VEP_EXPLAIN_LIST_OVERLAPS` | List selectable overlaps and exit (see also: vep corpus describe). |
| `--graph, --no_graph` | `false` | `ENSEMBL_VEP_EXPLAIN_GRAPH` | Draw the derivation as a diagram. |
| `--output_format {text,svg,dot,mermaid}` | – | `ENSEMBL_VEP_EXPLAIN_OUTPUT_FORMAT` | Report format (default: text, svg with --graph). |
| `--output_file OUTPUT_FILE` | – | `ENSEMBL_VEP_EXPLAIN_OUTPUT_FILE` | Write the explanation here. |
