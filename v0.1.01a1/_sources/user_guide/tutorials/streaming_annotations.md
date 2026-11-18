# Streaming annotations in Python

In this tutorial you use Ensembl VEP as a Python library. Instead of writing an output file and reading it back, you receive each annotated variant as soon as it is ready and process it in your own code. You will:

1. configure a run in Python;
2. iterate over the annotated variants;
3. collect your own statistics while the variants stream past;
4. convert the variants to Ensembl VEP JSON records;
5. write the stream to a file.

The tutorial uses the example data included in the repository. Follow the [source installation instructions](../install.md#source), then run the code from the root of the `ensembl-vep-py` repository.

## Configure the run

`Config.for_command` builds the configuration of a command from keyword arguments. The keywords are the option names you would give on the command line, without the leading dashes:

```python
from ensembl.vep.config import Config

config = Config.for_command(
    "run.default",
    input_file="examples/homo_sapiens/homo_sapiens_GRCh38.vcf",
    gtf="examples/homo_sapiens/genes.gtf.bgz",
    fasta="examples/homo_sapiens/unmasked.fa.bgz",
)
```

`"run.default"` is the `vep run default` command, which `vep run` runs. Every option of that command is accepted, so `processes=4` or `distance=2000` work here as they do on the command line. An option the command does not have raises a `ValueError` naming it.

## Iterate over the annotated variants

`Runner.annotate()` starts the annotation and returns the annotated variants as a stream. Use the stream in a `with` block, so the annotation stops and its resources are released even if you leave the loop early:

```python
from ensembl.vep.runner.run import Runner

annotations = Runner(config).annotate()
with annotations:
    for variant in annotations:
        print(variant.id, variant.seq_region_name, variant.start, variant.allele_string, variant.most_severe_consequence)
        break
```

```text
rs7289170 22 17181903 A/G synonymous_variant
```

Variants arrive in the order of the input file, whatever the number of worker processes. Each returned object is an `AnnotatedVariant`, a snapshot of the variant and its annotations. Its attributes, with the values of the first example variant:

| Attribute | Example | Description |
|---|---|---|
| `id` | `rs7289170` | The variant identifier from the input, or `None` when it has none. |
| `seq_region_name` | `22` | The sequence region name, spelled as in the input. |
| `start` | `17181903` | The genomic start position, inclusive. |
| `end` | `17181903` | The genomic end position, inclusive. |
| `strand` | `1` | The genomic strand of the variant. |
| `allele_string` | `A/G` | The reference and alternate alleles, separated by `/`. |
| `input_columns` | `('22', '17181903', 'rs7289170', ...)` | The original columns of the input line. |
| `annotations` | 25 dictionaries | One dictionary per allele and feature, with the same fields as the `CSQ` field of the [VCF output](../../reference/data_formats.md#vcf-output). |
| `most_severe_consequence` | `synonymous_variant` | The most severe consequence term over all the annotations. |
| `variant_class` | `SNV` | The Sequence Ontology class of the variant. |
| `assembly_name` | `None` | The assembly name, when it is known. |
| `skipped` | `False` | Whether the variant could not be annotated. |

A stream can be read only once. Call `Runner(config).annotate()` again to read the variants a second time.

## Collect your own statistics

Because each variant arrives as soon as it is annotated, you can aggregate results without keeping them all in memory. This counts the transcripts with a missense consequence for each gene:

```python
from collections import Counter

missense = Counter()
annotations = Runner(config).annotate()
with annotations:
    for variant in annotations:
        for record in variant.annotations:
            if "missense_variant" in record["Consequence"]:
                missense[record["Gene"]] += 1

print(missense.most_common(3))
```

```text
[('ENSG00000183963', 392), ('ENSG00000099889', 154), ('ENSG00000128242', 120)]
```

Once the stream has been read to the end, it also carries the run summary that `vep run` writes beside its output:

```python
summary = annotations.summary
print(summary.variants)
print(summary.most_severe[:3])
```

```text
136
(('missense_variant', 75), ('synonymous_variant', 59), ('splice_polypyrimidine_tract_variant', 1))
```

The summary is only available after the last variant: reading it earlier raises a `RuntimeError`.

## Convert the variants to JSON records

`JSONStreamFactory` turns each annotated variant into a dictionary in the Ensembl VEP JSON format, which is easier to pass to other tools than the `CSQ` field layout:

```python
from ensembl.vep.output.factory.stream import JSONStreamFactory

annotations = Runner(config).annotate()
with annotations:
    for record in JSONStreamFactory().records(annotations):
        consequences = record.get("transcript_consequences", [])
        print(record["id"], record["most_severe_consequence"], len(consequences))
        break
```

```text
rs7289170 synonymous_variant 25
```

The factory converts the variants one at a time as you read them, so the stream stays a stream.

## Write the stream to a file

The writers `vep run` uses are available too. `get_output_writer` returns the writer for an output format; it turns each variant into one line of text:

```python
from ensembl.vep.output.factory.write import get_output_writer

writer = get_output_writer("json")
annotations = Runner(config).annotate()
with annotations, open("annotated.json", "w", encoding="utf-8") as output:
    for line in writer.headers():
        output.write(line + "\n")
    for variant in annotations:
        output.write(writer.line(variant) + "\n")
```

`annotated.json` now holds one JSON object per line, the same content `vep run --output_format json` writes.

The VCF writer also passes the header lines of the input VCF through to the output, so give them to `writer.headers()`. The runner's input source reads them:

```python
runner = Runner(config)
writer = get_output_writer("vcf")
annotations = runner.annotate()
with annotations, open("annotated.vcf", "w", encoding="utf-8") as output:
    for line in writer.headers(runner.input_source.header_lines()):
        output.write(line + "\n")
    for variant in annotations:
        output.write(writer.line(variant) + "\n")
```

## Next steps

- [Configuring a run](../how_to/configuring_a_run.md) lists the other ways to set options.
- [Improving performance](../how_to/performance.md) explains the options that make a run faster.
