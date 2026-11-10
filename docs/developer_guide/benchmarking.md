# Benchmarking

`vep bench` measures annotation speed. Its two main commands answer different questions: `vep bench matrix` asks which annotator is fastest, and `vep bench endtoend` asks how much of that difference reaches a real run. They often disagree, and the gap between them is usually the finding.

Both run over a [test corpus](../concepts/corpus.md), selected with `--kind`. The [`vep bench` reference](../auto/reference/cli/bench.md) lists every option.

## The conformance gate

Before anything is timed, every configuration is run once and its output compared with the reference strategy's, variant by variant, as sets of consequence terms. A configuration that gives different terms is not timed: being fast because it is wrong is not a result. `--no_gate` skips this check when the code is already known to conform.

## `vep bench matrix`: the annotator alone

`vep bench matrix` times the annotators over overlaps that are already built, in the current process. No VCF is read, no transcript is loaded from a file and no output is written:

```bash
vep bench matrix --kind synthetic.heavy --strategies consequence.linear,consequence.columnar-staged --caches none,zonemap
```

It crosses each strategy in `--strategies` with each cache set in `--caches`, building fresh caches for every combination, because the caches are independent of the strategy. A cache set is `none`, or `zonemap`, `codon` and `mapper` joined with `+`, such as `zonemap+codon`. Alongside wall-clock time it reports the weighted number of fact computations per overlap, which depends only on the strategy and not on the language or machine it runs on.

Use it to compare strategies, tune a cache, or check that a rule change did not make annotation slower.

## `vep bench endtoend`: the whole command line

`vep bench endtoend` writes the corpus to disk as FASTA, GFF3 and VCF files, indexes them, and runs the real command line once per configuration: reading the VCF, loading transcripts through the indexed GFF3, fetching sequence, annotating and writing the output.

```bash
vep bench endtoend --kind synthetic.heavy --strategies consequence.linear,consequence.columnar --caches none,zonemap
```

`--processes` adds worker counts to the configurations. Use it to decide on a default, to confirm that an annotator's in-process gain reaches a real run, or when the change is not in the annotator at all: transcript loading, sequence access, batching or worker counts.

What it measures:

- **Every run is a separate process,** so no configuration benefits from another's warmed state, and start-up, imports and worker start are counted, as a real run pays them.
- **The run uses this working copy.** The child process gets the directory of the imported `ensembl.vep` first on its `PYTHONPATH`, so an older installed copy cannot be measured by mistake.
- **Peak memory is that of the largest single process,** not the sum across workers, so a run with several processes understates its total memory.
- **On sparse input, annotation is a small part of a run.** Most of the time can go to loading transcripts, so a strategy that wins several times over in `vep bench matrix` can move the end-to-end time by a few percent. That is the measurement working, not failing.

### Files

The corpus, its indexes and each run's output are written to `--work_dir`, to `--corpus_dir` when that is set, or otherwise to a temporary directory. `--clean`, the default, removes exactly the files the command created and leaves anything already in the directory, so pointing `--corpus_dir` at an existing corpus reuses it. `--no_clean` keeps everything, which helps when a result needs investigating. Each configuration writes its own output file, so a large corpus across many configurations needs disk space.

### Requirements

The corpus is annotated through a `bgzip` compressed, `tabix` indexed GFF3 file, because an unindexed file takes a different path through the reader and would answer a different question. `bgzip` and `tabix` must be on the `PATH`; the command stops with an error naming them otherwise. The FASTA index is written directly, so `samtools` is not needed.

## Other commands

- `vep bench sweep` runs one numeric annotation setting across a list of values, checking and timing each value.
- `vep bench oracle` compares the consequence strategies with the consequence code of Perl VEP.

## Recording and comparing results

`--repeats` and `--warmup` set how many timed and untimed passes each configuration gets, and `--output_format` the format of the results table. `--results_file` writes the results as JSON; a later run given that file as `--baseline_file` reports every configuration that is slower by more than `--regression_threshold`. A baseline measured on a different corpus is refused, because the comparison would mean nothing.
