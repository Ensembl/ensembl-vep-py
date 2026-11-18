# Improving performance

This guide lists the settings that change how fast a run is and how much memory it uses. None of them change the annotations a run writes. Measure your own data before and after a change: the best values depend on the number of variants, how densely they cover the genome, and your machine.

## Index the input files

Index the transcript file and the reference sequence, as described in [Choosing annotation sources](annotation_sources.md#indexing-the-files). Without an index, every batch of variants rescans the files, which is much slower than reading the indexed regions.

## Use several processes

`--processes` (or `--fork`) runs the annotation in several worker processes:

```bash
vep run --config vep.toml --processes 4
```

The input is split into buffers of `--buffer_size` records, and each buffer is annotated by one worker. The output is still written in input order. Each worker holds its own copy of the transcripts and caches it needs, so memory use grows with the number of processes. A run with fewer records than `--buffer_size` fills only one buffer and cannot use a second process; lower `--buffer_size` to share a small input between workers.

`--start_method` selects how worker processes are started (`spawn`, `fork` or `forkserver`). The default, `auto`, uses `fork` where the platform supports it, because a forked worker starts faster than a spawned one, and the platform's default start method when `--threads` is more than one.

## Tune the buffers

- `--buffer_size` sets the largest number of input records one worker annotates at a time. A buffer never holds records from two sequence regions.
- `--annotation_window` sets how many variants within a buffer are annotated together. Smaller windows hold less in memory at once.
- `--max_in_flight` limits how many buffers are handed to the workers at once, which bounds memory when the input is read faster than it is annotated. By default it matches the number of processes.
- `--transcript_pool_bytes` sets the memory budget for transcripts each worker keeps between buffers, so neighbouring buffers do not rebuild the same transcripts. `0` disables the pool.

## Optional caches

The annotation options include caches that help some inputs and not others:

- `--zone_map` helps inputs with many variants per transcript, such as a densely genotyped region.
- `--codon_cache` reuses codon calculations shared between variants.
- `--mapper_cache` reuses coordinate mappings between transcript and genome.

Each cache has a size or memory budget option. See the {ref}`annotation options <cli-section-annotation>` for all of them and their defaults.

## Threads

`--threads` lets a worker run independent annotation tracks in parallel. A run currently has only one annotation track, transcript consequences, so more threads do not make it faster.
