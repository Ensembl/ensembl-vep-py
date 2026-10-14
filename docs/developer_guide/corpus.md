# Corpus architecture

The corpus component supplies deterministic inputs for consequence tests, benchmarks, explanation,
and profiling. A consumer asks for a `Corpus`; it does not need to know whether the data was
generated, loaded from files, or obtained from another backend.

## Component boundary

`Corpus` is the consumer-facing contract. It exposes:

- the `CorpusSpec` used to construct it;
- every referenced sequence region;
- named variation alleles;
- transcripts;
- named transcript overlaps; and
- resource lifetime through `close()` and the context-manager protocol.

The contract is intentionally transcript-oriented because the current consequence annotators
consume `TranscriptVariationAllele`. A database-backed transcript corpus fits this contract. A
future regulatory or motif corpus needs an explicit extension of the biological feature and
overlap types; changing only the storage backend is not sufficient.

`CorpusSpec` is an open construction value. It contains a registered corpus name and scalar
factory parameters. The registry resolves the name and leaves interpretation and validation of
parameters to the selected factory. Built-in synthetic corpora use the canonical registered names
`synthetic.default`, `synthetic.heavy`, and `synthetic.all`.

```text
CorpusSpec ──► corpus registry ──► Corpus
                                  ├── SyntheticCorpus
                                  ├── LoadedCorpus
                                  ├── GoldenCorpus
                                  └── future backend corpus

Corpus ──► describer
       ├──► annotator / benchmark / tests
       └──► exporter registry ──► file bundle or future format
```

## Package structure

```text
corpus/
  base.py                 Corpus contract and named value aliases
  spec.py                 registered name plus factory-owned parameters
  registry.py             corpus factory registration and construction
  combine.py              combine compatible corpora
  validation.py           cross-collection identity and naming invariants
  overlap.py              construct transcript overlaps
  describe.py             descriptions and their builders

  synthetic/
    base.py               materialised SyntheticCorpus
    default.py            default synthetic recipe and registration
    heavy.py              dense shared-locus recipe and registration
    components/
      sequence.py         reusable sequence generation
      transcript.py       transcript construction and CDS placement
      variation.py        reusable variation generation

  loaded/
    base.py               LoadedCorpus and owned-resource lifetime
    files.py              FASTA + GFF3/GTF + VCF loader registration

  golden/
    base.py               GoldenCorpus paired with pinned expected terms
    resolve.py            load packaged or configured golden data
    terms.py              serialised expected annotation terms
    data/<dataset>/       packaged corpus files, manifest, and terms

  exporter/
    base.py               CorpusExporter protocol
    registry.py           exporter registration and lookup
    manifest.py           exported corpus metadata
    formats/file/         coordinated FASTA, GFF3, and VCF exporter
```

Tests mirror this structure under `src/python/tests/corpus/`, so ownership remains visible as the
component grows.

## Corpus implementations

### Synthetic corpus

A synthetic module defines one recipe. It builds sequence regions, transcripts, and variants,
then calls the shared overlap builder. Recipe modules own their parameters and register one
factory name, for example `synthetic.default`.

Reusable biological construction belongs in `synthetic/components/`; recipe-specific coordinates
and editing remain in the recipe module. Adding another synthetic corpus normally requires one
new recipe module, its registration, and its matching test module. Existing consumers require no
changes.

### Loaded corpus

`LoadedCorpus` holds materialised corpus values and any resources that must remain alive, such as
sequence sources or database connections. It closes owned resources once.

`loaded.files` is one factory, not the definition of a loaded corpus. It coordinates FASTA,
GFF3/GTF, and VCF sources, reconciles their sequence-region objects, and builds overlaps. A future
database loader should construct the same model objects, retain its connection as a corpus
resource when necessary, and register another name such as `loaded.database`.

### Golden corpus

`GoldenCorpus` decorates another corpus with pinned expected consequence terms. `golden.file`
loads an existing packaged or configured dataset through `loaded.files`. Missing or malformed data
is an error; the corpus layer never silently replaces a reference.

Recording terms and generating a missing dataset with `consequence.linear` are annotator concerns.
That integration should depend on this corpus package and call its public loading and terms APIs.
The corpus package must not import the annotator component. This dependency direction keeps file
loading and benchmark corpora usable before or without a particular annotation implementation.

## Exporters

Export is independent of corpus construction. Consumers resolve an exporter with
`get_exporter(name)` and call `export(corpus, destination, force=...)`. The `file` exporter writes
a coordinated FASTA, GFF3, VCF, and manifest bundle.

A database dump or another representation belongs in a new directory under
`exporter/formats/<format>/` and registers its exporter name. Corpus implementations do not gain
format-specific methods.

## Invariants and lifetime

- Sequence-region names are unique within a corpus.
- Named variants are unique and overlap names refer to those variants.
- Overlap feature and allele objects are the same instances held by the corpus collections.
- Corpus constructors accept materialised tuples and expose tuples through the common contract.
  This fixes membership and order without freezing the biological model objects themselves.
- A combined corpus preserves stable input order and rejects conflicting names.
- A corpus that owns external resources is used as a context manager or explicitly closed.
- Exporters render from `Corpus`; they do not depend on a synthetic or loaded subtype.
- Existing data that cannot be loaded is an error; fallback policy belongs to the consumer.

## Extension checklist

For a new synthetic recipe:

1. Add `synthetic/<name>.py` and a matching test module.
2. Reuse sequence, transcript, variation, and overlap component builders.
3. Register `synthetic.<name>` with a factory that validates its parameters.
4. Assert deterministic identity and representative biological boundaries.

For a new loading backend:

1. Convert backend records to the existing sequence, transcript, and variation models.
2. Build or load named overlaps in stable order.
3. Construct `LoadedCorpus`, passing owned connections or sources as resources.
4. Register a backend-specific corpus name and test multiple sequence regions and cleanup.

For a new export format:

1. Implement `CorpusExporter` under `exporter/formats/<format>/`.
2. Register the exporter name.
3. Add matching format-conformance or round-trip tests.
