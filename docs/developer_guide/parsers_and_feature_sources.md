# Parsers and Feature Sources

This page covers the two packages responsible for turning raw input files (and, eventually, databases) into domain `Feature` objects: `ensembl.vep.parsers` and `ensembl.vep.sources`. Per-method detail is in the docstrings, surfaced via the auto-generated [code reference](../reference/ensembl.vep.rst) — this page only covers the shape of the architecture and the reasoning behind it.

## Parser interface

`ensembl.vep.parsers` defines a small `Parser[T]` contract that every format-specific parser implements:

```python
open(file_path, mode="auto") -> Parser   # opens a file; mode is "auto" | "query" | "sequential"
close() -> None
seek(region) -> bool                     # restricts subsequent next() calls to a region
next() -> T | None                       # one raw record at a time, or None at EOF
records(limit=None) -> list[T]           # drains the remaining records from the current position
```

Plus two static, one-shot helpers that don't require `open()`/`close()` at all:

```python
parse_text(data, limit=None) -> list[T]
parse_file(file_path, limit=None) -> list[T]
```

`GFF3Parser` is the worked example: a thin Python `attrs` class wrapping the compiled Rust `ensembl_vep_core` extension. The facade exists so callers depend on a stable, documented Python contract rather than the Rust binding directly — `open()` returns the parser itself and supports the context-manager protocol (`with parser.open(path) as p: ...`), so the underlying file handle is always closed even if an exception is raised mid-read.

## FeatureSource hierarchy

`ensembl.vep.sources` mirrors the parser layer one level up — instead of raw records, it produces domain `Feature` objects (`Transcript`, `Variation`, etc.):

```
FeatureSource (abstract root)
│   get_all_features()            -> list[Feature]
│   iter_all_features()           -> Iterator[Feature]
│   get_features_by_region(...)   -> list[Feature]
│   iter_features_by_region(...)  -> Iterator[Feature]
│   close()
│   __enter__ / __exit__ (delegates to close())
│
├── FileSource[P, F]   (file-backed; P = parser type, F = feature type)
│   └── concrete sources implement `_records_to_features(records) -> list[F]`
│       e.g. GFF3TranscriptSource[GFF3Parser, Transcript]
│
└── DBSource[D, F]     (DB-backed; D = DB adaptor type, F = feature type)
    └── concrete sources implement `_rows_to_features(rows) -> list[F]`
```

A concrete source only ever implements that one conversion hook — `FileSource`/`DBSource` already provide the four read methods, `close()`, and context-manager support on top of it.

## Streaming vs. eager

`FileSource`'s default `iter_all_features`/`iter_features_by_region` convert and yield each record as soon as it comes off `next()`, rather than buffering the whole file first:

```python
def iter_all_features(self) -> Iterator[F]:
    with self.adaptor.open(self.file_path, mode="sequential") as adaptor:
        while (rec := adaptor.next()) is not None:
            yield from self._records_to_features([rec])
```

This default is correct whenever one raw record maps to one (or more) complete features — e.g. one VCF line is one Variation. It is **not** correct for formats where a feature is assembled from several records, which is exactly the GFF3 case: a `Transcript` needs its own record plus all of its child `exon`/`CDS` records, and possibly a parent `gene` record.

`GFF3TranscriptSource` therefore overrides both iterator methods. It reads records one at a time via `next()`, buffers them, and flushes a complete group as soon as the next "gene" record (or a second top-level transcript record) signals the previous group is finished — emitting `Transcript` objects without ever holding the whole file in memory at once. This relies on the GFF3 file being sorted in the conventional way (a gene record precedes its transcripts, and each transcript's exon/CDS children immediately follow it); the source logs this assumption once per call. Files that don't follow this ordering will still produce *some* result via the eager path (`get_all_features`/`get_features_by_region`, which read the whole file/region up front and don't depend on ordering), just not via the streaming path.
