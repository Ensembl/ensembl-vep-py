# Position in a Sequence Region

## Coordinate Systems and Slices

We use coordinate systems and slices to represent biological sequences, following the conventions established by Ensembl.

A **coordinate system** defines the reference frame for a biological sequence — for example, a chromosome, clone, or contig. Ensembl uses **1-based indexing**, meaning the first base of a sequence is position 1. Coordinate systems themselves carry no concept of strandedness.

A **slice** represents a defined region within a coordinate system. Slice coordinates are relative to the start of the underlying sequence region, and by convention the start of a slice is always less than or equal to its end, regardless of strandedness. Strandedness is introduced at the slice level — the strand of a slice represents its orientation relative to the default orientation of the sequence region (which should always be forward).

Features positioned on a slice have coordinates and strand values expressed relative to that slice. To obtain a feature's absolute position on the underlying sequence region, the `seq_region_start`, `seq_region_end`, and `seq_region_strand` methods should be used instead.

> **Note:** The concept of a "sequence region" is referenced in various parts of the codebase but has no dedicated class. A slice spanning the full length of a biological sequence with positive strandedness can be considered a functional equivalent of a sequence region.

## Coordinate Conventions

All positions on a slice or sequence feature use `start` and `end` coordinates, where the condition `start <= end` always holds true.

**Exception:** 
1. Insertion variants described in Ensembl's default format can have `end = start - 1`.
2. A circular slice on a bacterial or mitochondrial DNA can have `end < start`.

We use **1-based inclusive (closed)** coordinates for all positions. This means that the first base of a sequence is at position 1, and both `start` and `end` positions are included in the feature span.

## Insertion coordinates
As convention in Ensembl, for insertion variants, the `start` and `end` coordinates are defined such that `end = start - 1`. This is because in 1-based inclusive coordinates, length is represented as `end - start + 1`. An insertion has a length of 0, furthermore the adjacent start and end coordinates clearly indicate the position of the insertion between two bases.

This convention is also convenient for calculating pre and post sequence for shifting algorithm. Nonetheless, it can be a source of confusion and errors, so it is important to be aware of this special case when working with insertion variants.

## Removal of Reference Bases

When `Variation` objects are created any anchor bases are removed and the coordinates adjusted accordingly. For example, a deletion represented in VCF as `1   10  CAT C` becomes `AT/-` with `start=11` and `end=12`. An insertion such as `1   10  G   GAT` becomes `-/AT` with `start=11` and `end=10`.
The same is true for symbolic alleles. e.g. This variant `1   100 G <DEL>`, becomes `-/<DEL>` and `start` is updated to 101.

