// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

//! Batched fact kernels whose FFI surface is **flat integer arrays and nothing else**.
//!
//! Every function here accepts `Vec<i64>` (or a pair of them describing intervals) and returns
//! a packed bitmap. No struct crosses the boundary and no Python object is borrowed, nothing here
//! knows what a transcript is or looks like.
//!
//! The operations are the smallest units the fact layer is built from, they are already fully
//! specified by their signature, and there is no plausible future in which "do these two closed
//!  intervals overlap" grows a new field.
//!
//! **Why batched, and why bitmaps.** A per-call boundary crossing costs far more than any of
//! these operations, so a scalar `genomic_overlap(a, b, c, d)` FFI call would be pure loss. The
//! unit is therefore one column: given N contexts, compute the fact for all N and return the
//! answer as a little-endian bitmap the caller turns into one Python int with
//! `int.from_bytes(buf, "little")`.
//!
//! Bit `i` of the returned bitmap is the fact's value for input index `i`. Every function
//! returns `ceil(n / 8)` bytes.

use pyo3::prelude::*;

type VariantShapeBitmaps = (Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>);

/// Pack a per-index boolean predicate into a little-endian bitmap.
///
/// The one place the bit layout is written down, so a kernel cannot disagree with its siblings
/// about which bit is which context.
fn pack<F: FnMut(usize) -> bool>(count: usize, mut predicate: F) -> Vec<u8> {
    let mut bitmap = vec![0u8; count.div_ceil(8)];
    for index in 0..count {
        if predicate(index) {
            bitmap[index / 8] |= 1u8 << (index % 8);
        }
    }
    bitmap
}

/// Return whether two closed intervals share at least one position.
///
/// The scalar primitive `ruleset.facts.positions.genomic_overlap` ports to.
#[inline]
fn intervals_overlap(range_start: i64, range_end: i64, other_start: i64, other_end: i64) -> bool {
    range_start <= other_end && range_end >= other_start
}

/// Four variant shape facts from two length arrays, in one crossing.
///
/// Returns `(is_snv, is_short_insertion, is_short_deletion, length_change_inframe)`.
#[pyfunction]
pub fn batch_variant_shapes(
    reference_lengths: Vec<i64>,
    alternate_lengths: Vec<i64>,
) -> PyResult<VariantShapeBitmaps> {
    let count = reference_lengths.len();
    if alternate_lengths.len() != count {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "reference_lengths and alternate_lengths must be the same length",
        ));
    }
    Ok((
        pack(count, |i| {
            reference_lengths[i] == alternate_lengths[i] && reference_lengths[i] > 0
        }),
        pack(count, |i| alternate_lengths[i] > reference_lengths[i]),
        pack(count, |i| reference_lengths[i] > alternate_lengths[i]),
        // Python's `%` on a negative left operand returns a non-negative result, and Rust's
        // returns a negative one -- so a plain `% 3 == 0` would disagree for a deletion. Compare
        // against zero via rem_euclid to match `allele_length_difference(...) % 3 == 0` exactly.
        pack(count, |i| {
            (alternate_lengths[i] - reference_lengths[i]).rem_euclid(3) == 0
        }),
    ))
}

/// `overlap_feature` for a whole batch: does each variant span touch its own feature's span?
#[pyfunction]
pub fn batch_overlap_span(
    variant_mins: Vec<i64>,
    variant_maxs: Vec<i64>,
    feature_starts: Vec<i64>,
    feature_ends: Vec<i64>,
) -> PyResult<Vec<u8>> {
    let count = variant_mins.len();
    if variant_maxs.len() != count || feature_starts.len() != count || feature_ends.len() != count {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "every input array must be the same length",
        ));
    }
    Ok(pack(count, |i| {
        intervals_overlap(
            variant_mins[i],
            variant_maxs[i],
            feature_starts[i],
            feature_ends[i],
        )
    }))
}

/// `contain_feature` for a whole batch: does each variant span swallow its feature whole?
#[pyfunction]
pub fn batch_contain_span(
    variant_mins: Vec<i64>,
    variant_maxs: Vec<i64>,
    feature_starts: Vec<i64>,
    feature_ends: Vec<i64>,
) -> PyResult<Vec<u8>> {
    let count = variant_mins.len();
    if variant_maxs.len() != count || feature_starts.len() != count || feature_ends.len() != count {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "every input array must be the same length",
        ));
    }
    Ok(pack(count, |i| {
        variant_mins[i] <= feature_starts[i] && variant_maxs[i] >= feature_ends[i]
    }))
}

/// Does each variant span overlap *any* interval in its own context's interval list?
///
/// The kernel behind `overlap_exon`/`overlap_intron`. Interval lists are passed CSR-style -- one
/// flat `(starts, ends)` pair plus an `offsets` array of length `count + 1`, where context `i`
/// owns intervals `offsets[i]..offsets[i + 1]`, so a batch spanning many transcripts is still
/// exactly three arrays and no nesting. Contexts sharing a transcript can point at the same
/// slice, which is the common case for the heavy-locus corpus.
#[pyfunction]
pub fn batch_overlap_any_interval(
    variant_mins: Vec<i64>,
    variant_maxs: Vec<i64>,
    interval_starts: Vec<i64>,
    interval_ends: Vec<i64>,
    offsets: Vec<usize>,
) -> PyResult<Vec<u8>> {
    let count = variant_mins.len();
    if variant_maxs.len() != count || offsets.len() != count + 1 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "variant arrays must match and offsets must have count + 1 entries",
        ));
    }
    if interval_ends.len() != interval_starts.len() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "interval_starts and interval_ends must be the same length",
        ));
    }
    if offsets[count] > interval_starts.len() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "offsets run past the end of the interval arrays",
        ));
    }
    Ok(pack(count, |i| {
        let (variant_min, variant_max) = (variant_mins[i], variant_maxs[i]);
        (offsets[i]..offsets[i + 1]).any(|slot| {
            intervals_overlap(
                variant_min,
                variant_max,
                interval_starts[slot],
                interval_ends[slot],
            )
        })
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Read bit `index` out of a packed bitmap.
    fn bit(bitmap: &[u8], index: usize) -> bool {
        bitmap[index / 8] & (1u8 << (index % 8)) != 0
    }

    #[test]
    fn pack_sets_the_expected_bits_across_a_byte_boundary() {
        let bitmap = pack(10, |index| index == 0 || index == 7 || index == 9);
        assert_eq!(bitmap.len(), 2);
        for index in 0..10 {
            assert_eq!(
                bit(&bitmap, index),
                matches!(index, 0 | 7 | 9),
                "index {index}"
            );
        }
    }

    #[test]
    fn variant_shapes_match_the_python_definitions() {
        // (ref, alt): SNV, insertion, deletion, and an out-of-frame deletion.
        let reference_lengths = vec![1, 0, 3, 2];
        let alternate_lengths = vec![1, 4, 0, 0];
        let (snv, increases, decreases, inframe) =
            batch_variant_shapes(reference_lengths, alternate_lengths).unwrap();

        assert!(bit(&snv, 0) && !bit(&snv, 1) && !bit(&snv, 2));
        assert!(!bit(&increases, 0) && bit(&increases, 1) && !bit(&increases, 2));
        assert!(!bit(&decreases, 0) && !bit(&decreases, 1) && bit(&decreases, 2));
        // 0 diff and -3 diff are in frame; +4 and -2 are not.
        assert!(bit(&inframe, 0) && !bit(&inframe, 1) && bit(&inframe, 2) && !bit(&inframe, 3));
    }

    #[test]
    fn inframe_matches_python_modulo_for_negative_differences() {
        // Rust's `%` would give -0/-1/-2 here; Python's gives 0/1/2. Only rem_euclid agrees.
        let (_, _, _, inframe) = batch_variant_shapes(vec![3, 4, 5], vec![0, 0, 0]).unwrap();
        assert!(bit(&inframe, 0), "-3 is in frame");
        assert!(!bit(&inframe, 1), "-4 is not");
        assert!(!bit(&inframe, 2), "-5 is not");
    }

    #[test]
    fn zero_length_alleles_are_not_snvs() {
        let (snv, _, _, _) = batch_variant_shapes(vec![0], vec![0]).unwrap();
        assert!(!bit(&snv, 0), "an empty-to-empty edit is not an SNV");
    }

    #[test]
    fn overlap_span_is_inclusive_at_both_ends() {
        let mask = batch_overlap_span(
            vec![100, 99, 200, 201],
            vec![100, 99, 200, 201],
            vec![100, 100, 100, 100],
            vec![200, 200, 200, 200],
        )
        .unwrap();
        assert!(bit(&mask, 0), "first base of the feature is inside it");
        assert!(!bit(&mask, 1), "one before is outside");
        assert!(bit(&mask, 2), "last base of the feature is inside it");
        assert!(!bit(&mask, 3), "one after is outside");
    }

    #[test]
    fn overlap_span_uses_each_contexts_own_feature_bounds() {
        // Same variant position, different features: the kernel must not share bounds.
        let mask = batch_overlap_span(
            vec![150, 150],
            vec![150, 150],
            vec![100, 300],
            vec![200, 400],
        )
        .unwrap();
        assert!(bit(&mask, 0));
        assert!(!bit(&mask, 1));
    }

    #[test]
    fn contain_span_requires_covering_the_whole_feature() {
        let mask = batch_contain_span(
            vec![100, 101, 99],
            vec![200, 200, 201],
            vec![100, 100, 100],
            vec![200, 200, 200],
        )
        .unwrap();
        assert!(bit(&mask, 0), "exactly covering counts as containing");
        assert!(!bit(&mask, 1), "one base short at the start does not");
        assert!(bit(&mask, 2), "overhanging both ends does");
    }

    #[test]
    fn overlap_any_interval_respects_per_context_slices() {
        // Context 0 owns intervals [0,2), context 1 owns [2,3).
        let mask = batch_overlap_any_interval(
            vec![150, 150],
            vec![150, 150],
            vec![100, 300, 400],
            vec![120, 350, 450],
            vec![0, 2, 3],
        )
        .unwrap();
        assert!(!bit(&mask, 0), "150 misses both of context 0's intervals");
        assert!(!bit(&mask, 1), "150 misses context 1's interval too");

        let mask = batch_overlap_any_interval(
            vec![110, 410],
            vec![110, 410],
            vec![100, 300, 400],
            vec![120, 350, 450],
            vec![0, 2, 3],
        )
        .unwrap();
        assert!(bit(&mask, 0), "110 hits context 0's first interval");
        assert!(bit(&mask, 1), "410 hits context 1's interval");
    }

    #[test]
    fn overlap_any_interval_handles_an_empty_slice() {
        let mask =
            batch_overlap_any_interval(vec![150], vec![150], vec![], vec![], vec![0, 0]).unwrap();
        assert!(
            !bit(&mask, 0),
            "a context owning no intervals overlaps nothing"
        );
    }

    #[test]
    fn mismatched_array_lengths_are_rejected() {
        assert!(batch_variant_shapes(vec![1, 2], vec![1]).is_err());
        assert!(batch_overlap_span(vec![1], vec![1], vec![1], vec![]).is_err());
        assert!(batch_overlap_any_interval(vec![1], vec![1], vec![], vec![], vec![0]).is_err());
        assert!(batch_overlap_any_interval(vec![1], vec![1], vec![1], vec![], vec![0, 1]).is_err());
    }
}
