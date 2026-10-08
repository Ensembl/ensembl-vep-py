// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_gff::feature::record::{Phase, Strand};
use pyo3::{PyResult, exceptions::PyValueError};
use std::{
    fs::File,
    io::{self, BufRead, BufReader, ErrorKind, Read},
    path::{Path, PathBuf},
};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum ReaderMode {
    Auto,
    Sequential,
    Query,
}

impl ReaderMode {
    pub(crate) fn parse(value: Option<&str>) -> PyResult<Self> {
        match value.unwrap_or("auto") {
            "auto" => Ok(Self::Auto),
            "sequential" => Ok(Self::Sequential),
            "query" => Ok(Self::Query),
            value => Err(PyValueError::new_err(format!(
                "unsupported reader mode {value:?}: expected 'auto', 'query' or 'sequential'"
            ))),
        }
    }

    pub(crate) fn as_str(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Sequential => "sequential",
            Self::Query => "query",
        }
    }
}

pub(crate) fn is_bgzf_file(file_path: &Path) -> io::Result<bool> {
    let mut file = File::open(file_path)?;
    is_bgzf_reader(&mut file)
}

fn is_bgzf_reader<R>(reader: &mut R) -> io::Result<bool>
where
    R: Read,
{
    const GZIP_ID1: u8 = 0x1f;
    const GZIP_ID2: u8 = 0x8b;
    const DEFLATE: u8 = 0x08;
    const FEXTRA: u8 = 0x04;
    const BGZF_SUBFIELD_ID: [u8; 2] = *b"BC";
    const BGZF_SUBFIELD_LEN: usize = 2;
    const GZIP_HEADER_LEN: usize = 12;

    let mut header = [0; GZIP_HEADER_LEN];
    match reader.read_exact(&mut header) {
        Ok(()) => {}
        Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(false),
        Err(e) => return Err(e),
    }

    let is_gzip = header[0] == GZIP_ID1 && header[1] == GZIP_ID2 && header[2] == DEFLATE;
    if !is_gzip || header[3] & FEXTRA == 0 {
        return Ok(false);
    }

    let extra_len = u16::from_le_bytes([header[10], header[11]]) as usize;
    let mut extra = vec![0; extra_len];
    match reader.read_exact(&mut extra) {
        Ok(()) => {}
        Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(false),
        Err(e) => return Err(e),
    }

    let mut i = 0;
    while i + 4 <= extra.len() {
        let subfield_id = [extra[i], extra[i + 1]];
        let subfield_len = u16::from_le_bytes([extra[i + 2], extra[i + 3]]) as usize;
        i += 4;

        if i + subfield_len > extra.len() {
            return Ok(false);
        }

        if subfield_id == BGZF_SUBFIELD_ID && subfield_len == BGZF_SUBFIELD_LEN {
            return Ok(true);
        }

        i += subfield_len;
    }

    Ok(false)
}

pub fn open_sequential_reader(file_path: &Path) -> io::Result<Box<dyn BufRead + Send + Sync>> {
    let file = File::open(file_path)?;

    if is_bgzf_file(file_path)? {
        Ok(Box::new(bgzf::io::Reader::new(file)))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

pub fn associated_index_path(src: &Path, extension: &str) -> PathBuf {
    PathBuf::from(format!("{}.{}", src.to_string_lossy(), extension))
}

pub(crate) fn noodles_strand_to_string(strand: Strand) -> String {
    match strand {
        Strand::None => ".",
        Strand::Forward => "+",
        Strand::Reverse => "-",
        Strand::Unknown => "?",
    }
    .to_string()
}

pub(crate) fn noodles_phase_to_u8(phase: Phase) -> u8 {
    match phase {
        Phase::Zero => 0,
        Phase::One => 1,
        Phase::Two => 2,
    }
}

pub(crate) fn parse_region(region: &str) -> io::Result<Region> {
    region.parse().map_err(|e| {
        io::Error::new(
            ErrorKind::InvalidInput,
            format!("invalid tabix region {region:?}: {e}"),
        )
    })
}

/// Case-insensitive column-3 (feature type) filter, applied to a raw tabix/GFF/GTF line before it
/// is parsed into a record or converted into a Python object.
///
/// `records(feature_types=...)` on the GFF3 and GTF parsers uses this to skip non-matching lines
/// at the cost of one tab-split byte scan, instead of allocating a full record only to discard it.
#[derive(Clone, Debug)]
pub(crate) struct FeatureTypeFilter {
    lower_case_types: Vec<String>,
}

impl FeatureTypeFilter {
    pub(crate) fn new(feature_types: &[String]) -> Self {
        Self {
            lower_case_types: feature_types
                .iter()
                .map(|t| t.to_ascii_lowercase())
                .collect(),
        }
    }

    /// Return whether the raw line's column-3 feature type is in this filter's set.
    ///
    /// Compares case-insensitively without allocating: the field is sliced out of `line` by byte
    /// offset, and each candidate is compared with `eq_ignore_ascii_case`.
    pub(crate) fn matches(&self, line: &[u8]) -> bool {
        // A line with fewer than 3 tab-separated fields has no column 3 to compare: treated as a
        // non-match (silently skipped) rather than falling through to the full parser, which
        // would raise on it. Only reachable when a filter is active -- unfiltered parsing still
        // surfaces such a line as the parse error it always was.
        let Some(feature_type) = nth_tab_field(line, 2) else {
            return false;
        };
        self.lower_case_types
            .iter()
            .any(|candidate| feature_type.eq_ignore_ascii_case(candidate.as_bytes()))
    }
}

/// Return the `index`-th tab-delimited field (0-based) of a raw line, if present.
fn nth_tab_field(line: &[u8], index: usize) -> Option<&[u8]> {
    let mut field_index = 0;
    let mut start = 0;

    for (i, &byte) in line.iter().enumerate() {
        if byte == b'\t' {
            if field_index == index {
                return Some(&line[start..i]);
            }
            field_index += 1;
            start = i + 1;
        }
    }

    if field_index == index {
        Some(&line[start..])
    } else {
        None
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct SequenceRegion {
    pub(crate) seq_region_name: String,
    pub(crate) start: Option<usize>,
    pub(crate) end: Option<usize>,
}

pub(crate) fn parse_sequence_region(region: &str) -> io::Result<SequenceRegion> {
    let invalid_region = |message: String| io::Error::new(ErrorKind::InvalidInput, message);

    if region.is_empty() {
        return Err(invalid_region(
            "invalid sequence region: empty region".to_string(),
        ));
    }

    let Some((seq_region_name, rest)) = region.split_once(':') else {
        return Ok(SequenceRegion {
            seq_region_name: region.to_string(),
            start: None,
            end: None,
        });
    };

    if seq_region_name.is_empty() || rest.is_empty() {
        return Err(invalid_region(format!(
            "invalid sequence region {region:?}: missing sequence name or coordinates"
        )));
    }

    let parse_position = |value: &str| -> io::Result<usize> {
        let position = value
            .parse::<usize>()
            .map_err(|err| invalid_region(format!("invalid sequence region {region:?}: {err}")))?;

        if position == 0 {
            return Err(invalid_region(format!(
                "invalid sequence region {region:?}: coordinates must be 1-based"
            )));
        }

        Ok(position)
    };

    let Some((start, end)) = rest.split_once('-') else {
        let position = parse_position(rest)?;
        return Ok(SequenceRegion {
            seq_region_name: seq_region_name.to_string(),
            start: Some(position),
            end: None,
        });
    };

    let start = parse_position(start)?;
    let end = parse_position(end)?;

    if start > end {
        return Err(invalid_region(format!(
            "invalid sequence region {region:?}: start must be <= end"
        )));
    }

    Ok(SequenceRegion {
        seq_region_name: seq_region_name.to_string(),
        start: Some(start),
        end: Some(end),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        parsers::test_utils::write_bgzf_file,
        test_utils::{create_temp_dir, write_sample_file},
    };
    use noodles_core::Position;
    use pyo3::Python;
    use std::{fs, ops::Bound};

    // ---------------------------------------------------------------------------
    // Helpers
    // ---------------------------------------------------------------------------

    const SAMPLE: &[u8] = b"\
##gff-version 3
chr1\tsource\tgene\t1\t10\t.\t+\t.\tID=gene1
";

    fn sample_lines() -> Vec<String> {
        String::from_utf8_lossy(SAMPLE)
            .lines()
            .map(str::to_owned)
            .collect()
    }

    // ---------------------------------------------------------------------------
    // parse_sequence_region
    // ---------------------------------------------------------------------------

    #[test]
    fn test_parse_sequence_region_seq_region_name_only() -> io::Result<()> {
        assert_eq!(
            parse_sequence_region("chr1")?,
            SequenceRegion {
                seq_region_name: "chr1".to_string(),
                start: None,
                end: None,
            }
        );
        Ok(())
    }

    #[test]
    fn test_parse_sequence_region_range() -> io::Result<()> {
        assert_eq!(
            parse_sequence_region("chr1:10-20")?,
            SequenceRegion {
                seq_region_name: "chr1".to_string(),
                start: Some(10),
                end: Some(20),
            }
        );
        Ok(())
    }

    #[test]
    fn test_parse_sequence_region_start_only() -> io::Result<()> {
        assert_eq!(
            parse_sequence_region("chr1:10")?,
            SequenceRegion {
                seq_region_name: "chr1".to_string(),
                start: Some(10),
                end: None,
            }
        );
        Ok(())
    }

    #[test]
    fn test_parse_sequence_region_start_only_with_trailing_hyphen() {
        let err = parse_sequence_region("chr1:10-").unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidInput);
        assert!(
            err.to_string()
                .contains("cannot parse integer from empty string")
        );
    }

    #[test]
    fn test_parse_sequence_region_just_hyphen() {
        let err = parse_sequence_region("chr1:-").unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidInput);
        assert!(err.to_string().contains("invalid sequence region"));
    }

    #[test]
    fn test_parse_sequence_region_rejects_missing_coordinates() {
        let err = parse_sequence_region("chr1:").unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidInput);
        assert!(err.to_string().contains("chr1:"));
    }

    #[test]
    fn test_parse_sequence_region_rejects_empty_region() {
        let err = parse_sequence_region("").unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidInput);
        assert!(err.to_string().contains("empty region"));
    }

    #[test]
    fn test_parse_sequence_region_rejects_zero_position() {
        let err = parse_sequence_region("chr1:0-10").unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidInput);
        assert!(err.to_string().contains("1-based"));
    }

    #[test]
    fn test_parse_sequence_region_rejects_inverted_range() {
        let err = parse_sequence_region("chr1:20-10").unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidInput);
        assert!(err.to_string().contains("start must be <= end"));
    }

    // ---------------------------------------------------------------------------
    // parse_region
    // ---------------------------------------------------------------------------

    #[test]
    fn test_parse_region_valid() {
        assert!(parse_region("chr1:1-1000").is_ok());
    }

    #[test]
    fn test_parse_region_start_only() -> io::Result<()> {
        let region = parse_region("chr1:10")?;
        assert_eq!(region.name(), "chr1".to_string());
        let start = Position::try_from(10).unwrap();
        assert_eq!(region.start(), Bound::Included(start));
        assert_eq!(region.end(), Bound::Unbounded);
        Ok(())
    }

    #[test]
    fn test_parse_region_start_only_with_trailing_hyphen() {
        let err = parse_region("chr1:10-").unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidInput);
        assert!(err.to_string().contains("invalid tabix region"));
    }

    #[test]
    fn test_parse_region_just_hyphen() {
        let err = parse_region("chr1:-").unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidInput);
        assert!(err.to_string().contains("invalid tabix region"));
    }

    #[test]
    fn test_parse_region_invalid() {
        let err = parse_region("not_a_region::??").unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidInput);
        assert!(err.to_string().contains("not_a_region::??"));
    }

    #[test]
    fn test_parse_region_rejects_empty_region() {
        let err = parse_region("").unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidInput);
        assert!(err.to_string().contains("invalid tabix region"));
    }

    // -----------------------------------------------------------------------
    // noodles_phase_to_u8
    // -----------------------------------------------------------------------

    #[test]
    fn test_noodles_phase_to_u8() {
        assert_eq!(noodles_phase_to_u8(Phase::Zero), 0);
        assert_eq!(noodles_phase_to_u8(Phase::One), 1);
        assert_eq!(noodles_phase_to_u8(Phase::Two), 2);
    }

    // -----------------------------------------------------------------------
    // noodles_strand_to_string
    // -----------------------------------------------------------------------

    #[test]
    fn test_noodles_strand_to_string() {
        assert_eq!(noodles_strand_to_string(Strand::None), ".");
        assert_eq!(noodles_strand_to_string(Strand::Forward), "+");
        assert_eq!(noodles_strand_to_string(Strand::Reverse), "-");
        assert_eq!(noodles_strand_to_string(Strand::Unknown), "?");
    }

    // -----------------------------------------------------------------------
    // associated_index_path
    // -----------------------------------------------------------------------

    #[test]
    fn test_associated_index_path_csi() {
        assert_eq!(
            associated_index_path(Path::new("annotations.gff3.gz"), "csi"),
            PathBuf::from("annotations.gff3.gz.csi")
        );
    }

    #[test]
    fn test_associated_index_path_tbi() {
        assert_eq!(
            associated_index_path(Path::new("sample.vcf.gz"), "tbi"),
            PathBuf::from("sample.vcf.gz.tbi")
        );
    }

    // ---------------------------------------------------------------------------
    // is_bgzf_reader
    // ---------------------------------------------------------------------------

    #[test]
    fn test_is_bgzf_file_detects_bgzf_magic() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-gff3-is-bgzf")?;
        let path = write_bgzf_file(&dir, sample_lines(), Some("sample.gff3.gz"))?;

        assert!(is_bgzf_file(&path)?);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_is_bgzf_file_rejects_plain_file_with_gz_extension() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-gff3-is-bgzf-plain")?;
        let path = write_sample_file(&dir, "sample.gff3.gz", SAMPLE)?;

        assert!(!is_bgzf_file(&path)?);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_is_bgzf_file_rejects_plain_file() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-gff3-is-bgzf-plain-no-ext")?;
        let path = write_sample_file(&dir, "sample.gff3", SAMPLE)?;

        assert!(!is_bgzf_file(&path)?);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_is_bgzf_file_rejects_nonexistent_file() {
        let err = is_bgzf_file(Path::new("nonexistent_file.gz")).unwrap_err();
        assert_eq!(err.kind(), ErrorKind::NotFound);
    }

    // -----------------------------------------------------------------------
    // nth_tab_field
    // -----------------------------------------------------------------------

    #[test]
    fn test_nth_tab_field_returns_middle_field() {
        assert_eq!(
            nth_tab_field(b"chr1\tensembl\tgene\t10\t20", 2),
            Some(b"gene".as_slice())
        );
    }

    #[test]
    fn test_nth_tab_field_returns_last_field_without_trailing_tab() {
        assert_eq!(nth_tab_field(b"a\tb\tc", 2), Some(b"c".as_slice()));
    }

    #[test]
    fn test_nth_tab_field_missing_field_returns_none() {
        assert_eq!(nth_tab_field(b"a\tb", 2), None);
    }

    // -----------------------------------------------------------------------
    // FeatureTypeFilter
    // -----------------------------------------------------------------------

    #[test]
    fn test_feature_type_filter_matches_case_insensitively() {
        let filter = FeatureTypeFilter::new(&["mRNA".to_string(), "transcript".to_string()]);

        assert!(filter.matches(b"chr1\tensembl\tmRNA\t10\t20\t.\t+\t."));
        assert!(filter.matches(b"chr1\tensembl\tmrna\t10\t20\t.\t+\t."));
        assert!(filter.matches(b"chr1\tensembl\tMRNA\t10\t20\t.\t+\t."));
        assert!(filter.matches(b"chr1\tensembl\tTranscript\t10\t20\t.\t+\t."));
    }

    #[test]
    fn test_feature_type_filter_rejects_non_matching_type() {
        let filter = FeatureTypeFilter::new(&["mRNA".to_string()]);

        assert!(!filter.matches(b"chr1\tensembl\texon\t10\t20\t.\t+\t."));
    }

    #[test]
    fn test_feature_type_filter_rejects_line_missing_type_field() {
        let filter = FeatureTypeFilter::new(&["mRNA".to_string()]);

        assert!(!filter.matches(b"chr1\tensembl"));
    }

    // -----------------------------------------------------------------------
    // ReaderMode
    // -----------------------------------------------------------------------

    #[test]
    fn test_reader_mode_parse_defaults_to_auto() {
        assert_eq!(ReaderMode::parse(None).unwrap(), ReaderMode::Auto);
    }

    #[test]
    fn test_reader_mode_parse_accepts_aliases() {
        assert_eq!(ReaderMode::parse(Some("query")).unwrap(), ReaderMode::Query);
        assert_eq!(
            ReaderMode::parse(Some("sequential")).unwrap(),
            ReaderMode::Sequential
        );
        assert_eq!(ReaderMode::parse(Some("auto")).unwrap(), ReaderMode::Auto);
    }

    #[test]
    fn test_reader_mode_parse_rejects_invalid_mode() {
        Python::initialize();

        let err = ReaderMode::parse(Some("invalid-mode")).unwrap_err();
        assert!(err.to_string().contains(
            "unsupported reader mode \"invalid-mode\": expected 'auto', 'query' or 'sequential'"
        ));
    }

    #[test]
    fn test_reader_mode_as_str() {
        assert_eq!(ReaderMode::Sequential.as_str(), "sequential");
        assert_eq!(ReaderMode::Query.as_str(), "query");
        assert_eq!(ReaderMode::Auto.as_str(), "auto");
    }
}
