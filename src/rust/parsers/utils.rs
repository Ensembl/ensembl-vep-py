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

use noodles_core::Region;
use noodles_bgzf as bgzf;
use noodles_gff::{
    feature::record::{Phase, Strand},
};
use pyo3::{PyResult, exceptions::PyValueError};
use std::{
    fs::File,
    io::{self, ErrorKind, Read, BufRead, BufReader},
    path::Path,
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        parsers::test_utils::write_bgzf_file,
        test_utils::{create_temp_dir, write_sample_file},
    };
    use std::fs;

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
    // parse_region
    // ---------------------------------------------------------------------------

    #[test]
    fn test_parse_region_valid() {
        assert!(parse_region("chr1:1-1000").is_ok());
    }

    #[test]
    fn test_parse_region_invalid() {
        let err = parse_region("not_a_region::??").unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidInput);
        assert!(err.to_string().contains("not_a_region::??"));
    }

    #[test]
    fn test_parse_region_empty() {
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
