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

use crate::{
    parsers::utils::{ReaderMode, associated_index_path, is_bgzf_file, open_sequential_reader},
    utils::to_py_error,
};
use noodles_bgzf as bgzf;
use pyo3::{exceptions::PyRuntimeError, prelude::*, types::PyDict};
use std::{
    collections::HashMap,
    fs::File,
    io::{self, BufRead, BufReader, Cursor, Read, Seek, SeekFrom},
    path::{Path, PathBuf},
};

// ---------------------------------------------------------------------------
// FASTA record
// ---------------------------------------------------------------------------

#[pyclass(name = "FastaRecord")]
#[derive(Clone, Debug, PartialEq)]
pub struct FastaRecord {
    #[pyo3(get)]
    pub seqid: String,
    #[pyo3(get)]
    pub start: usize,
    #[pyo3(get)]
    pub end: usize,
    #[pyo3(get)]
    pub sequence: String,
}

#[pymethods]
impl FastaRecord {
    fn __repr__(&self) -> String {
        format!(
            "FastaRecord(seqid={:?}, start={}, end={}, sequence_length={})",
            self.seqid,
            self.start,
            self.end,
            self.sequence.len()
        )
    }
}

// ---------------------------------------------------------------------------
// FASTAMetadata
// ---------------------------------------------------------------------------

#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct FastaMetadata {
    sequence_lengths: HashMap<String, usize>,
}

// ---------------------------------------------------------------------------
// FASTA index
// ---------------------------------------------------------------------------

#[derive(Clone, Debug, PartialEq, Eq)]
struct FaiEntry {
    length: usize,
    offset: u64,
    line_bases: usize,
    line_width: usize,
}

impl FaiEntry {
    fn byte_offset(&self, pos_1_based: usize) -> u64 {
        let zero_based = pos_1_based.saturating_sub(1);
        let lines = zero_based / self.line_bases;
        let within_line = zero_based % self.line_bases;
        self.offset + (lines * self.line_width) as u64 + within_line as u64
    }
}

type FaiIndex = HashMap<String, FaiEntry>;

// ---------------------------------------------------------------------------
// Reader
// ---------------------------------------------------------------------------

pub struct FastaReader {
    file_path: PathBuf,
    backend: Backend,
    mode: ReaderMode,
    metadata: Option<FastaMetadata>,
}

impl FastaReader {
    fn open(file_path: &Path, mode: ReaderMode) -> io::Result<Self> {
        let (backend, mode) = Backend::open(file_path, mode)?;
        Ok(Self {
            file_path: file_path.to_path_buf(),
            backend,
            mode,
            metadata: None,
        })
    }

    fn metadata(&mut self) -> io::Result<&HashMap<String, usize>> {
        if self.metadata.is_none() {
            self.metadata = Some(read_fasta_metadata(&self.file_path)?);
        }

        Ok(&self
            .metadata
            .as_ref()
            .expect("metadata should be cached")
            .sequence_lengths)
    }

    fn next(&mut self) -> io::Result<Option<FastaRecord>> {
        self.backend.next()
    }

    fn seek(&mut self, region: &str) -> io::Result<bool> {
        self.backend.seek(region)
    }
}

enum Backend {
    Sequential(SequentialBackend),
    Fai(FaiBackend),
    Bgzf(BgzfBackend),
}

impl Backend {
    fn open(file_path: &Path, mode: ReaderMode) -> io::Result<(Self, ReaderMode)> {
        match mode {
            ReaderMode::Auto => {
                if is_bgzf_file(file_path)? {
                    if has_bgzf_sequence_index(file_path) {
                        BgzfBackend::new(file_path)
                            .map(Self::Bgzf)
                            .map(|backend| (backend, ReaderMode::Query))
                    } else {
                        SequentialBackend::new(file_path)
                            .map(Self::Sequential)
                            .map(|backend| (backend, ReaderMode::Sequential))
                    }
                } else if has_fai_index(file_path) {
                    FaiBackend::new(file_path)
                        .map(Self::Fai)
                        .map(|backend| (backend, ReaderMode::Query))
                } else {
                    SequentialBackend::new(file_path)
                        .map(Self::Sequential)
                        .map(|backend| (backend, ReaderMode::Sequential))
                }
            }
            ReaderMode::Sequential => SequentialBackend::new(file_path)
                .map(Self::Sequential)
                .map(|backend| (backend, ReaderMode::Sequential)),
            ReaderMode::Query => {
                if is_bgzf_file(file_path)? {
                    if has_bgzf_sequence_index(file_path) {
                        BgzfBackend::new(file_path)
                            .map(Self::Bgzf)
                            .map(|backend| (backend, ReaderMode::Query))
                    } else {
                        Err(io::Error::new(
                            io::ErrorKind::NotFound,
                            format!(
                                "no FASTA index found for file {:?}; expected {:?}",
                                file_path,
                                expected_index_paths(file_path)
                            ),
                        ))
                    }
                } else if has_fai_index(file_path) {
                    FaiBackend::new(file_path)
                        .map(Self::Fai)
                        .map(|backend| (backend, ReaderMode::Query))
                } else {
                    Err(io::Error::new(
                        io::ErrorKind::NotFound,
                        format!(
                            "no FASTA index found for file {:?}; expected {:?}",
                            file_path,
                            expected_index_paths(file_path)
                        ),
                    ))
                }
            }
        }
    }

    fn next(&mut self) -> io::Result<Option<FastaRecord>> {
        match self {
            Self::Sequential(reader) => reader.next(),
            Self::Fai(reader) => reader.next(),
            Self::Bgzf(reader) => reader.next(),
        }
    }

    fn seek(&mut self, region: &str) -> io::Result<bool> {
        match self {
            Self::Sequential(reader) => reader.seek(region),
            Self::Fai(reader) => reader.seek(region),
            Self::Bgzf(reader) => reader.seek(region),
        }
    }
}

struct SequentialBackend {
    reader: Box<dyn BufRead + Send + Sync>,
    pending_header: Option<String>,
    query_result: Option<FastaRecord>,
}

impl SequentialBackend {
    fn new(file_path: &Path) -> io::Result<Self> {
        Ok(Self {
            reader: open_sequential_reader(file_path)?,
            pending_header: None,
            query_result: None,
        })
    }

    fn next(&mut self) -> io::Result<Option<FastaRecord>> {
        if let Some(record) = self.query_result.take() {
            return Ok(Some(record));
        }

        read_next_record(&mut self.reader, &mut self.pending_header)
    }

    fn seek(&mut self, _region: &str) -> io::Result<bool> {
        Err(io::Error::other(
            "FASTA reader is set to 'sequential'; set reader to 'query' before seek",
        ))
    }
}

struct FaiBackend {
    index: FaiIndex,
    reader: File,
    query_result: Option<FastaRecord>,
}

impl FaiBackend {
    fn new(file_path: &Path) -> io::Result<Self> {
        Ok(Self {
            index: read_fai_index(file_path)?,
            reader: File::open(file_path)?,
            query_result: None,
        })
    }

    fn next(&mut self) -> io::Result<Option<FastaRecord>> {
        Ok(self.query_result.take())
    }

    fn seek(&mut self, region: &str) -> io::Result<bool> {
        let parsed = parse_region(region)?;
        let Some(entry) = self.index.get(&parsed.seqid) else {
            return Ok(false);
        };

        let (start, end) = match (parsed.start, parsed.end) {
            (Some(start), Some(end)) => (start, end),
            _ => (1, entry.length),
        };

        self.query_result = Some(FastaRecord {
            seqid: parsed.seqid,
            start,
            end,
            sequence: read_indexed_sequence(&mut self.reader, entry, start, end)?,
        });

        Ok(true)
    }
}

struct BgzfBackend {
    index: FaiIndex,
    reader: bgzf::io::IndexedReader<File>,
    query_result: Option<FastaRecord>,
}

impl BgzfBackend {
    fn new(file_path: &Path) -> io::Result<Self> {
        let file = File::open(file_path)?;
        let gzi_index = read_gzi_index(file_path)?;

        Ok(Self {
            index: read_fai_index(file_path)?,
            reader: bgzf::io::IndexedReader::new(file, gzi_index),
            query_result: None,
        })
    }

    fn next(&mut self) -> io::Result<Option<FastaRecord>> {
        Ok(self.query_result.take())
    }

    fn seek(&mut self, region: &str) -> io::Result<bool> {
        let parsed = parse_region(region)?;
        let Some(entry) = self.index.get(&parsed.seqid) else {
            return Ok(false);
        };

        let (start, end) = match (parsed.start, parsed.end) {
            (Some(start), Some(end)) => (start, end),
            _ => (1, entry.length),
        };

        self.query_result = Some(FastaRecord {
            seqid: parsed.seqid,
            start,
            end,
            sequence: read_indexed_sequence(&mut self.reader, entry, start, end)?,
        });

        Ok(true)
    }
}

// ---------------------------------------------------------------------------
// PyO3 wrapper
// ---------------------------------------------------------------------------

#[pyclass(name = "FastaParser")]
pub struct FastaParser {
    file_path: Option<PathBuf>,
    reader: Option<FastaReader>,
    reader_mode: ReaderMode,
}

#[pymethods]
impl FastaParser {
    #[new]
    #[pyo3(signature = (file_path=None, mode=None))]
    fn new(file_path: Option<String>, mode: Option<String>) -> PyResult<Self> {
        let reader_mode = ReaderMode::parse(mode.as_deref())?;

        let mut parser = Self {
            file_path: None,
            reader: None,
            reader_mode,
        };

        if let Some(path) = file_path {
            parser.open(&path, mode)?;
        }

        Ok(parser)
    }

    #[getter]
    fn format(&self) -> &'static str {
        "fasta"
    }

    #[getter]
    fn metadata<'py>(&mut self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);

        let Some(reader) = self.reader.as_mut() else {
            return Ok(dict);
        };

        let metadata = reader.metadata().map_err(to_py_error)?;
        if !metadata.is_empty() {
            let lengths = PyDict::new(py);
            for (seqid, length) in metadata {
                lengths.set_item(seqid, length)?;
            }
            dict.set_item("sequence_lengths", lengths)?;
        }

        Ok(dict)
    }

    #[getter]
    fn reader_mode(&self) -> &'static str {
        self.reader_mode.as_str()
    }

    fn open(&mut self, file_path: &str, mode: Option<String>) -> PyResult<()> {
        let path = PathBuf::from(file_path);
        let reader_mode = ReaderMode::parse(mode.as_deref())?;

        self.reader = Some(FastaReader::open(&path, reader_mode).map_err(to_py_error)?);
        self.reader_mode = self
            .reader
            .as_ref()
            .map(|reader| reader.mode)
            .unwrap_or(self.reader_mode);
        self.file_path = Some(path);

        Ok(())
    }

    fn close(&mut self) {
        self.reader = None;
    }

    fn seek(&mut self, region: &str) -> PyResult<bool> {
        self.reader()?.seek(region).map_err(to_py_error)
    }

    #[allow(clippy::should_implement_trait)]
    fn next(&mut self) -> PyResult<Option<FastaRecord>> {
        self.reader()?.next().map_err(to_py_error)
    }

    #[pyo3(signature = (limit=10000))]
    fn records(&mut self, limit: Option<usize>) -> PyResult<Vec<FastaRecord>> {
        let mut records = Vec::new();

        while limit.is_none_or(|n| records.len() < n) {
            match self.next()? {
                Some(record) => records.push(record),
                None => break,
            }
        }

        Ok(records)
    }

    #[staticmethod]
    #[pyo3(signature = (data, limit=None))]
    fn parse_text(data: &str, limit: Option<usize>) -> PyResult<Vec<FastaRecord>> {
        parse_fasta_bytes(data.as_bytes(), limit).map_err(to_py_error)
    }

    #[staticmethod]
    #[pyo3(signature = (file_path, limit=None))]
    fn parse_file(file_path: &str, limit: Option<usize>) -> PyResult<Vec<FastaRecord>> {
        let mut reader = open_sequential_reader(Path::new(file_path)).map_err(to_py_error)?;
        read_records(&mut reader, limit).map_err(to_py_error)
    }
}

impl FastaParser {
    fn reader(&mut self) -> PyResult<&mut FastaReader> {
        self.reader
            .as_mut()
            .ok_or_else(|| PyRuntimeError::new_err("FASTA parser is not open"))
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

#[derive(Clone, Debug, PartialEq, Eq)]
struct ParsedRegion {
    seqid: String,
    start: Option<usize>,
    end: Option<usize>,
}

fn parse_region(region: &str) -> io::Result<ParsedRegion> {
    let Some((seqid, rest)) = region.split_once(':') else {
        return Ok(ParsedRegion {
            seqid: region.to_string(),
            start: None,
            end: None,
        });
    };

    let Some((start, end)) = rest.split_once('-') else {
        let value = rest.parse::<usize>().map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("invalid FASTA region {region:?}: {err}"),
            )
        })?;

        return Ok(ParsedRegion {
            seqid: seqid.to_string(),
            start: Some(value),
            end: Some(value),
        });
    };

    let start = start.parse::<usize>().map_err(|err| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid FASTA region {region:?}: {err}"),
        )
    })?;
    let end = end.parse::<usize>().map_err(|err| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid FASTA region {region:?}: {err}"),
        )
    })?;

    Ok(ParsedRegion {
        seqid: seqid.to_string(),
        start: Some(start),
        end: Some(end),
    })
}

fn has_bgzf_sequence_index(file_path: &Path) -> bool {
    has_fai_index(file_path) && has_gzi_index(file_path)
}

fn expected_index_paths(file_path: &Path) -> Vec<PathBuf> {
    if is_bgzf_file(file_path).unwrap_or(false) {
        vec![
            associated_index_path(file_path, "fai"),
            associated_index_path(file_path, "gzi"),
        ]
    } else {
        vec![associated_index_path(file_path, "fai")]
    }
}

fn has_fai_index(file_path: &Path) -> bool {
    associated_index_path(file_path, "fai").exists()
}

fn has_gzi_index(file_path: &Path) -> bool {
    associated_index_path(file_path, "gzi").exists()
}

fn read_fai_index(file_path: &Path) -> io::Result<FaiIndex> {
    let fai_path = associated_index_path(file_path, "fai");
    let file = File::open(&fai_path)?;
    let mut reader = BufReader::new(file);
    let mut index = HashMap::new();
    let mut line = String::new();

    while reader.read_line(&mut line)? > 0 {
        let trimmed = line.trim_end();
        if trimmed.is_empty() {
            line.clear();
            continue;
        }

        let mut fields = trimmed.split('\t');
        let seqid = fields
            .next()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing seqid"))?;
        let length = fields
            .next()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing length"))?
            .parse::<usize>()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let offset = fields
            .next()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing offset"))?
            .parse::<u64>()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let line_bases = fields
            .next()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing line_bases"))?
            .parse::<usize>()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let line_width = fields
            .next()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing line_width"))?
            .parse::<usize>()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        index.insert(
            seqid.to_string(),
            FaiEntry {
                length,
                offset,
                line_bases,
                line_width,
            },
        );

        line.clear();
    }

    Ok(index)
}

fn read_gzi_index(file_path: &Path) -> io::Result<bgzf::gzi::Index> {
    bgzf::gzi::fs::read(associated_index_path(file_path, "gzi"))
}

fn read_indexed_sequence<R>(
    reader: &mut R,
    entry: &FaiEntry,
    start: usize,
    end: usize,
) -> io::Result<String>
where
    R: Read + Seek,
{
    if start < 1 || end < start {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "FASTA regions must be 1-based and start must be <= end",
        ));
    }

    if start > entry.length {
        return Ok("N".repeat(end - start + 1));
    }

    let actual_end = end.min(entry.length);
    let mut sequence = String::with_capacity(end - start + 1);

    for position in start..=actual_end {
        let byte_offset = entry.byte_offset(position);
        let mut byte = [0u8; 1];
        reader.seek(SeekFrom::Start(byte_offset))?;
        reader.read_exact(&mut byte)?;
        sequence.push(byte[0] as char);
    }

    if end > entry.length {
        sequence.push_str(&"N".repeat(end - entry.length));
    }

    Ok(sequence)
}

fn read_fasta_metadata(file_path: &Path) -> io::Result<FastaMetadata> {
    let mut reader = open_sequential_reader(file_path)?;
    read_fasta_metadata_from_reader(&mut reader)
}

fn read_fasta_metadata_from_reader<R>(reader: &mut R) -> io::Result<FastaMetadata>
where
    R: BufRead,
{
    let mut metadata = FastaMetadata::default();
    let mut pending_header = None;

    while let Some(record) = read_next_record(reader, &mut pending_header)? {
        metadata
            .sequence_lengths
            .insert(record.seqid, record.sequence.len());
    }

    Ok(metadata)
}

fn read_next_record<R>(
    reader: &mut R,
    pending_header: &mut Option<String>,
) -> io::Result<Option<FastaRecord>>
where
    R: BufRead,
{
    let mut line = String::new();
    let header = if let Some(header) = pending_header.take() {
        header
    } else {
        loop {
            line.clear();
            if reader.read_line(&mut line)? == 0 {
                return Ok(None);
            }

            let trimmed = line.trim_end();
            if trimmed.is_empty() {
                continue;
            }

            if let Some(header) = trimmed.strip_prefix('>') {
                break header.to_string();
            }
        }
    };

    let mut sequence = String::new();

    loop {
        line.clear();
        let bytes_read = reader.read_line(&mut line)?;
        if bytes_read == 0 {
            break;
        }

        let trimmed = line.trim_end();
        if trimmed.is_empty() {
            continue;
        }

        if let Some(next_header) = trimmed.strip_prefix('>') {
            *pending_header = Some(next_header.to_string());
            break;
        }

        sequence.push_str(&trimmed.to_uppercase());
    }

    let (seqid, _description) = parse_header(&header);
    let end = sequence.len();

    Ok(Some(FastaRecord {
        seqid,
        start: 1,
        end,
        sequence,
    }))
}

fn read_records<R>(reader: &mut R, limit: Option<usize>) -> io::Result<Vec<FastaRecord>>
where
    R: BufRead,
{
    let mut records = Vec::new();
    let mut pending_header = None;

    while limit.is_none_or(|n| records.len() < n) {
        let Some(record) = read_next_record(reader, &mut pending_header)? else {
            break;
        };
        records.push(record);
    }

    Ok(records)
}

fn parse_fasta_bytes(data: &[u8], limit: Option<usize>) -> io::Result<Vec<FastaRecord>> {
    let mut reader = Cursor::new(data);
    read_records(&mut reader, limit)
}

fn parse_header(header: &str) -> (String, Option<String>) {
    let mut parts = header.splitn(2, char::is_whitespace);
    let seqid = parts.next().unwrap_or_default().to_string();
    let description = parts.next().map(str::to_string);
    (seqid, description)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{parsers::test_utils::write_bgzf_file, test_utils::create_temp_dir};
    use noodles_bgzf::gzi;
    use std::fs;

    const SAMPLE: &[u8] = b"\
>chr1
ACGTTGCA
>chr2
AAAA
";

    #[test]
    fn test_parse_text() -> io::Result<()> {
        let records = parse_fasta_bytes(SAMPLE, None)?;
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].seqid, "chr1");
        assert_eq!(records[0].sequence, "ACGTTGCA");
        Ok(())
    }

    #[test]
    fn test_parse_fasta_bytes_honours_limit() -> io::Result<()> {
        let records = parse_fasta_bytes(SAMPLE, Some(1))?;
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].seqid, "chr1");
        Ok(())
    }

    #[test]
    fn test_read_fasta_metadata_from_reader_collects_lengths() -> io::Result<()> {
        let mut reader = Cursor::new(SAMPLE);
        let metadata = read_fasta_metadata_from_reader(&mut reader)?;
        assert_eq!(metadata.sequence_lengths["chr1"], 8);
        assert_eq!(metadata.sequence_lengths["chr2"], 4);
        Ok(())
    }

    #[test]
    fn test_fai_index() -> io::Result<()> {
        let dir = create_temp_dir("fasta-index-test")?;
        let path = dir.join("sample.fa");
        fs::write(&path, SAMPLE)?;
        fs::write(
            format!("{}.fai", path.to_string_lossy()),
            "chr1\t8\t6\t8\t9\nchr2\t4\t23\t4\t5\n",
        )?;

        let index = read_fai_index(&path)?;
        assert_eq!(index["chr1"].length, 8);
        assert_eq!(index["chr2"].length, 4);
        Ok(())
    }

    #[test]
    fn test_bgzf_backend_reads_indexed_sequence_with_gzi() -> io::Result<()> {
        let dir = create_temp_dir("fasta-bgzf-index-test")?;
        let path = write_bgzf_file(
            &dir,
            vec![">chr1".to_string(), "ACGTTGCA".to_string()],
            Some("sample.fa.gz"),
        )?;
        gzi::fs::write(associated_index_path(&path, "gzi"), &gzi::Index::default())?;
        fs::write(
            format!("{}.fai", path.to_string_lossy()),
            "chr1\t8\t6\t8\t9\n",
        )?;

        let mut backend = BgzfBackend::new(&path)?;
        assert!(backend.seek("chr1:2-5")?);
        let record = backend.next()?.expect("indexed record should be available");
        assert_eq!(record.sequence, "CGTT");
        Ok(())
    }
}
