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
    parsers::utils::{
        ReaderMode, associated_index_path, is_bgzf_file, open_sequential_reader,
        parse_sequence_region,
    },
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

/// Represents a sequence record in a FASTA file.
#[pyclass(name = "FastaRecord")]
#[derive(Clone, Debug, PartialEq)]
pub struct FastaRecord {
    /// str: Name of the seq region
    #[pyo3(get)]
    pub seqid: String,
    /// int: Start coordinate of the sequence
    #[pyo3(get)]
    pub start: usize,
    /// int: End coordinate of the sequence
    #[pyo3(get)]
    pub end: usize,
    /// str: Sequence of the seq region
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
// FastaMetadata
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
    fn byte_offset(&self, pos_1_based: usize) -> io::Result<u64> {
        if self.line_bases == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "FASTA index line_bases must be greater than 0",
            ));
        }

        let zero_based = pos_1_based.saturating_sub(1);
        let lines = zero_based / self.line_bases;
        let within_line = zero_based % self.line_bases;
        Ok(self.offset + (lines * self.line_width) as u64 + within_line as u64)
    }
}

type FaiIndex = HashMap<String, FaiEntry>;

// ---------------------------------------------------------------------------
// Reader
// ---------------------------------------------------------------------------

pub struct FastaReader {
    backend: Backend,
    mode: ReaderMode,
    metadata: Option<FastaMetadata>,
}

impl FastaReader {
    fn open(file_path: &Path, mode: ReaderMode) -> io::Result<Self> {
        let (backend, mode) = Backend::open(file_path, mode)?;
        Ok(Self {
            backend,
            mode,
            metadata: None,
        })
    }

    fn metadata(&mut self) -> io::Result<&FastaMetadata> {
        if self.metadata.is_none() {
            self.metadata = Some(self.backend.metadata()?);
        }

        Ok(self.metadata.as_ref().expect("metadata should be cached"))
    }

    fn next(&mut self) -> io::Result<Option<FastaRecord>> {
        self.backend.next()
    }

    fn seek(&mut self, region: &str) -> io::Result<bool> {
        self.backend.seek(region)
    }
}

// ---------------------------------------------------------------------------
// Backend
// ---------------------------------------------------------------------------

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

    /// Sequence lengths from the `.fai` index, or from a one-off scan when there is none.
    fn metadata(&mut self) -> io::Result<FastaMetadata> {
        let sequence_lengths = match self {
            Self::Sequential(reader) => reader
                .scanned_sequences()?
                .iter()
                .map(|(seqid, sequence)| (seqid.clone(), sequence.length))
                .collect(),
            Self::Fai(reader) => fai_sequence_lengths(&reader.index),
            Self::Bgzf(reader) => fai_sequence_lengths(&reader.index),
        };

        Ok(FastaMetadata { sequence_lengths })
    }

    fn seek(&mut self, region: &str) -> io::Result<bool> {
        match self {
            Self::Sequential(reader) => reader.seek(region),
            Self::Fai(reader) => reader.seek(region),
            Self::Bgzf(reader) => reader.seek(region),
        }
    }
}

fn fai_sequence_lengths(index: &FaiIndex) -> HashMap<String, usize> {
    index
        .iter()
        .map(|(seqid, entry)| (seqid.clone(), entry.length))
        .collect()
}

// ---------------------------------------------------------------------------
// SequentialBackend
// ---------------------------------------------------------------------------

/// Reads an unindexed FASTA.
///
/// Before any `seek`, `next` streams whole records. The first `seek` (or `metadata`) scans the
/// file once and records every sequence's length and, where its lines are regular, the same byte
/// layout a `.fai` would hold. Later seeks read those sequences directly; the rest, and every
/// sequence in a BGZF file (whose offsets are only meaningful with a `.gzi`), are rescanned from
/// the start of the file. The first record wins when a name is repeated.
struct SequentialBackend {
    file_path: PathBuf,
    is_bgzf: bool,
    reader: Box<dyn BufRead + Send + Sync>,
    pending_header: Option<String>,
    scanned_sequences: Option<HashMap<String, ScannedSequence>>,
    /// Handle for direct reads of regular sequences, opened on first use.
    file: Option<File>,
    query_result: Option<FastaRecord>,
    /// Whether `seek` has been called, which switches `next` from streaming to query results.
    query_active: bool,
}

impl SequentialBackend {
    fn new(file_path: &Path) -> io::Result<Self> {
        Ok(Self {
            file_path: file_path.to_path_buf(),
            is_bgzf: is_bgzf_file(file_path)?,
            reader: open_sequential_reader(file_path)?,
            pending_header: None,
            scanned_sequences: None,
            file: None,
            query_result: None,
            query_active: false,
        })
    }

    fn next(&mut self) -> io::Result<Option<FastaRecord>> {
        if self.query_active {
            return Ok(self.query_result.take());
        }
        read_next_record(&mut self.reader, &mut self.pending_header)
    }

    fn seek(&mut self, region: &str) -> io::Result<bool> {
        self.query_active = true;
        self.query_result = None;

        let parsed = parse_sequence_region(region)?;
        let Some(sequence) = self.scanned_sequences()?.get(&parsed.seq_region_name) else {
            return Ok(false);
        };

        let start = parsed.start.unwrap_or(1);
        // WARNING: With no end coordinate, this may return a very long sequence.
        let end = parsed.end.unwrap_or(sequence.length);
        let fai_entry = sequence.fai_entry.clone();

        let bases = match fai_entry {
            Some(entry) => {
                if self.file.is_none() {
                    self.file = Some(File::open(&self.file_path)?);
                }
                let file = self.file.as_mut().expect("file should be open");
                read_indexed_sequence(file, &entry, start, end)?
            }
            None => {
                let mut reader = open_sequential_reader(&self.file_path)?;
                scan_sequence_chunk(&mut reader, &parsed.seq_region_name, start, end)?
            }
        };

        self.query_result = Some(FastaRecord {
            seqid: parsed.seq_region_name,
            start,
            end,
            sequence: bases,
        });

        Ok(true)
    }

    /// Scan the file once for sequence lengths and byte layouts.
    fn scanned_sequences(&mut self) -> io::Result<&HashMap<String, ScannedSequence>> {
        if self.scanned_sequences.is_none() {
            let mut reader = open_sequential_reader(&self.file_path)?;
            self.scanned_sequences = Some(scan_sequences(&mut reader, !self.is_bgzf)?);
        }

        Ok(self
            .scanned_sequences
            .as_ref()
            .expect("scanned sequences should be cached"))
    }
}

// ---------------------------------------------------------------------------
// FaiBackend
// ---------------------------------------------------------------------------

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
        let parsed = parse_sequence_region(region)?;
        let Some(entry) = self.index.get(&parsed.seq_region_name) else {
            self.query_result = None;
            return Ok(false);
        };

        let start = parsed.start.unwrap_or(1);
        // WARNING: With no end coordinate, this may return a very long sequence.
        let end = parsed.end.unwrap_or(entry.length);

        self.query_result = Some(FastaRecord {
            seqid: parsed.seq_region_name,
            start,
            end,
            sequence: read_indexed_sequence(&mut self.reader, entry, start, end)?,
        });

        Ok(true)
    }
}

// ---------------------------------------------------------------------------
// BgzfBackend
// ---------------------------------------------------------------------------

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
        let parsed = parse_sequence_region(region)?;
        let Some(entry) = self.index.get(&parsed.seq_region_name) else {
            self.query_result = None;
            return Ok(false);
        };

        let start = parsed.start.unwrap_or(1);
        // WARNING: With no end coordinate, this may return a very long sequence.
        let end = parsed.end.unwrap_or(entry.length);

        self.query_result = Some(FastaRecord {
            seqid: parsed.seq_region_name,
            start,
            end,
            sequence: read_indexed_sequence(&mut self.reader, entry, start, end)?,
        });

        Ok(true)
    }
}

// ---------------------------------------------------------------------------
// SequenceScan — one-pass sequential scan and `.fai` layout detection
// ---------------------------------------------------------------------------

/// One sequence recorded by a sequential scan.
#[derive(Clone, Debug, PartialEq, Eq)]
struct ScannedSequence {
    length: usize,
    /// Byte layout for direct reads, or `None` when the lines are irregular or offsets were not
    /// tracked.
    fai_entry: Option<FaiEntry>,
}

/// Running state for the sequence currently being scanned.
struct SequenceScan {
    seq_region_name: String,
    length: usize,
    first_base_offset: Option<u64>,
    line_bases: usize,
    line_width: usize,
    /// A short or blank line has been seen, so any further bases make the layout irregular.
    ended: bool,
    regular: bool,
}

impl SequenceScan {
    fn new(seq_region_name: String) -> Self {
        Self {
            seq_region_name,
            length: 0,
            first_base_offset: None,
            line_bases: 0,
            line_width: 0,
            ended: false,
            regular: true,
        }
    }

    /// Count a sequence line's bases and check it keeps the layout `.fai`-compatible.
    fn add_line(&mut self, line: &[u8], line_offset: u64) {
        let bases = line
            .iter()
            .filter(|byte| !byte.is_ascii_whitespace())
            .count();
        self.length += bases;
        if !self.regular {
            return;
        }
        if bases == 0 {
            self.ended = true;
            return;
        }

        let has_newline = line.last() == Some(&b'\n');
        let content_length = line_content(line).len();
        if self.ended || content_length != bases {
            self.regular = false;
            return;
        }

        if self.first_base_offset.is_none() {
            self.first_base_offset = Some(line_offset);
            self.line_bases = bases;
            self.line_width = line.len();
        } else if bases > self.line_bases
            || (bases == self.line_bases && has_newline && line.len() != self.line_width)
        {
            self.regular = false;
            return;
        }

        if bases < self.line_bases || !has_newline {
            self.ended = true;
        }
    }

    fn finish(self, track_offsets: bool) -> (String, ScannedSequence) {
        let fai_entry = (track_offsets && self.regular).then(|| FaiEntry {
            length: self.length,
            offset: self.first_base_offset.unwrap_or(0),
            line_bases: self.line_bases,
            line_width: self.line_width,
        });

        (
            self.seq_region_name,
            ScannedSequence {
                length: self.length,
                fai_entry,
            },
        )
    }
}

/// Scan a FASTA once, recording each sequence's length and, with `track_offsets`, its byte layout.
fn scan_sequences<R>(
    reader: &mut R,
    track_offsets: bool,
) -> io::Result<HashMap<String, ScannedSequence>>
where
    R: BufRead,
{
    let mut sequences = HashMap::new();
    let mut current: Option<SequenceScan> = None;
    let mut line = Vec::new();
    let mut offset = 0u64;

    loop {
        line.clear();
        let bytes_read = reader.read_until(b'\n', &mut line)?;
        if bytes_read == 0 {
            break;
        }
        let line_offset = offset;
        offset += bytes_read as u64;

        if let Some(header) = line.strip_prefix(b">") {
            if let Some(scan) = current.take() {
                let (seq_region_name, sequence) = scan.finish(track_offsets);
                sequences.insert(seq_region_name, sequence);
            }
            let seq_region_name = fasta_header_name(&String::from_utf8_lossy(header)).to_string();
            if !sequences.contains_key(&seq_region_name) {
                current = Some(SequenceScan::new(seq_region_name));
            }
            continue;
        }

        if let Some(scan) = current.as_mut() {
            scan.add_line(&line, line_offset);
        }
    }

    if let Some(scan) = current.take() {
        let (seq_region_name, sequence) = scan.finish(track_offsets);
        sequences.insert(seq_region_name, sequence);
    }

    Ok(sequences)
}

/// Read `[start, end]` of the first sequence named `seq_region_name` by scanning from the current
/// position, padding with `N` past the sequence end.
fn scan_sequence_chunk<R>(
    reader: &mut R,
    seq_region_name: &str,
    start: usize,
    end: usize,
) -> io::Result<String>
where
    R: BufRead,
{
    if start < 1 || end < start {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "FASTA regions must be 1-based and start must be <= end",
        ));
    }

    let mut chunk = Vec::with_capacity(end - start + 1);
    let mut line = Vec::new();
    let mut line_bases = Vec::new();
    let mut position = 0usize;
    let mut in_requested_region = false;

    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }

        if let Some(header) = line.strip_prefix(b">") {
            if in_requested_region {
                break;
            }
            in_requested_region =
                fasta_header_name(&String::from_utf8_lossy(header)) == seq_region_name;
            position = 0;
            continue;
        }
        if !in_requested_region {
            continue;
        }

        line_bases.clear();
        line_bases.extend(
            line.iter()
                .filter(|byte| !byte.is_ascii_whitespace())
                .map(u8::to_ascii_uppercase),
        );
        let line_start = position + 1;
        let line_end = position + line_bases.len();
        let overlap_start = start.max(line_start);
        let overlap_end = end.min(line_end);
        if overlap_start <= overlap_end {
            chunk.extend_from_slice(
                &line_bases[overlap_start - line_start..=overlap_end - line_start],
            );
        }
        position = line_end;
        if position >= end {
            break;
        }
    }

    chunk.resize(end - start + 1, b'N');
    String::from_utf8(chunk).map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))
}

/// Strip a trailing `\n` or `\r\n` from a raw line.
fn line_content(line: &[u8]) -> &[u8] {
    let line = line.strip_suffix(b"\n").unwrap_or(line);
    line.strip_suffix(b"\r").unwrap_or(line)
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

        if !metadata.sequence_lengths.is_empty() {
            let lengths = PyDict::new(py);
            for (seqid, length) in &metadata.sequence_lengths {
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
    let expected_bases = actual_end - start + 1;
    let start_offset = entry.byte_offset(start)?;
    let end_offset = entry.byte_offset(actual_end)?;
    let raw_length = end_offset
        .checked_sub(start_offset)
        .and_then(|length| length.checked_add(1))
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid FASTA index offsets"))?;
    let mut remaining = usize::try_from(raw_length)
        .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
    let mut sequence = Vec::with_capacity(end - start + 1);
    let mut buffer = [0u8; 64 * 1024];

    reader.seek(SeekFrom::Start(start_offset))?;
    while remaining > 0 {
        let chunk_length = remaining.min(buffer.len());
        reader.read_exact(&mut buffer[..chunk_length])?;
        for byte in &buffer[..chunk_length] {
            if !byte.is_ascii_whitespace() {
                sequence.push(byte.to_ascii_uppercase());
            }
        }
        remaining -= chunk_length;
    }

    if sequence.len() != expected_bases {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "FASTA index described {expected_bases} bases but the indexed range contained {}",
                sequence.len()
            ),
        ));
    }

    if end > entry.length {
        sequence.extend(std::iter::repeat_n(b'N', end - entry.length));
    }

    String::from_utf8(sequence).map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))
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

    let end = sequence.len();

    Ok(Some(FastaRecord {
        seqid: fasta_header_name(&header).to_string(),
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

/// Return the sequence-region name from a FASTA header (without its `>`).
fn fasta_header_name(header: &str) -> &str {
    header.split_whitespace().next().unwrap_or_default()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        parsers::test_utils::write_bgzf_file, test_utils::create_temp_dir, utils::from_py_error,
    };
    use noodles_bgzf::gzi;
    use std::{fs, io::Cursor, str};

    const SAMPLE: &[u8] = b"\
>chr1
ACGTTGCA
>chr2
AAAA
";

    const SAMPLE_FAI: &str = "\
chr1\t8\t6\t8\t9
chr2\t4\t21\t4\t5
";

    const DESCRIBED_SAMPLE: &[u8] = b"\
>chr1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF
ACGT
TGCA
>chr2 alternate sequence
AAAA
";

    // -----------------------------------------------------------------------
    // Helpers
    // -----------------------------------------------------------------------

    fn write_sample_fasta_with_fai(prefix: &str) -> io::Result<(PathBuf, PathBuf)> {
        let dir = create_temp_dir(prefix)?;
        let path = dir.join("sample.fa");
        fs::write(&path, SAMPLE)?;
        fs::write(associated_index_path(&path, "fai"), SAMPLE_FAI)?;
        Ok((dir, path))
    }

    struct CountingReader {
        inner: Cursor<Vec<u8>>,
        reads: usize,
        seeks: usize,
    }

    impl CountingReader {
        fn new(data: &[u8]) -> Self {
            Self {
                inner: Cursor::new(data.to_vec()),
                reads: 0,
                seeks: 0,
            }
        }
    }

    impl Read for CountingReader {
        fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
            self.reads += 1;
            self.inner.read(buffer)
        }
    }

    impl Seek for CountingReader {
        fn seek(&mut self, position: SeekFrom) -> io::Result<u64> {
            self.seeks += 1;
            self.inner.seek(position)
        }
    }

    fn write_sample_bgzf_with_indexes(prefix: &str) -> io::Result<(PathBuf, PathBuf)> {
        let dir = create_temp_dir(prefix)?;
        let path = write_bgzf_file(&dir, sample_lines(), Some("sample.fa.gz"))?;
        gzi::fs::write(associated_index_path(&path, "gzi"), &gzi::Index::default())?;
        fs::write(associated_index_path(&path, "fai"), SAMPLE_FAI)?;
        Ok((dir, path))
    }

    fn sample_lines() -> Vec<String> {
        str::from_utf8(SAMPLE)
            .expect("sample should be valid UTF-8")
            .lines()
            .map(str::to_owned)
            .collect()
    }

    fn assert_first_record(record: &FastaRecord) {
        assert_eq!(record.seqid, "chr1");
        assert_eq!(record.start, 1);
        assert_eq!(record.end, 8);
        assert_eq!(record.sequence, "ACGTTGCA");
    }

    fn assert_second_record(record: &FastaRecord) {
        assert_eq!(record.seqid, "chr2");
        assert_eq!(record.start, 1);
        assert_eq!(record.end, 4);
        assert_eq!(record.sequence, "AAAA");
    }

    // -----------------------------------------------------------------------
    // fasta_header_name
    // -----------------------------------------------------------------------

    #[test]
    fn test_fasta_header_name_drops_description() {
        assert_eq!(
            fasta_header_name("chr1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF"),
            "chr1"
        );
    }

    #[test]
    fn test_fasta_header_name_skips_leading_whitespace() {
        assert_eq!(fasta_header_name("  chr1\tdescription"), "chr1");
        assert_eq!(fasta_header_name("chr1"), "chr1");
        assert_eq!(fasta_header_name(""), "");
    }

    // -----------------------------------------------------------------------
    // parse_fasta_bytes
    // -----------------------------------------------------------------------

    #[test]
    fn test_parse_fasta_bytes_decodes_records() -> io::Result<()> {
        let records = parse_fasta_bytes(SAMPLE, None)?;
        assert_eq!(records.len(), 2);
        assert_first_record(&records[0]);
        assert_second_record(&records[1]);
        Ok(())
    }

    #[test]
    fn test_parse_fasta_bytes_honours_limit() -> io::Result<()> {
        let records = parse_fasta_bytes(SAMPLE, Some(1))?;
        assert_eq!(records.len(), 1);
        assert_first_record(&records[0]);
        Ok(())
    }

    #[test]
    fn test_parse_fasta_bytes_honours_zero_limit() -> io::Result<()> {
        let records = parse_fasta_bytes(SAMPLE, Some(0))?;
        assert!(records.is_empty());
        Ok(())
    }

    // -----------------------------------------------------------------------
    // read_records
    // -----------------------------------------------------------------------

    #[test]
    fn test_read_records_reads_multiple_records() -> io::Result<()> {
        let mut reader = Cursor::new(SAMPLE);
        let records = read_records(&mut reader, None)?;
        assert_eq!(records.len(), 2);
        assert_first_record(&records[0]);
        assert_second_record(&records[1]);
        Ok(())
    }

    #[test]
    fn test_read_records_honours_limit() -> io::Result<()> {
        let mut reader = Cursor::new(SAMPLE);
        let records = read_records(&mut reader, Some(1))?;
        assert_eq!(records.len(), 1);
        assert_first_record(&records[0]);
        Ok(())
    }

    #[test]
    fn test_read_records_honours_zero_limit() -> io::Result<()> {
        let mut reader = Cursor::new(SAMPLE);
        let records = read_records(&mut reader, Some(0))?;
        assert!(records.is_empty());
        Ok(())
    }

    // -----------------------------------------------------------------------
    // read_next_record
    // -----------------------------------------------------------------------

    #[test]
    fn test_read_next_record_skips_blank_lines_and_preserves_pending_header() -> io::Result<()> {
        let mut reader = Cursor::new(b"\n>chr1\nACGT\n\n>chr2 desc\nTT\n");
        let mut pending_header = None;

        let first = read_next_record(&mut reader, &mut pending_header)?
            .expect("first record should be available");
        let second = read_next_record(&mut reader, &mut pending_header)?
            .expect("second record should be available");

        assert_eq!(first.seqid, "chr1");
        assert_eq!(first.sequence, "ACGT");
        assert_eq!(second.seqid, "chr2");
        assert_eq!(second.sequence, "TT");
        assert!(read_next_record(&mut reader, &mut pending_header)?.is_none());
        Ok(())
    }

    #[test]
    fn test_read_next_record_uppercases_sequence() -> io::Result<()> {
        let mut reader = Cursor::new(b">chr1\nacgt\n");
        let mut pending_header = None;
        let record = read_next_record(&mut reader, &mut pending_header)?.expect("expected record");

        assert_eq!(record.sequence, "ACGT");
        Ok(())
    }

    #[test]
    fn test_read_next_record_returns_none_without_header() -> io::Result<()> {
        let mut reader = Cursor::new(b"\nACGT\n");
        let mut pending_header = None;
        assert!(read_next_record(&mut reader, &mut pending_header)?.is_none());
        Ok(())
    }

    // -----------------------------------------------------------------------
    // read_indexed_sequence
    // -----------------------------------------------------------------------

    #[test]
    fn test_fai_entry_byte_offset_rejects_zero_line_bases() {
        let entry = FaiEntry {
            length: 8,
            offset: 0,
            line_bases: 0,
            line_width: 8,
        };

        let err = entry.byte_offset(1).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("line_bases"));
    }

    #[test]
    fn test_read_indexed_sequence_reads_range() -> io::Result<()> {
        let mut reader = Cursor::new(b"ACGTTGCA".to_vec());
        let entry = FaiEntry {
            length: 8,
            offset: 0,
            line_bases: 8,
            line_width: 8,
        };

        assert_eq!(read_indexed_sequence(&mut reader, &entry, 2, 5)?, "CGTT");
        Ok(())
    }

    #[test]
    fn test_read_indexed_sequence_buffers_across_line_breaks() -> io::Result<()> {
        let mut reader = CountingReader::new(b"ACGT\nTGCA\n");
        let entry = FaiEntry {
            length: 8,
            offset: 0,
            line_bases: 4,
            line_width: 5,
        };

        assert_eq!(read_indexed_sequence(&mut reader, &entry, 3, 7)?, "GTTGC");
        assert_eq!(reader.seeks, 1);
        assert_eq!(reader.reads, 1);
        Ok(())
    }

    #[test]
    fn test_read_indexed_sequence_handles_crlf_line_breaks() -> io::Result<()> {
        let mut reader = Cursor::new(b"ACGT\r\nTGCA\r\n".to_vec());
        let entry = FaiEntry {
            length: 8,
            offset: 0,
            line_bases: 4,
            line_width: 6,
        };

        assert_eq!(read_indexed_sequence(&mut reader, &entry, 3, 7)?, "GTTGC");
        Ok(())
    }

    #[test]
    fn test_read_indexed_sequence_returns_uppercase() -> io::Result<()> {
        let mut reader = Cursor::new(b"ACgtTGCA".to_vec());
        let entry = FaiEntry {
            length: 8,
            offset: 0,
            line_bases: 8,
            line_width: 8,
        };

        assert_eq!(read_indexed_sequence(&mut reader, &entry, 2, 5)?, "CGTT");
        Ok(())
    }

    #[test]
    fn test_read_indexed_sequence_pads_beyond_sequence_end() -> io::Result<()> {
        let mut reader = Cursor::new(b"ACGTTGCA".to_vec());
        let entry = FaiEntry {
            length: 8,
            offset: 0,
            line_bases: 8,
            line_width: 8,
        };

        assert_eq!(read_indexed_sequence(&mut reader, &entry, 7, 10)?, "CANN");
        Ok(())
    }

    #[test]
    fn test_read_indexed_sequence_pads_when_start_is_after_sequence_end() -> io::Result<()> {
        let mut reader = Cursor::new(b"ACGTTGCA".to_vec());
        let entry = FaiEntry {
            length: 8,
            offset: 0,
            line_bases: 8,
            line_width: 8,
        };

        assert_eq!(read_indexed_sequence(&mut reader, &entry, 9, 11)?, "NNN");
        Ok(())
    }

    #[test]
    fn test_read_indexed_sequence_rejects_invalid_range() {
        let mut reader = Cursor::new(b"ACGTTGCA".to_vec());
        let entry = FaiEntry {
            length: 8,
            offset: 0,
            line_bases: 8,
            line_width: 8,
        };

        let err = read_indexed_sequence(&mut reader, &entry, 0, 1).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);

        let err = read_indexed_sequence(&mut reader, &entry, 5, 2).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
    }

    // -----------------------------------------------------------------------
    // read_gzi_index
    // -----------------------------------------------------------------------

    #[test]
    fn test_read_gzi_index_reads_associated_index() -> io::Result<()> {
        let dir = create_temp_dir("fasta-gzi-index-test")?;
        let path = dir.join("sample.fa.gz");
        fs::write(&path, SAMPLE)?;
        gzi::fs::write(associated_index_path(&path, "gzi"), &gzi::Index::default())?;

        read_gzi_index(&path)?;
        Ok(())
    }

    #[test]
    fn test_read_gzi_index_missing_file_errors() -> io::Result<()> {
        let dir = create_temp_dir("fasta-gzi-missing-test")?;
        let path = dir.join("sample.fa.gz");
        fs::write(&path, SAMPLE)?;

        let err = read_gzi_index(&path).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::NotFound);
        Ok(())
    }

    // -----------------------------------------------------------------------
    // read_fai_index
    // -----------------------------------------------------------------------

    #[test]
    fn test_read_fai_index_reads_associated_index() -> io::Result<()> {
        let (_dir, path) = write_sample_fasta_with_fai("fasta-index-test")?;

        let index = read_fai_index(&path)?;
        assert_eq!(index["chr1"].length, 8);
        assert_eq!(index["chr2"].length, 4);
        Ok(())
    }

    #[test]
    fn test_read_fai_index_skips_empty_lines() -> io::Result<()> {
        let dir = create_temp_dir("fasta-index-empty-lines-test")?;
        let path = dir.join("sample.fa");
        fs::write(&path, SAMPLE)?;
        fs::write(
            associated_index_path(&path, "fai"),
            "\nchr1\t8\t6\t8\t9\n\nchr2\t4\t21\t4\t5\n",
        )?;

        let index = read_fai_index(&path)?;
        assert_eq!(index.len(), 2);
        Ok(())
    }

    #[test]
    fn test_read_fai_index_missing_file_errors() -> io::Result<()> {
        let dir = create_temp_dir("fasta-index-missing-test")?;
        let path = dir.join("sample.fa");
        fs::write(&path, SAMPLE)?;

        let err = read_fai_index(&path).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::NotFound);
        Ok(())
    }

    #[test]
    fn test_read_fai_index_rejects_missing_fields() -> io::Result<()> {
        let dir = create_temp_dir("fasta-index-missing-fields-test")?;
        let path = dir.join("sample.fa");
        fs::write(&path, SAMPLE)?;
        fs::write(associated_index_path(&path, "fai"), "chr1\t8\t6\t8\n")?;

        let err = read_fai_index(&path).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("missing line_width"));
        Ok(())
    }

    #[test]
    fn test_read_fai_index_rejects_invalid_numeric_field() -> io::Result<()> {
        let dir = create_temp_dir("fasta-index-invalid-number-test")?;
        let path = dir.join("sample.fa");
        fs::write(&path, SAMPLE)?;
        fs::write(associated_index_path(&path, "fai"), "chr1\tbad\t6\t8\t9\n")?;

        let err = read_fai_index(&path).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        Ok(())
    }

    // -----------------------------------------------------------------------
    // has_gzi_index and has_fai_index
    // -----------------------------------------------------------------------

    #[test]
    fn test_has_fai_index_and_has_gzi_index() -> io::Result<()> {
        let dir = create_temp_dir("fasta-index-presence-test")?;
        let path = dir.join("sample.fa.gz");
        fs::write(&path, SAMPLE)?;

        assert!(!has_fai_index(&path));
        assert!(!has_gzi_index(&path));

        fs::write(associated_index_path(&path, "fai"), SAMPLE_FAI)?;
        gzi::fs::write(associated_index_path(&path, "gzi"), &gzi::Index::default())?;

        assert!(has_fai_index(&path));
        assert!(has_gzi_index(&path));
        Ok(())
    }

    // -----------------------------------------------------------------------
    // has_bgzf_sequence_index
    // -----------------------------------------------------------------------

    #[test]
    fn test_has_bgzf_sequence_index_requires_fai_and_gzi() -> io::Result<()> {
        let dir = create_temp_dir("fasta-bgzf-index-presence-test")?;
        let path = dir.join("sample.fa.gz");
        fs::write(&path, SAMPLE)?;

        assert!(!has_bgzf_sequence_index(&path));

        fs::write(associated_index_path(&path, "fai"), SAMPLE_FAI)?;
        assert!(!has_bgzf_sequence_index(&path));

        gzi::fs::write(associated_index_path(&path, "gzi"), &gzi::Index::default())?;
        assert!(has_bgzf_sequence_index(&path));
        Ok(())
    }

    // -----------------------------------------------------------------------
    // expected_index_paths
    // -----------------------------------------------------------------------

    #[test]
    fn test_expected_index_paths_for_plain_fasta() -> io::Result<()> {
        let dir = create_temp_dir("fasta-expected-index-plain-test")?;
        let path = dir.join("sample.fa");
        fs::write(&path, SAMPLE)?;

        assert_eq!(
            expected_index_paths(&path),
            vec![associated_index_path(&path, "fai")]
        );
        Ok(())
    }

    #[test]
    fn test_expected_index_paths_for_bgzf_fasta() -> io::Result<()> {
        let dir = create_temp_dir("fasta-expected-index-bgzf-test")?;
        let path = write_bgzf_file(&dir, sample_lines(), Some("sample.fa.gz"))?;

        assert_eq!(
            expected_index_paths(&path),
            vec![
                associated_index_path(&path, "fai"),
                associated_index_path(&path, "gzi"),
            ]
        );
        Ok(())
    }

    // -----------------------------------------------------------------------
    // FastaParser:metadata
    // -----------------------------------------------------------------------

    #[test]
    fn test_fasta_parser_metadata_without_open_reader() -> io::Result<()> {
        let mut parser = FastaParser::new(None, None).map_err(from_py_error)?;

        Python::initialize();
        Python::attach(|py| {
            assert!(parser.metadata(py).unwrap().is_empty());
        });

        Ok(())
    }

    #[test]
    fn test_fasta_parser_metadata_sequential_reports_lengths_only() -> io::Result<()> {
        let dir = create_temp_dir("fasta-parser-metadata-seq")?;
        let path = dir.join("sample.fa");
        fs::write(&path, DESCRIBED_SAMPLE)?;

        let mut parser =
            FastaParser::new(None, Some("sequential".to_string())).map_err(from_py_error)?;
        parser
            .open(&path.to_string_lossy(), None)
            .map_err(from_py_error)?;

        Python::initialize();
        Python::attach(|py| {
            let metadata = parser.metadata(py).unwrap();
            let keys: Vec<String> = metadata.keys().extract().unwrap();
            assert_eq!(keys, ["sequence_lengths"]);
            let lengths: HashMap<String, usize> = metadata
                .get_item("sequence_lengths")
                .unwrap()
                .unwrap()
                .extract()
                .unwrap();
            assert_eq!(
                lengths,
                HashMap::from([("chr1".to_string(), 8), ("chr2".to_string(), 4)])
            );
        });

        Ok(())
    }

    #[test]
    fn test_fasta_parser_metadata_does_not_consume_sequential_reader() -> io::Result<()> {
        let dir = create_temp_dir("fasta-parser-metadata-preserve-seq")?;
        let path = dir.join("sample.fa");
        fs::write(&path, DESCRIBED_SAMPLE)?;

        let mut parser =
            FastaParser::new(None, Some("sequential".to_string())).map_err(from_py_error)?;
        parser
            .open(&path.to_string_lossy(), None)
            .map_err(from_py_error)?;

        let first = parser
            .next()
            .map_err(from_py_error)?
            .expect("expected first record");
        assert_first_record(&first);

        Python::initialize();
        Python::attach(|py| {
            assert!(!parser.metadata(py).unwrap().is_empty());
        });

        let second = parser
            .next()
            .map_err(from_py_error)?
            .expect("expected second record");
        assert_second_record(&second);
        Ok(())
    }

    // -----------------------------------------------------------------------
    // FastaParser:open
    // -----------------------------------------------------------------------

    #[test]
    fn test_fasta_parser_without_open_reader() -> io::Result<()> {
        let mut parser = FastaParser::new(None, None).map_err(from_py_error)?;

        assert_eq!(parser.format(), "fasta");
        assert_eq!(parser.reader_mode(), "auto");
        assert!(parser.next().unwrap_err().to_string().contains("not open"));
        assert!(
            parser
                .seek("chr1:10-20")
                .unwrap_err()
                .to_string()
                .contains("not open")
        );
        assert!(
            parser
                .records(None)
                .unwrap_err()
                .to_string()
                .contains("not open")
        );

        parser.close();
        Ok(())
    }

    #[test]
    fn test_fasta_parser_open_missing_file_errors() -> io::Result<()> {
        Python::initialize();
        let mut parser = FastaParser::new(None, None).map_err(from_py_error)?;
        let err = parser
            .open("/definitely/missing/sample.fa", None)
            .unwrap_err();

        assert!(
            err.to_string().contains("No such file") || err.to_string().contains("os error 2"),
            "unexpected error: {err}"
        );
        Ok(())
    }

    #[test]
    fn test_fasta_parser_open_override_reader_mode() -> io::Result<()> {
        let (_dir, path) = write_sample_fasta_with_fai("fasta-parser-open-override")?;

        let mut parser =
            FastaParser::new(None, Some("sequential".to_string())).map_err(from_py_error)?;
        assert_eq!(parser.reader_mode(), "sequential");

        parser
            .open(&path.to_string_lossy(), Some("query".to_string()))
            .map_err(from_py_error)?;
        assert_eq!(parser.reader_mode(), "query");
        Ok(())
    }

    #[test]
    fn test_fasta_parser_open_file_auto_mode_sequential() -> io::Result<()> {
        let dir = create_temp_dir("fasta-parser-open-auto-seq")?;
        let path = dir.join("sample.fa");
        fs::write(&path, SAMPLE)?;

        let mut parser = FastaParser::new(None, None).map_err(from_py_error)?;
        parser
            .open(&path.to_string_lossy(), None)
            .map_err(from_py_error)?;

        assert_eq!(parser.reader_mode(), "sequential");
        let records = parser.records(None).map_err(from_py_error)?;
        assert_eq!(records.len(), 2);
        Ok(())
    }

    #[test]
    fn test_fasta_parser_open_file_auto_mode_query() -> io::Result<()> {
        let (_dir, path) = write_sample_fasta_with_fai("fasta-parser-open-auto-query")?;

        let mut parser = FastaParser::new(None, None).map_err(from_py_error)?;
        parser
            .open(&path.to_string_lossy(), None)
            .map_err(from_py_error)?;

        assert_eq!(parser.reader_mode(), "query");
        assert!(parser.seek("chr1:2-5").map_err(from_py_error)?);
        let record = parser
            .next()
            .map_err(from_py_error)?
            .expect("expected query result");
        assert_eq!(record.sequence, "CGTT");
        Ok(())
    }

    #[test]
    fn test_fasta_parser_open_file_auto_mode_query_then_sequential() -> io::Result<()> {
        let (_dir, indexed_path) = write_sample_fasta_with_fai("fasta-parser-open-query-then-seq")?;
        let plain_dir = create_temp_dir("fasta-parser-open-query-then-seq-plain")?;
        let plain_path = plain_dir.join("plain.fa");
        fs::write(&plain_path, SAMPLE)?;

        let mut parser = FastaParser::new(None, None).map_err(from_py_error)?;
        parser
            .open(&indexed_path.to_string_lossy(), None)
            .map_err(from_py_error)?;
        assert_eq!(parser.reader_mode(), "query");

        parser
            .open(&plain_path.to_string_lossy(), None)
            .map_err(from_py_error)?;
        assert_eq!(parser.reader_mode(), "sequential");
        Ok(())
    }

    #[test]
    fn test_fasta_parser_open_file_query_mode_without_index_errors() -> io::Result<()> {
        let dir = create_temp_dir("fasta-parser-open-query-no-index")?;
        let path = dir.join("sample.fa");
        fs::write(&path, SAMPLE)?;

        Python::initialize();
        let mut parser =
            FastaParser::new(None, Some("query".to_string())).map_err(from_py_error)?;
        let err = parser
            .open(&path.to_string_lossy(), Some("query".to_string()))
            .unwrap_err();
        assert!(
            err.to_string().contains("no FASTA index found"),
            "unexpected error: {err}"
        );
        Ok(())
    }

    #[test]
    fn test_fasta_parser_open_file_query_mode_with_index() -> io::Result<()> {
        let (_dir, path) = write_sample_fasta_with_fai("fasta-parser-open-query-index")?;

        let mut parser =
            FastaParser::new(None, Some("query".to_string())).map_err(from_py_error)?;
        parser
            .open(&path.to_string_lossy(), None)
            .map_err(from_py_error)?;

        assert_eq!(parser.reader_mode(), "query");
        assert!(parser.seek("chr1:2-5").map_err(from_py_error)?);
        let records = parser.records(None).map_err(from_py_error)?;
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].sequence, "CGTT");
        Ok(())
    }

    #[test]
    fn test_fasta_parser_open_file_sequential_mode_with_index() -> io::Result<()> {
        let (_dir, path) = write_sample_fasta_with_fai("fasta-parser-open-seq-index")?;

        let mut parser =
            FastaParser::new(None, Some("sequential".to_string())).map_err(from_py_error)?;
        parser
            .open(&path.to_string_lossy(), Some("sequential".to_string()))
            .map_err(from_py_error)?;

        assert_eq!(parser.reader_mode(), "sequential");
        let records = parser.records(None).map_err(from_py_error)?;
        assert_eq!(records.len(), 2);
        Ok(())
    }

    #[test]
    fn test_fasta_parser_rejects_invalid_mode() {
        Python::initialize();
        let err = match FastaParser::new(None, Some("random".to_string())) {
            Ok(_) => panic!("expected invalid reader mode to error"),
            Err(err) => err,
        };

        assert!(
            err.to_string()
                .contains("unsupported reader mode \"random\""),
            "unexpected error: {err}"
        );
    }

    // -----------------------------------------------------------------------
    // FastaParser:close
    // -----------------------------------------------------------------------

    #[test]
    fn test_fasta_parser_close() -> io::Result<()> {
        let dir = create_temp_dir("fasta-parser-close")?;
        let path = dir.join("sample.fa");
        fs::write(&path, SAMPLE)?;

        Python::initialize();
        let mut parser = FastaParser::new(None, None).map_err(from_py_error)?;
        parser
            .open(&path.to_string_lossy(), None)
            .map_err(from_py_error)?;

        assert!(parser.reader.is_some());
        parser.close();
        assert!(parser.reader.is_none());
        assert!(parser.next().unwrap_err().to_string().contains("not open"));
        assert!(
            parser
                .seek("chr1:1-4")
                .unwrap_err()
                .to_string()
                .contains("not open")
        );
        assert!(
            parser
                .records(None)
                .unwrap_err()
                .to_string()
                .contains("not open")
        );
        Ok(())
    }

    // -----------------------------------------------------------------------
    // FastaParser:parse_text and parse_file
    // -----------------------------------------------------------------------

    #[test]
    fn test_fasta_parser_parse_text() -> io::Result<()> {
        let records = FastaParser::parse_text(str::from_utf8(SAMPLE).unwrap(), Some(1))
            .map_err(from_py_error)?;
        assert_eq!(records.len(), 1);
        assert_first_record(&records[0]);
        Ok(())
    }

    #[test]
    fn test_fasta_parser_parse_file() -> io::Result<()> {
        let dir = create_temp_dir("fasta-parser-parse-file")?;
        let path = dir.join("sample.fa");
        fs::write(&path, SAMPLE)?;

        let records =
            FastaParser::parse_file(&path.to_string_lossy(), None).map_err(from_py_error)?;
        assert_eq!(records.len(), 2);
        assert_first_record(&records[0]);
        assert_second_record(&records[1]);
        Ok(())
    }

    #[test]
    fn test_fasta_parser_parse_file_bgzf() -> io::Result<()> {
        let dir = create_temp_dir("fasta-parser-parse-file-bgzf")?;
        let path = write_bgzf_file(&dir, sample_lines(), Some("sample.fa.gz"))?;

        let records =
            FastaParser::parse_file(&path.to_string_lossy(), None).map_err(from_py_error)?;
        assert_eq!(records.len(), 2);
        assert_first_record(&records[0]);
        assert_second_record(&records[1]);
        Ok(())
    }

    #[test]
    fn test_fasta_parser_parse_file_missing_path_errors() {
        Python::initialize();
        let err = FastaParser::parse_file("/definitely/missing/sample.fa", None).unwrap_err();
        assert!(
            err.to_string().contains("No such file") || err.to_string().contains("os error 2"),
            "unexpected error: {err}"
        );
    }

    // -----------------------------------------------------------------------
    // BgzfBackend
    // -----------------------------------------------------------------------

    #[test]
    fn test_bgzf_backend_next_without_query_returns_none() -> io::Result<()> {
        let (_dir, path) = write_sample_bgzf_with_indexes("fasta-bgzf-next-empty-test")?;
        let mut backend = BgzfBackend::new(&path)?;
        assert!(backend.next()?.is_none());
        Ok(())
    }

    #[test]
    fn test_bgzf_backend_seek_returns_requested_range() -> io::Result<()> {
        let (_dir, path) = write_sample_bgzf_with_indexes("fasta-bgzf-range-test")?;

        let mut backend = BgzfBackend::new(&path)?;
        assert!(backend.seek("chr1:2-5")?);
        let record = backend.next()?.expect("indexed record should be available");
        assert_eq!(record.sequence, "CGTT");
        Ok(())
    }

    #[test]
    fn test_bgzf_backend_seek_uses_sequence_length_when_end_is_omitted() -> io::Result<()> {
        let (_dir, path) = write_sample_bgzf_with_indexes("fasta-bgzf-open-ended-test")?;

        let mut backend = BgzfBackend::new(&path)?;
        assert!(backend.seek("chr1:5")?);
        let record = backend.next()?.expect("indexed record should be available");
        assert_eq!(record.start, 5);
        assert_eq!(record.end, 8);
        assert_eq!(record.sequence, "TGCA");
        Ok(())
    }

    #[test]
    fn test_bgzf_backend_seek_returns_whole_sequence_for_seq_region_only() -> io::Result<()> {
        let (_dir, path) = write_sample_bgzf_with_indexes("fasta-bgzf-whole-test")?;

        let mut backend = BgzfBackend::new(&path)?;
        assert!(backend.seek("chr2")?);
        let record = backend.next()?.expect("indexed record should be available");
        assert_second_record(&record);
        Ok(())
    }

    #[test]
    fn test_bgzf_backend_seek_returns_false_for_missing_sequence() -> io::Result<()> {
        let (_dir, path) = write_sample_bgzf_with_indexes("fasta-bgzf-missing-test")?;

        let mut backend = BgzfBackend::new(&path)?;
        assert!(!backend.seek("chr3")?);
        assert!(backend.next()?.is_none());
        Ok(())
    }

    #[test]
    fn test_bgzf_backend_seek_followed_by_failed_seek_returns_none() -> io::Result<()> {
        let (_dir, path) = write_sample_bgzf_with_indexes("fasta-bgzf-seek-failed-seek-test")?;

        let mut backend = BgzfBackend::new(&path)?;
        assert!(backend.seek("chr1:2-5")?);
        assert!(!backend.seek("chr3")?);
        assert!(backend.next()?.is_none());
        Ok(())
    }

    // -----------------------------------------------------------------------
    // FaiBackend
    // -----------------------------------------------------------------------

    #[test]
    fn test_fai_backend_next_without_query_returns_none() -> io::Result<()> {
        let (_dir, path) = write_sample_fasta_with_fai("fasta-fai-next-empty-test")?;
        let mut backend = FaiBackend::new(&path)?;
        assert!(backend.next()?.is_none());
        Ok(())
    }

    #[test]
    fn test_fai_backend_seek_returns_requested_range() -> io::Result<()> {
        let (_dir, path) = write_sample_fasta_with_fai("fasta-fai-range-test")?;

        let mut backend = FaiBackend::new(&path)?;
        assert!(backend.seek("chr1:2-5")?);
        let record = backend.next()?.expect("indexed record should be available");
        assert_eq!(record.seqid, "chr1");
        assert_eq!(record.start, 2);
        assert_eq!(record.end, 5);
        assert_eq!(record.sequence, "CGTT");
        assert!(backend.next()?.is_none());
        Ok(())
    }

    #[test]
    fn test_fai_backend_seek_uses_sequence_length_when_end_is_omitted() -> io::Result<()> {
        let (_dir, path) = write_sample_fasta_with_fai("fasta-fai-open-ended-test")?;

        let mut backend = FaiBackend::new(&path)?;
        assert!(backend.seek("chr1:5")?);
        let record = backend.next()?.expect("indexed record should be available");
        assert_eq!(record.start, 5);
        assert_eq!(record.end, 8);
        assert_eq!(record.sequence, "TGCA");
        Ok(())
    }

    #[test]
    fn test_fai_backend_seek_returns_whole_sequence_for_seq_region_only() -> io::Result<()> {
        let (_dir, path) = write_sample_fasta_with_fai("fasta-fai-whole-test")?;

        let mut backend = FaiBackend::new(&path)?;
        assert!(backend.seek("chr2")?);
        let record = backend.next()?.expect("indexed record should be available");
        assert_second_record(&record);
        Ok(())
    }

    #[test]
    fn test_fai_backend_seek_pads_beyond_sequence_end() -> io::Result<()> {
        let (_dir, path) = write_sample_fasta_with_fai("fasta-fai-padding-test")?;

        let mut backend = FaiBackend::new(&path)?;
        assert!(backend.seek("chr1:7-10")?);
        let record = backend.next()?.expect("indexed record should be available");
        assert_eq!(record.sequence, "CANN");
        Ok(())
    }

    #[test]
    fn test_fai_backend_seek_returns_false_for_missing_sequence() -> io::Result<()> {
        let (_dir, path) = write_sample_fasta_with_fai("fasta-fai-missing-test")?;

        let mut backend = FaiBackend::new(&path)?;
        assert!(!backend.seek("chr3")?);
        assert!(backend.next()?.is_none());
        Ok(())
    }

    #[test]
    fn test_fai_backend_seek_followed_by_failed_seek_returns_none() -> io::Result<()> {
        let (_dir, path) = write_sample_fasta_with_fai("fasta-fai-seek-failed-seek-test")?;

        let mut backend = FaiBackend::new(&path)?;
        assert!(backend.seek("chr1:2-5")?);
        assert!(!backend.seek("chr3")?);
        assert!(backend.next()?.is_none());
        Ok(())
    }

    // -----------------------------------------------------------------------
    // SequentialBackend
    // -----------------------------------------------------------------------

    #[test]
    fn test_sequential_backend_next_reads_until_eof() -> io::Result<()> {
        let dir = create_temp_dir("fasta-sequential-backend-test")?;
        let path = dir.join("sample.fa");
        fs::write(&path, SAMPLE)?;

        let mut backend = SequentialBackend::new(&path)?;
        let first = backend.next()?.expect("first record should be available");
        let second = backend.next()?.expect("second record should be available");

        assert_first_record(&first);
        assert_second_record(&second);
        assert!(backend.next()?.is_none());
        Ok(())
    }

    #[test]
    fn test_sequential_backend_next_reads_bgzf_until_eof() -> io::Result<()> {
        let dir = create_temp_dir("fasta-sequential-bgzf-test")?;
        let path = write_bgzf_file(&dir, sample_lines(), Some("sample.fa.gz"))?;

        let mut backend = SequentialBackend::new(&path)?;
        let first = backend.next()?.expect("first record should be available");
        let second = backend.next()?.expect("second record should be available");

        assert_first_record(&first);
        assert_second_record(&second);
        assert!(backend.next()?.is_none());
        Ok(())
    }

    fn write_sequential_fasta(prefix: &str, data: &[u8]) -> io::Result<PathBuf> {
        let dir = create_temp_dir(prefix)?;
        let path = dir.join("sample.fa");
        fs::write(&path, data)?;
        Ok(path)
    }

    fn sequential_seek(
        backend: &mut SequentialBackend,
        region: &str,
    ) -> io::Result<Option<String>> {
        if !backend.seek(region)? {
            return Ok(None);
        }
        Ok(backend.next()?.map(|record| record.sequence))
    }

    const MULTILINE_SAMPLE: &[u8] = b"\
>chr1 first
ACGT
tgca
AC
>chr2
GGGG
CC
";

    #[test]
    fn test_sequential_backend_seek_reads_regular_sequence_directly() -> io::Result<()> {
        let path = write_sequential_fasta("fasta-sequential-seek-direct", MULTILINE_SAMPLE)?;
        let mut backend = SequentialBackend::new(&path)?;

        assert_eq!(
            sequential_seek(&mut backend, "chr1:3-7")?.as_deref(),
            Some("GTTGC")
        );
        assert!(backend.scanned_sequences()?["chr1"].fai_entry.is_some());
        assert!(backend.file.is_some());
        assert_eq!(
            sequential_seek(&mut backend, "chr2:5-8")?.as_deref(),
            Some("CCNN")
        );
        assert_eq!(
            sequential_seek(&mut backend, "chr1")?.as_deref(),
            Some("ACGTTGCAAC")
        );
        assert_eq!(backend.next()?, None);
        Ok(())
    }

    #[test]
    fn test_sequential_backend_seek_returns_false_for_missing_sequence() -> io::Result<()> {
        let path = write_sequential_fasta("fasta-sequential-seek-missing", SAMPLE)?;
        let mut backend = SequentialBackend::new(&path)?;

        assert_eq!(
            sequential_seek(&mut backend, "chr1:1-4")?.as_deref(),
            Some("ACGT")
        );
        assert!(!backend.seek("chr3:1-4")?);
        assert_eq!(backend.next()?, None);
        Ok(())
    }

    #[test]
    fn test_sequential_backend_seek_rejects_start_after_end() -> io::Result<()> {
        let path = write_sequential_fasta("fasta-sequential-seek-invalid", SAMPLE)?;
        let mut backend = SequentialBackend::new(&path)?;

        assert!(backend.seek("chr1:5-4").is_err());
        Ok(())
    }

    #[test]
    fn test_sequential_backend_seek_rescans_irregular_sequence() -> io::Result<()> {
        let path = write_sequential_fasta(
            "fasta-sequential-seek-irregular",
            b">chr1\nAC\nGTTG\nCA\n>chr2\nAA AA\n",
        )?;
        let mut backend = SequentialBackend::new(&path)?;

        assert_eq!(
            sequential_seek(&mut backend, "chr1:2-7")?.as_deref(),
            Some("CGTTGC")
        );
        assert_eq!(
            sequential_seek(&mut backend, "chr2")?.as_deref(),
            Some("AAAA")
        );
        assert!(backend.scanned_sequences()?["chr1"].fai_entry.is_none());
        assert!(backend.scanned_sequences()?["chr2"].fai_entry.is_none());
        assert!(backend.file.is_none());
        Ok(())
    }

    #[test]
    fn test_sequential_backend_seek_rescans_bgzf() -> io::Result<()> {
        let dir = create_temp_dir("fasta-sequential-seek-bgzf")?;
        let path = write_bgzf_file(&dir, sample_lines(), Some("sample.fa.gz"))?;
        let mut backend = SequentialBackend::new(&path)?;

        assert_eq!(
            sequential_seek(&mut backend, "chr1:3-6")?.as_deref(),
            Some("GTTG")
        );
        assert_eq!(
            sequential_seek(&mut backend, "chr2:3-6")?.as_deref(),
            Some("AANN")
        );
        assert!(backend.scanned_sequences()?["chr1"].fai_entry.is_none());
        Ok(())
    }

    #[test]
    fn test_sequential_backend_seek_uses_first_duplicate_record() -> io::Result<()> {
        let path = write_sequential_fasta(
            "fasta-sequential-seek-duplicate",
            b">chr1\nAC\n>chr1\nGGGG\n",
        )?;
        let mut backend = SequentialBackend::new(&path)?;

        assert_eq!(
            sequential_seek(&mut backend, "chr1")?.as_deref(),
            Some("AC")
        );
        assert_eq!(backend.scanned_sequences()?["chr1"].length, 2);
        Ok(())
    }

    #[test]
    fn test_sequential_backend_seek_matches_fai_backend() -> io::Result<()> {
        let data = b">chr1\nACGTACGTAC\nGTACGTACGT\nACG\n>chr2 described\nTTTTTGGGGG\nCC\n";
        let path = write_sequential_fasta("fasta-sequential-seek-parity", data)?;
        fs::write(
            associated_index_path(&path, "fai"),
            "chr1\t23\t6\t10\t11\nchr2\t12\t48\t10\t11\n",
        )?;
        let mut sequential_backend = SequentialBackend::new(&path)?;
        let mut fai_backend = FaiBackend::new(&path)?;

        for (seq_region_name, length) in [("chr1", 23), ("chr2", 12)] {
            for start in 1..=length + 2 {
                for end in start..=length + 2 {
                    let region = format!("{seq_region_name}:{start}-{end}");
                    assert!(fai_backend.seek(&region)?);
                    let expected = fai_backend.next()?.map(|record| record.sequence);
                    assert_eq!(
                        sequential_seek(&mut sequential_backend, &region)?,
                        expected,
                        "{region}"
                    );
                }
            }
        }
        Ok(())
    }

    #[test]
    fn test_sequential_backend_next_streams_until_first_seek() -> io::Result<()> {
        let path = write_sequential_fasta("fasta-sequential-stream-then-seek", SAMPLE)?;
        let mut backend = SequentialBackend::new(&path)?;

        assert_first_record(&backend.next()?.expect("first record should be available"));
        assert_eq!(
            sequential_seek(&mut backend, "chr1:1-2")?.as_deref(),
            Some("AC")
        );
        assert_eq!(backend.next()?, None);
        Ok(())
    }

    // -----------------------------------------------------------------------
    // scan_sequences / scan_sequence_chunk
    // -----------------------------------------------------------------------

    fn scanned(data: &[u8]) -> io::Result<HashMap<String, ScannedSequence>> {
        scan_sequences(&mut Cursor::new(data), true)
    }

    #[test]
    fn test_scan_sequences_records_fai_layout() -> io::Result<()> {
        let sequences = scanned(MULTILINE_SAMPLE)?;

        assert_eq!(
            sequences["chr1"].fai_entry,
            Some(FaiEntry {
                length: 10,
                offset: 12,
                line_bases: 4,
                line_width: 5,
            })
        );
        assert_eq!(sequences["chr2"].length, 6);
        Ok(())
    }

    #[test]
    fn test_scan_sequences_accepts_crlf_and_missing_final_newline() -> io::Result<()> {
        let sequences = scanned(b">chr1\r\nACGT\r\nAC")?;

        assert_eq!(
            sequences["chr1"].fai_entry,
            Some(FaiEntry {
                length: 6,
                offset: 7,
                line_bases: 4,
                line_width: 6,
            })
        );
        Ok(())
    }

    #[test]
    fn test_scan_sequences_allows_trailing_blank_lines() -> io::Result<()> {
        let sequences = scanned(b">chr1\nACGT\nAC\n\n>chr2\nA\n")?;

        assert!(sequences["chr1"].fai_entry.is_some());
        assert_eq!(sequences["chr1"].length, 6);
        Ok(())
    }

    #[test]
    fn test_scan_sequences_marks_irregular_layouts() -> io::Result<()> {
        let sequences = scanned(
            b">long_middle\nAC\nACGT\n>short_middle\nACGT\nAC\nACGT\n>blank_middle\nAC\n\nAC\n\
              >inner_space\nAC GT\n>leading_space\n ACGT\n",
        )?;

        for seq_region_name in [
            "long_middle",
            "short_middle",
            "blank_middle",
            "inner_space",
            "leading_space",
        ] {
            assert!(
                sequences[seq_region_name].fai_entry.is_none(),
                "{seq_region_name}"
            );
        }
        assert_eq!(sequences["inner_space"].length, 4);
        assert_eq!(sequences["leading_space"].length, 4);
        Ok(())
    }

    #[test]
    fn test_scan_sequences_without_offsets_has_no_layout() -> io::Result<()> {
        let sequences = scan_sequences(&mut Cursor::new(SAMPLE), false)?;

        assert_eq!(sequences["chr1"].length, 8);
        assert!(sequences["chr1"].fai_entry.is_none());
        Ok(())
    }

    #[test]
    fn test_scan_sequences_names_and_ignores_sequence_before_header() -> io::Result<()> {
        let sequences = scanned(b"ACGT\n>  chr1 description\nAC\n>chr1\nGGGG\n>empty\n")?;

        assert_eq!(sequences.len(), 2);
        assert_eq!(sequences["chr1"].length, 2);
        assert_eq!(sequences["empty"].length, 0);
        Ok(())
    }

    #[test]
    fn test_scan_sequence_chunk_reads_across_lines_and_pads() -> io::Result<()> {
        let chunk = |start, end| {
            scan_sequence_chunk(&mut Cursor::new(MULTILINE_SAMPLE), "chr1", start, end)
        };

        assert_eq!(chunk(1, 10)?, "ACGTTGCAAC");
        assert_eq!(chunk(4, 5)?, "TT");
        assert_eq!(chunk(9, 12)?, "ACNN");
        assert_eq!(chunk(20, 21)?, "NN");
        Ok(())
    }

    #[test]
    fn test_scan_sequence_chunk_missing_sequence_is_all_n() -> io::Result<()> {
        let chunk = scan_sequence_chunk(&mut Cursor::new(SAMPLE), "chr3", 1, 3)?;

        assert_eq!(chunk, "NNN");
        Ok(())
    }

    #[test]
    fn test_scan_sequence_chunk_rejects_invalid_range() {
        assert!(scan_sequence_chunk(&mut Cursor::new(SAMPLE), "chr1", 0, 3).is_err());
        assert!(scan_sequence_chunk(&mut Cursor::new(SAMPLE), "chr1", 4, 3).is_err());
    }

    // -----------------------------------------------------------------------
    // Backend
    // -----------------------------------------------------------------------

    #[test]
    fn test_backend_auto_uses_fai_for_plain_indexed_fasta() -> io::Result<()> {
        let (_dir, path) = write_sample_fasta_with_fai("fasta-backend-auto-fai-test")?;

        let (backend, mode) = Backend::open(&path, ReaderMode::Auto)?;
        assert!(matches!(backend, Backend::Fai(_)));
        assert_eq!(mode, ReaderMode::Query);
        Ok(())
    }

    #[test]
    fn test_backend_auto_uses_sequential_for_plain_unindexed_fasta() -> io::Result<()> {
        let dir = create_temp_dir("fasta-backend-auto-sequential-test")?;
        let path = dir.join("sample.fa");
        fs::write(&path, SAMPLE)?;

        let (backend, mode) = Backend::open(&path, ReaderMode::Auto)?;
        assert!(matches!(backend, Backend::Sequential(_)));
        assert_eq!(mode, ReaderMode::Sequential);
        Ok(())
    }

    #[test]
    fn test_backend_auto_uses_bgzf_when_fai_and_gzi_are_available() -> io::Result<()> {
        let (_dir, path) = write_sample_bgzf_with_indexes("fasta-backend-auto-bgzf-test")?;

        let (backend, mode) = Backend::open(&path, ReaderMode::Auto)?;
        assert!(matches!(backend, Backend::Bgzf(_)));
        assert_eq!(mode, ReaderMode::Query);
        Ok(())
    }

    #[test]
    fn test_backend_auto_uses_sequential_for_bgzf_without_indexes() -> io::Result<()> {
        let dir = create_temp_dir("fasta-backend-auto-bgzf-seq-test")?;
        let path = write_bgzf_file(&dir, sample_lines(), Some("sample.fa.gz"))?;

        let (backend, mode) = Backend::open(&path, ReaderMode::Auto)?;
        assert!(matches!(backend, Backend::Sequential(_)));
        assert_eq!(mode, ReaderMode::Sequential);
        Ok(())
    }

    #[test]
    fn test_backend_query_uses_fai_for_plain_indexed_fasta() -> io::Result<()> {
        let (_dir, path) = write_sample_fasta_with_fai("fasta-backend-query-fai-test")?;

        let (backend, mode) = Backend::open(&path, ReaderMode::Query)?;
        assert!(matches!(backend, Backend::Fai(_)));
        assert_eq!(mode, ReaderMode::Query);
        Ok(())
    }

    #[test]
    fn test_backend_query_uses_bgzf_when_fai_and_gzi_are_available() -> io::Result<()> {
        let (_dir, path) = write_sample_bgzf_with_indexes("fasta-backend-query-bgzf-test")?;

        let (backend, mode) = Backend::open(&path, ReaderMode::Query)?;
        assert!(matches!(backend, Backend::Bgzf(_)));
        assert_eq!(mode, ReaderMode::Query);
        Ok(())
    }

    #[test]
    fn test_backend_query_requires_plain_fasta_index() -> io::Result<()> {
        let dir = create_temp_dir("fasta-backend-query-missing-index-test")?;
        let path = dir.join("sample.fa");
        fs::write(&path, SAMPLE)?;

        let err = match Backend::open(&path, ReaderMode::Query) {
            Ok(_) => panic!("query mode should require a FASTA index"),
            Err(err) => err,
        };
        assert_eq!(err.kind(), io::ErrorKind::NotFound);
        assert!(err.to_string().contains("no FASTA index found"));
        Ok(())
    }

    #[test]
    fn test_backend_query_requires_bgzf_fai_and_gzi_indexes() -> io::Result<()> {
        let dir = create_temp_dir("fasta-backend-query-bgzf-missing-index-test")?;
        let path = write_bgzf_file(&dir, sample_lines(), Some("sample.fa.gz"))?;
        fs::write(associated_index_path(&path, "fai"), SAMPLE_FAI)?;

        let err = match Backend::open(&path, ReaderMode::Query) {
            Ok(_) => panic!("query mode should require FASTA and GZI indexes"),
            Err(err) => err,
        };
        assert_eq!(err.kind(), io::ErrorKind::NotFound);
        assert!(err.to_string().contains("no FASTA index found"));
        assert!(err.to_string().contains("gzi"));
        Ok(())
    }

    #[test]
    fn test_backend_sequential_mode_uses_sequential_backend_with_index() -> io::Result<()> {
        let (_dir, path) = write_sample_fasta_with_fai("fasta-backend-seq-index-test")?;

        let (backend, mode) = Backend::open(&path, ReaderMode::Sequential)?;
        assert!(matches!(backend, Backend::Sequential(_)));
        assert_eq!(mode, ReaderMode::Sequential);
        Ok(())
    }

    #[test]
    fn test_backend_next_and_seek_dispatch_to_fai_backend() -> io::Result<()> {
        let (_dir, path) = write_sample_fasta_with_fai("fasta-backend-dispatch-fai-test")?;
        let mut backend = Backend::Fai(FaiBackend::new(&path)?);

        assert!(backend.seek("chr1:2-5")?);
        let record = backend.next()?.expect("expected dispatched query result");
        assert_eq!(record.sequence, "CGTT");
        Ok(())
    }

    #[test]
    fn test_backend_next_and_seek_dispatch_to_bgzf_backend() -> io::Result<()> {
        let (_dir, path) = write_sample_bgzf_with_indexes("fasta-backend-dispatch-bgzf-test")?;
        let mut backend = Backend::Bgzf(BgzfBackend::new(&path)?);

        assert!(backend.seek("chr1:2-5")?);
        let record = backend.next()?.expect("expected dispatched query result");
        assert_eq!(record.sequence, "CGTT");
        Ok(())
    }

    #[test]
    fn test_backend_next_and_seek_dispatch_to_sequential_backend() -> io::Result<()> {
        let dir = create_temp_dir("fasta-backend-dispatch-seq-test")?;
        let path = dir.join("sample.fa");
        fs::write(&path, SAMPLE)?;
        let mut backend = Backend::Sequential(SequentialBackend::new(&path)?);

        let record = backend.next()?.expect("expected sequential record");
        assert_first_record(&record);
        assert!(backend.seek("chr1:1-4")?);
        assert_eq!(
            backend.next()?.map(|record| record.sequence).as_deref(),
            Some("ACGT")
        );
        assert_eq!(
            backend.metadata()?.sequence_lengths,
            HashMap::from([("chr1".to_string(), 8), ("chr2".to_string(), 4)])
        );
        Ok(())
    }

    // -----------------------------------------------------------------------
    // FastaReader
    // -----------------------------------------------------------------------

    #[test]
    fn test_fasta_reader_query_metadata_is_cached_from_fai() -> io::Result<()> {
        let (_dir, path) = write_sample_fasta_with_fai("fasta-reader-metadata-test")?;
        let mut reader = FastaReader::open(&path, ReaderMode::Query)?;
        let metadata = reader.metadata()?;
        assert_eq!(metadata.sequence_lengths["chr1"], 8);
        assert_eq!(metadata.sequence_lengths["chr2"], 4);
        assert!(reader.metadata.is_some());
        Ok(())
    }

    #[test]
    fn test_fasta_reader_next_and_seek_delegate_to_backend() -> io::Result<()> {
        let (_dir, path) = write_sample_fasta_with_fai("fasta-reader-query-test")?;
        let mut reader = FastaReader::open(&path, ReaderMode::Query)?;

        assert!(reader.seek("chr1:2-5")?);
        let record = reader.next()?.expect("expected query result");
        assert_eq!(record.sequence, "CGTT");
        assert!(reader.next()?.is_none());
        Ok(())
    }

    // -----------------------------------------------------------------------
    // FastaRecord
    // -----------------------------------------------------------------------

    #[test]
    fn test_fasta_record_repr() {
        let record = FastaRecord {
            seqid: "chr1".to_string(),
            start: 2,
            end: 5,
            sequence: "CGTT".to_string(),
        };

        assert_eq!(
            record.__repr__(),
            "FastaRecord(seqid=\"chr1\", start=2, end=5, sequence_length=4)"
        );
    }
}
