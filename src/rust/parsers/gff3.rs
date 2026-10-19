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

use crate::{
    parsers::{
        tabix::{self, has_tabix_index},
        utils::{
            FeatureTypeFilter, ReaderMode, associated_index_path, noodles_phase_to_u8,
            noodles_strand_to_string, open_sequential_reader,
        },
    },
    utils::{bytes_to_string, to_py_error},
};
use noodles_gff::{self as gff3, record::attributes::field::Value as AttributeValue};
use pyo3::{
    exceptions::PyRuntimeError,
    prelude::*,
    types::{PyDict, PyTuple},
};
use std::{
    io::{self, BufRead},
    path::{Path, PathBuf},
};

type AttributePairs = Vec<(String, Vec<String>)>;

// ---------------------------------------------------------------------------
// GFF3Record - represents a GFF3 record
// ---------------------------------------------------------------------------

/// Represents a record from one line of a GFF3 file.
#[pyclass(name = "GFF3Record")]
#[derive(Clone, Debug, PartialEq)]
pub struct Gff3Record {
    /// str: Name of the seq region
    #[pyo3(get)]
    pub seqid: String,
    /// str: Source of the annotation
    #[pyo3(get)]
    pub source: String,
    /// str: Class/type of feature (e.g. gene, transcript, exon)
    #[pyo3(get)]
    pub feature_type: String,
    /// int: Start coordinate of the feature
    #[pyo3(get)]
    pub start: usize,
    /// int: End coordinate of the feature
    #[pyo3(get)]
    pub end: usize,
    /// float | None: Score given to the feature
    #[pyo3(get)]
    pub score: Option<f64>,
    /// str: Symbol representing the strand the feature is on (+, -, . or ?)
    #[pyo3(get)]
    pub strand_symbol: String,
    /// int | None: Value for the phase of the feature (0, 1 or 2)
    #[pyo3(get)]
    pub phase_value: Option<u8>,
    attributes: AttributePairs,
    #[pyo3(get)]
    pub raw_attributes: String,
}

#[pymethods]
impl Gff3Record {
    /// dict[list]: The information from the attributes column
    #[getter]
    fn attributes<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);

        for (key, values) in &self.attributes {
            dict.set_item(key, PyTuple::new(py, values)?)?;
        }

        Ok(dict)
    }

    /// Return the first value for an attribute key, or `None` when the key is absent.
    ///
    /// Equivalent to `attributes[key][0]`, but a linear scan of the stored pairs instead of
    /// building a whole `dict` per call.
    ///
    /// Unlike GTF, GFF3 lets a key repeat as separate `key=value` pairs rather than one key
    /// with a comma-separated array, and noodles-gff keeps each occurrence as its own pair in
    /// `self.attributes`. The `attributes` getter above builds a `dict` via `set_item`, so a
    /// repeated key resolves there to its *last* occurrence -- scanning in reverse here matches
    /// that, instead of silently returning a different value for the same key.
    fn first_attribute(&self, key: &str) -> Option<&str> {
        self.attributes
            .iter()
            .rev()
            .find(|(pair_key, _)| pair_key == key)
            .and_then(|(_, values)| values.first())
            .map(String::as_str)
    }

    fn __repr__(&self) -> String {
        format!(
            "GFF3Record(seqid={:?}, feature_type={:?}, start={}, end={}, strand_symbol={:?}, phase_value={:?})",
            self.seqid,
            self.feature_type,
            self.start,
            self.end,
            self.strand_symbol,
            self.phase_value
        )
    }
}

// ---------------------------------------------------------------------------
// Gff3Metadata - holds GFF3 metadata (directives and header lines)
// ---------------------------------------------------------------------------

#[derive(Clone, Debug, Default, PartialEq)]
struct Gff3Metadata {
    header: Vec<String>,
    directives: Vec<(String, String)>,
}

// ---------------------------------------------------------------------------
// Gff3Reader - public entry point
// ---------------------------------------------------------------------------

// A reusable GFF3 reader that supports both sequential and indexed (tabix) reading modes.

// The reader by default uses indexed reading if a tabix index is available,
// but can be forced into sequential mode via the `mode` parameter on `Gff3Parser`.
pub struct Gff3Reader {
    file_path: PathBuf,
    backend: Backend,
    mode: ReaderMode,
    metadata: Option<Gff3Metadata>,
}

impl Gff3Reader {
    fn open(file_path: &Path, tabix_path: &str, mode: ReaderMode) -> io::Result<Self> {
        let (backend, mode) = Backend::open(file_path, tabix_path, mode)?;
        Ok(Self {
            file_path: file_path.to_path_buf(),
            backend,
            mode,
            metadata: None,
        })
    }

    fn next(&mut self) -> io::Result<Option<Gff3Record>> {
        self.backend.next()
    }

    /// Read the next record whose column-3 feature type is in `filter`, skipping non-matching
    /// lines before they are parsed into a record.
    fn next_filtered(&mut self, filter: &FeatureTypeFilter) -> io::Result<Option<Gff3Record>> {
        self.backend.next_filtered(filter)
    }

    fn seek(&mut self, region: &str) -> io::Result<bool> {
        self.backend.seek(region)
    }

    fn index_sequence_names(&mut self) -> io::Result<Vec<String>> {
        self.backend.index_sequence_names()
    }

    fn metadata(&mut self) -> io::Result<&Gff3Metadata> {
        if self.metadata.is_none() {
            self.metadata = Some(read_gff3_metadata(&self.file_path)?);
        }

        Ok(self.metadata.as_ref().expect("metadata should be cached"))
    }
}

// ---------------------------------------------------------------------------
// Backend
// ---------------------------------------------------------------------------

enum Backend {
    Sequential(SequentialBackend),
    Tabix(Box<TabixBackend>),
}

impl Backend {
    fn open(
        file_path: &Path,
        tabix_path: &str,
        mode: ReaderMode,
    ) -> io::Result<(Self, ReaderMode)> {
        match mode {
            ReaderMode::Auto => {
                if has_tabix_index(file_path) {
                    TabixBackend::new(file_path, tabix_path)
                        .map(Box::new)
                        .map(|backend| (Self::Tabix(backend), ReaderMode::Query))
                } else {
                    SequentialBackend::new(file_path)
                        .map(|backend| (Self::Sequential(backend), ReaderMode::Sequential))
                }
            }
            ReaderMode::Sequential => SequentialBackend::new(file_path)
                .map(|backend| (Self::Sequential(backend), ReaderMode::Sequential)),
            ReaderMode::Query => {
                if has_tabix_index(file_path) {
                    TabixBackend::new(file_path, tabix_path)
                        .map(Box::new)
                        .map(|backend| (Self::Tabix(backend), ReaderMode::Query))
                } else {
                    io::Result::Err(io::Error::new(
                        io::ErrorKind::NotFound,
                        format!(
                            "no tabix index found for file {:?}; expected {:?} or {:?}",
                            file_path,
                            associated_index_path(file_path, "tbi"),
                            associated_index_path(file_path, "csi")
                        ),
                    ))
                }
            }
        }
    }

    fn next(&mut self) -> io::Result<Option<Gff3Record>> {
        match self {
            Self::Sequential(reader) => reader.next(),
            Self::Tabix(reader) => reader.next(),
        }
    }

    fn next_filtered(&mut self, filter: &FeatureTypeFilter) -> io::Result<Option<Gff3Record>> {
        match self {
            Self::Sequential(reader) => reader.next_filtered(filter),
            Self::Tabix(reader) => reader.next_filtered(filter),
        }
    }

    fn seek(&mut self, region: &str) -> io::Result<bool> {
        match self {
            Self::Sequential(reader) => reader.seek(region),
            Self::Tabix(reader) => reader.seek(region),
        }
    }

    fn index_sequence_names(&mut self) -> io::Result<Vec<String>> {
        match self {
            Self::Sequential(_) => Ok(Vec::new()),
            Self::Tabix(reader) => reader.index_sequence_names(),
        }
    }
}

// ---------------------------------------------------------------------------
// SequentialBackend - simple line-by-line reader
// ---------------------------------------------------------------------------

struct SequentialBackend {
    reader: gff3::io::Reader<Box<dyn BufRead + Send + Sync>>,
    /// Reused across calls so each record read does not allocate a fresh `Line` buffer.
    line: gff3::Line,
}

impl SequentialBackend {
    fn new<P>(file_path: P) -> io::Result<Self>
    where
        P: AsRef<Path>,
    {
        Ok(Self {
            reader: open_sequential_gff3_reader(file_path.as_ref())?,
            line: gff3::Line::default(),
        })
    }

    fn seek(&mut self, _region: &str) -> io::Result<bool> {
        Err(io::Error::other(
            "GFF3 reader is set to 'sequential'; set reader to 'query' before seek",
        ))
    }

    fn next(&mut self) -> io::Result<Option<Gff3Record>> {
        read_next_record(&mut self.reader, &mut self.line, None)
    }

    fn next_filtered(&mut self, filter: &FeatureTypeFilter) -> io::Result<Option<Gff3Record>> {
        read_next_record(&mut self.reader, &mut self.line, Some(filter))
    }
}

// ---------------------------------------------------------------------------
// TabixBackend - uses tabix for indexed reading
// ---------------------------------------------------------------------------

struct TabixBackend {
    reader: Option<tabix::TabixReader>,
    query_active: bool,
    /// Reused across calls so each queried line does not allocate a fresh `Line` buffer.
    line: gff3::Line,
}

impl TabixBackend {
    fn new<P>(file_path: P, tabix_path: &str) -> io::Result<Self>
    where
        P: AsRef<Path>,
    {
        Ok(Self {
            reader: Some(tabix::TabixReader::open(file_path.as_ref(), tabix_path)?),
            query_active: false,
            line: gff3::Line::default(),
        })
    }

    fn seek(&mut self, region: &str) -> io::Result<bool> {
        let has_records = self
            .reader
            .as_mut()
            .expect("tabix reader should be initialized")
            .seek(region)?;

        self.query_active = has_records;
        Ok(has_records)
    }

    fn next(&mut self) -> io::Result<Option<Gff3Record>> {
        if self.query_active {
            let Some(reader) = self.reader.as_mut() else {
                return Err(io::Error::other("tabix reader is not initialized"));
            };

            loop {
                let Some(line) = reader.next()? else {
                    self.query_active = false;
                    return Ok(None);
                };

                if let Some(record) = parse_gff3_line_into(&line, &mut self.line, None)? {
                    return Ok(Some(record));
                }
            }
        }

        Ok(None)
    }

    fn next_filtered(&mut self, filter: &FeatureTypeFilter) -> io::Result<Option<Gff3Record>> {
        if self.query_active {
            let Some(reader) = self.reader.as_mut() else {
                return Err(io::Error::other("tabix reader is not initialized"));
            };

            loop {
                let Some(line) = reader.next()? else {
                    self.query_active = false;
                    return Ok(None);
                };

                if let Some(record) = parse_gff3_line_into(&line, &mut self.line, Some(filter))? {
                    return Ok(Some(record));
                }
            }
        }

        Ok(None)
    }

    /// Reference sequence names declared by the tabix index, or empty when unopened.
    fn index_sequence_names(&mut self) -> io::Result<Vec<String>> {
        match self.reader.as_mut() {
            Some(reader) => reader.reference_names(),
            None => Ok(Vec::new()),
        }
    }
}

// ---------------------------------------------------------------------------
// PyO3 wrapper
// ---------------------------------------------------------------------------

#[pyclass(name = "GFF3Parser")]
pub struct Gff3Parser {
    file_path: Option<PathBuf>,
    reader: Option<Gff3Reader>,
    tabix_path: String,
    reader_mode: ReaderMode,
}

#[pymethods]
impl Gff3Parser {
    #[new]
    #[pyo3(signature = (file_path=None, tabix_path=None, mode=None))]
    fn new(
        file_path: Option<String>,
        tabix_path: Option<String>,
        mode: Option<String>,
    ) -> PyResult<Self> {
        let reader_mode = ReaderMode::parse(mode.as_deref())?;

        let mut parser = Self {
            file_path: None,
            reader: None,
            tabix_path: tabix_path.unwrap_or_else(|| "tabix".to_string()),
            reader_mode,
        };

        if let Some(path) = file_path {
            parser.open(&path, mode)?;
        }

        Ok(parser)
    }

    #[getter]
    fn format(&self) -> &'static str {
        "gff3"
    }

    #[getter]
    fn metadata<'py>(&mut self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);

        let Some(reader) = self.reader.as_mut() else {
            return Ok(dict);
        };

        let metadata = reader.metadata().map_err(to_py_error)?;

        if !metadata.header.is_empty() {
            let header: Vec<_> = metadata.header.iter().map(String::as_str).collect();
            dict.set_item("header", header)?;
        }

        for (key, value) in &metadata.directives {
            dict.set_item(key, value)?;
        }

        Ok(dict)
    }

    #[getter]
    fn reader_mode(&self) -> &'static str {
        self.reader_mode.as_str()
    }

    /// Reference sequence names declared by the tabix index, in index order.
    ///
    /// Empty when the file is unopened or read sequentially. Callers can use these names to
    /// resolve a query contig to the exact spelling the index declares before seeking, instead of
    /// retrying alternate spellings on a miss.
    #[getter]
    fn index_sequence_names(&mut self) -> PyResult<Vec<String>> {
        match self.reader.as_mut() {
            Some(reader) => reader.index_sequence_names().map_err(to_py_error),
            None => Ok(Vec::new()),
        }
    }

    fn open(&mut self, file_path: &str, mode: Option<String>) -> PyResult<()> {
        let path = PathBuf::from(file_path);
        let reader_mode = ReaderMode::parse(mode.as_deref())?;

        self.reader =
            Some(Gff3Reader::open(&path, &self.tabix_path, reader_mode).map_err(to_py_error)?);
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
    fn next(&mut self) -> PyResult<Option<Gff3Record>> {
        self.reader()?.next().map_err(to_py_error)
    }

    /// Return up to `limit` records, optionally keeping only rows whose column 3 (feature type)
    /// is in `feature_types` (compared case-insensitively).
    ///
    /// The filter is applied on the Rust side before a matching line is parsed into a record or
    /// converted into a Python object: a non-matching line costs one tab-split byte scan, not a
    /// full record build.
    #[pyo3(signature = (limit=10000, feature_types=None))]
    fn records(
        &mut self,
        limit: Option<usize>,
        feature_types: Option<Vec<String>>,
    ) -> PyResult<Vec<Gff3Record>> {
        let filter = feature_types.as_deref().map(FeatureTypeFilter::new);
        let mut records = Vec::new();

        while limit.is_none_or(|n| records.len() < n) {
            let next = match filter.as_ref() {
                Some(filter) => self.reader()?.next_filtered(filter),
                None => self.reader()?.next(),
            }
            .map_err(to_py_error)?;

            match next {
                Some(record) => records.push(record),
                None => break,
            }
        }

        Ok(records)
    }

    #[staticmethod]
    #[pyo3(signature = (data, limit=None))]
    fn parse_text(data: &str, limit: Option<usize>) -> PyResult<Vec<Gff3Record>> {
        parse_gff3_bytes(data.as_bytes(), limit).map_err(to_py_error)
    }

    #[staticmethod]
    #[pyo3(signature = (file_path, limit=None))]
    fn parse_file(file_path: &str, limit: Option<usize>) -> PyResult<Vec<Gff3Record>> {
        let mut reader = open_sequential_gff3_reader(Path::new(file_path)).map_err(to_py_error)?;

        read_records(&mut reader, limit).map_err(to_py_error)
    }
}

impl Gff3Parser {
    fn reader(&mut self) -> PyResult<&mut Gff3Reader> {
        self.reader
            .as_mut()
            .ok_or_else(|| PyRuntimeError::new_err("GFF3 parser is not open"))
    }
}

// ---------------------------------------------------------------------------
// Helpers: Noodles GFF3 decoder
// ---------------------------------------------------------------------------

fn open_sequential_gff3_reader(
    file_path: &Path,
) -> io::Result<gff3::io::Reader<Box<dyn BufRead + Send + Sync>>> {
    let reader = open_sequential_reader(file_path)?;
    Ok(gff3::io::Reader::new(reader))
}

fn read_gff3_metadata(file_path: &Path) -> io::Result<Gff3Metadata> {
    let mut reader = open_sequential_gff3_reader(file_path)?;
    read_gff3_metadata_from_reader(&mut reader)
}

fn read_gff3_metadata_from_reader<R>(reader: &mut gff3::io::Reader<R>) -> io::Result<Gff3Metadata>
where
    R: BufRead,
{
    let mut metadata = Gff3Metadata::default();
    let mut line = gff3::Line::default();

    loop {
        if reader.read_line(&mut line)? == 0 {
            break;
        }

        if line.as_ref().is_empty() {
            continue;
        }

        if let Some(directive) = line.as_directive() {
            metadata
                .header
                .push(bytes_to_string(line.as_ref().as_ref()));

            if let Some(value) = directive.value() {
                metadata.directives.push((
                    bytes_to_string(directive.key().as_ref()),
                    bytes_to_string(value.as_ref()),
                ));
            }

            continue;
        }

        if line.as_comment().is_some() {
            metadata
                .header
                .push(bytes_to_string(line.as_ref().as_ref()));
            continue;
        }

        break;
    }

    Ok(metadata)
}

/// Read the next record, reusing a caller-owned `Line` buffer across calls.
///
/// Hoisting the buffer out lets a hot loop (one `Line` per `TabixBackend` or
/// `SequentialBackend`, rather than a fresh one per record) avoid a per-record allocation.
///
/// `filter`, when given, is checked against the raw line's column 3 before `line.as_record()` is
/// called: a non-matching line is skipped without being parsed into a record or converted into a
/// Python object. `None` keeps today's unfiltered behaviour exactly.
fn read_next_record<R>(
    reader: &mut gff3::io::Reader<R>,
    line: &mut gff3::Line,
    filter: Option<&FeatureTypeFilter>,
) -> io::Result<Option<Gff3Record>>
where
    R: BufRead,
{
    loop {
        if reader.read_line(line)? == 0 {
            return Ok(None);
        }

        if line.as_ref().is_empty() {
            continue;
        }

        if let Some(directive) = line.as_directive() {
            if directive.key() == "FASTA" {
                return Ok(None);
            }

            continue;
        }

        if let Some(filter) = filter
            && !filter.matches(line.as_ref().as_ref())
        {
            continue;
        }

        let Some(result) = line.as_record() else {
            continue;
        };

        let raw_attributes = raw_attributes(line);
        let record = result?;

        return record_to_py(record, raw_attributes).map(Some);
    }
}

fn read_records<R>(
    reader: &mut gff3::io::Reader<R>,
    limit: Option<usize>,
) -> io::Result<Vec<Gff3Record>>
where
    R: BufRead,
{
    let mut records = Vec::new();
    let mut line = gff3::Line::default();

    while limit.is_none_or(|n| records.len() < n) {
        match read_next_record(reader, &mut line, None)? {
            Some(record) => records.push(record),
            None => break,
        }
    }

    Ok(records)
}

fn parse_gff3_bytes(data: &[u8], limit: Option<usize>) -> io::Result<Vec<Gff3Record>> {
    let mut reader = gff3::io::Reader::new(data);
    read_records(&mut reader, limit)
}

/// Parse one already-dechomped line, reusing a caller-owned `Line` buffer.
///
/// `&[u8]` implements `BufRead` directly, so this reads straight from `line`'s bytes: no
/// intermediate `Vec<u8>` copy and no new `Line` allocation per call, the two costs the
/// `TabixBackend` hot loop previously paid for every queried row.
///
/// `line` is not re-terminated with `\n` before parsing (the earlier `Vec<u8>` copy used to add
/// one back). `read_next_record`'s inner `read_line` only strips a trailing `\r` when it finds a
/// trailing `\n` to strip first, so if `line` still carries a residual `\r` -- a doubled `\r`
/// before the newline `TabixReader` already trimmed, or a final `\r` with no newline at all on
/// the last line of a file -- that `\r` now passes through into the record, matching what the
/// sequential reader (`SequentialBackend`, reading `\r`-terminated bytes with no `\n` to strip
/// against) already does. `test_parse_gff3_line_into_residual_cr_matches_sequential_reader` pins
/// this parity.
fn parse_gff3_line_into(
    line: &str,
    buffer: &mut gff3::Line,
    filter: Option<&FeatureTypeFilter>,
) -> io::Result<Option<Gff3Record>> {
    let mut reader = gff3::io::Reader::new(line.as_bytes());
    read_next_record(&mut reader, buffer, filter)
}

/// Test-only convenience over [`parse_gff3_line_into`]: production code always has a `TabixBackend`
/// to own the reused `Line` buffer, so only tests need a one-shot, buffer-allocating variant.
#[cfg(test)]
fn parse_gff3_line(line: &str) -> io::Result<Option<Gff3Record>> {
    parse_gff3_line_into(line, &mut gff3::Line::default(), None)
}

fn record_to_py(record: gff3::Record<'_>, raw_attributes: String) -> io::Result<Gff3Record> {
    Ok(Gff3Record {
        seqid: bytes_to_string(record.reference_sequence_name()),
        source: bytes_to_string(record.source()),
        feature_type: bytes_to_string(record.ty()),
        start: record.start()?.get(),
        end: record.end()?.get(),
        score: record.score().transpose()?.map(f64::from),
        strand_symbol: noodles_strand_to_string(record.strand()?),
        phase_value: record.phase().transpose()?.map(noodles_phase_to_u8),
        attributes: parse_attributes(record.attributes())?,
        raw_attributes,
    })
}

fn parse_attributes(attributes: gff3::record::Attributes<'_>) -> io::Result<AttributePairs> {
    let mut pairs = Vec::new();

    for result in attributes.iter() {
        let (tag, value) = result?;
        let values = match value {
            AttributeValue::String(value) => vec![bytes_to_string(value.as_ref())],
            AttributeValue::Array(values) => values
                .iter()
                .map(|value| bytes_to_string(value.as_ref()))
                .collect(),
        };

        pairs.push((bytes_to_string(tag.as_ref()), values));
    }

    Ok(pairs)
}

fn raw_attributes(line: &gff3::Line) -> String {
    const FIELD_DELIMITER: u8 = b'\t';
    const ATTRIBUTE_FIELD_INDEX: usize = 8;

    let bytes: &[u8] = line.as_ref().as_ref();
    let mut field_count = 0;

    for (i, b) in bytes.iter().enumerate() {
        if *b == FIELD_DELIMITER {
            field_count += 1;

            if field_count == ATTRIBUTE_FIELD_INDEX {
                return bytes_to_string(&bytes[i + 1..]);
            }
        }
    }

    String::new()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        parsers::test_utils::{write_bgzf_file, write_gff3_with_tbi},
        test_utils::{create_temp_dir, write_sample_file},
        utils::{from_py_error, trim_newline},
    };
    use std::{fs, str};

    const SAMPLE: &[u8] = b"\
##gff-version 3
chr1\tensembl\tgene\t10\t20\t.\t+\t.\tID=gene:1;Name=BRCA2;Alias=a,b
chr1\tensembl\tCDS\t12\t18\t4.5\t-\t2\tID=cds:1;Parent=gene:1
##FASTA
>chr1
ACGT
";

    const SKIPPED_LINE_SAMPLE: &[u8] = b"\
##gff-version 3

# a regular comment
chr1\tensembl\tgene\t10\t20\t.\t.\t.\tID=gene:1
";

    // -----------------------------------------------------------------------
    // Helpers
    // -----------------------------------------------------------------------

    fn to_gff3_line(s: &str) -> gff3::Line {
        let mut reader = gff3::io::Reader::new(s.as_bytes());
        let mut line = gff3::Line::default();
        reader.read_line(&mut line).expect("failed to read line");
        line
    }

    fn sample_lines() -> Vec<String> {
        str::from_utf8(SAMPLE)
            .expect("sample should be valid UTF-8")
            .lines()
            .map(str::to_owned)
            .collect()
    }

    fn assert_first_sample_record(record: &Gff3Record) {
        assert_eq!(record.seqid, "chr1");
        assert_eq!(record.source, "ensembl");
        assert_eq!(record.feature_type, "gene");
        assert_eq!(record.start, 10);
        assert_eq!(record.end, 20);
        assert_eq!(record.score, None);
        assert_eq!(record.strand_symbol, "+");
        assert_eq!(record.phase_value, None);
        assert_eq!(record.raw_attributes, "ID=gene:1;Name=BRCA2;Alias=a,b");
        assert_eq!(
            record.attributes,
            vec![
                ("ID".to_string(), vec!["gene:1".to_string()]),
                ("Name".to_string(), vec!["BRCA2".to_string()]),
                ("Alias".to_string(), vec!["a".to_string(), "b".to_string()]),
            ]
        );
    }

    fn assert_second_sample_record(record: &Gff3Record) {
        assert_eq!(record.seqid, "chr1");
        assert_eq!(record.feature_type, "CDS");
        assert_eq!(record.start, 12);
        assert_eq!(record.end, 18);
        assert_eq!(record.score, Some(4.5));
        assert_eq!(record.strand_symbol, "-");
        assert_eq!(record.phase_value, Some(2));
        assert_eq!(
            record.attributes,
            vec![
                ("ID".to_string(), vec!["cds:1".to_string()]),
                ("Parent".to_string(), vec!["gene:1".to_string()]),
            ]
        );
    }

    // -----------------------------------------------------------------------
    // raw_attributes
    // -----------------------------------------------------------------------

    #[test]
    fn test_raw_attributes_with_attributes() -> io::Result<()> {
        let line =
            to_gff3_line("chr1\tensembl\tgene\t10\t20\t.\t+\t.\tID=gene:1;Name=BRCA2;Alias=a,b");
        assert_eq!(raw_attributes(&line), "ID=gene:1;Name=BRCA2;Alias=a,b");
        Ok(())
    }

    #[test]
    fn test_raw_attributes_missing_attribute_field_returns_empty() -> io::Result<()> {
        let line = to_gff3_line("chr1\tensembl\tgene\t10\t20\t.\t+\t.");
        assert_eq!(raw_attributes(&line), "");
        Ok(())
    }

    // -----------------------------------------------------------------------
    // parse_attributes
    // -----------------------------------------------------------------------

    #[test]
    fn test_parse_attributes_single_valued() -> io::Result<()> {
        let line = to_gff3_line("chr1\tensembl\tgene\t10\t20\t.\t+\t.\tID=gene:1;Name=BRCA2");
        let record = line.as_record().unwrap().unwrap();
        let attributes = record.attributes();
        let pairs = parse_attributes(attributes)?;

        assert_eq!(pairs.len(), 2);
        assert_eq!(pairs[0].0, "ID");
        assert_eq!(pairs[0].1, vec!["gene:1".to_string()]);
        assert_eq!(pairs[1].0, "Name");
        assert_eq!(pairs[1].1, vec!["BRCA2".to_string()]);
        Ok(())
    }

    #[test]
    fn test_parse_attributes_multi_valued() -> io::Result<()> {
        let line = to_gff3_line("chr1\tensembl\tgene\t10\t20\t.\t+\t.\tAlias=a,b,c");
        let record = line.as_record().unwrap().unwrap();
        let attributes = record.attributes();
        let pairs = parse_attributes(attributes)?;

        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0].0, "Alias");
        assert_eq!(
            pairs[0].1,
            vec!["a".to_string(), "b".to_string(), "c".to_string()]
        );
        Ok(())
    }

    #[test]
    fn test_parse_attributes_percent_encoded() -> io::Result<()> {
        let line = to_gff3_line("chr1\tensembl\tgene\t10\t20\t.\t+\t.\tName=BRCA%202");
        let record = line.as_record().unwrap().unwrap();
        let attributes = record.attributes();
        let pairs = parse_attributes(attributes)?;

        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0].0, "Name");
        assert_eq!(pairs[0].1, vec!["BRCA 2".to_string()]);
        Ok(())
    }

    #[test]
    fn test_parse_attributes_invalid() -> io::Result<()> {
        let line = to_gff3_line("chr1\tensembl\tgene\t10\t20\t.\t+\t.\tID=gene:1;Invalid");
        let record = line.as_record().unwrap().unwrap();
        let attributes = record.attributes();
        let err = parse_attributes(attributes).unwrap_err();

        assert!(
            err.to_string().contains("invalid tag"),
            "unexpected error: {}",
            err
        );
        Ok(())
    }

    #[test]
    fn test_parse_attributes_empty() -> io::Result<()> {
        let line = to_gff3_line("chr1\tensembl\tgene\t10\t20\t.\t+\t.\t");
        let record = line.as_record().unwrap().unwrap();
        let attributes = record.attributes();
        let pairs = parse_attributes(attributes)?;

        assert!(pairs.is_empty());
        Ok(())
    }

    // -----------------------------------------------------------------------
    // record_to_py
    // -----------------------------------------------------------------------

    #[test]
    fn test_record_to_py() -> io::Result<()> {
        let line =
            to_gff3_line("chr1\tensembl\tgene\t10\t20\t.\t+\t.\tID=gene:1;Name=BRCA2;Alias=a,b");
        let record = line.as_record().unwrap()?;
        let raw_attributes = raw_attributes(&line);
        let py_record = record_to_py(record, raw_attributes)?;

        assert_eq!(py_record.seqid, "chr1");
        assert_eq!(py_record.source, "ensembl");
        assert_eq!(py_record.feature_type, "gene");
        assert_eq!(py_record.start, 10);
        assert_eq!(py_record.end, 20);
        assert_eq!(py_record.score, None);
        assert_eq!(py_record.strand_symbol, "+");
        assert_eq!(py_record.phase_value, None);
        assert_eq!(py_record.raw_attributes, "ID=gene:1;Name=BRCA2;Alias=a,b");
        assert_eq!(
            py_record.attributes,
            vec![
                ("ID".to_string(), vec!["gene:1".to_string()]),
                ("Name".to_string(), vec!["BRCA2".to_string()]),
                ("Alias".to_string(), vec!["a".to_string(), "b".to_string()]),
            ]
        );
        Ok(())
    }

    // -----------------------------------------------------------------------
    // parse_gff3_line
    // -----------------------------------------------------------------------

    #[test]
    fn test_parse_gff3_line_decodes_record() -> io::Result<()> {
        let record = parse_gff3_line(
            "chr1\tensembl\tgene\t10\t20\t.\t?\t.\tID=gene:1;Dbxref=GeneID:1,HGNC:1",
        )?
        .expect("expected record");

        assert_eq!(record.seqid, "chr1");
        assert_eq!(record.source, "ensembl");
        assert_eq!(record.feature_type, "gene");
        assert_eq!(record.start, 10);
        assert_eq!(record.end, 20);
        assert_eq!(record.score, None);
        assert_eq!(record.strand_symbol, "?");
        assert_eq!(record.phase_value, None);
        assert_eq!(
            record.attributes,
            vec![
                ("ID".to_string(), vec!["gene:1".to_string()]),
                (
                    "Dbxref".to_string(),
                    vec!["GeneID:1".to_string(), "HGNC:1".to_string()]
                ),
            ]
        );

        Ok(())
    }

    #[test]
    fn test_parse_gff3_line_invalid_line_bad_fields() -> io::Result<()> {
        let err =
            parse_gff3_line("chr1\tensembl\tgene\tbad_start\t20\t.\t+\t.\tID=gene:1").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);

        let err =
            parse_gff3_line("chr1\tensembl\tgene\t10\tbad_end\t.\t+\t.\tID=gene:1").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);

        let err =
            parse_gff3_line("chr1\tensembl\tgene\t10\t20\tbad_score\t+\t.\tID=gene:1").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);

        let err = parse_gff3_line("chr1\tensembl\tgene\t10\t20\t.\tbad_strand\t.\tID=gene:1")
            .unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);

        let err =
            parse_gff3_line("chr1\tensembl\tgene\t10\t20\t.\t+\tbad_phase\tID=gene:1").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);

        let err =
            parse_gff3_line("chr1\tensembl\tgene\t10\t20\t.\t+\t.\tInvalidAttribute").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        Ok(())
    }

    #[test]
    fn test_parse_gff3_line_invalid_line_missing_fields() -> io::Result<()> {
        let err = parse_gff3_line("chr1\tensembl\tgene\t10\t20\t.\t+\t.").unwrap_err();

        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
        Ok(())
    }

    #[test]
    fn test_parse_gff3_line_ignores_non_record_lines() -> io::Result<()> {
        assert_eq!(parse_gff3_line("")?, None);
        assert_eq!(parse_gff3_line("##gff-version 3")?, None);
        assert_eq!(parse_gff3_line("# a regular comment")?, None);

        Ok(())
    }

    // -----------------------------------------------------------------------
    // parse_gff3_bytes
    // -----------------------------------------------------------------------

    #[test]
    fn test_parse_gff3_bytes_decodes_records() -> io::Result<()> {
        let records = parse_gff3_bytes(SAMPLE, None)?;

        assert_eq!(records.len(), 2);
        assert_first_sample_record(&records[0]);
        assert_second_sample_record(&records[1]);

        Ok(())
    }

    #[test]
    fn test_parse_gff3_bytes_honours_limit() -> io::Result<()> {
        let records = parse_gff3_bytes(SAMPLE, Some(1))?;

        assert_eq!(records.len(), 1);
        Ok(())
    }

    #[test]
    fn test_parse_gff3_bytes_honours_zero_limit() -> io::Result<()> {
        let records = parse_gff3_bytes(SAMPLE, Some(0))?;

        assert!(records.is_empty());
        Ok(())
    }

    #[test]
    fn test_parse_gff3_bytes_skips_empty_directive_and_comment_lines() -> io::Result<()> {
        let records = parse_gff3_bytes(SKIPPED_LINE_SAMPLE, None)?;

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].strand_symbol, ".");
        assert_eq!(records[0].raw_attributes, "ID=gene:1");

        Ok(())
    }

    #[test]
    fn test_parse_gff3_bytes_stops_at_fasta() -> io::Result<()> {
        let records = parse_gff3_bytes(SAMPLE, Some(10))?;

        assert_eq!(records.len(), 2);
        Ok(())
    }

    #[test]
    fn test_parse_gff3_bytes_invalid_record_errors() {
        let data = b"chr1\tensembl\tgene\tbad_start\t20\t.\t+\t.\tID=gene:1\n";
        let err = parse_gff3_bytes(data, None).unwrap_err();

        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    }

    // -----------------------------------------------------------------------
    // read_records
    // -----------------------------------------------------------------------

    #[test]
    fn test_read_records_reads_multiple_records() -> io::Result<()> {
        let mut reader = gff3::io::Reader::new(SAMPLE);
        let records = read_records(&mut reader, None)?;

        assert_eq!(records.len(), 2);
        assert_first_sample_record(&records[0]);
        assert_second_sample_record(&records[1]);
        Ok(())
    }

    #[test]
    fn test_read_records_honours_limit() -> io::Result<()> {
        let mut reader = gff3::io::Reader::new(SAMPLE);
        let records = read_records(&mut reader, Some(1))?;

        assert_eq!(records.len(), 1);
        Ok(())
    }

    #[test]
    fn test_read_records_honours_zero_limit() -> io::Result<()> {
        let mut reader = gff3::io::Reader::new(SAMPLE);
        let records = read_records(&mut reader, Some(0))?;

        assert!(records.is_empty());
        Ok(())
    }

    // -----------------------------------------------------------------------
    // read_next_record
    // -----------------------------------------------------------------------

    #[test]
    fn test_read_next_record_reads_records_one_by_one() -> io::Result<()> {
        let mut reader = gff3::io::Reader::new(SAMPLE);
        let mut line = gff3::Line::default();

        let first = read_next_record(&mut reader, &mut line, None)?.expect("expected first record");
        assert_first_sample_record(&first);

        let second =
            read_next_record(&mut reader, &mut line, None)?.expect("expected second record");
        assert_second_sample_record(&second);
        Ok(())
    }

    #[test]
    fn test_read_next_record_returns_none_at_eof() -> io::Result<()> {
        let mut reader = gff3::io::Reader::new(SAMPLE);
        let mut line = gff3::Line::default();
        read_next_record(&mut reader, &mut line, None)?; // first record
        read_next_record(&mut reader, &mut line, None)?; // second record

        assert_eq!(read_next_record(&mut reader, &mut line, None)?, None);
        Ok(())
    }

    #[test]
    fn test_read_next_record_no_records() -> io::Result<()> {
        let mut reader =
            gff3::io::Reader::new(b"##gff-version 3\n##FASTA\n>chr1\nACGT\n".as_slice());
        let mut line = gff3::Line::default();

        assert_eq!(read_next_record(&mut reader, &mut line, None)?, None);
        Ok(())
    }

    /// Reusing one `Line` buffer across calls must not leak a previous record's bytes into the
    /// next result — the earlier record has already been converted into an owned `Gff3Record`
    /// before the buffer is reused, so there is nothing to leak, but this pins that invariant.
    #[test]
    fn test_read_next_record_reused_buffer_does_not_leak_previous_line() -> io::Result<()> {
        let mut reader = gff3::io::Reader::new(SAMPLE);
        let mut line = gff3::Line::default();

        let first = read_next_record(&mut reader, &mut line, None)?.expect("expected first record");
        let second =
            read_next_record(&mut reader, &mut line, None)?.expect("expected second record");

        assert_ne!(first.feature_type, second.feature_type);
        assert_first_sample_record(&first);
        assert_second_sample_record(&second);
        Ok(())
    }

    // -----------------------------------------------------------------------
    // read_gff3_metadata
    // -----------------------------------------------------------------------

    #[test]
    fn test_read_gff3_metadata_from_reader_collects_header_and_directives() -> io::Result<()> {
        let mut reader = gff3::io::Reader::new(SKIPPED_LINE_SAMPLE);
        let metadata = read_gff3_metadata_from_reader(&mut reader)?;

        assert_eq!(
            metadata.header,
            vec![
                "##gff-version 3".to_string(),
                "# a regular comment".to_string()
            ]
        );
        assert_eq!(
            metadata.directives,
            vec![("gff-version".to_string(), "3".to_string())]
        );
        Ok(())
    }

    #[test]
    fn test_read_gff3_metadata_with_no_header() -> io::Result<()> {
        let mut reader =
            gff3::io::Reader::new(b"chr1\tensembl\tgene\t10\t20\t.\t+\t.\tID=gene:1\n".as_slice());
        let metadata = read_gff3_metadata_from_reader(&mut reader)?;

        assert!(metadata.header.is_empty());
        assert!(metadata.directives.is_empty());
        Ok(())
    }

    // -----------------------------------------------------------------------
    // open_sequential_gff3_reader
    // -----------------------------------------------------------------------

    #[test]
    fn test_open_sequential_gff3_reader_reads_plain_file() -> io::Result<()> {
        let dir = create_temp_dir("gff3-open-sequential-plain")?;
        let path = write_sample_file(&dir, "sample.gff3", SAMPLE)?;
        let mut reader = open_sequential_gff3_reader(&path)?;
        let records = read_records(&mut reader, None)?;

        assert_eq!(records.len(), 2);
        assert_first_sample_record(&records[0]);
        assert_second_sample_record(&records[1]);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_open_sequential_gff3_reader_reads_bgzf_file() -> io::Result<()> {
        let dir = create_temp_dir("gff3-open-sequential-bgzf")?;
        let path = write_bgzf_file(&dir, sample_lines(), Some("sample.gff3.gz"))?;
        let mut reader = open_sequential_gff3_reader(&path)?;
        let records = read_records(&mut reader, None)?;

        assert_eq!(records.len(), 2);
        assert_first_sample_record(&records[0]);
        assert_second_sample_record(&records[1]);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // Gff3Parser:metadata
    // -----------------------------------------------------------------------

    #[test]
    fn test_gff3_parser_metadata_sequential() -> io::Result<()> {
        let dir = create_temp_dir("gff3-parser-metadata-seq")?;
        let path = write_sample_file(&dir, "sample.gff3", SKIPPED_LINE_SAMPLE)?;
        let file_path = path.to_string_lossy().into_owned();

        let mut parser =
            Gff3Parser::new(None, None, Some("sequential".to_string())).map_err(from_py_error)?;
        parser.open(&file_path, None).map_err(from_py_error)?;

        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let metadata = parser.metadata(py)?;
            let header: Vec<String> = metadata
                .get_item("header")?
                .expect("metadata should contain header")
                .extract()?;

            assert_eq!(
                header,
                vec![
                    "##gff-version 3".to_string(),
                    "# a regular comment".to_string()
                ]
            );
            assert_eq!(
                metadata
                    .get_item("gff-version")?
                    .unwrap()
                    .extract::<String>()?,
                "3"
            );
            Ok(())
        })
        .map_err(from_py_error)?;

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_gff3_parser_metadata_does_not_consume_sequential_reader() -> io::Result<()> {
        let dir = create_temp_dir("gff3-parser-metadata-preserve-seq")?;
        let path = write_sample_file(&dir, "sample.gff3", SAMPLE)?;
        let file_path = path.to_string_lossy().into_owned();

        let mut parser =
            Gff3Parser::new(None, None, Some("sequential".to_string())).map_err(from_py_error)?;
        parser.open(&file_path, None).map_err(from_py_error)?;

        let first = parser
            .next()
            .map_err(from_py_error)?
            .expect("expected first record");
        assert_first_sample_record(&first);

        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let metadata = parser.metadata(py)?;
            assert_eq!(
                metadata
                    .get_item("gff-version")?
                    .unwrap()
                    .extract::<String>()?,
                "3"
            );
            Ok(())
        })
        .map_err(from_py_error)?;

        let second = parser
            .next()
            .map_err(from_py_error)?
            .expect("expected second record");
        assert_second_sample_record(&second);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_gff3_parser_metadata_query_bgzf() -> io::Result<()> {
        let dir = create_temp_dir("gff3-parser-metadata-query")?;
        let (gz_path, _tbi_path) = write_gff3_with_tbi(
            &dir,
            str::from_utf8(SKIPPED_LINE_SAMPLE)
                .ok()
                .map(|s| s.lines().map(str::to_owned).collect()),
        )?;
        let file_path = gz_path.to_string_lossy().into_owned();

        let mut parser = Gff3Parser::new(None, None, None).map_err(from_py_error)?;
        parser.open(&file_path, None).map_err(from_py_error)?;

        assert_eq!(parser.reader_mode(), "query");

        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let metadata = parser.metadata(py)?;
            let header: Vec<String> = metadata
                .get_item("header")?
                .expect("metadata should contain header")
                .extract()?;

            assert_eq!(
                header,
                vec![
                    "##gff-version 3".to_string(),
                    "# a regular comment".to_string()
                ]
            );
            assert_eq!(
                metadata
                    .get_item("gff-version")?
                    .unwrap()
                    .extract::<String>()?,
                "3"
            );
            Ok(())
        })
        .map_err(from_py_error)?;

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_gff3_parser_open_missing_file_errors() -> io::Result<()> {
        Python::initialize();
        let mut parser = Gff3Parser::new(None, None, None).map_err(from_py_error)?;
        let err = parser
            .open("/definitely/missing/sample.gff3", None)
            .unwrap_err();

        assert!(
            err.to_string().contains("No such file") || err.to_string().contains("os error 2"),
            "unexpected error: {err}"
        );
        Ok(())
    }

    // -----------------------------------------------------------------------
    // Gff3Parser:open
    // -----------------------------------------------------------------------

    #[test]
    fn test_gff3_parser_without_open_reader() -> io::Result<()> {
        let mut parser = Gff3Parser::new(None, None, None).map_err(from_py_error)?;

        assert_eq!(parser.format(), "gff3");
        Python::initialize();
        Python::attach(|py| {
            assert!(parser.metadata(py).unwrap().is_empty());
        });
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
                .records(None, None)
                .unwrap_err()
                .to_string()
                .contains("not open")
        );

        parser.close();
        Ok(())
    }

    #[test]
    fn test_gff3_parser_open_override_reader_mode() -> io::Result<()> {
        let dir = create_temp_dir("gff3-parser-open-override")?;
        let (gz_path, _tbi_path) = write_gff3_with_tbi(
            &dir,
            str::from_utf8(SAMPLE)
                .ok()
                .map(|s| s.lines().map(str::to_owned).collect()),
        )?;
        let file_path = gz_path.to_string_lossy().into_owned();

        let mut parser =
            Gff3Parser::new(None, None, Some("sequential".to_string())).map_err(from_py_error)?;
        assert_eq!(parser.reader_mode(), "sequential");

        parser
            .open(&file_path, Some("query".to_string()))
            .map_err(from_py_error)?;
        assert_eq!(parser.reader_mode(), "query");

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // Gff3Parser:open - auto reader modes
    // -----------------------------------------------------------------------

    #[test]
    fn test_gff3_parser_open_file_auto_mode_sequential() -> io::Result<()> {
        let dir = create_temp_dir("gff3-parser-open-auto")?;
        let path = write_sample_file(&dir, "sample.gff3", SAMPLE)?;
        let file_path = path.to_string_lossy().into_owned();

        let mut parser = Gff3Parser::new(None, None, None).map_err(from_py_error)?;
        parser.open(&file_path, None).map_err(from_py_error)?;

        assert_eq!(parser.reader_mode(), "sequential");

        let records = parser.records(None, None).map_err(from_py_error)?;
        assert_eq!(records.len(), 2);
        assert_first_sample_record(&records[0]);
        assert_second_sample_record(&records[1]);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_gff3_parser_open_file_auto_mode_query() -> io::Result<()> {
        let dir = create_temp_dir("gff3-parser-open-auto-query")?;
        let (gz_path, _tbi_path) = write_gff3_with_tbi(
            &dir,
            str::from_utf8(SAMPLE)
                .ok()
                .map(|s| s.lines().map(str::to_owned).collect()),
        )?;
        let file_path = gz_path.to_string_lossy().into_owned();

        let mut parser = Gff3Parser::new(None, None, None).map_err(from_py_error)?;
        parser.open(&file_path, None).map_err(from_py_error)?;

        assert_eq!(parser.reader_mode(), "query");
        parser.seek("chr1:10-20").map_err(from_py_error)?;

        let records = parser.records(None, None).map_err(from_py_error)?;
        assert_eq!(records.len(), 2);
        assert_first_sample_record(&records[0]);
        assert_second_sample_record(&records[1]);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // Gff3Parser:records - feature_types filter
    // -----------------------------------------------------------------------

    const FILTER_SAMPLE: &[u8] = b"\
##gff-version 3
chr1\tensembl\tgene\t1\t100\t.\t+\t.\tID=gene:1
chr1\tensembl\tmRNA\t1\t100\t.\t+\t.\tID=transcript:1;Parent=gene:1
chr1\tensembl\texon\t1\t50\t.\t+\t.\tID=exon:1;Parent=transcript:1
# a comment between records
chr1\tensembl\tCDS\t1\t50\t.\t+\t0\tParent=transcript:1
chr1\tensembl\tmrna\t60\t100\t.\t+\t.\tID=transcript:2;Parent=gene:1
";

    #[test]
    fn test_gff3_parser_records_feature_types_filter_sequential() -> io::Result<()> {
        let dir = create_temp_dir("gff3-parser-records-filter-seq")?;
        let path = write_sample_file(&dir, "sample.gff3", FILTER_SAMPLE)?;
        let file_path = path.to_string_lossy().into_owned();

        let mut parser =
            Gff3Parser::new(None, None, Some("sequential".to_string())).map_err(from_py_error)?;
        parser.open(&file_path, None).map_err(from_py_error)?;

        // Mixed case in both the file ("mRNA"/"mrna") and the filter ("MRNA") must still match.
        let records = parser
            .records(None, Some(vec!["MRNA".to_string()]))
            .map_err(from_py_error)?;

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].feature_type, "mRNA");
        assert_eq!(records[1].feature_type, "mrna");

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_gff3_parser_records_feature_types_filter_query() -> io::Result<()> {
        let dir = create_temp_dir("gff3-parser-records-filter-query")?;
        let (gz_path, _tbi_path) = write_gff3_with_tbi(
            &dir,
            str::from_utf8(FILTER_SAMPLE)
                .ok()
                .map(|s| s.lines().map(str::to_owned).collect()),
        )?;
        let file_path = gz_path.to_string_lossy().into_owned();

        let mut parser = Gff3Parser::new(None, None, None).map_err(from_py_error)?;
        parser.open(&file_path, None).map_err(from_py_error)?;
        assert_eq!(parser.reader_mode(), "query");
        parser.seek("chr1:1-100").map_err(from_py_error)?;

        let records = parser
            .records(None, Some(vec!["mrna".to_string()]))
            .map_err(from_py_error)?;

        assert_eq!(records.len(), 2);
        assert!(
            records
                .iter()
                .all(|record| record.feature_type.eq_ignore_ascii_case("mrna"))
        );

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_gff3_parser_records_no_filter_matches_unfiltered_records() -> io::Result<()> {
        let dir = create_temp_dir("gff3-parser-records-no-filter")?;
        let path = write_sample_file(&dir, "sample.gff3", FILTER_SAMPLE)?;
        let file_path = path.to_string_lossy().into_owned();

        let mut unfiltered_parser = Gff3Parser::new(None, None, None).map_err(from_py_error)?;
        unfiltered_parser
            .open(&file_path, None)
            .map_err(from_py_error)?;
        let unfiltered = unfiltered_parser
            .records(None, None)
            .map_err(from_py_error)?;

        let mut none_filter_parser = Gff3Parser::new(None, None, None).map_err(from_py_error)?;
        none_filter_parser
            .open(&file_path, None)
            .map_err(from_py_error)?;
        let explicit_none = none_filter_parser
            .records(None, None)
            .map_err(from_py_error)?;

        assert_eq!(unfiltered, explicit_none);
        assert_eq!(unfiltered.len(), 5);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_gff3_parser_records_feature_types_filter_skips_comments() -> io::Result<()> {
        let dir = create_temp_dir("gff3-parser-records-filter-comments")?;
        let path = write_sample_file(&dir, "sample.gff3", FILTER_SAMPLE)?;
        let file_path = path.to_string_lossy().into_owned();

        let mut parser = Gff3Parser::new(None, None, None).map_err(from_py_error)?;
        parser.open(&file_path, None).map_err(from_py_error)?;

        let records = parser
            .records(None, Some(vec!["CDS".to_string()]))
            .map_err(from_py_error)?;

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].feature_type, "CDS");

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    /// A line with fewer than 3 tab-separated fields has no column 3 to compare against the
    /// filter, so it is silently skipped rather than reaching the full parser (which would raise
    /// on it, as the unfiltered case below confirms it still does without a filter active).
    #[test]
    fn test_gff3_parser_records_feature_types_filter_skips_short_lines_without_erroring()
    -> io::Result<()> {
        let dir = create_temp_dir("gff3-parser-records-filter-short-line")?;
        let short_line_sample = b"\
##gff-version 3
chr1\tensembl
chr1\tensembl\tCDS\t1\t50\t.\t+\t0\tParent=transcript:1
";
        let path = write_sample_file(&dir, "sample.gff3", short_line_sample)?;
        let file_path = path.to_string_lossy().into_owned();

        let mut parser = Gff3Parser::new(None, None, None).map_err(from_py_error)?;
        parser.open(&file_path, None).map_err(from_py_error)?;

        let records = parser
            .records(None, Some(vec!["CDS".to_string()]))
            .map_err(from_py_error)?;

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].feature_type, "CDS");

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    /// The same short line, unfiltered, still reaches the full parser and errors -- the filter
    /// changes what happens to a short line only when it is actually applied.
    #[test]
    fn test_gff3_parser_records_no_filter_still_errors_on_short_line() -> io::Result<()> {
        let dir = create_temp_dir("gff3-parser-records-no-filter-short-line")?;
        let short_line_sample = b"chr1\tensembl\n";
        let path = write_sample_file(&dir, "sample.gff3", short_line_sample)?;
        let file_path = path.to_string_lossy().into_owned();

        Python::initialize();
        let mut parser = Gff3Parser::new(None, None, None).map_err(from_py_error)?;
        parser.open(&file_path, None).map_err(from_py_error)?;

        let err = parser.records(None, None).unwrap_err();
        assert!(
            err.to_string().contains("unexpected end of file"),
            "unexpected error: {err}"
        );

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_gff3_parser_records_feature_types_filter_honours_limit() -> io::Result<()> {
        let dir = create_temp_dir("gff3-parser-records-filter-limit")?;
        let path = write_sample_file(&dir, "sample.gff3", FILTER_SAMPLE)?;
        let file_path = path.to_string_lossy().into_owned();

        let mut parser = Gff3Parser::new(None, None, None).map_err(from_py_error)?;
        parser.open(&file_path, None).map_err(from_py_error)?;

        let records = parser
            .records(Some(1), Some(vec!["mrna".to_string()]))
            .map_err(from_py_error)?;

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].feature_type, "mRNA");

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_gff3_parser_open_file_auto_mode_query_then_sequential() -> io::Result<()> {
        let dir = create_temp_dir("gff3-parser-open-auto-query-sequential")?;
        let (gz_path, _tbi_path) = write_gff3_with_tbi(
            &dir,
            str::from_utf8(SAMPLE)
                .ok()
                .map(|s| s.lines().map(str::to_owned).collect()),
        )?;
        let first_file_path = gz_path.to_string_lossy().into_owned();

        let second_path = write_sample_file(&dir, "sample.gff3", SAMPLE)?;
        let second_file_path = second_path.to_string_lossy().into_owned();

        let mut parser = Gff3Parser::new(None, None, None).map_err(from_py_error)?;
        parser.open(&first_file_path, None).map_err(from_py_error)?;

        assert_eq!(parser.reader_mode(), "query");
        parser.seek("chr1:10-20").map_err(from_py_error)?;

        parser
            .open(&second_file_path, None)
            .map_err(from_py_error)?;
        assert!(parser.reader_mode() == "sequential");

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // Gff3Parser:open - explicit query mode
    // -----------------------------------------------------------------------

    #[test]
    fn test_gff3_parser_open_file_query_mode_without_index_errors() -> io::Result<()> {
        let dir = create_temp_dir("gff3-parser-open-query-no-index")?;
        let path = write_sample_file(&dir, "sample.gff3", SAMPLE)?;
        let file_path = path.to_string_lossy().into_owned();

        Python::initialize();
        let mut parser =
            Gff3Parser::new(None, None, Some("query".to_string())).map_err(from_py_error)?;

        let err = parser
            .open(&file_path, Some("query".to_string()))
            .unwrap_err();
        assert!(
            err.to_string().contains("no tabix index found"),
            "unexpected error: {err}"
        );

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_gff3_parser_open_file_query_mode_with_index() -> io::Result<()> {
        let dir = create_temp_dir("gff3-parser-open-query-index")?;
        let (gz_path, _tbi_path) = write_gff3_with_tbi(
            &dir,
            str::from_utf8(SAMPLE)
                .ok()
                .map(|s| s.lines().map(str::to_owned).collect()),
        )?;
        let file_path = gz_path.to_string_lossy().into_owned();

        let mut parser =
            Gff3Parser::new(None, None, Some("query".to_string())).map_err(from_py_error)?;
        parser.open(&file_path, None).map_err(from_py_error)?;

        assert_eq!(parser.reader_mode(), "query");
        parser.seek("chr1:10-20").map_err(from_py_error)?;
        let records = parser.records(None, None).map_err(from_py_error)?;
        assert_eq!(records.len(), 2);
        assert_first_sample_record(&records[0]);
        assert_second_sample_record(&records[1]);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_gff3_parser_open_file_query_mode_invalid_tabix_path() -> io::Result<()> {
        let dir = create_temp_dir("gff3-parser-open-query-invalid-tabix")?;
        let path = write_sample_file(&dir, "sample.gff3.gz", SAMPLE)?;
        let file_path = path.to_string_lossy().into_owned();

        Python::initialize();
        let mut parser = Gff3Parser::new(
            None,
            Some(dir.join("missing-tabix").to_string_lossy().into_owned()),
            Some("query".to_string()),
        )
        .map_err(from_py_error)?;

        let err = parser
            .open(&file_path, Some("query".to_string()))
            .unwrap_err();
        assert!(
            err.to_string().contains("no tabix index found"),
            "unexpected error: {err}"
        );

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // Gff3Parser:open - explicit sequential mode
    // -----------------------------------------------------------------------

    #[test]
    fn test_gff3_parser_open_file_sequential_mode() -> io::Result<()> {
        let dir = create_temp_dir("gff3-parser-open-seq")?;
        let path = write_sample_file(&dir, "sample.gff3", SAMPLE)?;
        let file_path = path.to_string_lossy().into_owned();

        let mut parser = Gff3Parser::new(
            None,
            Some(dir.join("missing-tabix").to_string_lossy().into_owned()),
            Some("sequential".to_string()),
        )
        .map_err(from_py_error)?;

        parser.open(&file_path, None).map_err(from_py_error)?;
        assert_eq!(parser.reader_mode(), "sequential");

        let records = parser.records(None, None).map_err(from_py_error)?;
        assert_eq!(records.len(), 2);
        assert_first_sample_record(&records[0]);
        assert_second_sample_record(&records[1]);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_gff3_parser_open_file_sequential_mode_with_index() -> io::Result<()> {
        let dir = create_temp_dir("gff3-parser-open-seq-index")?;
        let (gz_path, _tbi_path) = write_gff3_with_tbi(
            &dir,
            str::from_utf8(SAMPLE)
                .ok()
                .map(|s| s.lines().map(str::to_owned).collect()),
        )?;
        let file_path = gz_path.to_string_lossy().into_owned();

        Python::initialize();
        let mut parser =
            Gff3Parser::new(None, None, Some("sequential".to_string())).map_err(from_py_error)?;

        parser
            .open(&file_path, Some("sequential".to_string()))
            .map_err(from_py_error)?;
        assert_eq!(parser.reader_mode(), "sequential");

        let records = parser.records(None, None).map_err(from_py_error)?;
        assert_eq!(records.len(), 2);
        assert_first_sample_record(&records[0]);
        assert_second_sample_record(&records[1]);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_sequential_mode_rejects_seek() -> io::Result<()> {
        let dir = create_temp_dir("gff3-seq")?;
        let path = write_sample_file(&dir, "sample.gff3", SAMPLE)?;

        let mut parser = Gff3Parser::new(
            Some(path.to_string_lossy().into_owned()),
            None,
            Some("sequential".to_string()),
        )
        .map_err(from_py_error)?;

        assert_eq!(parser.reader_mode(), "sequential");
        assert!(
            parser
                .seek("chr1:10-20")
                .unwrap_err()
                .to_string()
                .contains("sequential")
        );

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // Gff3Parser:open - invalid mode
    // -----------------------------------------------------------------------

    #[test]
    fn test_gff3_parser_rejects_invalid_mode() {
        Python::initialize();
        let err = match Gff3Parser::new(None, None, Some("random".to_string())) {
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
    // Gff3Parser:close
    // -----------------------------------------------------------------------

    #[test]
    fn test_gff3_parser_close() -> io::Result<()> {
        let dir = create_temp_dir("gff3-parser-close")?;
        let path = write_sample_file(&dir, "sample.gff3", SAMPLE)?;

        Python::initialize();
        let mut parser = Gff3Parser::new(None, None, None).map_err(from_py_error)?;
        parser
            .open(&path.to_string_lossy(), None)
            .map_err(from_py_error)?;

        assert!(parser.reader.is_some());
        parser.close();
        assert!(parser.reader.is_none());
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
                .records(None, None)
                .unwrap_err()
                .to_string()
                .contains("not open")
        );

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // Gff3Parser:parse_text and parse_file
    // -----------------------------------------------------------------------

    #[test]
    fn test_gff3_parser_parse_text() -> io::Result<()> {
        let records = Gff3Parser::parse_text(str::from_utf8(SAMPLE).unwrap(), Some(1))
            .map_err(from_py_error)?;
        assert_eq!(records.len(), 1);
        assert_first_sample_record(&records[0]);

        Ok(())
    }

    #[test]
    fn test_gff3_parser_parse_file() -> io::Result<()> {
        let dir = create_temp_dir("gff3-parse-file")?;
        let path = write_sample_file(&dir, "sample.gff3", SAMPLE)?;

        let records =
            Gff3Parser::parse_file(&path.to_string_lossy(), None).map_err(from_py_error)?;
        assert_eq!(records.len(), 2);
        assert_first_sample_record(&records[0]);
        assert_second_sample_record(&records[1]);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_gff3_parser_parse_file_bgzf() -> io::Result<()> {
        let dir = create_temp_dir("gff3-parse-file-bgzf")?;
        let path = write_bgzf_file(&dir, sample_lines(), Some("sample.gff3.gz"))?;

        let records =
            Gff3Parser::parse_file(&path.to_string_lossy(), None).map_err(from_py_error)?;
        assert_eq!(records.len(), 2);
        assert_first_sample_record(&records[0]);
        assert_second_sample_record(&records[1]);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_gff3_parser_parse_file_missing_path_errors() {
        Python::initialize();
        let err = Gff3Parser::parse_file("/definitely/missing/sample.gff3", None).unwrap_err();
        assert!(
            err.to_string().contains("No such file") || err.to_string().contains("os error 2"),
            "unexpected error: {err}"
        );
    }

    // -----------------------------------------------------------------------
    // SequentialBackend
    // -----------------------------------------------------------------------

    #[test]
    fn test_sequential_backend_next_reads_until_eof() -> io::Result<()> {
        let dir = create_temp_dir("gff3-sequential-next")?;
        let path = write_sample_file(&dir, "sample.gff3", SAMPLE)?;
        let mut backend = SequentialBackend::new(&path)?;

        let first = backend.next()?.expect("expected first record");
        let second = backend.next()?.expect("expected second record");

        assert_first_sample_record(&first);
        assert_second_sample_record(&second);
        assert_eq!(backend.next()?, None);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_sequential_backend_next_reads_bgzf_until_eof() -> io::Result<()> {
        let dir = create_temp_dir("gff3-sequential-bgzf-next")?;
        let path = write_bgzf_file(&dir, sample_lines(), Some("sample.gff3.gz"))?;
        let mut backend = SequentialBackend::new(&path)?;

        let first = backend.next()?.expect("expected first record");
        let second = backend.next()?.expect("expected second record");

        assert_first_sample_record(&first);
        assert_second_sample_record(&second);
        assert_eq!(backend.next()?, None);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_sequential_backend_rejects_seek() -> io::Result<()> {
        let dir = create_temp_dir("gff3-sequential-seek")?;
        let path = write_sample_file(&dir, "sample.gff3", SAMPLE)?;
        let mut backend = SequentialBackend::new(&path)?;

        let err = backend.seek("chr1:10-20").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::Other);
        assert!(err.to_string().contains("sequential"));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // TabixBackend
    // -----------------------------------------------------------------------

    #[test]
    fn test_tabix_backend_seek_reads_records() -> io::Result<()> {
        let dir = create_temp_dir("gff3-tabix-seek")?;
        let (gz_path, _tbi_path) = write_gff3_with_tbi(&dir, None)?;
        let mut backend = TabixBackend::new(&gz_path, "tabix-nonexistent-binary-for-test")?;

        assert!(backend.seek("chr1:150-250")?);

        let first = backend.next()?.expect("expected first record");
        let second = backend.next()?.expect("expected second record");

        assert_eq!(first.raw_attributes, "ID=gene2");
        assert_eq!(second.raw_attributes, "ID=gene3");
        assert_eq!(backend.next()?, None);
        assert!(!backend.query_active);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_tabix_backend_seek_empty_region() -> io::Result<()> {
        let dir = create_temp_dir("gff3-tabix-empty-region")?;
        let (gz_path, _tbi_path) = write_gff3_with_tbi(&dir, None)?;
        let mut backend = TabixBackend::new(&gz_path, "tabix-nonexistent-binary-for-test")?;

        assert!(!backend.seek("chr1:600-700")?);
        assert!(!backend.query_active);
        assert_eq!(backend.next()?, None);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_tabix_backend_seek_invalid_region_errors() -> io::Result<()> {
        let dir = create_temp_dir("gff3-tabix-invalid-region")?;
        let (gz_path, _tbi_path) = write_gff3_with_tbi(&dir, None)?;
        let mut backend = TabixBackend::new(&gz_path, "tabix-nonexistent-binary-for-test")?;

        // "chr1:10-" fails region *syntax* parsing (trailing hyphen, no numeric end) -- distinct
        // from a syntactically valid but unindexed name, which now returns `false` instead.
        let err = backend.seek("chr1:10-").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_tabix_backend_next_without_query_returns_none() -> io::Result<()> {
        let dir = create_temp_dir("gff3-tabix-no-query")?;
        let (gz_path, _tbi_path) = write_gff3_with_tbi(&dir, None)?;
        let mut backend = TabixBackend::new(&gz_path, "tabix-nonexistent-binary-for-test")?;

        assert_eq!(backend.next()?, None);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // Gff3Reader::open / Backend::open
    // -----------------------------------------------------------------------

    #[test]
    fn test_gff3_reader_query_mode_with_index_uses_tabix_backend() -> io::Result<()> {
        let dir = create_temp_dir("gff3-reader-tabix")?;
        let (gz_path, _tbi_path) = write_gff3_with_tbi(&dir, None)?;

        let reader = Gff3Reader::open(
            &gz_path,
            "tabix-nonexistent-binary-for-test",
            ReaderMode::Query,
        )?;

        match reader.backend {
            Backend::Tabix(backend) => {
                assert!(backend.reader.is_some());
                assert!(!backend.query_active);
            }
            Backend::Sequential(_) => panic!("expected tabix backend"),
        }

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_gff3_reader_query_mode_without_index_errors() -> io::Result<()> {
        let dir = create_temp_dir("gff3-reader-no-index")?;
        let path = write_sample_file(&dir, "sample.gff3", SAMPLE)?;

        let err = match Gff3Reader::open(
            &path,
            "tabix-nonexistent-binary-for-test",
            ReaderMode::Query,
        ) {
            Ok(_) => panic!("expected error"),
            Err(err) => err,
        };

        assert_eq!(err.kind(), io::ErrorKind::NotFound);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_gff3_reader_sequential_mode_uses_sequential_backend_with_index() -> io::Result<()> {
        let dir = create_temp_dir("gff3-reader-sequential")?;
        let (gz_path, _tbi_path) = write_gff3_with_tbi(&dir, None)?;

        let reader = Gff3Reader::open(
            &gz_path,
            "tabix-nonexistent-binary-for-test",
            ReaderMode::Sequential,
        )?;
        assert!(matches!(reader.backend, Backend::Sequential(_)));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // GFF3Record
    // -----------------------------------------------------------------------

    #[test]
    fn test_gff3_record_attributes_getter_and_repr() -> io::Result<()> {
        let records = parse_gff3_bytes(SAMPLE, Some(1))?;
        let record = records.first().expect("expected record");

        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let attributes = record.attributes(py)?;

            assert_eq!(attributes.len(), 3);
            assert_eq!(
                attributes
                    .get_item("ID")?
                    .expect("expected ID")
                    .extract::<Vec<String>>()?,
                vec!["gene:1".to_string()]
            );
            assert_eq!(
                attributes
                    .get_item("Alias")?
                    .expect("expected Alias")
                    .extract::<Vec<String>>()?,
                vec!["a".to_string(), "b".to_string()]
            );

            Ok(())
        })
        .map_err(from_py_error)?;

        assert_eq!(
            record.__repr__(),
            "GFF3Record(seqid=\"chr1\", feature_type=\"gene\", start=10, end=20, strand_symbol=\"+\", phase_value=None)"
        );

        Ok(())
    }

    // -----------------------------------------------------------------------
    // GFF3Record::first_attribute
    // -----------------------------------------------------------------------

    #[test]
    fn test_first_attribute_present_single_valued() -> io::Result<()> {
        let records = parse_gff3_bytes(SAMPLE, Some(1))?;
        let record = records.first().expect("expected record");

        assert_eq!(record.first_attribute("ID"), Some("gene:1"));
        assert_eq!(record.first_attribute("Name"), Some("BRCA2"));
        Ok(())
    }

    #[test]
    fn test_first_attribute_present_multi_valued_returns_first() -> io::Result<()> {
        let records = parse_gff3_bytes(SAMPLE, Some(1))?;
        let record = records.first().expect("expected record");

        // "Alias=a,b" — first_attribute returns only the first value, matching
        // attributes[key][0].
        assert_eq!(record.first_attribute("Alias"), Some("a"));
        Ok(())
    }

    #[test]
    fn test_first_attribute_missing_key_returns_none() -> io::Result<()> {
        let records = parse_gff3_bytes(SAMPLE, Some(1))?;
        let record = records.first().expect("expected record");

        assert_eq!(record.first_attribute("Parent"), None);
        Ok(())
    }

    /// Assert `first_attribute` agrees with `attributes[key][0]` from the `dict` getter, for
    /// every key the record has plus the key under test. This is the oracle the parity tests
    /// below check against, so it never hardcodes an expected string itself.
    fn assert_first_attribute_matches_dict(record: &Gff3Record, key: &str) {
        Python::initialize();
        Python::attach(|py| {
            let attributes = record.attributes(py).expect("attributes should build");
            let expected = attributes
                .get_item(key)
                .expect("dict lookup should not error")
                .map(|values| {
                    values
                        .extract::<Vec<String>>()
                        .expect("attribute values should extract as strings")
                        .first()
                        .cloned()
                        .expect("a present key always has at least one value")
                });

            assert_eq!(
                record.first_attribute(key),
                expected.as_deref(),
                "first_attribute({key:?}) disagreed with attributes[{key:?}][0]"
            );
        });
    }

    #[test]
    fn test_first_attribute_duplicate_key_matches_last_occurrence() -> io::Result<()> {
        // GFF3 lets a key repeat as separate pairs (unlike GTF, which merges repeats into one
        // array). The `attributes` dict getter overwrites on `set_item`, so the last occurrence
        // wins there; `first_attribute` must agree rather than returning the first occurrence.
        let records = parse_gff3_bytes(
            b"chr1\tensembl\tgene\t10\t20\t.\t+\t.\tID=a;Name=x%3By,z;ID=b\n",
            Some(1),
        )?;
        let record = records.first().expect("expected record");

        assert_eq!(record.first_attribute("ID"), Some("b"));
        assert_first_attribute_matches_dict(record, "ID");
        Ok(())
    }

    #[test]
    fn test_first_attribute_percent_encoded_value_matches_dict() -> io::Result<()> {
        let records = parse_gff3_bytes(
            b"chr1\tensembl\tgene\t10\t20\t.\t+\t.\tName=BRCA%202\n",
            Some(1),
        )?;
        let record = records.first().expect("expected record");

        assert_eq!(record.first_attribute("Name"), Some("BRCA 2"));
        assert_first_attribute_matches_dict(record, "Name");
        Ok(())
    }

    #[test]
    fn test_first_attribute_empty_value_matches_dict() -> io::Result<()> {
        let records = parse_gff3_bytes(b"chr1\tensembl\tgene\t10\t20\t.\t+\t.\tNote=\n", Some(1))?;
        let record = records.first().expect("expected record");

        assert_first_attribute_matches_dict(record, "Note");
        Ok(())
    }

    #[test]
    fn test_first_attribute_empty_array_element_matches_dict() -> io::Result<()> {
        let records = parse_gff3_bytes(b"chr1\tensembl\tgene\t10\t20\t.\t+\t.\tX=,\n", Some(1))?;
        let record = records.first().expect("expected record");

        assert_first_attribute_matches_dict(record, "X");
        Ok(())
    }

    #[test]
    fn test_first_attribute_comma_array_matches_dict() -> io::Result<()> {
        let records = parse_gff3_bytes(SAMPLE, Some(1))?;
        let record = records.first().expect("expected record");

        // "Alias=a,b" from SAMPLE.
        assert_first_attribute_matches_dict(record, "Alias");
        Ok(())
    }

    #[test]
    fn test_first_attribute_crlf_line_matches_dict() -> io::Result<()> {
        let records = parse_gff3_bytes(
            b"chr1\tensembl\tgene\t10\t20\t.\t+\t.\tID=gene:1\r\n",
            Some(1),
        )?;
        let record = records.first().expect("expected record");

        assert_eq!(record.first_attribute("ID"), Some("gene:1"));
        assert_first_attribute_matches_dict(record, "ID");
        Ok(())
    }

    #[test]
    fn test_first_attribute_after_skipped_comment_and_directive_lines_matches_dict()
    -> io::Result<()> {
        // Comment and directive lines are skipped by `read_next_record` without producing a
        // record; this pins that skipping them first does not leave stale bytes in the reused
        // `Line` buffer that could corrupt the next real record's attributes.
        let records = parse_gff3_bytes(SKIPPED_LINE_SAMPLE, None)?;
        let record = records.first().expect("expected record");

        assert_eq!(record.first_attribute("ID"), Some("gene:1"));
        assert_first_attribute_matches_dict(record, "ID");
        Ok(())
    }

    #[test]
    fn test_first_attribute_before_fasta_directive_matches_dict() -> io::Result<()> {
        // SAMPLE ends with a ##FASTA directive followed by sequence lines; parsing must stop
        // there without disturbing the attributes of records read beforehand.
        let records = parse_gff3_bytes(SAMPLE, None)?;
        let record = records.last().expect("expected record");

        assert_eq!(record.first_attribute("ID"), Some("cds:1"));
        assert_first_attribute_matches_dict(record, "ID");
        Ok(())
    }

    // -----------------------------------------------------------------------
    // parse_gff3_line_into — reused buffer parity with parse_gff3_line
    // -----------------------------------------------------------------------

    #[test]
    fn test_parse_gff3_line_into_reused_buffer_matches_parse_gff3_line() -> io::Result<()> {
        let mut buffer = gff3::Line::default();

        let first = parse_gff3_line_into(
            "chr1\tensembl\tgene\t10\t20\t.\t+\t.\tID=gene:1;Name=BRCA2;Alias=a,b",
            &mut buffer,
            None,
        )?
        .expect("expected first record");
        assert_first_sample_record(&first);

        let second = parse_gff3_line_into(
            "chr1\tensembl\tCDS\t12\t18\t4.5\t-\t2\tID=cds:1;Parent=gene:1",
            &mut buffer,
            None,
        )?
        .expect("expected second record");
        assert_second_sample_record(&second);

        Ok(())
    }

    /// A residual `\r` left behind by `TabixReader`'s single-pass `trim_newline` (a doubled `\r`
    /// before the real line terminator, or a bare trailing `\r` on a file's last, newline-less
    /// line) must parse the same way whether the line arrived through the tabix/query path or
    /// the sequential reader. See the note on `parse_gff3_line_into`.
    #[test]
    fn test_parse_gff3_line_into_residual_cr_matches_sequential_reader() -> io::Result<()> {
        let raw_line = "chr1\tensembl\tgene\t10\t20\t.\t+\t.\tID=gene:1\r\r\n";

        // Tabix/query path: TabixReader::next() would already have run trim_newline over the
        // raw bytes before handing this string to parse_gff3_line_into.
        let mut trimmed = raw_line.to_string();
        trim_newline(&mut trimmed);
        let mut buffer = gff3::Line::default();
        let tabix_path_record =
            parse_gff3_line_into(&trimmed, &mut buffer, None)?.expect("expected record");

        // Sequential path: SequentialBackend reads the same raw bytes directly.
        let sequential_records = parse_gff3_bytes(raw_line.as_bytes(), None)?;
        let sequential_record = sequential_records.first().expect("expected record");

        assert_eq!(
            tabix_path_record.raw_attributes,
            sequential_record.raw_attributes
        );
        assert_eq!(
            tabix_path_record.first_attribute("ID"),
            sequential_record.first_attribute("ID")
        );
        // Pins the actual (slightly surprising) shared behaviour: with only one `\r` stripped by
        // either path, the other survives into the attribute value on both.
        assert_eq!(tabix_path_record.first_attribute("ID"), Some("gene:1\r"));

        Ok(())
    }

    /// The same parity, for a file's last line ending in a bare `\r` with no `\n` at all.
    #[test]
    fn test_parse_gff3_line_into_trailing_cr_without_newline_matches_sequential_reader()
    -> io::Result<()> {
        let raw_line = "chr1\tensembl\tgene\t10\t20\t.\t+\t.\tID=gene:1\r";

        let mut trimmed = raw_line.to_string();
        trim_newline(&mut trimmed);
        let mut buffer = gff3::Line::default();
        let tabix_path_record =
            parse_gff3_line_into(&trimmed, &mut buffer, None)?.expect("expected record");

        let sequential_records = parse_gff3_bytes(raw_line.as_bytes(), None)?;
        let sequential_record = sequential_records.first().expect("expected record");

        assert_eq!(
            tabix_path_record.first_attribute("ID"),
            sequential_record.first_attribute("ID")
        );
        assert_eq!(tabix_path_record.first_attribute("ID"), Some("gene:1\r"));

        Ok(())
    }

    // -----------------------------------------------------------------------
    // TabixBackend::index_sequence_names / Gff3Parser::index_sequence_names
    // -----------------------------------------------------------------------

    #[test]
    fn test_tabix_backend_index_sequence_names() -> io::Result<()> {
        let dir = create_temp_dir("gff3-tabix-index-names")?;
        let (gz_path, _tbi_path) = write_gff3_with_tbi(&dir, None)?;
        let mut backend = TabixBackend::new(&gz_path, "tabix-nonexistent-binary-for-test")?;

        assert_eq!(
            backend.index_sequence_names()?,
            vec!["chr1".to_string(), "chr2".to_string()]
        );

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_sequential_backend_index_sequence_names_is_empty() -> io::Result<()> {
        let dir = create_temp_dir("gff3-sequential-index-names")?;
        let path = write_sample_file(&dir, "sample.gff3", SAMPLE)?;

        let mut reader = Gff3Reader::open(
            &path,
            "tabix-nonexistent-binary-for-test",
            ReaderMode::Sequential,
        )?;
        assert!(reader.index_sequence_names()?.is_empty());

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_gff3_parser_index_sequence_names() -> io::Result<()> {
        let dir = create_temp_dir("gff3-parser-index-names")?;
        let (gz_path, _tbi_path) = write_gff3_with_tbi(&dir, None)?;
        let file_path = gz_path.to_string_lossy().into_owned();

        let mut parser = Gff3Parser::new(None, None, None).map_err(from_py_error)?;
        assert!(
            parser
                .index_sequence_names()
                .map_err(from_py_error)?
                .is_empty()
        );

        parser.open(&file_path, None).map_err(from_py_error)?;
        assert_eq!(
            parser.index_sequence_names().map_err(from_py_error)?,
            vec!["chr1".to_string(), "chr2".to_string()]
        );

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_gff3_parser_seek_unknown_contig_returns_false() -> io::Result<()> {
        let dir = create_temp_dir("gff3-parser-unknown-contig")?;
        let (gz_path, _tbi_path) = write_gff3_with_tbi(&dir, None)?;
        let file_path = gz_path.to_string_lossy().into_owned();

        let mut parser = Gff3Parser::new(None, None, None).map_err(from_py_error)?;
        parser.open(&file_path, None).map_err(from_py_error)?;

        assert!(!parser.seek("chrZ:1-100").map_err(from_py_error)?);
        assert_eq!(parser.next().map_err(from_py_error)?, None);

        fs::remove_dir_all(dir)?;
        Ok(())
    }
}
