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
    parsers::tabix::{self, associated_index_path},
    parsers::utils::{is_bgzf_file, ReaderMode},
    utils::to_py_error,
};
use noodles_bgzf as bgzf;
use noodles_gff::{
    self as gff,
    feature::record::{Phase, Strand},
    record::attributes::field::Value as AttributeValue,
};
use pyo3::{
    exceptions::PyRuntimeError,
    prelude::*,
    types::{PyDict, PyTuple},
};
use std::{
    fs::File,
    io::{self, BufRead, BufReader},
    path::{Path, PathBuf},
};

type AttributePairs = Vec<(String, Vec<String>)>;

// ---------------------------------------------------------------------------
// GFFRecord - represents a GFF record
// ---------------------------------------------------------------------------

#[pyclass(name = "GFFRecord")]
#[derive(Clone, Debug, PartialEq)]
pub struct GffRecord {
    #[pyo3(get)]
    pub seqid: String,
    #[pyo3(get)]
    pub source: String,
    #[pyo3(get)]
    pub feature_type: String,
    #[pyo3(get)]
    pub start: usize,
    #[pyo3(get)]
    pub end: usize,
    #[pyo3(get)]
    pub score: Option<f64>,
    #[pyo3(get)]
    pub strand_symbol: String,
    #[pyo3(get)]
    pub phase_value: Option<u8>,
    attributes: AttributePairs,
    #[pyo3(get)]
    pub raw_attributes: String,
}

#[pymethods]
impl GffRecord {
    #[getter]
    fn attributes<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);

        for (key, values) in &self.attributes {
            dict.set_item(key, PyTuple::new(py, values)?)?;
        }

        Ok(dict)
    }

    fn __repr__(&self) -> String {
        format!(
            "GFFRecord(seqid={:?}, feature_type={:?}, start={}, end={}, strand_symbol={:?}, phase_value={:?})",
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
// GffReader - public entry point
// ---------------------------------------------------------------------------

// A reusable GFF reader that supports both sequential and indexed (tabix) reading modes.

// The reader by default uses indexed reading if a tabix index is available,
// but can be forced into sequential mode via the `mode` parameter on `GffParser`.
pub struct GffReader {
    file_path: PathBuf,
    backend: Backend,
    mode: ReaderMode,
    metadata: Option<GffMetadata>,
}

impl GffReader {
    fn open(file_path: &Path, tabix_path: &str, mode: ReaderMode) -> io::Result<Self> {
        let (backend, mode) = Backend::open(file_path, tabix_path, mode)?;
        Ok(Self {
            file_path: file_path.to_path_buf(),
            backend,
            mode,
            metadata: None,
        })
    }

    fn next(&mut self) -> io::Result<Option<GffRecord>> {
        self.backend.next()
    }

    fn seek(&mut self, region: &str) -> io::Result<bool> {
        self.backend.seek(region)
    }

    fn metadata(&mut self) -> io::Result<&GffMetadata> {
        if self.metadata.is_none() {
            self.metadata = Some(read_gff_metadata(&self.file_path)?);
        }

        Ok(self.metadata.as_ref().expect("metadata should be cached"))
    }
}

// ---------------------------------------------------------------------------
// GffMetadata - holds GFF metadata (directives and header lines)
// ---------------------------------------------------------------------------

#[derive(Clone, Debug, Default, PartialEq)]
struct GffMetadata {
    header: Vec<String>,
    directives: Vec<(String, String)>,
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

    fn next(&mut self) -> io::Result<Option<GffRecord>> {
        match self {
            Self::Sequential(reader) => reader.next(),
            Self::Tabix(reader) => reader.next(),
        }
    }

    fn seek(&mut self, region: &str) -> io::Result<bool> {
        match self {
            Self::Sequential(reader) => reader.seek(region),
            Self::Tabix(reader) => reader.seek(region),
        }
    }
}

// ---------------------------------------------------------------------------
// SequentialBackend - simple line-by-line reader
// ---------------------------------------------------------------------------

struct SequentialBackend {
    reader: gff::io::Reader<Box<dyn BufRead + Send + Sync>>,
}

impl SequentialBackend {
    fn new<P>(file_path: P) -> io::Result<Self>
    where
        P: AsRef<Path>,
    {
        Ok(Self {
            reader: open_sequential_gff_reader(file_path.as_ref())?,
        })
    }

    fn seek(&mut self, _region: &str) -> io::Result<bool> {
        Err(io::Error::other(
            "GFF reader is set to 'sequential'; set reader to 'query' before seek",
        ))
    }

    fn next(&mut self) -> io::Result<Option<GffRecord>> {
        read_next_record(&mut self.reader)
    }
}

// ---------------------------------------------------------------------------
// TabixBackend - uses tabix for indexed reading
// ---------------------------------------------------------------------------

struct TabixBackend {
    file_path: PathBuf,
    tabix_path: String,
    reader: Option<tabix::TabixReader>,
    query_active: bool,
}

impl TabixBackend {
    fn new<P>(file_path: P, tabix_path: &str) -> io::Result<Self>
    where
        P: AsRef<Path>,
    {
        Ok(Self {
            file_path: file_path.as_ref().to_path_buf(),
            tabix_path: tabix_path.to_string(),
            reader: None,
            query_active: false,
        })
    }

    fn seek(&mut self, region: &str) -> io::Result<bool> {
        if self.reader.is_none() {
            self.reader = Some(tabix::TabixReader::open(&self.file_path, &self.tabix_path)?);
        }

        let has_records = self
            .reader
            .as_mut()
            .expect("tabix reader should be initialized")
            .seek(region)?;

        self.query_active = has_records;
        Ok(has_records)
    }

    fn next(&mut self) -> io::Result<Option<GffRecord>> {
        if self.query_active {
            let Some(reader) = self.reader.as_mut() else {
                return Err(io::Error::other("tabix reader is not initialized"));
            };

            loop {
                let Some(line) = reader.next()? else {
                    self.query_active = false;
                    return Ok(None);
                };

                if let Some(record) = parse_gff_line(&line)? {
                    return Ok(Some(record));
                }
            }
        }

        Ok(None)
    }
}

// ---------------------------------------------------------------------------
// PyO3 wrapper
// ---------------------------------------------------------------------------

#[pyclass(name = "GFFParser")]
pub struct GffParser {
    file_path: Option<PathBuf>,
    reader: Option<GffReader>,
    tabix_path: String,
    reader_mode: ReaderMode,
}

#[pymethods]
impl GffParser {
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
            parser.open(&path)?;
        }

        Ok(parser)
    }

    #[getter]
    fn format(&self) -> &'static str {
        "gff"
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

    fn open(&mut self, file_path: &str) -> PyResult<()> {
        let path = PathBuf::from(file_path);

        self.reader =
            Some(GffReader::open(&path, &self.tabix_path, self.reader_mode).map_err(to_py_error)?);
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
    fn next(&mut self) -> PyResult<Option<GffRecord>> {
        self.reader()?.next().map_err(to_py_error)
    }

    #[pyo3(signature = (limit=10000))]
    fn records(&mut self, limit: Option<usize>) -> PyResult<Vec<GffRecord>> {
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
    fn parse_text(data: &str, limit: Option<usize>) -> PyResult<Vec<GffRecord>> {
        parse_gff_bytes(data.as_bytes(), limit).map_err(to_py_error)
    }

    #[staticmethod]
    #[pyo3(signature = (file_path, limit=None))]
    fn parse_file(file_path: &str, limit: Option<usize>) -> PyResult<Vec<GffRecord>> {
        let mut reader = open_sequential_gff_reader(Path::new(file_path)).map_err(to_py_error)?;

        read_records(&mut reader, limit).map_err(to_py_error)
    }
}

impl GffParser {
    fn reader(&mut self) -> PyResult<&mut GffReader> {
        self.reader
            .as_mut()
            .ok_or_else(|| PyRuntimeError::new_err("GFF parser is not open"))
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn has_tabix_index(file_path: &Path) -> bool {
    associated_index_path(file_path, "tbi").exists()
        || associated_index_path(file_path, "csi").exists()
}

fn open_sequential_gff_reader(
    file_path: &Path,
) -> io::Result<gff::io::Reader<Box<dyn BufRead + Send + Sync>>> {
    let file = File::open(file_path)?;
    let reader: Box<dyn BufRead + Send + Sync> = if is_bgzf_file(file_path)? {
        Box::new(bgzf::io::Reader::new(file))
    } else {
        Box::new(BufReader::new(file))
    };

    Ok(gff::io::Reader::new(reader))
}

fn read_gff_metadata(file_path: &Path) -> io::Result<GffMetadata> {
    let mut reader = open_sequential_gff_reader(file_path)?;
    read_gff_metadata_from_reader(&mut reader)
}

fn read_gff_metadata_from_reader<R>(reader: &mut gff::io::Reader<R>) -> io::Result<GffMetadata>
where
    R: BufRead,
{
    let mut metadata = GffMetadata::default();
    let mut line = gff::Line::default();

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

fn read_next_record<R>(reader: &mut gff::io::Reader<R>) -> io::Result<Option<GffRecord>>
where
    R: BufRead,
{
    let mut line = gff::Line::default();

    loop {
        if reader.read_line(&mut line)? == 0 {
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

        let Some(result) = line.as_record() else {
            continue;
        };

        let raw_attributes = raw_attributes(&line);
        let record = result?;

        return record_to_py(record, raw_attributes).map(Some);
    }
}

fn read_records<R>(
    reader: &mut gff::io::Reader<R>,
    limit: Option<usize>,
) -> io::Result<Vec<GffRecord>>
where
    R: BufRead,
{
    let mut records = Vec::new();

    while limit.is_none_or(|n| records.len() < n) {
        match read_next_record(reader)? {
            Some(record) => records.push(record),
            None => break,
        }
    }

    Ok(records)
}

fn parse_gff_bytes(data: &[u8], limit: Option<usize>) -> io::Result<Vec<GffRecord>> {
    let mut reader = gff::io::Reader::new(data);
    read_records(&mut reader, limit)
}

fn parse_gff_line(line: &str) -> io::Result<Option<GffRecord>> {
    let mut data = Vec::with_capacity(line.len() + 1);
    data.extend_from_slice(line.as_bytes());
    data.push(b'\n');

    Ok(parse_gff_bytes(&data, Some(1))?.pop())
}

fn record_to_py(record: gff::Record<'_>, raw_attributes: String) -> io::Result<GffRecord> {
    Ok(GffRecord {
        seqid: bytes_to_string(record.reference_sequence_name()),
        source: bytes_to_string(record.source()),
        feature_type: bytes_to_string(record.ty()),
        start: record.start()?.get(),
        end: record.end()?.get(),
        score: record.score().transpose()?.map(f64::from),
        strand_symbol: strand_to_string(record.strand()?),
        phase_value: record.phase().transpose()?.map(phase_to_u8),
        attributes: parse_attributes(record.attributes())?,
        raw_attributes,
    })
}

fn parse_attributes(attributes: gff::record::Attributes<'_>) -> io::Result<AttributePairs> {
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

fn raw_attributes(line: &gff::Line) -> String {
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

fn bytes_to_string(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

fn strand_to_string(strand: Strand) -> String {
    match strand {
        Strand::None => ".",
        Strand::Forward => "+",
        Strand::Reverse => "-",
        Strand::Unknown => "?",
    }
    .to_string()
}

fn phase_to_u8(phase: Phase) -> u8 {
    match phase {
        Phase::Zero => 0,
        Phase::One => 1,
        Phase::Two => 2,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        parsers::test_utils::{write_bgzf_file, write_gff3_with_tbi},
        test_utils::{create_temp_dir, write_sample_file},
        utils::from_py_error,
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

    fn to_gff_line(s: &str) -> gff::Line {
        let mut reader = gff::io::Reader::new(s.as_bytes());
        let mut line = gff::Line::default();
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

    fn assert_first_sample_record(record: &GffRecord) {
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

    fn assert_second_sample_record(record: &GffRecord) {
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
    // phase_to_u8
    // -----------------------------------------------------------------------

    #[test]
    fn test_phase_to_u8() {
        assert_eq!(phase_to_u8(Phase::Zero), 0);
        assert_eq!(phase_to_u8(Phase::One), 1);
        assert_eq!(phase_to_u8(Phase::Two), 2);
    }

    // -----------------------------------------------------------------------
    // strand_to_string
    // -----------------------------------------------------------------------

    #[test]
    fn test_strand_to_string() {
        assert_eq!(strand_to_string(Strand::None), ".");
        assert_eq!(strand_to_string(Strand::Forward), "+");
        assert_eq!(strand_to_string(Strand::Reverse), "-");
        assert_eq!(strand_to_string(Strand::Unknown), "?");
    }

    // -----------------------------------------------------------------------
    // bytes_to_string
    // -----------------------------------------------------------------------

    #[test]
    fn test_bytes_to_string_valid_utf8() {
        assert_eq!(bytes_to_string(b"abc"), "abc");
    }

    #[test]
    fn test_bytes_to_string_replaces_invalid_utf8() {
        assert_eq!(
            bytes_to_string(b"a\xffb"),
            format!("a{}b", char::REPLACEMENT_CHARACTER)
        );
    }

    // -----------------------------------------------------------------------
    // raw_attributes
    // -----------------------------------------------------------------------

    #[test]
    fn test_raw_attributes_with_attributes() -> io::Result<()> {
        let line =
            to_gff_line("chr1\tensembl\tgene\t10\t20\t.\t+\t.\tID=gene:1;Name=BRCA2;Alias=a,b");
        assert_eq!(raw_attributes(&line), "ID=gene:1;Name=BRCA2;Alias=a,b");
        Ok(())
    }

    #[test]
    fn test_raw_attributes_missing_attribute_field_returns_empty() -> io::Result<()> {
        let line = to_gff_line("chr1\tensembl\tgene\t10\t20\t.\t+\t.");
        assert_eq!(raw_attributes(&line), "");
        Ok(())
    }

    // -----------------------------------------------------------------------
    // parse_attributes
    // -----------------------------------------------------------------------

    #[test]
    fn test_parse_attributes_single_valued() -> io::Result<()> {
        let line = to_gff_line("chr1\tensembl\tgene\t10\t20\t.\t+\t.\tID=gene:1;Name=BRCA2");
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
        let line = to_gff_line("chr1\tensembl\tgene\t10\t20\t.\t+\t.\tAlias=a,b,c");
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
        let line = to_gff_line("chr1\tensembl\tgene\t10\t20\t.\t+\t.\tName=BRCA%202");
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
        let line = to_gff_line("chr1\tensembl\tgene\t10\t20\t.\t+\t.\tID=gene:1;Invalid");
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
        let line = to_gff_line("chr1\tensembl\tgene\t10\t20\t.\t+\t.\t");
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
            to_gff_line("chr1\tensembl\tgene\t10\t20\t.\t+\t.\tID=gene:1;Name=BRCA2;Alias=a,b");
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
    // parse_gff_line
    // -----------------------------------------------------------------------

    #[test]
    fn test_parse_gff_line_decodes_record() -> io::Result<()> {
        let record = parse_gff_line(
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
    fn test_parse_gff_line_invalid_line_bad_fields() -> io::Result<()> {
        let err =
            parse_gff_line("chr1\tensembl\tgene\tbad_start\t20\t.\t+\t.\tID=gene:1").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);

        let err =
            parse_gff_line("chr1\tensembl\tgene\t10\tbad_end\t.\t+\t.\tID=gene:1").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);

        let err =
            parse_gff_line("chr1\tensembl\tgene\t10\t20\tbad_score\t+\t.\tID=gene:1").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);

        let err =
            parse_gff_line("chr1\tensembl\tgene\t10\t20\t.\tbad_strand\t.\tID=gene:1").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);

        let err =
            parse_gff_line("chr1\tensembl\tgene\t10\t20\t.\t+\tbad_phase\tID=gene:1").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);

        let err =
            parse_gff_line("chr1\tensembl\tgene\t10\t20\t.\t+\t.\tInvalidAttribute").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        Ok(())
    }

    #[test]
    fn test_parse_gff_line_invalid_line_missing_fields() -> io::Result<()> {
        let err = parse_gff_line("chr1\tensembl\tgene\t10\t20\t.\t+\t.").unwrap_err();

        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
        Ok(())
    }

    #[test]
    fn test_parse_gff_line_ignores_non_record_lines() -> io::Result<()> {
        assert_eq!(parse_gff_line("")?, None);
        assert_eq!(parse_gff_line("##gff-version 3")?, None);
        assert_eq!(parse_gff_line("# a regular comment")?, None);

        Ok(())
    }

    // -----------------------------------------------------------------------
    // parse_gff_bytes
    // -----------------------------------------------------------------------

    #[test]
    fn test_parse_gff_bytes_decodes_records() -> io::Result<()> {
        let records = parse_gff_bytes(SAMPLE, None)?;

        assert_eq!(records.len(), 2);
        assert_first_sample_record(&records[0]);
        assert_second_sample_record(&records[1]);

        Ok(())
    }

    #[test]
    fn test_parse_gff_bytes_honours_limit() -> io::Result<()> {
        let records = parse_gff_bytes(SAMPLE, Some(1))?;

        assert_eq!(records.len(), 1);
        Ok(())
    }

    #[test]
    fn test_parse_gff_bytes_honours_zero_limit() -> io::Result<()> {
        let records = parse_gff_bytes(SAMPLE, Some(0))?;

        assert!(records.is_empty());
        Ok(())
    }

    #[test]
    fn test_parse_gff_bytes_skips_empty_directive_and_comment_lines() -> io::Result<()> {
        let records = parse_gff_bytes(SKIPPED_LINE_SAMPLE, None)?;

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].strand_symbol, ".");
        assert_eq!(records[0].raw_attributes, "ID=gene:1");

        Ok(())
    }

    #[test]
    fn test_parse_gff_bytes_stops_at_fasta() -> io::Result<()> {
        let records = parse_gff_bytes(SAMPLE, Some(10))?;

        assert_eq!(records.len(), 2);
        Ok(())
    }

    #[test]
    fn test_parse_gff_bytes_invalid_record_errors() {
        let data = b"chr1\tensembl\tgene\tbad_start\t20\t.\t+\t.\tID=gene:1\n";
        let err = parse_gff_bytes(data, None).unwrap_err();

        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    }

    // -----------------------------------------------------------------------
    // read_records
    // -----------------------------------------------------------------------

    #[test]
    fn test_read_records_reads_multiple_records() -> io::Result<()> {
        let mut reader = gff::io::Reader::new(SAMPLE);
        let records = read_records(&mut reader, None)?;

        assert_eq!(records.len(), 2);
        assert_first_sample_record(&records[0]);
        assert_second_sample_record(&records[1]);
        Ok(())
    }

    #[test]
    fn test_read_records_honours_limit() -> io::Result<()> {
        let mut reader = gff::io::Reader::new(SAMPLE);
        let records = read_records(&mut reader, Some(1))?;

        assert_eq!(records.len(), 1);
        Ok(())
    }

    #[test]
    fn test_read_records_honours_zero_limit() -> io::Result<()> {
        let mut reader = gff::io::Reader::new(SAMPLE);
        let records = read_records(&mut reader, Some(0))?;

        assert!(records.is_empty());
        Ok(())
    }

    // -----------------------------------------------------------------------
    // read_next_record
    // -----------------------------------------------------------------------

    #[test]
    fn test_read_next_record_reads_records_one_by_one() -> io::Result<()> {
        let mut reader = gff::io::Reader::new(SAMPLE);

        let first = read_next_record(&mut reader)?.expect("expected first record");
        assert_first_sample_record(&first);

        let second = read_next_record(&mut reader)?.expect("expected second record");
        assert_second_sample_record(&second);
        Ok(())
    }

    #[test]
    fn test_read_next_record_returns_none_at_eof() -> io::Result<()> {
        let mut reader = gff::io::Reader::new(SAMPLE);
        read_next_record(&mut reader)?; // first record
        read_next_record(&mut reader)?; // second record

        assert_eq!(read_next_record(&mut reader)?, None);
        Ok(())
    }

    #[test]
    fn test_read_next_record_no_records() -> io::Result<()> {
        let mut reader =
            gff::io::Reader::new(b"##gff-version 3\n##FASTA\n>chr1\nACGT\n".as_slice());

        assert_eq!(read_next_record(&mut reader)?, None);
        Ok(())
    }

    // -----------------------------------------------------------------------
    // read_gff_metadata
    // -----------------------------------------------------------------------

    #[test]
    fn test_read_gff_metadata_from_reader_collects_header_and_directives() -> io::Result<()> {
        let mut reader = gff::io::Reader::new(SKIPPED_LINE_SAMPLE);
        let metadata = read_gff_metadata_from_reader(&mut reader)?;

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
    fn test_read_gff_metadata_with_no_header() -> io::Result<()> {
        let mut reader = gff::io::Reader::new(b"chr1\tensembl\tgene\t10\t20\t.\t+\t.\tID=gene:1\n".as_slice());
        let metadata = read_gff_metadata_from_reader(&mut reader)?;

        assert!(metadata.header.is_empty());
        assert!(metadata.directives.is_empty());
        Ok(())
    }

    // -----------------------------------------------------------------------
    // open_sequential_gff_reader
    // -----------------------------------------------------------------------

    #[test]
    fn test_open_sequential_gff_reader_reads_plain_file() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-gff-open-sequential-plain")?;
        let path = write_sample_file(&dir, "sample.gff3", SAMPLE)?;
        let mut reader = open_sequential_gff_reader(&path)?;
        let records = read_records(&mut reader, None)?;

        assert_eq!(records.len(), 2);
        assert_first_sample_record(&records[0]);
        assert_second_sample_record(&records[1]);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_open_sequential_gff_reader_reads_bgzf_file() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-gff-open-sequential-bgzf")?;
        let path = write_bgzf_file(&dir, sample_lines(), Some("sample.gff3.gz"))?;
        let mut reader = open_sequential_gff_reader(&path)?;
        let records = read_records(&mut reader, None)?;

        assert_eq!(records.len(), 2);
        assert_first_sample_record(&records[0]);
        assert_second_sample_record(&records[1]);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // has_tabix_index
    // -----------------------------------------------------------------------

    #[test]
    fn test_has_tabix_index_exists() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-gff-has-tabix")?;
        let (gz_path, _tbi_path) = write_gff3_with_tbi(
            &dir,
            str::from_utf8(SAMPLE)
                .ok()
                .map(|s| s.lines().map(str::to_owned).collect()),
        )?;

        assert!(has_tabix_index(&gz_path));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_has_tabix_index_missing_index() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-gff-has-tabix-missing")?;
        let path = write_sample_file(&dir, "sample.gff3", SAMPLE)?;

        assert!(!has_tabix_index(&path));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // GffParser:open
    // -----------------------------------------------------------------------

    #[test]
    fn test_gff_parser_without_open_reader() -> io::Result<()> {
        let mut parser = GffParser::new(None, None, None).map_err(from_py_error)?;

        assert_eq!(parser.format(), "gff");
        Python::initialize();
        Python::attach(|py| {
            assert!(parser.metadata(py).unwrap().is_empty());
        });
        assert_eq!(parser.reader_mode(), "auto");
        assert!(parser.next().unwrap_err().to_string().contains("not open"));
        assert!(parser
            .seek("chr1:10-20")
            .unwrap_err()
            .to_string()
            .contains("not open"));
        assert!(parser
            .records(None)
            .unwrap_err()
            .to_string()
            .contains("not open"));

        parser.close();
        Ok(())
    }

    #[test]
    fn test_gff_parser_metadata_sequential() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-gff-parser-metadata-seq")?;
        let path = write_sample_file(&dir, "sample.gff3", SKIPPED_LINE_SAMPLE)?;
        let file_path = path.to_string_lossy().into_owned();

        let mut parser =
            GffParser::new(None, None, Some("sequential".to_string())).map_err(from_py_error)?;
        parser.open(&file_path).map_err(from_py_error)?;

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
    fn test_gff_parser_metadata_does_not_consume_sequential_reader() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-gff-parser-metadata-preserve-seq")?;
        let path = write_sample_file(&dir, "sample.gff3", SAMPLE)?;
        let file_path = path.to_string_lossy().into_owned();

        let mut parser =
            GffParser::new(None, None, Some("sequential".to_string())).map_err(from_py_error)?;
        parser.open(&file_path).map_err(from_py_error)?;

        let first = parser
            .next()
            .map_err(from_py_error)?
            .expect("expected first record");
        assert_first_sample_record(&first);

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
    fn test_gff_parser_metadata_query_bgzf() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-gff-parser-metadata-query")?;
        let (gz_path, _tbi_path) = write_gff3_with_tbi(
            &dir,
            str::from_utf8(SKIPPED_LINE_SAMPLE)
                .ok()
                .map(|s| s.lines().map(str::to_owned).collect()),
        )?;
        let file_path = gz_path.to_string_lossy().into_owned();

        let mut parser = GffParser::new(None, None, None).map_err(from_py_error)?;
        parser.open(&file_path).map_err(from_py_error)?;

        assert_eq!(parser.reader_mode(), "query");

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
    fn test_gff_parser_open_missing_file_errors() -> io::Result<()> {
        Python::initialize();
        let mut parser = GffParser::new(None, None, None).map_err(from_py_error)?;
        let err = parser.open("/definitely/missing/sample.gff3").unwrap_err();

        assert!(
            err.to_string().contains("No such file") || err.to_string().contains("os error 2"),
            "unexpected error: {err}"
        );
        Ok(())
    }

    // -----------------------------------------------------------------------
    // GffParser:open - auto reader modes
    // -----------------------------------------------------------------------

    #[test]
    fn test_gff_parser_open_file_auto_mode_sequential() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-gff-parser-open-auto")?;
        let path = write_sample_file(&dir, "sample.gff3", SAMPLE)?;
        let file_path = path.to_string_lossy().into_owned();

        let mut parser = GffParser::new(None, None, None).map_err(from_py_error)?;
        parser.open(&file_path).map_err(from_py_error)?;

        assert_eq!(parser.reader_mode(), "sequential");

        let records = parser.records(None).map_err(from_py_error)?;
        assert_eq!(records.len(), 2);
        assert_first_sample_record(&records[0]);
        assert_second_sample_record(&records[1]);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_gff_parser_open_file_auto_mode_query() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-gff-parser-open-auto-query")?;
        let (gz_path, _tbi_path) = write_gff3_with_tbi(
            &dir,
            str::from_utf8(SAMPLE)
                .ok()
                .map(|s| s.lines().map(str::to_owned).collect()),
        )?;
        let file_path = gz_path.to_string_lossy().into_owned();

        let mut parser = GffParser::new(None, None, None).map_err(from_py_error)?;
        parser.open(&file_path).map_err(from_py_error)?;

        assert_eq!(parser.reader_mode(), "query");
        parser.seek("chr1:10-20").map_err(from_py_error)?;

        let records = parser.records(None).map_err(from_py_error)?;
        assert_eq!(records.len(), 2);
        assert_first_sample_record(&records[0]);
        assert_second_sample_record(&records[1]);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // GffParser:open - explicit query mode
    // -----------------------------------------------------------------------

    #[test]
    fn test_gff_parser_open_file_query_mode_without_index_errors() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-gff-parser-open-query-no-index")?;
        let path = write_sample_file(&dir, "sample.gff3", SAMPLE)?;
        let file_path = path.to_string_lossy().into_owned();

        Python::initialize();
        let mut parser =
            GffParser::new(None, None, Some("query".to_string())).map_err(from_py_error)?;

        let err = parser.open(&file_path).unwrap_err();
        assert!(
            err.to_string().contains("no tabix index found"),
            "unexpected error: {err}"
        );

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_gff_parser_open_file_query_mode_with_index() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-gff-parser-open-query-index")?;
        let (gz_path, _tbi_path) = write_gff3_with_tbi(
            &dir,
            str::from_utf8(SAMPLE)
                .ok()
                .map(|s| s.lines().map(str::to_owned).collect()),
        )?;
        let file_path = gz_path.to_string_lossy().into_owned();

        let mut parser =
            GffParser::new(None, None, Some("query".to_string())).map_err(from_py_error)?;
        parser.open(&file_path).map_err(from_py_error)?;

        assert_eq!(parser.reader_mode(), "query");
        parser.seek("chr1:10-20").map_err(from_py_error)?;
        let records = parser.records(None).map_err(from_py_error)?;
        assert_eq!(records.len(), 2);
        assert_first_sample_record(&records[0]);
        assert_second_sample_record(&records[1]);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_gff_parser_open_file_query_mode_invalid_tabix_path() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-gff-parser-open-query-invalid-tabix")?;
        let path = write_sample_file(&dir, "sample.gff3.gz", SAMPLE)?;
        let file_path = path.to_string_lossy().into_owned();

        Python::initialize();
        let mut parser = GffParser::new(
            None,
            Some(dir.join("missing-tabix").to_string_lossy().into_owned()),
            Some("query".to_string()),
        )
        .map_err(from_py_error)?;

        let err = parser.open(&file_path).unwrap_err();
        assert!(
            err.to_string().contains("no tabix index found"),
            "unexpected error: {err}"
        );

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // GffParser:open - explicit sequential mode
    // -----------------------------------------------------------------------

    #[test]
    fn test_gff_parser_open_file_sequential_mode() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-gff-parser-open-seq")?;
        let path = write_sample_file(&dir, "sample.gff3", SAMPLE)?;
        let file_path = path.to_string_lossy().into_owned();

        let mut parser = GffParser::new(
            None,
            Some(dir.join("missing-tabix").to_string_lossy().into_owned()),
            Some("sequential".to_string()),
        )
        .map_err(from_py_error)?;

        parser.open(&file_path).map_err(from_py_error)?;
        assert_eq!(parser.reader_mode(), "sequential");

        let records = parser.records(None).map_err(from_py_error)?;
        assert_eq!(records.len(), 2);
        assert_first_sample_record(&records[0]);
        assert_second_sample_record(&records[1]);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_gff_parser_open_file_sequential_mode_with_index() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-gff-parser-open-seq-index")?;
        let (gz_path, _tbi_path) = write_gff3_with_tbi(
            &dir,
            str::from_utf8(SAMPLE)
                .ok()
                .map(|s| s.lines().map(str::to_owned).collect()),
        )?;
        let file_path = gz_path.to_string_lossy().into_owned();

        Python::initialize();
        let mut parser =
            GffParser::new(None, None, Some("sequential".to_string())).map_err(from_py_error)?;

        parser.open(&file_path).map_err(from_py_error)?;
        assert_eq!(parser.reader_mode(), "sequential");

        let records = parser.records(None).map_err(from_py_error)?;
        assert_eq!(records.len(), 2);
        assert_first_sample_record(&records[0]);
        assert_second_sample_record(&records[1]);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_sequential_mode_rejects_seek() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-gff-seq")?;
        let path = write_sample_file(&dir, "sample.gff3", SAMPLE)?;

        let mut parser = GffParser::new(
            Some(path.to_string_lossy().into_owned()),
            None,
            Some("sequential".to_string()),
        )
        .map_err(from_py_error)?;

        assert_eq!(parser.reader_mode(), "sequential");
        assert!(parser
            .seek("chr1:10-20")
            .unwrap_err()
            .to_string()
            .contains("sequential"));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // GffParser:open - invalid mode
    // -----------------------------------------------------------------------

    #[test]
    fn test_gff_parser_rejects_invalid_mode() {
        Python::initialize();
        let err = match GffParser::new(None, None, Some("random".to_string())) {
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
    // GffParser:close
    // -----------------------------------------------------------------------

    #[test]
    fn test_gff_parser_close() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-gff-parser-close")?;
        let path = write_sample_file(&dir, "sample.gff3", SAMPLE)?;

        Python::initialize();
        let mut parser = GffParser::new(None, None, None).map_err(from_py_error)?;
        parser
            .open(&path.to_string_lossy())
            .map_err(from_py_error)?;

        assert!(parser.reader.is_some());
        parser.close();
        assert!(parser.reader.is_none());
        assert!(parser.next().unwrap_err().to_string().contains("not open"));
        assert!(parser
            .seek("chr1:10-20")
            .unwrap_err()
            .to_string()
            .contains("not open"));
        assert!(parser
            .records(None)
            .unwrap_err()
            .to_string()
            .contains("not open"));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // GffParser:parse_text and parse_file
    // -----------------------------------------------------------------------

    #[test]
    fn test_gff_parser_parse_text() -> io::Result<()> {
        let records = GffParser::parse_text(str::from_utf8(SAMPLE).unwrap(), Some(1))
            .map_err(from_py_error)?;
        assert_eq!(records.len(), 1);
        assert_first_sample_record(&records[0]);

        Ok(())
    }

    #[test]
    fn test_gff_parser_parse_file() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-gff-parse-file")?;
        let path = write_sample_file(&dir, "sample.gff3", SAMPLE)?;

        let records =
            GffParser::parse_file(&path.to_string_lossy(), None).map_err(from_py_error)?;
        assert_eq!(records.len(), 2);
        assert_first_sample_record(&records[0]);
        assert_second_sample_record(&records[1]);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_gff_parser_parse_file_bgzf() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-gff-parse-file-bgzf")?;
        let path = write_bgzf_file(&dir, sample_lines(), Some("sample.gff3.gz"))?;

        let records =
            GffParser::parse_file(&path.to_string_lossy(), None).map_err(from_py_error)?;
        assert_eq!(records.len(), 2);
        assert_first_sample_record(&records[0]);
        assert_second_sample_record(&records[1]);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_gff_parser_parse_file_missing_path_errors() {
        Python::initialize();
        let err = GffParser::parse_file("/definitely/missing/sample.gff3", None).unwrap_err();
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
        let dir = create_temp_dir("ensembl-vep-gff-sequential-next")?;
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
        let dir = create_temp_dir("ensembl-vep-gff-sequential-bgzf-next")?;
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
        let dir = create_temp_dir("ensembl-vep-gff-sequential-seek")?;
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
        let dir = create_temp_dir("ensembl-vep-gff-tabix-seek")?;
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
        let dir = create_temp_dir("ensembl-vep-gff-tabix-empty-region")?;
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
        let dir = create_temp_dir("ensembl-vep-gff-tabix-invalid-region")?;
        let (gz_path, _tbi_path) = write_gff3_with_tbi(&dir, None)?;
        let mut backend = TabixBackend::new(&gz_path, "tabix-nonexistent-binary-for-test")?;

        let err = backend.seek(":::invalid:::").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_tabix_backend_next_without_query_returns_none() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-gff-tabix-no-query")?;
        let (gz_path, _tbi_path) = write_gff3_with_tbi(&dir, None)?;
        let mut backend = TabixBackend::new(&gz_path, "tabix-nonexistent-binary-for-test")?;

        assert_eq!(backend.next()?, None);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_tabix_backend_next_errors_when_query_active_without_reader() -> io::Result<()> {
        let mut backend = TabixBackend {
            file_path: PathBuf::from("sample.gff3.gz"),
            tabix_path: "tabix-nonexistent-binary-for-test".to_string(),
            reader: None,
            query_active: true,
        };

        let err = backend.next().unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::Other);
        assert!(
            err.to_string().contains("not initialized"),
            "unexpected error: {err}"
        );

        Ok(())
    }

    // -----------------------------------------------------------------------
    // GffReader::open / Backend::open
    // -----------------------------------------------------------------------

    #[test]
    fn test_gff_reader_query_mode_with_index_uses_tabix_backend() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-gff-reader-tabix")?;
        let (gz_path, _tbi_path) = write_gff3_with_tbi(&dir, None)?;

        let reader = GffReader::open(
            &gz_path,
            "tabix-nonexistent-binary-for-test",
            ReaderMode::Query,
        )?;

        match reader.backend {
            Backend::Tabix(backend) => {
                assert!(backend.reader.is_none());
                assert!(!backend.query_active);
            }
            Backend::Sequential(_) => panic!("expected tabix backend"),
        }

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_gff_reader_query_mode_without_index_errors() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-gff-reader-no-index")?;
        let path = write_sample_file(&dir, "sample.gff3", SAMPLE)?;

        let err = match GffReader::open(
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
    fn test_gff_reader_sequential_mode_uses_sequential_backend_with_index() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-gff-reader-sequential")?;
        let (gz_path, _tbi_path) = write_gff3_with_tbi(&dir, None)?;

        let reader = GffReader::open(
            &gz_path,
            "tabix-nonexistent-binary-for-test",
            ReaderMode::Sequential,
        )?;
        assert!(matches!(reader.backend, Backend::Sequential(_)));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // GFFRecord
    // -----------------------------------------------------------------------

    #[test]
    fn test_gff_record_attributes_getter_and_repr() -> io::Result<()> {
        let records = parse_gff_bytes(SAMPLE, Some(1))?;
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
            "GFFRecord(seqid=\"chr1\", feature_type=\"gene\", start=10, end=20, strand_symbol=\"+\", phase_value=None)"
        );

        Ok(())
    }
}
