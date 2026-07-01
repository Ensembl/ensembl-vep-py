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
    parsers::tabix::{self, associated_index_path, has_tabix_index},
    parsers::utils::{
        ReaderMode, noodles_phase_to_u8, noodles_strand_to_string, open_sequential_reader,
    },
    utils::{bytes_to_string, to_py_error},
};
use noodles_gtf::{self as gtf, record::attributes::field::Value as AttributeValue};
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
// GTFRecord - represents a GTF record
// ---------------------------------------------------------------------------

#[pyclass(name = "GTFRecord")]
#[derive(Clone, Debug, PartialEq)]
pub struct GtfRecord {
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
impl GtfRecord {
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
            "GTFRecord(seqid={:?}, feature_type={:?}, start={}, end={}, strand_symbol={:?}, phase_value={:?})",
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
// GtfMetadata - holds GTF metadata (directives and header lines)
// ---------------------------------------------------------------------------

#[derive(Clone, Debug, Default, PartialEq)]
struct GtfMetadata {
    header: Vec<String>,
}

// ---------------------------------------------------------------------------
// GtfReader - public entry point
// ---------------------------------------------------------------------------

// A reusable GTF reader that supports both sequential and indexed (tabix) reading modes.

// The reader by default uses indexed reading if a tabix index is available,
// but can be forced into sequential mode via the `mode` parameter on `GtfParser`.
pub struct GtfReader {
    file_path: PathBuf,
    backend: Backend,
    mode: ReaderMode,
    metadata: Option<GtfMetadata>,
}

impl GtfReader {
    fn open(file_path: &Path, tabix_path: &str, mode: ReaderMode) -> io::Result<Self> {
        let (backend, mode) = Backend::open(file_path, tabix_path, mode)?;
        Ok(Self {
            file_path: file_path.to_path_buf(),
            backend,
            mode,
            metadata: None,
        })
    }

    fn next(&mut self) -> io::Result<Option<GtfRecord>> {
        self.backend.next()
    }

    fn seek(&mut self, region: &str) -> io::Result<bool> {
        self.backend.seek(region)
    }

    fn metadata(&mut self) -> io::Result<&GtfMetadata> {
        if self.metadata.is_none() {
            self.metadata = Some(read_gtf_metadata(&self.file_path)?);
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

    fn next(&mut self) -> io::Result<Option<GtfRecord>> {
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
    reader: gtf::io::Reader<Box<dyn BufRead + Send + Sync>>,
}

impl SequentialBackend {
    fn new<P>(file_path: P) -> io::Result<Self>
    where
        P: AsRef<Path>,
    {
        Ok(Self {
            reader: open_sequential_gtf_reader(file_path.as_ref())?,
        })
    }

    fn seek(&mut self, _region: &str) -> io::Result<bool> {
        Err(io::Error::other(
            "GTF reader is set to 'sequential'; set reader to 'query' before seek",
        ))
    }

    fn next(&mut self) -> io::Result<Option<GtfRecord>> {
        read_next_record(&mut self.reader)
    }
}

// ---------------------------------------------------------------------------
// TabixBackend - uses tabix for indexed reading
// ---------------------------------------------------------------------------

struct TabixBackend {
    reader: Option<tabix::TabixReader>,
    query_active: bool,
}

impl TabixBackend {
    fn new<P>(file_path: P, tabix_path: &str) -> io::Result<Self>
    where
        P: AsRef<Path>,
    {
        Ok(Self {
            reader: Some(tabix::TabixReader::open(file_path.as_ref(), tabix_path)?),
            query_active: false,
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

    fn next(&mut self) -> io::Result<Option<GtfRecord>> {
        if self.query_active {
            let Some(reader) = self.reader.as_mut() else {
                return Err(io::Error::other("tabix reader is not initialized"));
            };

            loop {
                let Some(line) = reader.next()? else {
                    self.query_active = false;
                    return Ok(None);
                };

                if let Some(record) = parse_gtf_line(&line)? {
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

#[pyclass(name = "GTFParser")]
pub struct GtfParser {
    file_path: Option<PathBuf>,
    reader: Option<GtfReader>,
    tabix_path: String,
    reader_mode: ReaderMode,
}

#[pymethods]
impl GtfParser {
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
        "gtf"
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

        Ok(dict)
    }

    #[getter]
    fn reader_mode(&self) -> &'static str {
        self.reader_mode.as_str()
    }

    fn open(&mut self, file_path: &str) -> PyResult<()> {
        let path = PathBuf::from(file_path);

        self.reader =
            Some(GtfReader::open(&path, &self.tabix_path, self.reader_mode).map_err(to_py_error)?);
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
    fn next(&mut self) -> PyResult<Option<GtfRecord>> {
        self.reader()?.next().map_err(to_py_error)
    }

    #[pyo3(signature = (limit=10000))]
    fn records(&mut self, limit: Option<usize>) -> PyResult<Vec<GtfRecord>> {
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
    fn parse_text(data: &str, limit: Option<usize>) -> PyResult<Vec<GtfRecord>> {
        parse_gtf_bytes(data.as_bytes(), limit).map_err(to_py_error)
    }

    #[staticmethod]
    #[pyo3(signature = (file_path, limit=None))]
    fn parse_file(file_path: &str, limit: Option<usize>) -> PyResult<Vec<GtfRecord>> {
        let mut reader = open_sequential_gtf_reader(Path::new(file_path)).map_err(to_py_error)?;

        read_records(&mut reader, limit).map_err(to_py_error)
    }
}

impl GtfParser {
    fn reader(&mut self) -> PyResult<&mut GtfReader> {
        self.reader
            .as_mut()
            .ok_or_else(|| PyRuntimeError::new_err("GTF parser is not open"))
    }
}

// ---------------------------------------------------------------------------
// Helpers: Noodles GTF decoder
// ---------------------------------------------------------------------------

fn open_sequential_gtf_reader(
    file_path: &Path,
) -> io::Result<gtf::io::Reader<Box<dyn BufRead + Send + Sync>>> {
    let reader = open_sequential_reader(file_path)?;
    Ok(gtf::io::Reader::new(reader))
}

fn read_gtf_metadata(file_path: &Path) -> io::Result<GtfMetadata> {
    let mut reader = open_sequential_gtf_reader(file_path)?;
    read_gtf_metadata_from_reader(&mut reader)
}

fn read_gtf_metadata_from_reader<R>(reader: &mut gtf::io::Reader<R>) -> io::Result<GtfMetadata>
where
    R: BufRead,
{
    let mut metadata = GtfMetadata::default();
    let mut line = gtf::Line::default();

    loop {
        if reader.read_line(&mut line)? == 0 {
            break;
        }

        if line.as_ref().is_empty() {
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

fn read_next_record<R>(reader: &mut gtf::io::Reader<R>) -> io::Result<Option<GtfRecord>>
where
    R: BufRead,
{
    let mut line = gtf::Line::default();

    loop {
        if reader.read_line(&mut line)? == 0 {
            return Ok(None);
        }

        if line.as_ref().is_empty() {
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
    reader: &mut gtf::io::Reader<R>,
    limit: Option<usize>,
) -> io::Result<Vec<GtfRecord>>
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

fn parse_gtf_bytes(data: &[u8], limit: Option<usize>) -> io::Result<Vec<GtfRecord>> {
    let mut reader = gtf::io::Reader::new(data);
    read_records(&mut reader, limit)
}

fn parse_gtf_line(line: &str) -> io::Result<Option<GtfRecord>> {
    let mut data = Vec::with_capacity(line.len() + 1);
    data.extend_from_slice(line.as_bytes());
    data.push(b'\n');

    Ok(parse_gtf_bytes(&data, Some(1))?.pop())
}

fn record_to_py(record: gtf::Record<'_>, raw_attributes: String) -> io::Result<GtfRecord> {
    Ok(GtfRecord {
        seqid: bytes_to_string(record.reference_sequence_name()),
        source: bytes_to_string(record.source()),
        feature_type: bytes_to_string(record.ty()),
        start: record.start()?.get(),
        end: record.end()?.get(),
        score: record.score().transpose()?.map(f64::from),
        strand_symbol: noodles_strand_to_string(record.strand()?),
        phase_value: record.phase().transpose()?.map(noodles_phase_to_u8),
        attributes: parse_attributes(record.attributes()?)?,
        raw_attributes,
    })
}

fn parse_attributes(attributes: gtf::record::Attributes<'_>) -> io::Result<AttributePairs> {
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

fn raw_attributes(line: &gtf::Line) -> String {
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
        parsers::test_utils::{write_bgzf_file, write_gtf_with_tbi},
        test_utils::{create_temp_dir, write_sample_file},
        utils::from_py_error,
    };
    use std::{fs, str};

    const SAMPLE: &[u8] = b"\
#!genome-version GRCh38
chr1\tensembl\tgene\t10\t20\t.\t+\t.\tgene_id \"gene:1\"; gene_name \"BRCA2\";
chr1\tensembl\ttranscript\t12\t18\t4.5\t-\t2\tgene_id \"gene:1\"; transcript_id \"transcript:1\"; transcript_name \"BRCA2-001\";
";

    const SKIPPED_LINE_SAMPLE: &[u8] = b"\
#!genome-version GRCh38

# a regular comment
chr1\tensembl\tgene\t10\t20\t.\t.\t.\tgene_id \"gene:1\"; gene_name \"BRCA2\";
";

    // -----------------------------------------------------------------------
    // Helpers
    // -----------------------------------------------------------------------

    fn to_gtf_line(s: &str) -> gtf::Line {
        let mut reader = gtf::io::Reader::new(s.as_bytes());
        let mut line = gtf::Line::default();
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

    fn assert_first_sample_record(record: &GtfRecord) {
        assert_eq!(record.seqid, "chr1");
        assert_eq!(record.source, "ensembl");
        assert_eq!(record.feature_type, "gene");
        assert_eq!(record.start, 10);
        assert_eq!(record.end, 20);
        assert_eq!(record.score, None);
        assert_eq!(record.strand_symbol, "+");
        assert_eq!(record.phase_value, None);
        assert_eq!(
            record.raw_attributes,
            "gene_id \"gene:1\"; gene_name \"BRCA2\";"
        );
        assert_eq!(
            record.attributes,
            vec![
                ("gene_id".to_string(), vec!["gene:1".to_string()]),
                ("gene_name".to_string(), vec!["BRCA2".to_string()]),
            ]
        );
    }

    fn assert_second_sample_record(record: &GtfRecord) {
        assert_eq!(record.seqid, "chr1");
        assert_eq!(record.feature_type, "transcript");
        assert_eq!(record.start, 12);
        assert_eq!(record.end, 18);
        assert_eq!(record.score, Some(4.5));
        assert_eq!(record.strand_symbol, "-");
        assert_eq!(record.phase_value, Some(2));
        assert_eq!(
            record.attributes,
            vec![
                ("gene_id".to_string(), vec!["gene:1".to_string()]),
                (
                    "transcript_id".to_string(),
                    vec!["transcript:1".to_string()]
                ),
                ("transcript_name".to_string(), vec!["BRCA2-001".to_string()]),
            ]
        );
    }

    // -----------------------------------------------------------------------
    // raw_attributes
    // -----------------------------------------------------------------------

    #[test]
    fn test_raw_attributes_with_attributes() -> io::Result<()> {
        let line = to_gtf_line(
            "chr1\tensembl\tgene\t10\t20\t.\t+\t.\tgene_id \"gene:1\"; gene_name \"BRCA2\";",
        );
        assert_eq!(
            raw_attributes(&line),
            "gene_id \"gene:1\"; gene_name \"BRCA2\";"
        );
        Ok(())
    }

    #[test]
    fn test_raw_attributes_missing_attribute_field_returns_empty() -> io::Result<()> {
        let line = to_gtf_line("chr1\tensembl\tgene\t10\t20\t.\t+\t.");
        assert_eq!(raw_attributes(&line), "");
        Ok(())
    }

    // -----------------------------------------------------------------------
    // parse_attributes
    // -----------------------------------------------------------------------

    #[test]
    fn test_parse_attributes_single_valued() -> io::Result<()> {
        let line = to_gtf_line(
            "chr1\tensembl\tgene\t10\t20\t.\t+\t.\tgene_id \"gene:1\"; gene_name \"BRCA2\";",
        );
        let record = line.as_record().unwrap().unwrap();
        let attributes = record.attributes();
        let pairs = parse_attributes(attributes?)?;

        assert_eq!(pairs.len(), 2);
        assert_eq!(pairs[0].0, "gene_id");
        assert_eq!(pairs[0].1, vec!["gene:1".to_string()]);
        assert_eq!(pairs[1].0, "gene_name");
        assert_eq!(pairs[1].1, vec!["BRCA2".to_string()]);
        Ok(())
    }

    #[test]
    fn test_parse_attributes_repeated_item() -> io::Result<()> {
        let line = to_gtf_line(
            "chr1\tensembl\tgene\t10\t20\t.\t+\t.\ttag \"CCDS\"; tag \"basic\"; tag \"Ensembl_canonical\";",
        );
        let record = line.as_record().unwrap().unwrap();
        let attributes = record.attributes();
        let pairs = parse_attributes(attributes?)?;

        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0].0, "tag");
        assert_eq!(
            pairs[0].1,
            vec![
                "CCDS".to_string(),
                "basic".to_string(),
                "Ensembl_canonical".to_string()
            ]
        );
        Ok(())
    }

    #[test]
    fn test_parse_attributes_multi_valued() -> io::Result<()> {
        let line = to_gtf_line(
            "chr1\tensembl\tgene\t10\t20\t.\t+\t.\tgene_id \"gene:1\"; db_xref \"GO:0006412,GO:0003735\";",
        );
        let record = line.as_record().unwrap().unwrap();
        let attributes = record.attributes();
        let pairs = parse_attributes(attributes?)?;

        assert_eq!(pairs.len(), 2);
        assert_eq!(pairs[0].0, "gene_id");
        assert_eq!(pairs[0].1, vec!["gene:1".to_string()]);
        assert_eq!(pairs[1].0, "db_xref");
        assert_eq!(pairs[1].1, vec!["GO:0006412,GO:0003735".to_string()]);
        Ok(())
    }

    #[test]
    fn test_parse_attributes_invalid() -> io::Result<()> {
        let line = to_gtf_line("chr1\tensembl\tgene\t10\t20\t.\t+\t.\tgene_id \"gene:1\"; Invalid");
        let record = line.as_record().unwrap().unwrap();
        let attributes = record.attributes();
        let err = attributes.unwrap_err();

        assert!(
            err.to_string().contains("invalid data"),
            "unexpected error: {}",
            err
        );
        Ok(())
    }

    #[test]
    fn test_parse_attributes_empty() -> io::Result<()> {
        let line = to_gtf_line("chr1\tensembl\tgene\t10\t20\t.\t+\t.\t");
        let record = line.as_record().unwrap().unwrap();
        let attributes = record.attributes();
        let pairs = parse_attributes(attributes?)?;

        assert!(pairs.is_empty());
        Ok(())
    }

    // -----------------------------------------------------------------------
    // record_to_py
    // -----------------------------------------------------------------------

    #[test]
    fn test_record_to_py() -> io::Result<()> {
        let line = to_gtf_line(
            "chr1\tensembl\tgene\t10\t20\t.\t+\t.\tgene_id \"gene:1\"; gene_name \"BRCA2\";",
        );
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
        assert_eq!(
            py_record.raw_attributes,
            "gene_id \"gene:1\"; gene_name \"BRCA2\";"
        );
        assert_eq!(
            py_record.attributes,
            vec![
                ("gene_id".to_string(), vec!["gene:1".to_string()]),
                ("gene_name".to_string(), vec!["BRCA2".to_string()]),
            ]
        );
        Ok(())
    }

    // -----------------------------------------------------------------------
    // parse_gtf_line
    // -----------------------------------------------------------------------

    #[test]
    fn test_parse_gtf_line_decodes_record() -> io::Result<()> {
        let record = parse_gtf_line(
            "chr1\tensembl\tgene\t10\t20\t.\t+\t.\tgene_id \"gene:1\"; gene_name \"BRCA2\";",
        )?
        .expect("expected record");

        assert_eq!(record.seqid, "chr1");
        assert_eq!(record.source, "ensembl");
        assert_eq!(record.feature_type, "gene");
        assert_eq!(record.start, 10);
        assert_eq!(record.end, 20);
        assert_eq!(record.score, None);
        assert_eq!(record.strand_symbol, "+");
        assert_eq!(record.phase_value, None);
        assert_eq!(
            record.attributes,
            vec![
                ("gene_id".to_string(), vec!["gene:1".to_string()]),
                ("gene_name".to_string(), vec!["BRCA2".to_string()]),
            ]
        );

        Ok(())
    }

    #[test]
    fn test_parse_gtf_line_invalid_line_bad_fields() -> io::Result<()> {
        let err =
            parse_gtf_line("chr1\tensembl\tgene\tbad_start\t20\t.\t+\t.\tID=gene:1").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);

        let err =
            parse_gtf_line("chr1\tensembl\tgene\t10\tbad_end\t.\t+\t.\tID=gene:1").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);

        let err =
            parse_gtf_line("chr1\tensembl\tgene\t10\t20\tbad_score\t+\t.\tID=gene:1").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);

        let err =
            parse_gtf_line("chr1\tensembl\tgene\t10\t20\t.\tbad_strand\t.\tID=gene:1").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);

        let err =
            parse_gtf_line("chr1\tensembl\tgene\t10\t20\t.\t+\tbad_phase\tID=gene:1").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);

        let err =
            parse_gtf_line("chr1\tensembl\tgene\t10\t20\t.\t+\t.\tInvalidAttribute").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        Ok(())
    }

    #[test]
    fn test_parse_gtf_line_invalid_line_missing_fields() -> io::Result<()> {
        let err = parse_gtf_line("chr1\tensembl\tgene\t10\t20\t.\t+\t.").unwrap_err();

        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
        Ok(())
    }

    #[test]
    fn test_parse_gtf_line_ignores_non_record_lines() -> io::Result<()> {
        assert_eq!(parse_gtf_line("")?, None);
        assert_eq!(parse_gtf_line("#!genome-version GRCh38")?, None);
        assert_eq!(parse_gtf_line("# a regular comment")?, None);

        Ok(())
    }

    // -----------------------------------------------------------------------
    // parse_gtf_bytes
    // -----------------------------------------------------------------------

    #[test]
    fn test_parse_gtf_bytes_decodes_records() -> io::Result<()> {
        let records = parse_gtf_bytes(SAMPLE, None)?;

        assert_eq!(records.len(), 2);
        assert_first_sample_record(&records[0]);
        assert_second_sample_record(&records[1]);

        Ok(())
    }

    #[test]
    fn test_parse_gtf_bytes_honours_limit() -> io::Result<()> {
        let records = parse_gtf_bytes(SAMPLE, Some(1))?;

        assert_eq!(records.len(), 1);
        Ok(())
    }

    #[test]
    fn test_parse_gtf_bytes_honours_zero_limit() -> io::Result<()> {
        let records = parse_gtf_bytes(SAMPLE, Some(0))?;

        assert!(records.is_empty());
        Ok(())
    }

    #[test]
    fn test_parse_gtf_bytes_skips_empty_and_comment_lines() -> io::Result<()> {
        let records = parse_gtf_bytes(SKIPPED_LINE_SAMPLE, None)?;

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].strand_symbol, ".");
        assert_eq!(
            records[0].raw_attributes,
            "gene_id \"gene:1\"; gene_name \"BRCA2\";"
        );

        Ok(())
    }

    #[test]
    fn test_parse_gtf_bytes_invalid_record_errors() {
        let data = b"chr1\tensembl\tgene\tbad_start\t20\t.\t+\t.\tID=gene:1\n";
        let err = parse_gtf_bytes(data, None).unwrap_err();

        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    }

    // -----------------------------------------------------------------------
    // read_records
    // -----------------------------------------------------------------------

    #[test]
    fn test_read_records_reads_multiple_records() -> io::Result<()> {
        let mut reader = gtf::io::Reader::new(SAMPLE);
        let records = read_records(&mut reader, None)?;

        assert_eq!(records.len(), 2);
        assert_first_sample_record(&records[0]);
        assert_second_sample_record(&records[1]);
        Ok(())
    }

    #[test]
    fn test_read_records_honours_limit() -> io::Result<()> {
        let mut reader = gtf::io::Reader::new(SAMPLE);
        let records = read_records(&mut reader, Some(1))?;

        assert_eq!(records.len(), 1);
        Ok(())
    }

    #[test]
    fn test_read_records_honours_zero_limit() -> io::Result<()> {
        let mut reader = gtf::io::Reader::new(SAMPLE);
        let records = read_records(&mut reader, Some(0))?;

        assert!(records.is_empty());
        Ok(())
    }

    // -----------------------------------------------------------------------
    // read_next_record
    // -----------------------------------------------------------------------

    #[test]
    fn test_read_next_record_reads_records_one_by_one() -> io::Result<()> {
        let mut reader = gtf::io::Reader::new(SAMPLE);

        let first = read_next_record(&mut reader)?.expect("expected first record");
        assert_first_sample_record(&first);

        let second = read_next_record(&mut reader)?.expect("expected second record");
        assert_second_sample_record(&second);
        Ok(())
    }

    #[test]
    fn test_read_next_record_returns_none_at_eof() -> io::Result<()> {
        let mut reader = gtf::io::Reader::new(SAMPLE);
        read_next_record(&mut reader)?; // first record
        read_next_record(&mut reader)?; // second record

        assert_eq!(read_next_record(&mut reader)?, None);
        Ok(())
    }

    #[test]
    fn test_read_next_record_no_records() -> io::Result<()> {
        let mut reader = gtf::io::Reader::new(b"#!genome-version GRCh38\n# a comment".as_slice());

        assert_eq!(read_next_record(&mut reader)?, None);
        Ok(())
    }

    // -----------------------------------------------------------------------
    // read_gtf_metadata
    // -----------------------------------------------------------------------

    #[test]
    fn test_read_gtf_metadata_from_reader_collects_header() -> io::Result<()> {
        let mut reader = gtf::io::Reader::new(SKIPPED_LINE_SAMPLE);
        let metadata = read_gtf_metadata_from_reader(&mut reader)?;

        assert_eq!(
            metadata.header,
            vec![
                "#!genome-version GRCh38".to_string(),
                "# a regular comment".to_string()
            ]
        );
        Ok(())
    }

    #[test]
    fn test_read_gtf_metadata_with_no_header() -> io::Result<()> {
        let mut reader =
            gtf::io::Reader::new(b"chr1\tensembl\tgene\t10\t20\t.\t+\t.\tID=gene:1\n".as_slice());
        let metadata = read_gtf_metadata_from_reader(&mut reader)?;

        assert!(metadata.header.is_empty());
        Ok(())
    }

    // -----------------------------------------------------------------------
    // open_sequential_gtf_reader
    // -----------------------------------------------------------------------

    #[test]
    fn test_open_sequential_gtf_reader_reads_plain_file() -> io::Result<()> {
        let dir = create_temp_dir("gtf-open-sequential-plain")?;
        let path = write_sample_file(&dir, "sample.gtf", SAMPLE)?;
        let mut reader = open_sequential_gtf_reader(&path)?;
        let records = read_records(&mut reader, None)?;

        assert_eq!(records.len(), 2);
        assert_first_sample_record(&records[0]);
        assert_second_sample_record(&records[1]);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_open_sequential_gtf_reader_reads_bgzf_file() -> io::Result<()> {
        let dir = create_temp_dir("gtf-open-sequential-bgzf")?;
        let path = write_bgzf_file(&dir, sample_lines(), Some("sample.gtf.gz"))?;
        let mut reader = open_sequential_gtf_reader(&path)?;
        let records = read_records(&mut reader, None)?;

        assert_eq!(records.len(), 2);
        assert_first_sample_record(&records[0]);
        assert_second_sample_record(&records[1]);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // GtfParser:metadata
    // -----------------------------------------------------------------------

    #[test]
    fn test_gtf_parser_metadata_sequential() -> io::Result<()> {
        let dir = create_temp_dir("gtf-parser-metadata-seq")?;
        let path = write_sample_file(&dir, "sample.gtf", SKIPPED_LINE_SAMPLE)?;
        let file_path = path.to_string_lossy().into_owned();

        let mut parser =
            GtfParser::new(None, None, Some("sequential".to_string())).map_err(from_py_error)?;
        parser.open(&file_path).map_err(from_py_error)?;

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
                    "#!genome-version GRCh38".to_string(),
                    "# a regular comment".to_string()
                ]
            );
            Ok(())
        })
        .map_err(from_py_error)?;

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_gtf_parser_metadata_does_not_consume_sequential_reader() -> io::Result<()> {
        let dir = create_temp_dir("gtf-parser-metadata-preserve-seq")?;
        let path = write_sample_file(&dir, "sample.gtf", SAMPLE)?;
        let file_path = path.to_string_lossy().into_owned();

        let mut parser =
            GtfParser::new(None, None, Some("sequential".to_string())).map_err(from_py_error)?;
        parser.open(&file_path).map_err(from_py_error)?;

        let first = parser
            .next()
            .map_err(from_py_error)?
            .expect("expected first record");
        assert_first_sample_record(&first);

        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let metadata = parser.metadata(py)?;
            let header: Vec<String> = metadata
                .get_item("header")?
                .expect("metadata should contain header")
                .extract()?;

            assert_eq!(header, vec!["#!genome-version GRCh38".to_string(),]);
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
    fn test_gtf_parser_metadata_query_bgzf() -> io::Result<()> {
        let dir = create_temp_dir("gtf-parser-metadata-query")?;
        let (gz_path, _tbi_path) = write_gtf_with_tbi(
            &dir,
            str::from_utf8(SKIPPED_LINE_SAMPLE)
                .ok()
                .map(|s| s.lines().map(str::to_owned).collect()),
        )?;
        let file_path = gz_path.to_string_lossy().into_owned();

        let mut parser = GtfParser::new(None, None, None).map_err(from_py_error)?;
        parser.open(&file_path).map_err(from_py_error)?;

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
                    "#!genome-version GRCh38".to_string(),
                    "# a regular comment".to_string()
                ]
            );
            Ok(())
        })
        .map_err(from_py_error)?;

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_gtf_parser_open_missing_file_errors() -> io::Result<()> {
        Python::initialize();
        let mut parser = GtfParser::new(None, None, None).map_err(from_py_error)?;
        let err = parser.open("/definitely/missing/sample.gtf").unwrap_err();

        assert!(
            err.to_string().contains("No such file") || err.to_string().contains("os error 2"),
            "unexpected error: {err}"
        );
        Ok(())
    }

    // -----------------------------------------------------------------------
    // GtfParser:open
    // -----------------------------------------------------------------------

    #[test]
    fn test_gtf_parser_without_open_reader() -> io::Result<()> {
        let mut parser = GtfParser::new(None, None, None).map_err(from_py_error)?;

        assert_eq!(parser.format(), "gtf");
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
                .records(None)
                .unwrap_err()
                .to_string()
                .contains("not open")
        );

        parser.close();
        Ok(())
    }

    // -----------------------------------------------------------------------
    // GtfParser:open - auto reader modes
    // -----------------------------------------------------------------------

    #[test]
    fn test_gtf_parser_open_file_auto_mode_sequential() -> io::Result<()> {
        let dir = create_temp_dir("gtf-parser-open-auto")?;
        let path = write_sample_file(&dir, "sample.gtf", SAMPLE)?;
        let file_path = path.to_string_lossy().into_owned();

        let mut parser = GtfParser::new(None, None, None).map_err(from_py_error)?;
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
    fn test_gtf_parser_open_file_auto_mode_query() -> io::Result<()> {
        let dir = create_temp_dir("gtf-parser-open-auto-query")?;
        let (gz_path, _tbi_path) = write_gtf_with_tbi(
            &dir,
            str::from_utf8(SAMPLE)
                .ok()
                .map(|s| s.lines().map(str::to_owned).collect()),
        )?;
        let file_path = gz_path.to_string_lossy().into_owned();

        let mut parser = GtfParser::new(None, None, None).map_err(from_py_error)?;
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
    // GtfParser:open - explicit query mode
    // -----------------------------------------------------------------------

    #[test]
    fn test_gtf_parser_open_file_query_mode_without_index_errors() -> io::Result<()> {
        let dir = create_temp_dir("gtf-parser-open-query-no-index")?;
        let path = write_sample_file(&dir, "sample.gtf", SAMPLE)?;
        let file_path = path.to_string_lossy().into_owned();

        Python::initialize();
        let mut parser =
            GtfParser::new(None, None, Some("query".to_string())).map_err(from_py_error)?;

        let err = parser.open(&file_path).unwrap_err();
        assert!(
            err.to_string().contains("no tabix index found"),
            "unexpected error: {err}"
        );

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_gtf_parser_open_file_query_mode_with_index() -> io::Result<()> {
        let dir = create_temp_dir("gtf-parser-open-query-index")?;
        let (gz_path, _tbi_path) = write_gtf_with_tbi(
            &dir,
            str::from_utf8(SAMPLE)
                .ok()
                .map(|s| s.lines().map(str::to_owned).collect()),
        )?;
        let file_path = gz_path.to_string_lossy().into_owned();

        let mut parser =
            GtfParser::new(None, None, Some("query".to_string())).map_err(from_py_error)?;
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
    fn test_gtf_parser_open_file_query_mode_invalid_tabix_path() -> io::Result<()> {
        let dir = create_temp_dir("gtf-parser-open-query-invalid-tabix")?;
        let path = write_sample_file(&dir, "sample.gtf.gz", SAMPLE)?;
        let file_path = path.to_string_lossy().into_owned();

        Python::initialize();
        let mut parser = GtfParser::new(
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
    // GtfParser:open - explicit sequential mode
    // -----------------------------------------------------------------------

    #[test]
    fn test_gtf_parser_open_file_sequential_mode() -> io::Result<()> {
        let dir = create_temp_dir("gtf-parser-open-seq")?;
        let path = write_sample_file(&dir, "sample.gtf", SAMPLE)?;
        let file_path = path.to_string_lossy().into_owned();

        let mut parser = GtfParser::new(
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
    fn test_gtf_parser_open_file_sequential_mode_with_index() -> io::Result<()> {
        let dir = create_temp_dir("gtf-parser-open-seq-index")?;
        let (gz_path, _tbi_path) = write_gtf_with_tbi(
            &dir,
            str::from_utf8(SAMPLE)
                .ok()
                .map(|s| s.lines().map(str::to_owned).collect()),
        )?;
        let file_path = gz_path.to_string_lossy().into_owned();

        Python::initialize();
        let mut parser =
            GtfParser::new(None, None, Some("sequential".to_string())).map_err(from_py_error)?;

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
        let dir = create_temp_dir("gtf-seq")?;
        let path = write_sample_file(&dir, "sample.gtf", SAMPLE)?;

        let mut parser = GtfParser::new(
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
    // GtfParser:open - invalid mode
    // -----------------------------------------------------------------------

    #[test]
    fn test_gtf_parser_rejects_invalid_mode() {
        Python::initialize();
        let err = match GtfParser::new(None, None, Some("random".to_string())) {
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
    // GtfParser:close
    // -----------------------------------------------------------------------

    #[test]
    fn test_gtf_parser_close() -> io::Result<()> {
        let dir = create_temp_dir("gtf-parser-close")?;
        let path = write_sample_file(&dir, "sample.gtf", SAMPLE)?;

        Python::initialize();
        let mut parser = GtfParser::new(None, None, None).map_err(from_py_error)?;
        parser
            .open(&path.to_string_lossy())
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
                .records(None)
                .unwrap_err()
                .to_string()
                .contains("not open")
        );

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // GtfParser:parse_text and parse_file
    // -----------------------------------------------------------------------

    #[test]
    fn test_gtf_parser_parse_text() -> io::Result<()> {
        let records = GtfParser::parse_text(str::from_utf8(SAMPLE).unwrap(), Some(1))
            .map_err(from_py_error)?;
        assert_eq!(records.len(), 1);
        assert_first_sample_record(&records[0]);

        Ok(())
    }

    #[test]
    fn test_gtf_parser_parse_file() -> io::Result<()> {
        let dir = create_temp_dir("gtf-parse-file")?;
        let path = write_sample_file(&dir, "sample.gtf", SAMPLE)?;

        let records =
            GtfParser::parse_file(&path.to_string_lossy(), None).map_err(from_py_error)?;
        assert_eq!(records.len(), 2);
        assert_first_sample_record(&records[0]);
        assert_second_sample_record(&records[1]);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_gtf_parser_parse_file_bgzf() -> io::Result<()> {
        let dir = create_temp_dir("gtf-parse-file-bgzf")?;
        let path = write_bgzf_file(&dir, sample_lines(), Some("sample.gtf.gz"))?;

        let records =
            GtfParser::parse_file(&path.to_string_lossy(), None).map_err(from_py_error)?;
        assert_eq!(records.len(), 2);
        assert_first_sample_record(&records[0]);
        assert_second_sample_record(&records[1]);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_gtf_parser_parse_file_missing_path_errors() {
        Python::initialize();
        let err = GtfParser::parse_file("/definitely/missing/sample.gtf", None).unwrap_err();
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
        let dir = create_temp_dir("gtf-sequential-next")?;
        let path = write_sample_file(&dir, "sample.gtf", SAMPLE)?;
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
        let dir = create_temp_dir("gtf-sequential-bgzf-next")?;
        let path = write_bgzf_file(&dir, sample_lines(), Some("sample.gtf.gz"))?;
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
        let dir = create_temp_dir("gtf-sequential-seek")?;
        let path = write_sample_file(&dir, "sample.gtf", SAMPLE)?;
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
        let dir = create_temp_dir("gtf-tabix-seek")?;
        let (gz_path, _tbi_path) = write_gtf_with_tbi(&dir, None)?;
        let mut backend = TabixBackend::new(&gz_path, "tabix-nonexistent-binary-for-test")?;

        assert!(backend.seek("chr1:150-250")?);

        let first = backend.next()?.expect("expected first record");
        let second = backend.next()?.expect("expected second record");

        assert_eq!(first.raw_attributes, "gene_id \"gene:2\";");
        assert_eq!(second.raw_attributes, "gene_id \"gene:3\";");
        assert_eq!(backend.next()?, None);
        assert!(!backend.query_active);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_tabix_backend_seek_empty_region() -> io::Result<()> {
        let dir = create_temp_dir("gtf-tabix-empty-region")?;
        let (gz_path, _tbi_path) = write_gtf_with_tbi(&dir, None)?;
        let mut backend = TabixBackend::new(&gz_path, "tabix-nonexistent-binary-for-test")?;

        assert!(!backend.seek("chr1:600-700")?);
        assert!(!backend.query_active);
        assert_eq!(backend.next()?, None);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_tabix_backend_seek_invalid_region_errors() -> io::Result<()> {
        let dir = create_temp_dir("gtf-tabix-invalid-region")?;
        let (gz_path, _tbi_path) = write_gtf_with_tbi(&dir, None)?;
        let mut backend = TabixBackend::new(&gz_path, "tabix-nonexistent-binary-for-test")?;

        let err = backend.seek(":::invalid:::").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_tabix_backend_next_without_query_returns_none() -> io::Result<()> {
        let dir = create_temp_dir("gtf-tabix-no-query")?;
        let (gz_path, _tbi_path) = write_gtf_with_tbi(&dir, None)?;
        let mut backend = TabixBackend::new(&gz_path, "tabix-nonexistent-binary-for-test")?;

        assert_eq!(backend.next()?, None);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // GtfReader::open / Backend::open
    // -----------------------------------------------------------------------

    #[test]
    fn test_gtf_reader_query_mode_with_index_uses_tabix_backend() -> io::Result<()> {
        let dir = create_temp_dir("gtf-reader-tabix")?;
        let (gz_path, _tbi_path) = write_gtf_with_tbi(&dir, None)?;

        let reader = GtfReader::open(
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
    fn test_gtf_reader_query_mode_without_index_errors() -> io::Result<()> {
        let dir = create_temp_dir("gtf-reader-no-index")?;
        let path = write_sample_file(&dir, "sample.gtf", SAMPLE)?;

        let err = match GtfReader::open(
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
    fn test_gtf_reader_sequential_mode_uses_sequential_backend_with_index() -> io::Result<()> {
        let dir = create_temp_dir("gtf-reader-sequential")?;
        let (gz_path, _tbi_path) = write_gtf_with_tbi(&dir, None)?;

        let reader = GtfReader::open(
            &gz_path,
            "tabix-nonexistent-binary-for-test",
            ReaderMode::Sequential,
        )?;
        assert!(matches!(reader.backend, Backend::Sequential(_)));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // GTFRecord
    // -----------------------------------------------------------------------

    #[test]
    fn test_gtf_record_attributes_getter_and_repr() -> io::Result<()> {
        let records = parse_gtf_bytes(SAMPLE, Some(1))?;
        let record = records.first().expect("expected record");

        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let attributes = record.attributes(py)?;

            assert_eq!(attributes.len(), 2);
            assert_eq!(
                attributes
                    .get_item("gene_id")?
                    .expect("expected gene_id")
                    .extract::<Vec<String>>()?,
                vec!["gene:1".to_string()]
            );
            assert_eq!(
                attributes
                    .get_item("gene_name")?
                    .expect("expected gene_name")
                    .extract::<Vec<String>>()?,
                vec!["BRCA2".to_string()]
            );

            Ok(())
        })
        .map_err(from_py_error)?;

        assert_eq!(
            record.__repr__(),
            "GTFRecord(seqid=\"chr1\", feature_type=\"gene\", start=10, end=20, strand_symbol=\"+\", phase_value=None)"
        );

        Ok(())
    }
}
