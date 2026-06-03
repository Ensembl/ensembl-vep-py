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

use std::{
    fs::File,
    io::{self, BufRead, BufReader},
    path::{Path, PathBuf},
};

use crate::{
    parsers::tabix::{self, associated_index_path, py_err},
    parsers::utils::ReaderMode,
    utils::py_error,
};
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

type AttributePairs = Vec<(String, Vec<String>)>;

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
    pub strand: String,
    #[pyo3(get)]
    pub phase: Option<u8>,
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
            "GFFRecord(seqid={:?}, feature_type={:?}, start={}, end={})",
            self.seqid, self.feature_type, self.start, self.end
        )
    }
}

pub struct GffReader {
    backend: Backend,
    mode: ReaderMode,
}

enum Backend {
    Sequential(SequentialBackend),
    Tabix(Box<TabixBackend>),
}

struct SequentialBackend {
    reader: gff::io::Reader<BufReader<File>>,
}

struct TabixBackend {
    file_path: PathBuf,
    tabix_path: String,
    reader: Option<tabix::TabixReader>,
    query_active: bool,
}

impl GffReader {
    fn open(file_path: &Path, tabix_path: &str, mode: ReaderMode) -> io::Result<Self> {
        Ok(Self {
            backend: Backend::open(file_path, tabix_path, mode)?,
            mode,
        })
    }

    fn next(&mut self) -> io::Result<Option<GffRecord>> {
        self.backend.next()
    }

    fn seek(&mut self, region: &str) -> io::Result<bool> {
        self.backend.seek(region)
    }
}

impl Backend {
    fn open(file_path: &Path, tabix_path: &str, mode: ReaderMode) -> io::Result<Self> {
        match mode {
            ReaderMode::Sequential => {
                SequentialBackend::new(file_path).map(Self::Sequential)
            }
            ReaderMode::Query => {
                if has_tabix_index(file_path) {
                    TabixBackend::new(file_path, tabix_path)
                        .map(Box::new)
                        .map(Self::Tabix)
                } else {
                    SequentialBackend::new(file_path).map(Self::Sequential)
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

impl SequentialBackend {
    fn new<P>(file_path: P) -> io::Result<Self>
    where
        P: AsRef<Path>,
    {
        let file = File::open(file_path)?;

        Ok(Self {
            reader: gff::io::Reader::new(BufReader::new(file)),
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

fn has_tabix_index(file_path: &Path) -> bool {
    associated_index_path(file_path, "tbi").exists()
        || associated_index_path(file_path, "csi").exists()
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
    fn metadata<'py>(&self, py: Python<'py>) -> Bound<'py, PyDict> {
        PyDict::new(py)
    }

    #[getter]
    fn reader_mode(&self) -> &'static str {
        self.reader_mode.as_str()
    }

    fn open(&mut self, file_path: &str) -> PyResult<()> {
        let path = PathBuf::from(file_path);

        self.reader = Some(GffReader::open(&path, &self.tabix_path, self.reader_mode).map_err(py_err)?);
        self.file_path = Some(path);

        Ok(())
    }

    fn close(&mut self) {
        self.reader = None;
    }

    fn seek(&mut self, region: &str) -> PyResult<bool> {
        self.reader()?.seek(region).map_err(py_err)
    }

    #[allow(clippy::should_implement_trait)]
    fn next(&mut self) -> PyResult<Option<GffRecord>> {
        self.reader()?.next().map_err(py_err)
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
        parse_gff_bytes(data.as_bytes(), limit).map_err(py_err)
    }

    #[staticmethod]
    #[pyo3(signature = (file_path, limit=None))]
    fn parse_file(file_path: &str, limit: Option<usize>) -> PyResult<Vec<GffRecord>> {
        let file = File::open(file_path).map_err(py_err)?;
        let mut reader = gff::io::Reader::new(BufReader::new(file));

        read_records(&mut reader, limit).map_err(py_err)
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
        strand: strand_to_string(record.strand()?),
        phase: record.phase().transpose()?.map(phase_to_u8),
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

    const SAMPLE: &[u8] = b"\
##gff-version 3
chr1\tensembl\tgene\t10\t20\t.\t+\t.\tID=gene:1;Name=BRCA2;Alias=a,b
chr1\tensembl\tCDS\t12\t18\t4.5\t-\t2\tID=cds:1;Parent=gene:1
##FASTA
>chr1
ACGT
";

    #[test]
    fn test_parse_gff_bytes_decodes_records() -> io::Result<()> {
        let records = parse_gff_bytes(SAMPLE, None)?;

        assert_eq!(records.len(), 2);

        let record = &records[0];
        assert_eq!(record.seqid, "chr1");
        assert_eq!(record.source, "ensembl");
        assert_eq!(record.feature_type, "gene");
        assert_eq!(record.start, 10);
        assert_eq!(record.end, 20);
        assert_eq!(record.score, None);
        assert_eq!(record.strand, "+");
        assert_eq!(record.phase, None);
        assert_eq!(record.raw_attributes, "ID=gene:1;Name=BRCA2;Alias=a,b");
        assert_eq!(
            record.attributes,
            vec![
                ("ID".to_string(), vec!["gene:1".to_string()]),
                ("Name".to_string(), vec!["BRCA2".to_string()]),
                ("Alias".to_string(), vec!["a".to_string(), "b".to_string()]),
            ]
        );

        let record = &records[1];
        assert_eq!(record.score, Some(4.5));
        assert_eq!(record.strand, "-");
        assert_eq!(record.phase, Some(2));

        Ok(())
    }

    #[test]
    fn test_parse_gff_bytes_honours_limit() -> io::Result<()> {
        let records = parse_gff_bytes(SAMPLE, Some(1))?;
        assert_eq!(records.len(), 1);
        Ok(())
    }

    #[test]
    fn test_parse_gff_bytes_stops_at_fasta() -> io::Result<()> {
        let records = parse_gff_bytes(SAMPLE, Some(10))?;
        assert_eq!(records.len(), 2);
        Ok(())
    }

    #[test]
    fn test_parser_next_does_not_require_tabix() -> PyResult<()> {
        let dir = std::env::temp_dir().join(format!("ensembl-vep-gff-test-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).map_err(py_err)?;

        let path = dir.join("sample.gff3");
        std::fs::write(&path, SAMPLE).map_err(py_err)?;

        let mut parser = GffParser::new(
            Some(path.to_string_lossy().into_owned()),
            Some(dir.join("missing-tabix").to_string_lossy().into_owned()),
            None,
        )?;

        assert_eq!(parser.reader_mode(), "query");

        let first = parser.next()?.expect("first record");
        let second = parser.next()?.expect("second record");

        assert_eq!(first.feature_type, "gene");
        assert_eq!(second.feature_type, "CDS");
        assert!(parser.next()?.is_none());

        std::fs::remove_dir_all(dir).map_err(py_err)?;

        Ok(())
    }

    #[test]
    fn test_query_mode_seek_without_index_errors() -> PyResult<()> {
        let dir = std::env::temp_dir().join(format!("ensembl-vep-gff-no-index-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).map_err(py_err)?;

        let path = dir.join("sample.gff3.gz");
        std::fs::write(&path, SAMPLE).map_err(py_err)?;

        let mut parser = GffParser::new(Some(path.to_string_lossy().into_owned()), None, None)?;
        let err = parser.seek("chr1:10-20").unwrap_err();
        assert!(err.to_string().contains("sequential"), "unexpected error: {err}");

        std::fs::remove_dir_all(dir).map_err(py_err)?;
        Ok(())
    }

    #[test]
    fn test_sequential_mode_rejects_seek() -> PyResult<()> {
        let dir = std::env::temp_dir().join(format!("ensembl-vep-gff-seq-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).map_err(py_err)?;

        let path = dir.join("sample.gff3");
        std::fs::write(&path, SAMPLE).map_err(py_err)?;

        let mut parser = GffParser::new(
            Some(path.to_string_lossy().into_owned()),
            None,
            Some("sequential".to_string()),
        )?;

        assert_eq!(parser.reader_mode(), "sequential");
        assert!(parser.seek("chr1:10-20").unwrap_err().to_string().contains("sequential"));

        std::fs::remove_dir_all(dir).map_err(py_err)?;
        Ok(())
    }
}
