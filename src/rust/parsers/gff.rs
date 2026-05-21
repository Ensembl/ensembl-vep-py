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

use crate::parsers::tabix;
use noodles_gff::{
    self as gff,
    feature::record::{Phase, Strand},
    record::attributes::field::Value as AttributeValue,
};
use pyo3::{
    exceptions::{PyIOError, PyRuntimeError, PyValueError},
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

#[pyclass(name = "GFFParser")]
pub struct GffParser {
    file_path: Option<PathBuf>,
    reader: Option<GffReaderKind>,
    tabix_path: String,
    reader_mode: GffReaderMode,
}

enum GffReaderKind {
    Sequential(SequentialGffReader),
    Tabix(Box<TabixGffReader>),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum GffReaderMode {
    Auto,
    Sequential,
    Tabix,
}

struct SequentialGffReader {
    reader: gff::io::Reader<BufReader<File>>,
}

struct TabixGffReader {
    reader: tabix::TabixReader,
}

impl GffReaderKind {
    fn next(&mut self) -> io::Result<Option<GffRecord>> {
        match self {
            Self::Sequential(reader) => reader.next(),
            Self::Tabix(reader) => reader.next(),
        }
    }

    fn error(&self, err: io::Error) -> PyErr {
        match self {
            Self::Sequential(_) => parse_error(err),
            Self::Tabix(_) => tabix_error(err),
        }
    }
}

impl SequentialGffReader {
    fn open<P>(file_path: P) -> io::Result<Self>
    where
        P: AsRef<Path>,
    {
        let file = File::open(file_path)?;

        Ok(Self {
            reader: gff::io::Reader::new(BufReader::new(file)),
        })
    }

    fn next(&mut self) -> io::Result<Option<GffRecord>> {
        read_next_record(&mut self.reader)
    }
}

impl TabixGffReader {
    fn open<P>(file_path: P, tabix_path: &str) -> io::Result<Self>
    where
        P: AsRef<Path>,
    {
        Ok(Self {
            reader: tabix::TabixReader::open(file_path, tabix_path)?,
        })
    }

    fn seek(&mut self, region: &str) -> io::Result<bool> {
        self.reader.seek(region)
    }

    fn next(&mut self) -> io::Result<Option<GffRecord>> {
        loop {
            let Some(line) = self.reader.next()? else {
                return Ok(None);
            };

            if let Some(record) = parse_gff_line(&line)? {
                return Ok(Some(record));
            }
        }
    }
}

impl GffReaderMode {
    fn parse(value: Option<&str>) -> PyResult<Self> {
        match value.unwrap_or("auto") {
            "auto" => Ok(Self::Auto),
            "sequential" => Ok(Self::Sequential),
            "tabix" => Ok(Self::Tabix),
            value => Err(PyValueError::new_err(format!(
                "unsupported GFF reader {value:?}: expected 'auto', 'sequential', or 'tabix'"
            ))),
        }
    }

    fn as_str(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Sequential => "sequential",
            Self::Tabix => "tabix",
        }
    }

    fn initial_reader(self, file_path: &Path, tabix_path: &str) -> io::Result<GffReaderKind> {
        match self {
            Self::Auto | Self::Sequential => {
                SequentialGffReader::open(file_path).map(GffReaderKind::Sequential)
            }
            Self::Tabix => TabixGffReader::open(file_path, tabix_path)
                .map(Box::new)
                .map(GffReaderKind::Tabix),
        }
    }

    fn open_error(self, err: io::Error) -> PyErr {
        match self {
            Self::Auto | Self::Sequential => io_error(err),
            Self::Tabix => tabix_error(err),
        }
    }
}

#[pymethods]
impl GffParser {
    #[new]
    #[pyo3(signature = (file_path=None, tabix_path=None, reader=None))]
    fn new(
        file_path: Option<String>,
        tabix_path: Option<String>,
        reader: Option<String>,
    ) -> PyResult<Self> {
        let mut parser = Self {
            file_path: None,
            reader: None,
            tabix_path: tabix_path.unwrap_or_else(|| "tabix".to_string()),
            reader_mode: GffReaderMode::parse(reader.as_deref())?,
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
    fn reader(&self) -> &'static str {
        self.reader_mode.as_str()
    }

    fn open(&mut self, file_path: &str) -> PyResult<()> {
        let path = PathBuf::from(file_path);

        self.reader = Some(
            self.reader_mode
                .initial_reader(&path, &self.tabix_path)
                .map_err(|err| self.reader_mode.open_error(err))?,
        );
        self.file_path = Some(path);

        Ok(())
    }

    fn close(&mut self) {
        self.reader = None;
    }

    fn set_reader(&mut self, reader: &str) -> PyResult<()> {
        self.reader_mode = GffReaderMode::parse(Some(reader))?;

        if let Some(path) = self.file_path.clone() {
            self.reader = Some(
                self.reader_mode
                    .initial_reader(&path, &self.tabix_path)
                    .map_err(|err| self.reader_mode.open_error(err))?,
            );
        }

        Ok(())
    }

    fn next(&mut self) -> PyResult<Option<GffRecord>> {
        let reader = self.reader_mut()?;

        match reader.next() {
            Ok(record) => Ok(record),
            Err(err) => Err(reader.error(err)),
        }
    }

    #[pyo3(signature = (limit=None))]
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

    fn seek(&mut self, region: &str) -> PyResult<bool> {
        if self.reader_mode == GffReaderMode::Sequential {
            return Err(PyRuntimeError::new_err(
                "GFF reader is set to 'sequential'; set reader to 'auto' or 'tabix' before seek",
            ));
        }

        let has_records = self.tabix_reader()?.seek(region).map_err(tabix_error)?;
        Ok(has_records)
    }

    #[staticmethod]
    #[pyo3(signature = (data, limit=None))]
    fn parse_text(data: &str, limit: Option<usize>) -> PyResult<Vec<GffRecord>> {
        parse_gff_bytes(data.as_bytes(), limit).map_err(parse_error)
    }

    #[staticmethod]
    #[pyo3(signature = (file_path, limit=None))]
    fn parse_file(file_path: &str, limit: Option<usize>) -> PyResult<Vec<GffRecord>> {
        let file = File::open(file_path).map_err(io_error)?;
        let mut reader = gff::io::Reader::new(BufReader::new(file));

        read_records(&mut reader, limit).map_err(parse_error)
    }
}

impl GffParser {
    fn reader_mut(&mut self) -> PyResult<&mut GffReaderKind> {
        self.reader
            .as_mut()
            .ok_or_else(|| PyRuntimeError::new_err("GFF parser is not open"))
    }

    fn tabix_reader(&mut self) -> PyResult<&mut TabixGffReader> {
        if !matches!(self.reader, Some(GffReaderKind::Tabix(_))) {
            let path = self
                .file_path
                .as_ref()
                .ok_or_else(|| PyRuntimeError::new_err("GFF parser is not open"))?;

            self.reader = Some(GffReaderKind::Tabix(Box::new(
                TabixGffReader::open(path, &self.tabix_path).map_err(tabix_error)?,
            )));
        }

        match self.reader.as_mut() {
            Some(GffReaderKind::Tabix(reader)) => Ok(reader.as_mut()),
            _ => Err(PyRuntimeError::new_err("failed to initialize tabix reader")),
        }
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

fn io_error(err: io::Error) -> PyErr {
    PyIOError::new_err(err.to_string())
}

fn parse_error(err: io::Error) -> PyErr {
    match err.kind() {
        io::ErrorKind::InvalidData | io::ErrorKind::InvalidInput => {
            PyValueError::new_err(err.to_string())
        }
        _ => io_error(err),
    }
}

fn tabix_error(err: io::Error) -> PyErr {
    match err.kind() {
        io::ErrorKind::InvalidData | io::ErrorKind::InvalidInput => {
            PyValueError::new_err(err.to_string())
        }
        io::ErrorKind::NotFound => PyIOError::new_err(err.to_string()),
        _ => PyRuntimeError::new_err(err.to_string()),
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
        std::fs::create_dir_all(&dir).map_err(io_error)?;

        let path = dir.join("sample.gff3");
        std::fs::write(&path, SAMPLE).map_err(io_error)?;

        let mut parser = GffParser::new(
            Some(path.to_string_lossy().into_owned()),
            Some(dir.join("missing-tabix").to_string_lossy().into_owned()),
            None,
        )?;

        let first = parser.next()?.expect("first record");
        let second = parser.next()?.expect("second record");

        assert_eq!(first.feature_type, "gene");
        assert_eq!(second.feature_type, "CDS");
        assert!(parser.next()?.is_none());

        std::fs::remove_dir_all(dir).map_err(io_error)?;

        Ok(())
    }
}
