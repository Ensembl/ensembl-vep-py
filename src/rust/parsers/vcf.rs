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
    parsers::utils::{ReaderMode, open_sequential_reader},
    utils::{bytes_to_string, to_py_error},
};
use noodles_vcf::{
    self as vcf,
    header::{
        FileFormat,
        record::value::{
            Collection, Map,
            map::{AlternativeAllele, Contig, Filter, Format, Info},
        },
    },
    variant::{
        self,
        record_buf::{
            info::field::Value as InfoFieldValue, info::field::value::Array as ValueArray,
            samples::sample::value::Array as SampleArray,
            samples::sample::value::Genotype as SampleGenotype,
            samples::sample::value::Value as SampleValue,
        },
    },
};
use pyo3::{
    exceptions::PyRuntimeError,
    prelude::*,
    types::{PyDict, PyList},
};
use std::{
    collections::HashMap,
    fmt::Display,
    io::{self, BufRead, Read},
    path::{Path, PathBuf},
};

// ---------------------------------------------------------------------------
// VCFRecord - represents a VCF record
// ---------------------------------------------------------------------------

#[pyclass(name = "VCFRecord")]
#[derive(Clone, Debug, PartialEq)]
pub struct VcfRecord {
    #[pyo3(get)]
    pub chrom: String,
    #[pyo3(get)]
    pub pos: Option<usize>,
    #[pyo3(get)]
    pub id: Vec<String>,
    #[pyo3(get)]
    pub reference: String,
    #[pyo3(get)]
    pub alt: Vec<String>,
    #[pyo3(get)]
    pub qual: Option<f32>,
    #[pyo3(get)]
    pub filter: Vec<String>,
    info: variant::record_buf::Info,
    samples: variant::record_buf::Samples,
}

#[pymethods]
impl VcfRecord {
    #[getter]
    fn info<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        let info = &self.info;

        for (key, value) in info.keys().cloned().zip(info.values()) {
            dict.set_item(key, info_value_to_string(value))?;
        }

        Ok(dict)
    }

    #[getter]
    fn samples<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let list = PyList::empty(py);

        for sample in self.samples.values() {
            let dict = PyDict::new(py);

            for (key, value) in sample.keys().as_ref().iter().cloned().zip(sample.values()) {
                dict.set_item(key, sample_value_to_string(value))?;
            }

            list.append(dict)?;
        }

        Ok(list)
    }

    fn __repr__(&self) -> String {
        format!(
            "VCFRecord(chrom={:?}, pos={:?}, id={:?}, reference={:?}, alt={:?}, qual={:?}, filter={:?})",
            self.chrom, self.pos, self.id, self.reference, self.alt, self.qual, self.filter
        )
    }
}

// ---------------------------------------------------------------------------
// VCFMetadata - holds VCF metadata
// ---------------------------------------------------------------------------

#[pyclass(name = "VCFMetadata", mapping)]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct VcfMetadata {
    #[pyo3(get)]
    pub file_format: String,
    #[pyo3(get)]
    pub info: HashMap<String, HashMap<String, String>>,
    #[pyo3(get)]
    pub filter: HashMap<String, HashMap<String, String>>,
    #[pyo3(get)]
    pub format: HashMap<String, HashMap<String, String>>,
    #[pyo3(get)]
    pub alt: HashMap<String, HashMap<String, String>>,
    #[pyo3(get)]
    pub contig: HashMap<String, HashMap<String, String>>,
    #[pyo3(get)]
    pub samples: Vec<String>,
    #[pyo3(get)]
    pub other: HashMap<String, Vec<String>>,
}

#[pymethods]
impl VcfMetadata {
    #[new]
    fn new() -> Self {
        Self::default()
    }

    fn __len__(&self) -> usize {
        self.keys().len()
    }

    fn keys(&self) -> Vec<String> {
        if self.is_empty() {
            return Vec::new();
        }

        vec![
            "file_format".to_string(),
            "info".to_string(),
            "filter".to_string(),
            "format".to_string(),
            "alt".to_string(),
            "contig".to_string(),
            "samples".to_string(),
            "other".to_string(),
        ]
    }

    fn __getitem__(&self, py: Python<'_>, key: String) -> PyResult<Py<PyAny>> {
        match key.as_str() {
            "file_format" => Ok(self
                .file_format
                .clone()
                .into_pyobject(py)?
                .into_any()
                .unbind()),
            "info" => Ok(self.info.clone().into_pyobject(py)?.into_any().unbind()),
            "filter" => Ok(self.filter.clone().into_pyobject(py)?.into_any().unbind()),
            "format" => Ok(self.format.clone().into_pyobject(py)?.into_any().unbind()),
            "alt" => Ok(self.alt.clone().into_pyobject(py)?.into_any().unbind()),
            "contig" => Ok(self.contig.clone().into_pyobject(py)?.into_any().unbind()),
            "samples" => Ok(self.samples.clone().into_pyobject(py)?.into_any().unbind()),
            "other" => Ok(self.other.clone().into_pyobject(py)?.into_any().unbind()),
            _ => Err(PyRuntimeError::new_err(format!(
                "unknown VCF metadata key {key:?}"
            ))),
        }
    }

    fn is_empty(&self) -> bool {
        self.file_format.is_empty()
            && self.info.is_empty()
            && self.filter.is_empty()
            && self.format.is_empty()
            && self.alt.is_empty()
            && self.contig.is_empty()
            && self.samples.is_empty()
            && self.other.is_empty()
    }
}

// ---------------------------------------------------------------------------
// VcfReader - public entry point
// ---------------------------------------------------------------------------

pub struct VcfReader {
    backend: Backend,
    mode: ReaderMode,
    metadata: Option<VcfMetadata>,
}

impl VcfReader {
    fn open(file_path: &Path, tabix_path: &str, mode: ReaderMode) -> io::Result<Self> {
        let (backend, mode) = Backend::open(file_path, tabix_path, mode)?;

        Ok(Self {
            backend,
            mode,
            metadata: None,
        })
    }

    fn next(&mut self) -> io::Result<Option<VcfRecord>> {
        self.backend.next()
    }

    fn seek(&mut self, region: &str) -> io::Result<bool> {
        self.backend.seek(region)
    }

    fn metadata(&mut self) -> io::Result<VcfMetadata> {
        if self.metadata.is_none() {
            let header: &vcf::Header = match &self.backend {
                Backend::Sequential(reader) => &reader.header,
                Backend::Tabix(reader) => &reader.header,
            };
            self.metadata = Some(header_to_metadata(header));
        }

        Ok(self
            .metadata
            .as_ref()
            .expect("metadata should be cached")
            .clone())
    }
}

// ---------------------------------------------------------------------------
// Backend
// ---------------------------------------------------------------------------

enum Backend {
    Sequential(Box<SequentialBackend>),
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
                        .map(Box::new)
                        .map(|backend| (Self::Sequential(backend), ReaderMode::Sequential))
                }
            }
            ReaderMode::Sequential => SequentialBackend::new(file_path)
                .map(Box::new)
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

    fn next(&mut self) -> io::Result<Option<VcfRecord>> {
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
// SequentialBackend - simple reader backed by noodles-vcf
// ---------------------------------------------------------------------------

struct SequentialBackend {
    reader: vcf::io::Reader<Box<dyn BufRead + Send + Sync>>,
    header: vcf::Header,
}

impl SequentialBackend {
    fn new<P>(file_path: P) -> io::Result<Self>
    where
        P: AsRef<Path>,
    {
        let (reader, header) = get_header_from_file(file_path.as_ref())?;

        Ok(Self { reader, header })
    }

    fn seek(&mut self, _region: &str) -> io::Result<bool> {
        Err(io::Error::other(
            "VCF reader is set to 'sequential'; set reader to 'query' before seek",
        ))
    }

    fn next(&mut self) -> io::Result<Option<VcfRecord>> {
        read_next_record(&mut self.reader, &self.header)
    }
}

// ---------------------------------------------------------------------------
// TabixBackend - uses tabix for indexed reading
// ---------------------------------------------------------------------------

struct TabixBackend {
    reader: Option<tabix::TabixReader>,
    header: vcf::Header,
    query_active: bool,
}

impl TabixBackend {
    fn new<P>(file_path: P, tabix_path: &str) -> io::Result<Self>
    where
        P: AsRef<Path>,
    {
        // read header first
        let (_reader, header) = get_header_from_file(file_path.as_ref())?;

        Ok(Self {
            reader: Some(tabix::TabixReader::open(file_path.as_ref(), tabix_path)?),
            header,
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

    fn next(&mut self) -> io::Result<Option<VcfRecord>> {
        if self.query_active {
            let Some(reader) = self.reader.as_mut() else {
                return Err(io::Error::other("tabix reader is not initialized"));
            };

            loop {
                let Some(line) = reader.next()? else {
                    self.query_active = false;
                    return Ok(None);
                };

                if let Some(record) = parse_vcf_line(&line, &self.header)? {
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

#[pyclass(name = "VCFParser")]
pub struct VcfParser {
    file_path: Option<PathBuf>,
    reader: Option<VcfReader>,
    tabix_path: String,
    reader_mode: ReaderMode,
}

#[pymethods]
impl VcfParser {
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
        "vcf"
    }

    #[getter]
    fn metadata(&mut self) -> PyResult<VcfMetadata> {
        let Some(reader) = self.reader.as_mut() else {
            return Ok(VcfMetadata::default());
        };

        reader.metadata().map_err(to_py_error)
    }

    #[getter]
    fn reader_mode(&self) -> &'static str {
        self.reader_mode.as_str()
    }

    fn open(&mut self, file_path: &str, mode: Option<String>) -> PyResult<()> {
        let path = PathBuf::from(file_path);
        let reader_mode = ReaderMode::parse(mode.as_deref())?;

        self.reader =
            Some(VcfReader::open(&path, &self.tabix_path, reader_mode).map_err(to_py_error)?);
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
    fn next(&mut self) -> PyResult<Option<VcfRecord>> {
        self.reader()?.next().map_err(to_py_error)
    }

    #[pyo3(signature = (limit=10000))]
    fn records(&mut self, limit: Option<usize>) -> PyResult<Vec<VcfRecord>> {
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
    fn parse_text(data: &str, limit: Option<usize>) -> PyResult<Vec<VcfRecord>> {
        parse_vcf_bytes(data.as_bytes(), limit).map_err(to_py_error)
    }

    #[staticmethod]
    #[pyo3(signature = (file_path, limit=None))]
    fn parse_file(file_path: &str, limit: Option<usize>) -> PyResult<Vec<VcfRecord>> {
        let (mut reader, header) =
            get_header_from_file(Path::new(file_path)).map_err(to_py_error)?;

        read_records(&mut reader, &header, limit).map_err(to_py_error)
    }
}

impl VcfParser {
    fn reader(&mut self) -> PyResult<&mut VcfReader> {
        self.reader
            .as_mut()
            .ok_or_else(|| PyRuntimeError::new_err("VCF parser is not open"))
    }
}

// ---------------------------------------------------------------------------
// Helpers: Noodles VCF decoder
// ---------------------------------------------------------------------------

fn read_next_record<R>(
    reader: &mut vcf::io::Reader<R>,
    header: &vcf::Header,
) -> io::Result<Option<VcfRecord>>
where
    R: BufRead,
{
    let mut record = variant::RecordBuf::default();

    match reader.read_record_buf(header, &mut record)? {
        0 => Ok(None),
        _ => record_to_py(record).map(Some),
    }
}

fn read_records<R>(
    reader: &mut vcf::io::Reader<R>,
    header: &vcf::Header,
    limit: Option<usize>,
) -> io::Result<Vec<VcfRecord>>
where
    R: BufRead,
{
    let mut records = Vec::new();

    while limit.is_none_or(|n| records.len() < n) {
        match read_next_record(reader, header)? {
            Some(record) => records.push(record),
            None => break,
        }
    }

    Ok(records)
}

fn parse_vcf_bytes(data: &[u8], limit: Option<usize>) -> io::Result<Vec<VcfRecord>> {
    let (mut reader, header) = get_header_from_bytes(data)?;
    read_records(&mut reader, &header, limit)
}

fn parse_vcf_line(line: &str, header: &vcf::Header) -> io::Result<Option<VcfRecord>> {
    let mut data = Vec::with_capacity(line.len() + 1);
    data.extend_from_slice(line.as_bytes());
    data.push(b'\n');

    let mut vcf_reader = vcf::io::Reader::new(data.as_slice());
    read_next_record(&mut vcf_reader, header)
}

fn has_header_prefix<R: BufRead>(reader: &mut R) -> io::Result<bool> {
    Ok(reader.fill_buf()?.first() == Some(&b'#'))
}

fn get_header_from_file(
    file_path: &Path,
) -> io::Result<(vcf::io::Reader<Box<dyn BufRead + Send + Sync>>, vcf::Header)> {
    let mut reader = open_sequential_reader(file_path)?;
    if has_header_prefix(&mut reader)? {
        let mut vcf_reader = vcf::io::Reader::new(reader);
        let header = vcf_reader.read_header()?;

        return Ok((vcf_reader, header));
    }

    let mut first_line = String::new();
    reader.read_line(&mut first_line)?;
    let header = header_from_first_data_line(&first_line);
    let prefixed_reader = io::Cursor::new(first_line.into_bytes()).chain(reader);

    Ok((vcf::io::Reader::new(Box::new(prefixed_reader)), header))
}

fn get_header_from_bytes(data: &[u8]) -> io::Result<(vcf::io::Reader<&[u8]>, vcf::Header)> {
    let mut reader = data;
    let mut vcf_reader;
    let header: vcf::Header;
    if has_header_prefix(&mut reader)? {
        vcf_reader = vcf::io::Reader::new(data);
        header = vcf_reader.read_header()?;
    } else {
        let first_line = data.split(|byte| *byte == b'\n').next().unwrap_or(&[]);
        header = header_from_first_data_line(
            std::str::from_utf8(first_line)
                .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?,
        );
        vcf_reader = vcf::io::Reader::new(data);
    }
    Ok((vcf_reader, header))
}

fn header_from_first_data_line(line: &str) -> vcf::Header {
    let fields: Vec<&str> = line.trim_end_matches(['\n', '\r']).split('\t').collect();
    let mut builder = vcf::Header::builder().set_file_format(FileFormat::default());

    if let Some(info_field) = fields.get(7) {
        for entry in info_field
            .split(';')
            .filter(|entry| !entry.is_empty() && *entry != ".")
        {
            let key = entry.split_once('=').map(|(key, _)| key).unwrap_or(entry);
            let info: Map<Info> = Map::<Info>::from(key);
            builder = builder.add_info(key, info);
        }
    }

    if let Some(format_field) = fields.get(8) {
        for key in format_field.split(':').filter(|key| !key.is_empty()) {
            let format: Map<Format> = Map::<Format>::from(key);
            builder = builder.add_format(key, format);
        }
    }

    for sample_index in 0..fields.len().saturating_sub(9) {
        builder = builder.add_sample_name(format!("sample{sample_index}"));
    }

    builder.build()
}

fn record_to_py(record: variant::RecordBuf) -> io::Result<VcfRecord> {
    let id = record.ids().as_ref().iter().cloned().collect();
    let alt = record.alternate_bases().as_ref().to_vec();
    let filter = record.filters().as_ref().iter().cloned().collect();

    Ok(VcfRecord {
        chrom: bytes_to_string(record.reference_sequence_name().as_bytes()),
        pos: record.variant_start().map(|position| position.get()),
        id,
        reference: record.reference_bases().to_string(),
        alt,
        qual: record.quality_score(),
        filter,
        info: record.info().clone(),
        samples: record.samples().clone(),
    })
}

fn info_value_to_string(value: Option<&InfoFieldValue>) -> Option<String> {
    match value {
        None => None,
        Some(InfoFieldValue::Flag) => None,
        Some(InfoFieldValue::Integer(n)) => Some(n.to_string()),
        Some(InfoFieldValue::Float(n)) => Some(n.to_string()),
        Some(InfoFieldValue::Character(c)) => Some(c.to_string()),
        Some(InfoFieldValue::String(s)) => Some(s.clone()),
        Some(InfoFieldValue::Array(values)) => Some(info_array_to_string(values)),
    }
}

fn info_array_to_string(values: &ValueArray) -> String {
    match values {
        ValueArray::Integer(values) => join_optional_values(values.iter()),
        ValueArray::Float(values) => join_optional_values(values.iter()),
        ValueArray::Character(values) => join_optional_values(values.iter()),
        ValueArray::String(values) => join_optional_values(values.iter()),
    }
}

fn sample_value_to_string(value: &Option<SampleValue>) -> Option<String> {
    match value {
        None => None,
        Some(SampleValue::Integer(n)) => Some(n.to_string()),
        Some(SampleValue::Float(n)) => Some(n.to_string()),
        Some(SampleValue::Character(c)) => Some(c.to_string()),
        Some(SampleValue::String(s)) => Some(s.clone()),
        Some(SampleValue::Genotype(genotype)) => Some(genotype_to_string(genotype)),
        Some(SampleValue::Array(values)) => Some(sample_array_to_string(values)),
    }
}

fn sample_array_to_string(values: &SampleArray) -> String {
    match values {
        SampleArray::Integer(values) => join_optional_values(values.iter()),
        SampleArray::Float(values) => join_optional_values(values.iter()),
        SampleArray::Character(values) => join_optional_values(values.iter()),
        SampleArray::String(values) => join_optional_values(values.iter()),
    }
}

fn genotype_to_string(genotype: &SampleGenotype) -> String {
    use variant::record::samples::series::value::genotype::Phasing;

    let mut rendered = String::new();

    for allele in genotype.as_ref() {
        if !rendered.is_empty() {
            rendered.push(match allele.phasing() {
                Phasing::Phased => '|',
                Phasing::Unphased => '/',
            });
        }

        match allele.position() {
            Some(position) => rendered.push_str(&position.to_string()),
            None => rendered.push('.'),
        }
    }

    rendered
}

fn join_optional_values<'a, T, I>(values: I) -> String
where
    I: IntoIterator<Item = &'a Option<T>>,
    T: ToString + 'a,
{
    values
        .into_iter()
        .map(|value| match value {
            Some(value) => value.to_string(),
            None => ".".to_string(),
        })
        .collect::<Vec<_>>()
        .join(",")
}

// ---------------------------------------------------------------------------
// Helpers: VCF metadata conversion
// ---------------------------------------------------------------------------

fn header_to_metadata(header: &vcf::Header) -> VcfMetadata {
    VcfMetadata {
        file_format: file_format_to_string(header.file_format()),
        info: header
            .infos()
            .iter()
            .map(|(id, map)| (id.clone(), info_map_to_map(map)))
            .collect(),
        filter: header
            .filters()
            .iter()
            .map(|(id, map)| (id.clone(), filter_map_to_map(map)))
            .collect(),
        format: header
            .formats()
            .iter()
            .map(|(id, map)| (id.clone(), format_map_to_map(map)))
            .collect(),
        alt: header
            .alternative_alleles()
            .iter()
            .map(|(id, map)| (id.clone(), alt_map_to_map(map)))
            .collect(),
        contig: header
            .contigs()
            .iter()
            .map(|(id, map)| (id.clone(), contig_map_to_map(map)))
            .collect(),
        samples: header.sample_names().iter().cloned().collect(),
        other: header
            .other_records()
            .iter()
            .map(|(id, collection)| (id.to_string(), collection_to_strings(collection)))
            .collect(),
    }
}

fn file_format_to_string(file_format: FileFormat) -> String {
    format!("VCFv{}.{}", file_format.major(), file_format.minor())
}

fn info_map_to_map(map: &Map<Info>) -> HashMap<String, String> {
    let mut result = HashMap::new();
    result.insert("number".to_string(), info_number_to_string(map.number()));
    result.insert("type".to_string(), map.ty().to_string());
    result.insert("description".to_string(), map.description().to_string());
    result.insert(
        "idx".to_string(),
        map.idx().map(|idx| idx.to_string()).unwrap_or_default(),
    );
    result.insert(
        "other_fields".to_string(),
        other_fields_to_string(map.other_fields().iter()),
    );
    result
}

fn filter_map_to_map(map: &Map<Filter>) -> HashMap<String, String> {
    let mut result = HashMap::new();
    result.insert("description".to_string(), map.description().to_string());
    result.insert(
        "idx".to_string(),
        map.idx().map(|idx| idx.to_string()).unwrap_or_default(),
    );
    result.insert(
        "other_fields".to_string(),
        other_fields_to_string(map.other_fields().iter()),
    );
    result
}

fn format_map_to_map(map: &Map<Format>) -> HashMap<String, String> {
    let mut result = HashMap::new();
    result.insert("number".to_string(), format_number_to_string(map.number()));
    result.insert("type".to_string(), map.ty().to_string());
    result.insert("description".to_string(), map.description().to_string());
    result.insert(
        "idx".to_string(),
        map.idx().map(|idx| idx.to_string()).unwrap_or_default(),
    );
    result.insert(
        "other_fields".to_string(),
        other_fields_to_string(map.other_fields().iter()),
    );
    result
}

fn alt_map_to_map(map: &Map<AlternativeAllele>) -> HashMap<String, String> {
    let mut result = HashMap::new();
    result.insert("description".to_string(), map.description().to_string());
    result.insert(
        "other_fields".to_string(),
        other_fields_to_string(map.other_fields().iter()),
    );
    result
}

fn contig_map_to_map(map: &Map<Contig>) -> HashMap<String, String> {
    let mut result = HashMap::new();
    result.insert(
        "length".to_string(),
        map.length()
            .map(|length| length.to_string())
            .unwrap_or_default(),
    );
    result.insert("md5".to_string(), map.md5().unwrap_or_default().to_string());
    result.insert("url".to_string(), map.url().unwrap_or_default().to_string());
    result.insert(
        "idx".to_string(),
        map.idx().map(|idx| idx.to_string()).unwrap_or_default(),
    );
    result
}

fn collection_to_strings(collection: &Collection) -> Vec<String> {
    match collection {
        Collection::Unstructured(values) => values.clone(),
        Collection::Structured(map) => map
            .iter()
            .map(|(id, value)| {
                format!(
                    "{id}:{}",
                    other_fields_to_string(value.other_fields().iter())
                )
            })
            .collect(),
    }
}

fn other_fields_to_string<'a, I, K>(fields: I) -> String
where
    I: IntoIterator<Item = (&'a K, &'a String)>,
    K: Display + 'a,
{
    fields
        .into_iter()
        .map(|(key, value)| format!("{key}={value}"))
        .collect::<Vec<_>>()
        .join(";")
}

fn info_number_to_string(number: vcf::header::record::value::map::info::Number) -> String {
    use vcf::header::record::value::map::info::Number as InfoNumber;

    match number {
        InfoNumber::Count(n) => n.to_string(),
        InfoNumber::AlternateBases => "A".to_string(),
        InfoNumber::ReferenceAlternateBases => "R".to_string(),
        InfoNumber::Samples => "G".to_string(),
        InfoNumber::Unknown => ".".to_string(),
    }
}

fn format_number_to_string(number: vcf::header::record::value::map::format::Number) -> String {
    use vcf::header::record::value::map::format::Number as FormatNumber;

    match number {
        FormatNumber::Count(n) => n.to_string(),
        FormatNumber::AlternateBases => "A".to_string(),
        FormatNumber::ReferenceAlternateBases => "R".to_string(),
        FormatNumber::Samples => "G".to_string(),
        FormatNumber::LocalAlternateBases => "LA".to_string(),
        FormatNumber::LocalReferenceAlternateBases => "LR".to_string(),
        FormatNumber::LocalSamples => "LG".to_string(),
        FormatNumber::Ploidy => "P".to_string(),
        FormatNumber::BaseModifications => "M".to_string(),
        FormatNumber::Unknown => ".".to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        test_utils::{create_temp_dir, write_sample_file},
        utils::from_py_error,
    };
    use noodles_vcf::variant::{
        record::info::field::key as InfoKey,
        record::samples::keys::key as SampleKey,
        record_buf::{
            Samples,
            info::field::value::Array as ValueArray,
            samples::{Keys, sample::Value},
        },
    };
    use std::{fs, path::PathBuf};

    const HEADER: &[u8] = b"\
##fileformat=VCFv4.3
##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">
##FILTER=<ID=q10,Description=\"Quality below 10\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample0\tsample1
";
    const QUERY_ROW_1: &str = "chr1\t10\trs1\tA\tC\t13\tPASS\tNS=2\tGT\t0/1\t1/1";
    const QUERY_ROW_2: &str = "chr1\t20\trs2\tG\tT\t.\tq10\tNS=2\tGT\t1/0\t0/0";
    const QUERY_ROW_3: &str = "chr2\t100\trs3\tC\tA\t.\tPASS\tNS=2\tGT\t1/1\t0/1";

    const SAMPLE: &[u8] = b"\
##fileformat=VCFv4.3
##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">
##FILTER=<ID=q10,Description=\"Quality below 10\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample0\tsample1
chr1\t10\trs1\tA\tC\t13\tPASS\tNS=2\tGT\t0/1\t1/1
chr1\t20\trs2\tG\tT\t.\tq10\tNS=2\tGT\t1/0\t0/0
chr2\t100\trs3\tC\tA\t.\tPASS\tNS=2\tGT\t1/1\t0/1
";
    const NO_HEADER_SAMPLE: &[u8] = b"\
chr1\t10\trs1\tA\tC\t13\tPASS\tNS=2\tGT\t0/1\t1/1
chr1\t20\trs2\tG\tT\t.\tq10\tNS=2\tGT\t1/0\t0/0
chr2\t100\trs3\tC\tA\t.\tPASS\tNS=2\tGT\t1/1\t0/1
";

    const SAMPLE_VCF_LINE: &str =
        "chr1\t10\trs1\tA\tC\t13\tPASS\tNS=2\tGT:GQ:DP\t0|0:21:4\t1|0:34:10";

    fn create_test_header() -> io::Result<vcf::Header> {
        let ns: &str = InfoKey::SAMPLES_WITH_DATA_COUNT;
        let info: Map<Info> = Map::<Info>::from(ns);
        let gt: &str = SampleKey::GENOTYPE;
        let format: Map<Format> = Map::<Format>::from(gt);
        let test_header_obj: vcf::Header = vcf::Header::builder()
            .set_file_format(FileFormat::default())
            .add_info(ns, info.clone())
            .add_format(gt, format.clone())
            .add_sample_name("sample0")
            .add_sample_name("sample1")
            .build();
        Ok(test_header_obj)
    }

    fn make_samples() -> Samples {
        let keys: Keys = ["GT", "GQ", "DP"].into_iter().map(String::from).collect();

        Samples::new(
            keys,
            vec![
                vec![
                    Some(Value::Genotype("0|0".parse().expect("valid genotype"))),
                    Some(Value::from(21)),
                    Some(Value::from(4)),
                ],
                vec![
                    Some(Value::Genotype("1|0".parse().expect("valid genotype"))),
                    Some(Value::from(34)),
                    Some(Value::from(10)),
                ],
            ],
        )
    }

    #[cfg(unix)]
    fn write_mock_tabix(dir: &std::path::Path, body: &str) -> io::Result<PathBuf> {
        use std::os::unix::fs::PermissionsExt;

        let tabix = dir.join("tabix");
        fs::write(&tabix, format!("#!/bin/sh\n{body}\n"))?;
        let mut permissions = fs::metadata(&tabix)?.permissions();
        permissions.set_mode(0o755);
        fs::set_permissions(&tabix, permissions)?;
        Ok(tabix)
    }

    #[cfg(unix)]
    fn write_query_tabix(dir: &std::path::Path) -> io::Result<PathBuf> {
        write_mock_tabix(
            dir,
            &format!(
                r#"
if [ "$1" = "--version" ]; then
    exit 0
fi

case "$2" in
    chr1)
        printf '%s\n' '{QUERY_ROW_1}'
        printf '%s\n' '{QUERY_ROW_2}'
        ;;
    chr1:10-20)
        printf '%s\n' '{QUERY_ROW_1}'
        printf '%s\n' '{QUERY_ROW_2}'
        ;;
    chr2:100-120)
        printf '%s\n' '{QUERY_ROW_3}'
        ;;
    chr3:1-100)
        ;;
esac
"#
            ),
        )
    }

    // -----------------------------------------------------------------------
    // Helpers
    // -----------------------------------------------------------------------

    #[test]
    fn test_format_number_to_string() {
        use vcf::header::record::value::map::format::Number::*;

        assert_eq!(format_number_to_string(Count(3)), "3");
        assert_eq!(format_number_to_string(AlternateBases), "A");
        assert_eq!(format_number_to_string(ReferenceAlternateBases), "R");
        assert_eq!(format_number_to_string(Samples), "G");
        assert_eq!(format_number_to_string(LocalAlternateBases), "LA");
        assert_eq!(format_number_to_string(LocalReferenceAlternateBases), "LR");
        assert_eq!(format_number_to_string(LocalSamples), "LG");
        assert_eq!(format_number_to_string(Ploidy), "P");
        assert_eq!(format_number_to_string(BaseModifications), "M");
        assert_eq!(format_number_to_string(Unknown), ".");
    }

    #[test]
    fn test_info_number_to_string() {
        use vcf::header::record::value::map::info::Number as InfoNumber;

        assert_eq!(info_number_to_string(InfoNumber::Count(3)), "3");
        assert_eq!(info_number_to_string(InfoNumber::AlternateBases), "A");
        assert_eq!(
            info_number_to_string(InfoNumber::ReferenceAlternateBases),
            "R"
        );
        assert_eq!(info_number_to_string(InfoNumber::Samples), "G");
        assert_eq!(info_number_to_string(InfoNumber::Unknown), ".");
    }

    #[test]
    fn test_collection_to_strings() {
        let unstructured = Collection::Unstructured(vec!["one".to_string(), "two".to_string()]);
        assert_eq!(
            collection_to_strings(&unstructured),
            vec!["one".to_string(), "two".to_string()]
        );

        let mut other = Map::<noodles_vcf::header::record::value::map::Other>::new();
        other
            .other_fields_mut()
            .insert("NOODLES".parse().expect("valid tag"), "vcf".to_string());
        let structured =
            Collection::Structured([("record".to_string(), other)].into_iter().collect());

        assert_eq!(
            collection_to_strings(&structured),
            vec!["record:NOODLES=vcf".to_string()]
        );
    }

    #[test]
    fn test_contig_map_to_map() {
        let mut map = Map::<Contig>::new();
        *map.length_mut() = Some(248_956_422);
        *map.md5_mut() = Some("d7eba311421bbc9d3ada44709dd61534".to_string());
        *map.url_mut() = Some("https://example.com/chr1.fa".to_string());
        *map.idx_mut() = Some(7);

        let result = contig_map_to_map(&map);

        assert_eq!(result.get("length"), Some(&"248956422".to_string()));
        assert_eq!(
            result.get("md5"),
            Some(&"d7eba311421bbc9d3ada44709dd61534".to_string())
        );
        assert_eq!(
            result.get("url"),
            Some(&"https://example.com/chr1.fa".to_string())
        );
        assert_eq!(result.get("idx"), Some(&"7".to_string()));
    }

    #[test]
    fn test_alt_map_to_map() {
        let map = Map::<AlternativeAllele>::new("Deletion");
        let result = alt_map_to_map(&map);

        assert_eq!(result.get("description"), Some(&"Deletion".to_string()));
        assert_eq!(result.get("other_fields"), Some(&"".to_string()));
    }

    #[test]
    fn test_header_to_metadata() -> io::Result<()> {
        let mut reader = vcf::io::Reader::new(HEADER);
        let header = reader.read_header()?;
        let metadata = header_to_metadata(&header);

        assert_eq!(metadata.file_format, "VCFv4.3");
        assert_eq!(
            metadata.info,
            HashMap::from([(
                "NS".to_string(),
                HashMap::from([
                    ("number".to_string(), "1".to_string()),
                    ("type".to_string(), "Integer".to_string()),
                    (
                        "description".to_string(),
                        "Number of samples with data".to_string()
                    ),
                    ("idx".to_string(), "".to_string()),
                    ("other_fields".to_string(), "".to_string()),
                ])
            ),])
        );
        assert_eq!(
            metadata.filter,
            HashMap::from([(
                "q10".to_string(),
                HashMap::from([
                    ("description".to_string(), "Quality below 10".to_string()),
                    ("idx".to_string(), "".to_string()),
                    ("other_fields".to_string(), "".to_string()),
                ])
            ),])
        );
        assert_eq!(
            metadata.format,
            HashMap::from([(
                "GT".to_string(),
                HashMap::from([
                    ("number".to_string(), "1".to_string()),
                    ("type".to_string(), "String".to_string()),
                    ("description".to_string(), "Genotype".to_string()),
                    ("idx".to_string(), "".to_string()),
                    ("other_fields".to_string(), "".to_string()),
                ])
            ),])
        );
        assert_eq!(
            metadata.samples,
            vec!["sample0".to_string(), "sample1".to_string()]
        );
        Ok(())
    }

    #[test]
    fn test_sample_array_to_string() {
        let values = variant::record_buf::samples::sample::value::Array::Integer(vec![
            Some(1),
            None,
            Some(3),
        ]);
        assert_eq!(sample_array_to_string(&values), "1,.,3");

        let values = variant::record_buf::samples::sample::value::Array::Float(vec![
            Some(1.5),
            None,
            Some(3.25),
        ]);
        assert_eq!(sample_array_to_string(&values), "1.5,.,3.25");

        let values = variant::record_buf::samples::sample::value::Array::Character(vec![
            Some('A'),
            None,
            Some('Z'),
        ]);
        assert_eq!(sample_array_to_string(&values), "A,.,Z");

        let values = variant::record_buf::samples::sample::value::Array::String(vec![
            Some("AB".to_string()),
            None,
            Some("XY".to_string()),
        ]);
        assert_eq!(sample_array_to_string(&values), "AB,.,XY");
    }

    #[test]
    fn test_genotype_to_string() -> io::Result<()> {
        let unphased: variant::record_buf::samples::sample::value::Genotype =
            "0/1".parse().expect("expected unphased genotype");
        assert_eq!(genotype_to_string(&unphased), "0/1");

        let phased: variant::record_buf::samples::sample::value::Genotype =
            "0|1".parse().expect("expected phased genotype");
        assert_eq!(genotype_to_string(&phased), "0|1");

        let missing: variant::record_buf::samples::sample::value::Genotype =
            "./.".parse().expect("expected missing genotype");
        assert_eq!(genotype_to_string(&missing), "./.");

        Ok(())
    }

    #[test]
    fn test_sample_value_to_string() -> io::Result<()> {
        let none: Option<SampleValue> = None;
        assert_eq!(sample_value_to_string(&none), None);

        let integer = Some(SampleValue::Integer(1));
        assert_eq!(sample_value_to_string(&integer), Some("1".to_string()));

        let float = Some(SampleValue::Float(1.4));
        assert_eq!(sample_value_to_string(&float), Some("1.4".to_string()));

        let character = Some(SampleValue::Character('A'));
        assert_eq!(sample_value_to_string(&character), Some("A".to_string()));

        let genotype = Some(SampleValue::Genotype(
            "0|1".parse().expect("expected valid genotype"),
        ));
        assert_eq!(sample_value_to_string(&genotype), Some("0|1".to_string()));

        let array = Some(SampleValue::Array(
            variant::record_buf::samples::sample::value::Array::String(vec![
                Some("GT".to_string()),
                Some("GQ".to_string()),
                Some("DP".to_string()),
            ]),
        ));
        assert_eq!(sample_value_to_string(&array), Some("GT,GQ,DP".to_string()));

        Ok(())
    }

    #[test]
    fn test_info_array_to_string() {
        let values = ValueArray::Integer(vec![Some(1), None, Some(3)]);
        assert_eq!(info_array_to_string(&values), "1,.,3");

        let values = ValueArray::Float(vec![Some(1.5), None, Some(3.2)]);
        assert_eq!(info_array_to_string(&values), "1.5,.,3.2");

        let values = ValueArray::Character(vec![Some('A'), None, Some('T')]);
        assert_eq!(info_array_to_string(&values), "A,.,T");

        let values = ValueArray::String(vec![
            Some("alpha".to_string()),
            None,
            Some("gamma".to_string()),
        ]);
        assert_eq!(info_array_to_string(&values), "alpha,.,gamma");
    }

    #[test]
    fn test_info_value_to_string() -> io::Result<()> {
        assert_eq!(info_value_to_string(None), None);
        assert_eq!(info_value_to_string(Some(&InfoFieldValue::Flag)), None);
        assert_eq!(
            info_value_to_string(Some(&InfoFieldValue::Integer(1))),
            Some("1".to_string())
        );
        assert_eq!(
            info_value_to_string(Some(&InfoFieldValue::Float(1.4))),
            Some("1.4".to_string())
        );
        assert_eq!(
            info_value_to_string(Some(&InfoFieldValue::Character('A'))),
            Some("A".to_string())
        );
        assert_eq!(
            info_value_to_string(Some(&InfoFieldValue::String("GT".to_string()))),
            Some("GT".to_string())
        );
        assert_eq!(
            info_value_to_string(Some(&InfoFieldValue::Array(ValueArray::String(vec![
                Some("GT".to_string()),
                Some("GQ".to_string()),
                Some("DP".to_string()),
            ])))),
            Some("GT,GQ,DP".to_string())
        );

        Ok(())
    }

    // -----------------------------------------------------------------------
    // Header
    // -----------------------------------------------------------------------

    #[test]
    fn test_minimal_vcf_header() -> io::Result<()> {
        let minimal_vcf: &[u8] = b"\
##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample0\tsample1
chr1\t10\trs1\tA\tC\t13\tPASS\tNS=2\tGT\t0/1\t1/1
";
        let records = parse_vcf_bytes(minimal_vcf, None);
        assert_eq!(records?.len(), 1);

        let no_header_lines_vcf: &[u8] = b"\
chr1\t10\trs1\tA\tC\t13\tPASS\tNS=2\tGT\t0/1\t1/1
";
        let records = parse_vcf_bytes(no_header_lines_vcf, None)?;
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].chrom, "chr1");

        let no_format_line_vcf: &[u8] = b"\
##FILTER=<ID=q10,Description=\"Quality below 10\">
chr1\t10\trs1\tA\tC\t13\tPASS\tNS=2\tGT\t0/1\t1/1
";
        let err = get_header_from_bytes(no_format_line_vcf).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert!(
            err.to_string().contains("missing fileformat"),
            "unexpected error: {}",
            err
        );

        let no_chrom_line_vcf: &[u8] = b"\
##fileformat=VCFv4.3
chr1\t10\trs1\tA\tC\t13\tPASS\tNS=2\tGT\t0/1\t1/1
";
        let err = get_header_from_bytes(no_chrom_line_vcf).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert!(
            err.to_string().contains("missing header"),
            "unexpected error: {}",
            err
        );

        let no_version_line_vcf: &[u8] = b"\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample0\tsample1
chr1\t10\trs1\tA\tC\t13\tPASS\tNS=2\tGT\t0/1\t1/1
";
        let err = get_header_from_bytes(no_version_line_vcf).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert!(
            err.to_string().contains("invalid record"),
            "unexpected error: {}",
            err
        );

        Ok(())
    }

    // -----------------------------------------------------------------------
    // parse_vcf_line
    // -----------------------------------------------------------------------

    #[test]
    fn test_parse_vcf_line() -> io::Result<()> {
        let header = create_test_header()?;
        let record = parse_vcf_line(SAMPLE_VCF_LINE, &header)?.expect("expected record");

        assert_eq!(record.chrom, "chr1");
        assert_eq!(record.pos, Some(10));
        assert_eq!(record.id, vec!["rs1".to_string()]);
        assert_eq!(record.reference, "A");
        assert_eq!(record.alt, vec!["C".to_string()]);
        assert_eq!(record.qual, Some(13.0));
        assert_eq!(record.filter, vec!["PASS".to_string()]);
        assert_eq!(
            record.info.get("NS"),
            Some(Some(&InfoFieldValue::Integer(2)))
        );
        let samples = make_samples();
        assert_eq!(record.samples, samples);

        Ok(())
    }

    #[test]
    fn test_parse_vcf_line_invalid_line_bad_fields() -> io::Result<()> {
        let header = create_test_header()?;

        // Allowed: Space in CHROM
        let _ = parse_vcf_line(
            "chr 1\t10\trs1\tA\tC\t13\tPASS\tNS=2\tGT\t0/1\t1/1",
            &header,
        );

        // Non-integer POS
        let err = parse_vcf_line(
            "chr1\tbad_pos\trs1\tA\tC\t13\tPASS\tNS=2\tGT\t0/1\t1/1",
            &header,
        )
        .unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);

        // Allowed: Int for REF and ALT
        let _ = parse_vcf_line("chr1\t10\trs1\t2\tC\t13\tPASS\tNS=2\tGT\t0/1\t1/1", &header);
        let _ = parse_vcf_line("chr1\t10\trs1\tA\t3\t13\tPASS\tNS=2\tGT\t0/1\t1/1", &header);

        // String for QUAL
        let err = parse_vcf_line(
            "chr1\t10\trs1\tA\tC\tONE\tPASS\tNS=2\tGT\t0/1\t1/1",
            &header,
        )
        .unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);

        // Incorrectly formatted INFO
        let err = parse_vcf_line(
            "chr1\t10\trs1\tA\tC\t13.5\tPASS\tNS=\tGT\t0/1\t1/1",
            &header,
        )
        .unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        Ok(())
    }

    // -----------------------------------------------------------------------
    // parse_vcf_bytes
    // -----------------------------------------------------------------------

    #[test]
    fn test_parse_vcf_bytes() -> io::Result<()> {
        let records = parse_vcf_bytes(SAMPLE, None)?;

        assert_eq!(records.len(), 3);
        assert_eq!(records[0].chrom, "chr1");
        assert_eq!(records[0].pos, Some(10));

        Ok(())
    }

    #[test]
    fn test_parse_vcf_bytes_without_header() -> io::Result<()> {
        let records = parse_vcf_bytes(NO_HEADER_SAMPLE, None)?;

        assert_eq!(records.len(), 3);
        assert_eq!(records[0].chrom, "chr1");
        assert_eq!(records[1].pos, Some(20));

        Ok(())
    }
    // -----------------------------------------------------------------------
    // VcfParser
    // -----------------------------------------------------------------------

    #[test]
    fn test_vcf_parser_metadata_empty_before_open() -> io::Result<()> {
        // create new VcfParser without file path
        let mut parser = VcfParser::new(None, None, Some("sequential".to_string()))?;

        assert_eq!(parser.metadata()?, VcfMetadata::default());

        Ok(())
    }

    #[test]
    fn test_vcf_parser_parse_text_returns_records() -> io::Result<()> {
        let records = VcfParser::parse_text(
            std::str::from_utf8(SAMPLE).expect("sample should be valid UTF-8"),
            None,
        )?;

        assert_eq!(records.len(), 3);
        assert_eq!(records[0].chrom, "chr1");
        assert_eq!(records[1].pos, Some(20));
        assert_eq!(records[2].chrom, "chr2");

        Ok(())
    }

    #[test]
    fn test_vcf_parser_parse_text_without_header_returns_records() -> io::Result<()> {
        let records = VcfParser::parse_text(
            std::str::from_utf8(NO_HEADER_SAMPLE).expect("sample should be valid UTF-8"),
            None,
        )?;

        assert_eq!(records.len(), 3);
        assert_eq!(records[0].chrom, "chr1");
        assert_eq!(records[1].pos, Some(20));
        assert_eq!(records[2].chrom, "chr2");

        Ok(())
    }

    #[test]
    fn test_vcf_parser_parse_text_respects_limit() -> io::Result<()> {
        let records = VcfParser::parse_text(
            std::str::from_utf8(SAMPLE).expect("sample should be valid UTF-8"),
            Some(2),
        )?;

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].pos, Some(10));
        assert_eq!(records[1].pos, Some(20));

        Ok(())
    }

    #[test]
    fn test_vcf_parser_parse_file_returns_records() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-vcf-parse-file")?;
        let path = write_sample_file(&dir, "sample.vcf", SAMPLE)?;

        let records = VcfParser::parse_file(path.to_str().expect("utf-8 path"), None)?;

        assert_eq!(records.len(), 3);
        assert_eq!(records[0].chrom, "chr1");
        assert_eq!(records[1].pos, Some(20));
        assert_eq!(records[2].chrom, "chr2");

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_vcf_parser_parse_file_without_header_returns_records() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-vcf-parse-file-no-header")?;
        let path = write_sample_file(&dir, "sample.vcf", NO_HEADER_SAMPLE)?;

        let records = VcfParser::parse_file(path.to_str().expect("utf-8 path"), None)?;

        assert_eq!(records.len(), 3);
        assert_eq!(records[0].chrom, "chr1");
        assert_eq!(records[1].pos, Some(20));
        assert_eq!(records[2].chrom, "chr2");

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_vcf_parser_parse_file_respects_limit() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-vcf-parse-file-limit")?;
        let path = write_sample_file(&dir, "sample.vcf", SAMPLE)?;

        let records = VcfParser::parse_file(path.to_str().expect("utf-8 path"), Some(1))?;

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].pos, Some(10));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_vcf_parser_open_close_and_records() -> io::Result<()> {
        Python::initialize();
        let dir = create_temp_dir("ensembl-vep-vcf-parser-sequential")?;
        let path = write_sample_file(&dir, "sample.vcf", SAMPLE)?;
        let mut parser = VcfParser::new(None, None, Some("sequential".to_string()))?;

        assert_eq!(parser.format(), "vcf");
        assert_eq!(parser.reader_mode(), "sequential");

        parser.open(
            path.to_str().expect("utf-8 path"),
            Some("sequential".to_string()),
        )?;
        assert_eq!(parser.metadata()?.file_format, "VCFv4.3");

        let record = parser.next()?.expect("expected record");
        assert_eq!(record.chrom, "chr1");
        assert_eq!(parser.records(Some(1))?.len(), 1);
        assert_eq!(parser.records(None)?.len(), 1);

        parser.close();
        assert!(parser.reader.is_none());

        let err = parser.next().unwrap_err();
        assert!(err.to_string().contains("not open"));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[cfg(unix)]
    #[test]
    fn test_vcf_parser_query_mode_seek_and_records() -> io::Result<()> {
        Python::initialize();
        let dir = create_temp_dir("ensembl-vep-vcf-parser-query")?;
        let path = write_sample_file(&dir, "sample.vcf.gz", SAMPLE)?;
        let tbi_path = PathBuf::from(format!("{}.tbi", path.display()));
        fs::write(&tbi_path, b"")?;
        let tabix = write_query_tabix(&dir)?;
        let mut parser = VcfParser::new(
            Some(path.to_str().expect("utf-8 path").to_string()),
            Some(tabix.to_str().expect("utf-8 path").to_string()),
            Some("auto".to_string()),
        )?;

        assert_eq!(parser.reader_mode(), "query");
        assert_eq!(parser.metadata()?.file_format, "VCFv4.3");

        assert!(parser.seek("chr1:10-20")?);
        let first = parser.next()?.expect("expected first record");
        assert_eq!(first.pos, Some(10));
        assert_eq!(parser.records(Some(1))?.len(), 1);
        assert_eq!(parser.next()?, None);

        assert!(!parser.seek("chr3:1-100")?);
        assert_eq!(parser.next()?, None);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // Backend
    // -----------------------------------------------------------------------

    #[test]
    fn test_backend_open_auto_uses_sequential_without_index() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-vcf-backend-auto-seq")?;
        let path = write_sample_file(&dir, "sample.vcf", SAMPLE)?;

        let (backend, mode) = Backend::open(&path, "tabix", ReaderMode::Auto)?;

        assert_eq!(mode, ReaderMode::Sequential);
        assert!(matches!(backend, Backend::Sequential(_)));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[cfg(unix)]
    #[test]
    fn test_backend_open_auto_uses_tabix_with_index() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-vcf-backend-auto-tabix")?;
        let path = write_sample_file(&dir, "sample.vcf.gz", SAMPLE)?;
        let tbi_path = PathBuf::from(format!("{}.tbi", path.display()));
        fs::write(&tbi_path, b"")?;
        let tabix = write_query_tabix(&dir)?;

        let (backend, mode) =
            Backend::open(&path, tabix.to_str().expect("utf-8 path"), ReaderMode::Auto)?;

        assert_eq!(mode, ReaderMode::Query);
        assert!(matches!(backend, Backend::Tabix(_)));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_backend_open_query_without_index_errors() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-vcf-backend-query-no-index")?;
        let path = write_sample_file(&dir, "sample.vcf", SAMPLE)?;

        let err = match Backend::open(&path, "tabix", ReaderMode::Query) {
            Ok(_) => panic!("expected missing index error"),
            Err(err) => err,
        };
        assert_eq!(err.kind(), io::ErrorKind::NotFound);
        assert!(err.to_string().contains("no tabix index found"));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[cfg(unix)]
    #[test]
    fn test_backend_open_query_uses_tabix_with_index() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-vcf-backend-query-with-index")?;
        let path = write_sample_file(&dir, "sample.vcf.gz", SAMPLE)?;
        let tbi_path = PathBuf::from(format!("{}.tbi", path.display()));
        fs::write(&tbi_path, b"")?;
        let tabix = write_query_tabix(&dir)?;

        let (backend, mode) = Backend::open(
            &path,
            tabix.to_str().expect("utf-8 path"),
            ReaderMode::Query,
        )?;

        assert_eq!(mode, ReaderMode::Query);
        assert!(matches!(backend, Backend::Tabix(_)));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_sequential_backend_next_reads_and_seek_rejects() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-vcf-sequential-backend")?;
        let path = write_sample_file(&dir, "sample.vcf", SAMPLE)?;
        let mut backend = SequentialBackend::new(&path)?;

        let record = backend.next()?.expect("expected record");
        assert_eq!(record.chrom, "chr1");
        assert!(backend.next()?.is_some());
        assert!(backend.next()?.is_some());
        assert_eq!(backend.next()?, None);

        let err = backend.seek("chr1:1-10").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::Other);
        assert!(err.to_string().contains("sequential"));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_backend_seek_wrapper_rejects_on_sequential_backend() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-vcf-backend-seek-wrapper")?;
        let path = write_sample_file(&dir, "sample.vcf", SAMPLE)?;
        let mut backend = Backend::Sequential(Box::new(SequentialBackend::new(&path)?));

        let err = backend.seek("chr1:1-10").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::Other);
        assert!(err.to_string().contains("sequential"));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[cfg(unix)]
    #[test]
    fn test_tabix_backend_seek_and_next_reads_records() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-vcf-tabix-backend")?;
        let path = write_sample_file(&dir, "sample.vcf.gz", SAMPLE)?;
        let tbi_path = PathBuf::from(format!("{}.tbi", path.display()));
        fs::write(&tbi_path, b"")?;
        let tabix = write_query_tabix(&dir)?;
        let mut backend = TabixBackend::new(&path, tabix.to_str().expect("utf-8 path"))?;

        assert!(backend.seek("chr1:10-20")?);

        let first = backend.next()?.expect("expected first record");
        let second = backend.next()?.expect("expected second record");
        assert_eq!(first.chrom, "chr1");
        assert_eq!(first.pos, Some(10));
        assert_eq!(second.chrom, "chr1");
        assert_eq!(second.pos, Some(20));
        assert_eq!(backend.next()?, None);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // VcfReader
    // -----------------------------------------------------------------------

    #[test]
    fn test_vcf_reader_sequential_reads_records_and_metadata() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-vcf-reader-sequential")?;
        let path = write_sample_file(&dir, "sample.vcf", SAMPLE)?;
        let mut reader = VcfReader::open(&path, "tabix", ReaderMode::Sequential)?;

        assert_eq!(reader.mode, ReaderMode::Sequential);
        assert_eq!(reader.metadata()?.file_format, "VCFv4.3");

        let first = reader.next()?.expect("expected record");
        assert_eq!(first.chrom, "chr1");
        let second = reader.next()?.expect("expected second record");
        assert_eq!(second.pos, Some(20));
        let third = reader.next()?.expect("expected third record");
        assert_eq!(third.pos, Some(100));
        assert_eq!(reader.next()?, None);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_vcf_reader_metadata_is_cached() -> io::Result<()> {
        let dir = create_temp_dir("ensembl-vep-vcf-reader-metadata-cache")?;
        let path = write_sample_file(&dir, "sample.vcf", SAMPLE)?;
        let mut reader = VcfReader::open(&path, "tabix", ReaderMode::Sequential)?;

        let first = reader.metadata()?;
        let second = reader.metadata()?;

        assert_eq!(first, second);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // VcfMetadata
    // -----------------------------------------------------------------------

    #[test]
    fn test_vcf_metadata_default_is_empty() {
        let metadata = VcfMetadata::default();

        assert!(metadata.is_empty());
        assert!(metadata.keys().is_empty());
        assert_eq!(metadata.__len__(), 0);

        let metadata = VcfMetadata::new();
        assert!(metadata.is_empty());
    }

    #[test]
    fn test_vcf_metadata_exposes_all_top_level_fields() {
        let metadata = VcfMetadata {
            file_format: "VCFv4.3".to_string(),
            info: HashMap::from([(
                "NS".to_string(),
                HashMap::from([
                    ("number".to_string(), "1".to_string()),
                    ("type".to_string(), "Integer".to_string()),
                    (
                        "description".to_string(),
                        "Number of samples with data".to_string(),
                    ),
                ]),
            )]),
            filter: HashMap::from([(
                "q10".to_string(),
                HashMap::from([("description".to_string(), "Quality below 10".to_string())]),
            )]),
            format: HashMap::from([(
                "GT".to_string(),
                HashMap::from([
                    ("number".to_string(), "1".to_string()),
                    ("type".to_string(), "String".to_string()),
                    ("description".to_string(), "Genotype".to_string()),
                ]),
            )]),
            alt: HashMap::from([(
                "DEL".to_string(),
                HashMap::from([("description".to_string(), "Deletion".to_string())]),
            )]),
            contig: HashMap::from([(
                "chr1".to_string(),
                HashMap::from([
                    ("length".to_string(), "248956422".to_string()),
                    ("md5".to_string(), "md5".to_string()),
                    ("url".to_string(), "https://example.com/chr1.fa".to_string()),
                ]),
            )]),
            samples: vec!["sample0".to_string(), "sample1".to_string()],
            other: HashMap::from([("fileDate".to_string(), vec!["20260714".to_string()])]),
        };

        assert!(!metadata.is_empty());
        assert_eq!(metadata.__len__(), 8);
        assert_eq!(
            metadata.keys(),
            vec![
                "file_format".to_string(),
                "info".to_string(),
                "filter".to_string(),
                "format".to_string(),
                "alt".to_string(),
                "contig".to_string(),
                "samples".to_string(),
                "other".to_string(),
            ]
        );
        assert_eq!(metadata.file_format, "VCFv4.3");
        assert_eq!(metadata.samples, vec!["sample0", "sample1"]);
        assert_eq!(
            metadata.other.get("fileDate"),
            Some(&vec!["20260714".to_string()])
        );
    }

    #[test]
    fn test_vcf_metadata_getitem_returns_top_level_values() -> io::Result<()> {
        let metadata = VcfMetadata {
            file_format: "VCFv4.3".to_string(),
            info: HashMap::from([(
                "NS".to_string(),
                HashMap::from([
                    ("number".to_string(), "1".to_string()),
                    ("type".to_string(), "Integer".to_string()),
                ]),
            )]),
            filter: HashMap::from([(
                "q10".to_string(),
                HashMap::from([("description".to_string(), "Quality below 10".to_string())]),
            )]),
            format: HashMap::from([(
                "GT".to_string(),
                HashMap::from([
                    ("number".to_string(), "1".to_string()),
                    ("type".to_string(), "String".to_string()),
                ]),
            )]),
            alt: HashMap::from([(
                "DEL".to_string(),
                HashMap::from([("description".to_string(), "Deletion".to_string())]),
            )]),
            contig: HashMap::from([(
                "chr1".to_string(),
                HashMap::from([("length".to_string(), "248956422".to_string())]),
            )]),
            samples: vec!["sample0".to_string(), "sample1".to_string()],
            other: HashMap::from([("fileDate".to_string(), vec!["20260714".to_string()])]),
        };

        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            assert_eq!(
                metadata
                    .__getitem__(py, "file_format".to_string())?
                    .bind(py)
                    .extract::<String>()?,
                "VCFv4.3"
            );
            assert_eq!(
                metadata
                    .__getitem__(py, "samples".to_string())?
                    .bind(py)
                    .extract::<Vec<String>>()?,
                vec!["sample0".to_string(), "sample1".to_string()]
            );
            assert_eq!(
                metadata
                    .__getitem__(py, "other".to_string())?
                    .bind(py)
                    .extract::<HashMap<String, Vec<String>>>()?,
                HashMap::from([("fileDate".to_string(), vec!["20260714".to_string()],)])
            );

            let info = metadata
                .__getitem__(py, "info".to_string())?
                .bind(py)
                .extract::<HashMap<String, HashMap<String, String>>>()?;
            assert_eq!(
                info.get("NS").and_then(|v| v.get("number")),
                Some(&"1".to_string())
            );
            assert_eq!(
                info.get("NS").and_then(|v| v.get("type")),
                Some(&"Integer".to_string())
            );

            Ok(())
        })
        .map_err(from_py_error)?;

        Ok(())
    }

    #[test]
    fn test_vcf_metadata_getitem_unknown_key_errors() {
        let metadata = VcfMetadata::default();

        Python::initialize();
        Python::attach(|py| {
            let err = metadata.__getitem__(py, "missing".to_string()).unwrap_err();
            assert!(err.to_string().contains("unknown VCF metadata key"));
        });
    }

    // -----------------------------------------------------------------------
    // VcfRecord
    // -----------------------------------------------------------------------

    #[test]
    fn test_vcf_record_info_getter_returns_python_dict() -> io::Result<()> {
        let header = create_test_header()?;
        let record = parse_vcf_line(SAMPLE_VCF_LINE, &header)?.expect("expected record");

        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let info: HashMap<String, Option<String>> = record.info(py)?.extract()?;
            assert_eq!(
                info,
                HashMap::from([("NS".to_string(), Some("2".to_string()))])
            );
            Ok(())
        })
        .map_err(from_py_error)?;

        Ok(())
    }

    #[test]
    fn test_vcf_record_samples_getter_returns_python_list_of_dicts() -> io::Result<()> {
        let header = create_test_header()?;
        let record = parse_vcf_line(SAMPLE_VCF_LINE, &header)?.expect("expected record");

        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let samples: Vec<HashMap<String, Option<String>>> = record.samples(py)?.extract()?;

            assert_eq!(samples.len(), 2);
            assert_eq!(
                samples[0],
                HashMap::from([
                    ("GT".to_string(), Some("0|0".to_string())),
                    ("GQ".to_string(), Some("21".to_string())),
                    ("DP".to_string(), Some("4".to_string())),
                ])
            );
            assert_eq!(
                samples[1],
                HashMap::from([
                    ("GT".to_string(), Some("1|0".to_string())),
                    ("GQ".to_string(), Some("34".to_string())),
                    ("DP".to_string(), Some("10".to_string())),
                ])
            );

            Ok(())
        })
        .map_err(from_py_error)?;

        Ok(())
    }

    #[test]
    fn test_vcf_record_repr() -> io::Result<()> {
        let header = create_test_header()?;
        let record = parse_vcf_line(SAMPLE_VCF_LINE, &header)?.expect("expected record");

        assert_eq!(
            record.__repr__(),
            "VCFRecord(chrom=\"chr1\", pos=Some(10), id=[\"rs1\"], reference=\"A\", alt=[\"C\"], qual=Some(13.0), filter=[\"PASS\"])"
        );

        Ok(())
    }
}
