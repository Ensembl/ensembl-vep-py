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

use std::{
    fs::File,
    io::{self, BufRead, BufReader, ErrorKind},
    path::{Path, PathBuf},
    process::{Child, ChildStdout, Command, Stdio},
};

use log::warn;
use noodles_bgzf::{self as bgzf, io::Seek as BgzfSeek};
use noodles_core::{Position, Region, region::Interval};
use noodles_csi::{
    self as csi, BinningIndex,
    binning_index::index::{
        Header, header::format::CoordinateSystem, reference_sequence::bin::Chunk,
    },
};
use noodles_tabix as tabix;
use pyo3::{
    exceptions::{PyIOError, PyRuntimeError, PyValueError},
    prelude::*,
    types::PyDict,
};

// ---------------------------------------------------------------------------
// Type aliases
// ---------------------------------------------------------------------------

type CsiReader = csi::io::IndexedReader<bgzf::io::Reader<File>, csi::Index>;
type TbiReader = csi::io::IndexedReader<bgzf::io::Reader<File>, tabix::Index>;

// ---------------------------------------------------------------------------
// TabixReader — public entry point
// ---------------------------------------------------------------------------

/// A reusable tabix reader with parser-like `seek` and `next` operations.
///
/// Prefers the system `tabix` executable for parity with htslib. Falls back
/// to noodles when the executable is unavailable, keeping the BGZF file and
/// index open for repeated region queries.
pub struct TabixReader {
    backend: Backend,
    buffered_line: Option<String>,
}

impl TabixReader {
    pub fn open<P: AsRef<Path>>(src: P, tabix_path: &str) -> io::Result<Self> {
        Ok(Self {
            backend: Backend::open(src.as_ref(), tabix_path)?,
            buffered_line: None,
        })
    }

    /// Seek to a region and prime the reader.
    ///
    /// Returns `true` when at least one row is available. The first row is
    /// buffered and returned by the next call to [`next`](Self::next).
    pub fn seek(&mut self, region: &str) -> io::Result<bool> {
        self.buffered_line = None;
        self.backend.seek(region)?;

        if let Some(line) = self.next()? {
            self.buffered_line = Some(line);
            Ok(true)
        } else {
            Ok(false)
        }
    }

    /// Return the next raw line from the current region query, or `None` when
    /// the region is exhausted.
    pub fn next(&mut self) -> io::Result<Option<String>> {
        if let Some(line) = self.buffered_line.take() {
            return Ok(Some(line));
        }
        self.backend.next()
    }
}

// ---------------------------------------------------------------------------
// Backend
// ---------------------------------------------------------------------------

enum Backend {
    System(SystemBackend),
    Csi(NoodlesBackend<csi::Index>),
    Tbi(NoodlesBackend<tabix::Index>),
}

impl Backend {
    fn open(src: &Path, tabix_path: &str) -> io::Result<Self> {
        if SystemBackend::is_available(tabix_path)? {
            return Ok(Self::System(SystemBackend::new(src, tabix_path)));
        }

        if associated_index_path(src, "csi").exists() {
            let index = csi::fs::read(associated_index_path(src, "csi"))?;
            let file = File::open(src)?;
            return Ok(Self::Csi(NoodlesBackend::new(
                csi::io::IndexedReader::new(file, index),
            )));
        }

        let reader = tabix::io::indexed_reader::Builder::default().build_from_path(src)?;
        Ok(Self::Tbi(NoodlesBackend::new(reader)))
    }

    fn seek(&mut self, region: &str) -> io::Result<()> {
        match self {
            Self::System(b) => b.seek(region),
            Self::Csi(b) => b.seek(region),
            Self::Tbi(b) => b.seek(region),
        }
    }

    fn next(&mut self) -> io::Result<Option<String>> {
        match self {
            Self::System(b) => b.next(),
            Self::Csi(b) => b.next(),
            Self::Tbi(b) => b.next(),
        }
    }
}

// ---------------------------------------------------------------------------
// SystemBackend
// ---------------------------------------------------------------------------

struct SystemBackend {
    src: PathBuf,
    tabix_path: String,
    query: Option<SystemQuery>,
}

struct SystemQuery {
    child: Child,
    stdout: BufReader<ChildStdout>,
}

impl SystemBackend {
    fn new(src: &Path, tabix_path: &str) -> Self {
        Self {
            src: src.to_path_buf(),
            tabix_path: tabix_path.to_string(),
            query: None,
        }
    }

    fn is_available(tabix_path: &str) -> io::Result<bool> {
        match Command::new(tabix_path).arg("--version").output() {
            Ok(_) => Ok(true),
            Err(e) if e.kind() == ErrorKind::NotFound => {
                warn!("tabix executable not found: {tabix_path}");
                Ok(false)
            }
            Err(e) => Err(e),
        }
    }

    fn seek(&mut self, region: &str) -> io::Result<()> {
        self.close_query();

        let mut child = Command::new(&self.tabix_path)
            .arg(&self.src)
            .arg(region)
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .spawn()?;

        let stdout = child
            .stdout
            .take()
            .ok_or_else(|| io::Error::other("failed to capture tabix stdout"))?;

        self.query = Some(SystemQuery {
            child,
            stdout: BufReader::new(stdout),
        });

        Ok(())
    }

    fn next(&mut self) -> io::Result<Option<String>> {
        let Some(query) = self.query.as_mut() else {
            return Ok(None);
        };

        let mut line = String::new();
        let len = query.stdout.read_line(&mut line)?;

        if len == 0 {
            let status = query.child.wait()?;
            self.query = None;

            return if status.success() {
                Ok(None)
            } else {
                Err(io::Error::other(format!(
                    "tabix executable failed with status {status}"
                )))
            };
        }

        trim_newline(&mut line);
        Ok(Some(line))
    }

    fn close_query(&mut self) {
        if let Some(mut q) = self.query.take() {
            let _ = q.child.kill();
            let _ = q.child.wait();
        }
    }
}

impl Drop for SystemBackend {
    fn drop(&mut self) {
        self.close_query();
    }
}

// ---------------------------------------------------------------------------
// NoodlesBackend
// ---------------------------------------------------------------------------

struct NoodlesBackend<I: BinningIndex> {
    reader: csi::io::IndexedReader<bgzf::io::Reader<File>, I>,
    query: Option<NoodlesQuery>,
}

struct NoodlesQuery {
    region: Region,
    header: QueryHeader,
    chunks: std::vec::IntoIter<Chunk>,
    current_chunk_end: Option<bgzf::VirtualPosition>,
    line_buf: String,
}

impl<I: BinningIndex> NoodlesBackend<I> {
    fn new(reader: csi::io::IndexedReader<bgzf::io::Reader<File>, I>) -> Self {
        Self {
            reader,
            query: None,
        }
    }

    fn seek(&mut self, region: &str) -> io::Result<()> {
        let region = parse_region(region)?;

        let header = self
            .reader
            .index()
            .header()
            .map(QueryHeader::from)
            .ok_or_else(|| io::Error::new(ErrorKind::InvalidInput, "missing index header"))?;

        let reference_sequence_id = header
            .reference_sequence_id(self.reader.index(), &region)?;

        let chunks = self
            .reader
            .index()
            .query(reference_sequence_id, region.interval())?;

        self.query = Some(NoodlesQuery {
            region,
            header,
            chunks: chunks.into_iter(),
            current_chunk_end: None,
            line_buf: String::new(),
        });

        Ok(())
    }

    fn next(&mut self) -> io::Result<Option<String>> {
        let Some(query) = self.query.as_mut() else {
            return Ok(None);
        };

        let bgzf_reader = self.reader.get_mut();

        loop {
            // Advance to next chunk if current one is exhausted
            if query
                .current_chunk_end
                .is_none_or(|end| bgzf_reader.virtual_position() >= end)
            {
                let Some(chunk) = query.chunks.next() else {
                    self.query = None;
                    return Ok(None);
                };
                bgzf_reader.seek_to_virtual_position(chunk.start())?;
                query.current_chunk_end = Some(chunk.end());
            }

            query.line_buf.clear();
            let len = bgzf_reader.read_line(&mut query.line_buf)?;

            if len == 0 {
                self.query = None;
                return Ok(None);
            }

            trim_newline(&mut query.line_buf);

            if query.header.line_intersects_region(&query.line_buf, &query.region)? {
                return Ok(Some(query.line_buf.clone()));
            }
        }
    }
}

// ---------------------------------------------------------------------------
// QueryHeader — index metadata + region filtering logic
// ---------------------------------------------------------------------------

#[derive(Clone, Copy)]
struct QueryHeader {
    line_comment_prefix: u8,
    reference_sequence_name_index: usize,
    start_position_index: usize,
    end_position_index: Option<usize>,
    coordinate_system: CoordinateSystem,
}

impl QueryHeader {
    fn reference_sequence_id<I: BinningIndex>(
        &self,
        index: &I,
        region: &Region,
    ) -> io::Result<usize> {
        index
            .header()
            .and_then(|h| h.reference_sequence_names().get_index_of(region.name()))
            .ok_or_else(|| {
                io::Error::new(
                    ErrorKind::InvalidInput,
                    format!(
                        "missing reference sequence name: {}",
                        String::from_utf8_lossy(region.name())
                    ),
                )
            })
    }

    fn line_intersects_region(&self, line: &str, region: &Region) -> io::Result<bool> {
        // Skip comment lines
        if line
            .as_bytes()
            .first()
            .is_some_and(|b| *b == self.line_comment_prefix)
        {
            return Ok(false);
        }

        let fields: Vec<_> = line.split('\t').collect();

        // Check reference sequence name
        let seq_name = fields
            .get(self.reference_sequence_name_index)
            .ok_or_else(|| {
                io::Error::new(ErrorKind::InvalidData, "missing reference sequence name")
            })?;

        if seq_name.as_bytes() != region.name() {
            return Ok(false);
        }

        // Parse start position
        let start = fields
            .get(self.start_position_index)
            .ok_or_else(|| io::Error::new(ErrorKind::InvalidData, "missing start position"))
            .and_then(|s| self.parse_start_position(s))?;

        // Parse end position (fall back to start for point features)
        let end = match self.end_position_index {
            Some(i) => fields
                .get(i)
                .ok_or_else(|| io::Error::new(ErrorKind::InvalidData, "missing end position"))
                .and_then(|s| parse_end_position(s))?,
            None => start,
        };

        Ok(Interval::from(start..=end).intersects(region.interval()))
    }

    fn parse_start_position(&self, s: &str) -> io::Result<Position> {
        let raw = s
            .parse::<usize>()
            .map_err(|e| io::Error::new(ErrorKind::InvalidData, e))?;

        let pos = match self.coordinate_system {
            CoordinateSystem::Gff => raw,
            CoordinateSystem::Bed => raw
                .checked_add(1)
                .ok_or_else(|| io::Error::new(ErrorKind::InvalidData, "position overflow"))?,
        };

        pos.try_into()
            .map_err(|e| io::Error::new(ErrorKind::InvalidData, e))
    }
}

impl From<&Header> for QueryHeader {
    fn from(header: &Header) -> Self {
        Self {
            line_comment_prefix: header.line_comment_prefix(),
            reference_sequence_name_index: header.reference_sequence_name_index(),
            start_position_index: header.start_position_index(),
            end_position_index: header.end_position_index(),
            coordinate_system: header.format().coordinate_system(),
        }
    }
}

// ---------------------------------------------------------------------------
// PyO3 wrapper
// ---------------------------------------------------------------------------

#[pyclass(name = "TabixParser")]
pub struct TabixParser {
    file_path: Option<PathBuf>,
    reader: Option<TabixReader>,
    tabix_path: String,
}

#[pymethods]
impl TabixParser {
    #[new]
    #[pyo3(signature = (file_path=None, tabix_path=None))]
    fn new(file_path: Option<String>, tabix_path: Option<String>) -> PyResult<Self> {
        let mut parser = Self {
            file_path: None,
            reader: None,
            tabix_path: tabix_path.unwrap_or_else(|| "tabix".to_string()),
        };

        if let Some(path) = file_path {
            parser.open(&path)?;
        }

        Ok(parser)
    }

    #[getter]
    fn format(&self) -> &'static str {
        "tabix"
    }

    #[getter]
    fn metadata<'py>(&self, py: Python<'py>) -> Bound<'py, PyDict> {
        PyDict::new(py)
    }

    fn open(&mut self, file_path: &str) -> PyResult<()> {
        let path = PathBuf::from(file_path);
        self.reader = Some(TabixReader::open(&path, &self.tabix_path).map_err(tabix_py_error)?);
        self.file_path = Some(path);
        Ok(())
    }

    fn close(&mut self) {
        self.reader = None;
    }

    fn seek(&mut self, region: &str) -> PyResult<bool> {
        self.reader()?.seek(region).map_err(tabix_py_error)
    }

    fn next(&mut self) -> PyResult<Option<String>> {
        self.reader()?.next().map_err(tabix_py_error)
    }

    #[pyo3(signature = (limit=None))]
    fn records(&mut self, limit: Option<usize>) -> PyResult<Vec<String>> {
        let mut records = Vec::new();

        while limit.is_none_or(|n| records.len() < n) {
            match self.next()? {
                Some(record) => records.push(record),
                None => break,
            }
        }

        Ok(records)
    }
}

impl TabixParser {
    fn reader(&mut self) -> PyResult<&mut TabixReader> {
        self.reader
            .as_mut()
            .ok_or_else(|| PyRuntimeError::new_err("Tabix parser is not open"))
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn parse_end_position(s: &str) -> io::Result<Position> {
    s.parse::<usize>()
        .map_err(|e| io::Error::new(ErrorKind::InvalidData, e))?
        .try_into()
        .map_err(|e| io::Error::new(ErrorKind::InvalidData, e))
}

fn parse_region(region: &str) -> io::Result<Region> {
    region.parse().map_err(|e| {
        io::Error::new(
            ErrorKind::InvalidInput,
            format!("invalid tabix region {region:?}: {e}"),
        )
    })
}

fn associated_index_path(src: &Path, extension: &str) -> PathBuf {
    let mut path = src.as_os_str().to_os_string();
    path.push(".");
    path.push(extension);
    PathBuf::from(path)
}

fn trim_newline(line: &mut String) {
    if line.ends_with('\n') {
        line.pop();
        if line.ends_with('\r') {
            line.pop();
        }
    }
}

fn tabix_py_error(err: io::Error) -> PyErr {
    match err.kind() {
        ErrorKind::InvalidData | ErrorKind::InvalidInput => PyValueError::new_err(err.to_string()),
        ErrorKind::NotFound => PyIOError::new_err(err.to_string()),
        _ => PyRuntimeError::new_err(err.to_string()),
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::{env, fs};

    #[test]
    fn test_associated_index_path() {
        assert_eq!(
            associated_index_path(Path::new("annotations.gff3.gz"), "csi"),
            PathBuf::from("annotations.gff3.gz.csi")
        );
    }

    #[test]
    fn test_line_intersects_region_for_gff_coordinates() -> io::Result<()> {
        let header = Header::builder()
            .set_reference_sequence_name_index(0)
            .set_start_position_index(3)
            .set_end_position_index(Some(4))
            .build();
        let qh = QueryHeader::from(&header);
        let region = parse_region("chr1:15-25")?;

        assert!(qh.line_intersects_region("chr1\tsrc\tgene\t10\t20\t.\t+\t.\tID=1", &region)?);
        assert!(!qh.line_intersects_region("chr2\tsrc\tgene\t10\t20\t.\t+\t.\tID=1", &region)?);
        assert!(!qh.line_intersects_region("#comment", &region)?);

        Ok(())
    }

    #[test]
    #[cfg(unix)]
    fn test_query_prefers_system_tabix() -> io::Result<()> {
        use std::os::unix::fs::PermissionsExt;

        let dir = env::temp_dir().join(format!("ensembl-vep-tabix-test-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir)?;

        let tabix_path = dir.join("tabix");
        fs::write(
            &tabix_path,
            "#!/bin/sh\nprintf '%s\\n' 'sq0\ttabix\tgene\t8\t13\t.\t+\t.\tID=1'\n",
        )?;

        let mut permissions = fs::metadata(&tabix_path)?.permissions();
        permissions.set_mode(0o755);
        fs::set_permissions(&tabix_path, permissions)?;

        let data_path = dir.join("sample.gff3.gz");
        fs::write(&data_path, [])?;

        let mut reader = TabixReader::open(&data_path, tabix_path.to_str().unwrap())?;
        assert!(reader.seek("sq0:8-13")?);
        assert_eq!(
            reader.next()?.as_deref(),
            Some("sq0\ttabix\tgene\t8\t13\t.\t+\t.\tID=1")
        );
        assert_eq!(reader.next()?, None);

        fs::remove_dir_all(dir)?;

        Ok(())
    }
}