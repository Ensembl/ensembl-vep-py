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
    ops::Bound as Bnd,
    path::{Path, PathBuf},
    process::{Child, ChildStdout, Command, Stdio},
};

use crate::{
    parsers::utils::{associated_index_path, parse_region},
    utils::{bytes_to_string, to_py_error, trim_newline},
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
use pyo3::{exceptions::PyRuntimeError, prelude::*, types::PyDict};

// ---------------------------------------------------------------------------
// TabixReader — public entry point
// ---------------------------------------------------------------------------

/// A reusable tabix reader with parser-like `seek` and `next` operations.
///
/// Uses noodles when a CSI or TBI index is present, keeping the BGZF file and
/// index open for repeated region queries. Falls back to the system `tabix`
/// executable only when no index is found.
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
    #[allow(clippy::should_implement_trait)]
    pub fn next(&mut self) -> io::Result<Option<String>> {
        if let Some(line) = self.buffered_line.take() {
            return Ok(Some(line));
        }
        self.backend.next()
    }

    /// Return all header lines from the file, lazily loaded on first call.
    pub fn header_lines(&mut self) -> io::Result<&[String]> {
        self.backend.header_lines()
    }

    /// Return the reference sequence names declared by the index, in index order.
    ///
    /// The system `tabix` fallback lists them once via `tabix -l` and caches the result.
    pub fn reference_names(&mut self) -> io::Result<Vec<String>> {
        self.backend.reference_names()
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
    /// Open the in-process noodles reader rather than the system `tabix` executable.
    ///
    /// Each system query spawns a `tabix` child that reloads the index, which costs
    /// milliseconds per region; annotation issues thousands of them. The executable is used
    /// only when no index is present or noodles cannot read it.
    fn open(src: &Path, tabix_path: &str) -> io::Result<Self> {
        let noodles = Self::open_noodles(src);
        if noodles.is_ok() {
            return noodles;
        }

        if SystemBackend::is_available(tabix_path)? {
            return Ok(Self::System(SystemBackend::new(src, tabix_path)));
        }

        noodles
    }

    fn open_noodles(src: &Path) -> io::Result<Self> {
        if associated_index_path(src, "csi").exists() {
            let index = csi::fs::read(associated_index_path(src, "csi"))?;
            let file = File::open(src)?;
            return Ok(Self::Csi(NoodlesBackend::new(csi::io::IndexedReader::new(
                file, index,
            ))));
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

    fn header_lines(&mut self) -> io::Result<&[String]> {
        match self {
            Self::System(b) => b.header_lines(),
            Self::Csi(b) => b.header_lines(),
            Self::Tbi(b) => b.header_lines(),
        }
    }

    fn reference_names(&mut self) -> io::Result<Vec<String>> {
        match self {
            Self::System(b) => b.reference_names().map(<[String]>::to_vec),
            Self::Csi(b) => Ok(b.reference_names()),
            Self::Tbi(b) => Ok(b.reference_names()),
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
    /// Lazily loaded header lines (`None` = not yet loaded).
    cached_headers: Option<Vec<String>>,
    /// Lazily loaded reference sequence names (`None` = not yet loaded).
    cached_reference_names: Option<Vec<String>>,
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
            cached_headers: None,
            cached_reference_names: None,
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

    /// Load header lines via `tabix -H <src>`, which prints only the header.
    fn header_lines(&mut self) -> io::Result<&[String]> {
        if self.cached_headers.is_none() {
            let output = Command::new(&self.tabix_path)
                .arg("-H")
                .arg(&self.src)
                .output()?;

            if !output.status.success() {
                return Err(io::Error::other(format!(
                    "tabix -H failed with status {}",
                    output.status
                )));
            }

            let lines = String::from_utf8(output.stdout)
                .map_err(|e| io::Error::new(ErrorKind::InvalidData, e))?
                .lines()
                .map(str::to_owned)
                .collect();

            self.cached_headers = Some(lines);
        }

        Ok(self.cached_headers.as_deref().unwrap())
    }

    /// Load reference sequence names via `tabix -l <src>`, which reads only the index.
    fn reference_names(&mut self) -> io::Result<&[String]> {
        if self.cached_reference_names.is_none() {
            let output = Command::new(&self.tabix_path)
                .arg("-l")
                .arg(&self.src)
                .output()?;

            if !output.status.success() {
                return Err(io::Error::other(format!(
                    "tabix -l failed with status {}",
                    output.status
                )));
            }

            let names = bytes_to_string(&output.stdout)
                .lines()
                .map(str::to_owned)
                .collect();

            self.cached_reference_names = Some(names);
        }

        Ok(self.cached_reference_names.as_deref().unwrap())
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
    /// Lazily loaded header lines (`None` = not yet loaded).
    cached_headers: Option<Vec<String>>,
    /// Lines read from BGZF for the active or most recently spent query, matched or not. Reset
    /// on each `seek`. Test-only: it exists so a test can prove the early-stop behaviour in
    /// `next` actually reduces the number of lines scanned, rather than merely returning the
    /// right rows -- a claim that happens to hold even without early stop whenever the past-end
    /// row is also the chunk's last line, as it is in the small fixtures most other tests share.
    /// It lives here rather than on `NoodlesQuery` because early stop now clears `query` (moving
    /// `PastEnd` handling from "abandon this chunk" to "end the whole query"), and a test needs
    /// the count to survive that.
    #[cfg(test)]
    scanned_lines: usize,
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
            cached_headers: None,
            #[cfg(test)]
            scanned_lines: 0,
        }
    }

    fn seek(&mut self, region: &str) -> io::Result<()> {
        // Clear any previous query up front. Every early return below (an inverted range, a
        // missing header, a reference sequence the index does not declare, or `index.query`
        // failing) must leave no active query rather than the previous seek's stale rows -- a
        // seek that produces no results must behave like `next` was never primed.
        self.query = None;
        #[cfg(test)]
        {
            self.scanned_lines = 0;
        }

        let region = parse_region(region)?;
        if let (Bnd::Included(start), Bnd::Included(end)) = (region.start(), region.end())
            && start > end
        {
            return Ok(());
        }

        let index = self.reader.index();

        let header = index
            .header()
            .map(QueryHeader::from)
            .ok_or_else(|| io::Error::new(ErrorKind::InvalidInput, "missing index header"))?;

        // A reference sequence absent from the index has no rows; treating that as "no data"
        // lets callers fall back cleanly instead of catching a broad exception.
        let Some(reference_sequence_id) = self.reference_sequence_id(&region) else {
            return Ok(());
        };

        let chunks = index.query(reference_sequence_id, region.interval())?;

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
            #[cfg(test)]
            {
                if len > 0 {
                    self.scanned_lines += 1;
                }
            }

            if len == 0 {
                self.query = None;
                return Ok(None);
            }

            trim_newline(&mut query.line_buf);

            match query.header.classify_line(&query.line_buf, &query.region)? {
                LineRegionMatch::Intersects => {
                    // Move the buffer out instead of cloning it. A clone's real cost is copying
                    // the line's bytes; moving avoids that copy, at the price of one fresh
                    // (empty) allocation to replace the buffer -- sized to this line's length, so
                    // reading a next line no longer than this one does not need to grow it again.
                    let capacity_hint = query.line_buf.len();
                    let line = std::mem::replace(
                        &mut query.line_buf,
                        String::with_capacity(capacity_hint),
                    );
                    return Ok(Some(line));
                }
                LineRegionMatch::PastEnd => {
                    // Safe to end the whole query here, not just this chunk. `Index::query`
                    // (noodles-csi 0.55, `binning_index/index.rs`) scopes every returned chunk to
                    // the one queried `reference_sequence_id` and runs them through
                    // `optimize_chunks`, which sorts by ascending virtual position and merges
                    // overlaps -- so `query.chunks` is already a non-overlapping, file-ordered
                    // sequence for this contig alone. The file is coordinate-sorted per contig,
                    // so virtual-position order here is genomic order: once a same-contig row
                    // starts past the region end, every row still to come -- in this chunk or a
                    // later one -- is further along the file and so past the end too.
                    self.query = None;
                    return Ok(None);
                }
                LineRegionMatch::NoMatch => {}
            }
        }
    }

    /// Read header lines by seeking to the start of the BGZF file and
    /// collecting all comment lines before the first data line.
    ///
    /// The BGZF position is restored to where it was before the call so that
    /// any active query is not disturbed. In practice `seek()` always resets
    /// the position via `seek_to_virtual_position`, so this is a belt-and-
    /// braces precaution.
    fn header_lines(&mut self) -> io::Result<&[String]> {
        if self.cached_headers.is_none() {
            let comment_prefix = self
                .reader
                .index()
                .header()
                .map(|h| h.line_comment_prefix())
                .unwrap_or(b'#');

            let bgzf_reader = self.reader.get_mut();

            // Save position so we can restore it afterwards.
            let saved_pos = bgzf_reader.virtual_position();

            // Seek to the very start of the BGZF stream.
            bgzf_reader.seek_to_virtual_position(bgzf::VirtualPosition::default())?;

            let mut headers = Vec::new();
            let mut line_buf = String::new();

            loop {
                line_buf.clear();
                let len = bgzf_reader.read_line(&mut line_buf)?;

                if len == 0 {
                    break;
                }

                trim_newline(&mut line_buf);

                if line_buf
                    .as_bytes()
                    .first()
                    .is_some_and(|b| *b == comment_prefix)
                {
                    headers.push(line_buf.clone());
                } else {
                    // First non-comment line — stop scanning.
                    break;
                }
            }

            // Restore the reader position.
            bgzf_reader.seek_to_virtual_position(saved_pos)?;

            self.cached_headers = Some(headers);
        }

        Ok(self.cached_headers.as_deref().unwrap())
    }

    /// Return the index of `region`'s reference sequence, or `None` when the index does not
    /// declare it. A missing name is a normal "no data for this contig" outcome, not an error.
    fn reference_sequence_id(&self, region: &Region) -> Option<usize> {
        self.reader
            .index()
            .header()?
            .reference_sequence_names()
            .get_index_of(region.name())
    }

    /// Return every reference sequence name declared by the index header, in index order.
    fn reference_names(&self) -> Vec<String> {
        self.reader
            .index()
            .header()
            .map(|header| {
                header
                    .reference_sequence_names()
                    .iter()
                    .map(|name| bytes_to_string(name.as_ref()))
                    .collect()
            })
            .unwrap_or_default()
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

/// The relationship between one tabix data line and the region being queried.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum LineRegionMatch {
    /// The line's interval intersects the queried region.
    Intersects,
    /// The line is on a different reference sequence, or falls before or after the region
    /// without ruling out further matches later in the file (for example an unbounded region).
    NoMatch,
    /// The line is on the queried reference sequence and starts after the region end. Since the
    /// file is coordinate-sorted per contig, nothing later in this chunk can match either.
    PastEnd,
}

impl QueryHeader {
    /// Classify one line against `region` in a single left-to-right pass over its tab-separated
    /// fields, without collecting them into a `Vec` first.
    fn classify_line(&self, line: &str, region: &Region) -> io::Result<LineRegionMatch> {
        if line
            .as_bytes()
            .first()
            .is_some_and(|b| *b == self.line_comment_prefix)
        {
            return Ok(LineRegionMatch::NoMatch);
        }

        // The three needed columns (reference sequence name, start, end) can appear in any
        // order: GFF/BED/VCF all place the name first, but a generic tabix index (`tabix -s -b
        // -e`) can point the name column anywhere, including after start/end. A single forward
        // pass still avoids collecting a `Vec`, but it can only stop once *all three* configured
        // indices have actually been seen -- stopping as soon as start/end are read, regardless
        // of whether the name column has been reached yet, would wrongly report the name as
        // missing whenever it comes later in the line.
        let mut seq_name_found = false;
        let mut seq_matches = false;
        let mut start_field = None;
        let mut end_field = None;

        for (index, field) in line.split('\t').enumerate() {
            if index == self.reference_sequence_name_index {
                seq_name_found = true;
                seq_matches = field.as_bytes() == region.name();
                if !seq_matches {
                    break;
                }
            }
            if index == self.start_position_index {
                start_field = Some(field);
            }
            if self.end_position_index == Some(index) {
                end_field = Some(field);
            }

            let scanned_every_needed_field = seq_name_found
                && index >= self.start_position_index
                && self
                    .end_position_index
                    .is_none_or(|end_index| index >= end_index);
            if scanned_every_needed_field {
                break;
            }
        }

        if !seq_name_found {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "missing reference sequence name",
            ));
        }
        if !seq_matches {
            return Ok(LineRegionMatch::NoMatch);
        }

        let start = start_field
            .ok_or_else(|| io::Error::new(ErrorKind::InvalidData, "missing start position"))
            .and_then(|s| self.parse_start_position(s))?;
        let end = match self.end_position_index {
            Some(_) => end_field
                .ok_or_else(|| io::Error::new(ErrorKind::InvalidData, "missing end position"))
                .and_then(|s| self.parse_end_position(s))?,
            None => start,
        };

        if Interval::from(start..=end).intersects(region.interval()) {
            return Ok(LineRegionMatch::Intersects);
        }

        let past_end = matches!(region.end(), Bnd::Included(region_end) if start > region_end);
        Ok(if past_end {
            LineRegionMatch::PastEnd
        } else {
            LineRegionMatch::NoMatch
        })
    }

    /// Test-only convenience over [`Self::classify_line`]: production code (`NoodlesBackend::next`)
    /// needs the richer `LineRegionMatch` to drive the early-stop decision, so only tests need the
    /// plain intersects/doesn't-intersect question.
    #[cfg(test)]
    fn line_intersects_region(&self, line: &str, region: &Region) -> io::Result<bool> {
        Ok(matches!(
            self.classify_line(line, region)?,
            LineRegionMatch::Intersects
        ))
    }

    fn parse_start_position(&self, start: &str) -> io::Result<Position> {
        let raw = start
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

    fn parse_end_position(&self, end: &str) -> io::Result<Position> {
        let raw = end
            .parse::<usize>()
            .map_err(|e| io::Error::new(ErrorKind::InvalidData, e))?;

        let pos = match self.coordinate_system {
            CoordinateSystem::Gff => raw,
            CoordinateSystem::Bed => raw.max(1),
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

// Note: we use a separate wrapper struct instead of implementing PyO3 traits directly on TabixReader
// because TabixReader's API is not designed for direct Python use and we want to keep it ergonomic for Rust usage.
// The wrapper delegates to TabixReader and translates errors into Python exceptions.
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
    pub fn new(file_path: Option<String>, tabix_path: Option<String>) -> PyResult<Self> {
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
    pub fn format(&self) -> &'static str {
        "tabix"
    }

    /// Return a metadata dict. The `"header"` key is lazily populated with
    /// header lines from the file on first access.
    #[getter]
    pub fn metadata<'py>(&mut self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);

        let lines = self.header_lines()?;
        if !lines.is_empty() {
            let py_lines: Vec<_> = lines.iter().map(|l| l.as_str()).collect();
            dict.set_item("header", py_lines)?;
        }

        Ok(dict)
    }

    pub fn open(&mut self, file_path: &str) -> PyResult<()> {
        let path = PathBuf::from(file_path);
        self.reader = Some(TabixReader::open(&path, &self.tabix_path).map_err(to_py_error)?);
        self.file_path = Some(path);
        Ok(())
    }

    pub fn close(&mut self) {
        self.reader = None;
    }

    pub fn seek(&mut self, region: &str) -> PyResult<bool> {
        self.reader()?.seek(region).map_err(to_py_error)
    }

    #[allow(clippy::should_implement_trait)]
    pub fn next(&mut self) -> PyResult<Option<String>> {
        self.reader()?.next().map_err(to_py_error)
    }

    #[pyo3(signature = (limit=None))]
    pub fn records(&mut self, limit: Option<usize>) -> PyResult<Vec<String>> {
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
    pub fn header_lines(&mut self) -> PyResult<Vec<String>> {
        if let Some(reader) = self.reader.as_mut() {
            return Ok(reader.header_lines().map_err(to_py_error)?.to_vec());
        }

        Ok(Vec::new())
    }

    fn reader(&mut self) -> PyResult<&mut TabixReader> {
        self.reader
            .as_mut()
            .ok_or_else(|| PyRuntimeError::new_err("Tabix parser is not open"))
    }
}

// ---------------------------------------------------------------------------
pub fn has_tabix_index(file_path: &Path) -> bool {
    associated_index_path(file_path, "tbi").exists()
        || associated_index_path(file_path, "csi").exists()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parsers::test_utils::{write_bed_with_tbi, write_bgzf_file, write_gff3_with_tbi};
    #[cfg(unix)]
    use crate::test_utils::{assert_process_is_gone, create_temp_dir, process_exists};
    use noodles_csi::binning_index::index::header::Format;
    use std::fs;

    // -----------------------------------------------------------------------
    // Helpers
    // -----------------------------------------------------------------------

    fn gff_query_header() -> QueryHeader {
        let header = Header::builder()
            .set_reference_sequence_name_index(0)
            .set_start_position_index(3)
            .set_end_position_index(Some(4))
            .build();
        QueryHeader::from(&header)
    }

    fn bed_query_header() -> QueryHeader {
        let header = Header::builder()
            .set_format(Format::Generic(CoordinateSystem::Bed))
            .set_reference_sequence_name_index(0)
            .set_start_position_index(1)
            .set_end_position_index(Some(2))
            .build();
        QueryHeader::from(&header)
    }

    /// A generic tabix index whose reference-name column comes *after* start/end, as produced by
    /// `tabix -s 3 -b 1 -e 2` (1-based CLI columns => 0-based indices seq=2, start=0, end=1).
    fn generic_late_contig_header() -> QueryHeader {
        let header = Header::builder()
            .set_reference_sequence_name_index(2)
            .set_start_position_index(0)
            .set_end_position_index(Some(1))
            .build();
        QueryHeader::from(&header)
    }

    #[cfg(unix)]
    fn get_system_query_child_id(reader: &TabixReader) -> u32 {
        match &reader.backend {
            Backend::System(backend) => backend
                .query
                .as_ref()
                .expect("expected active system query")
                .child
                .id(),
            _ => panic!("expected system backend"),
        }
    }

    // -----------------------------------------------------------------------
    // has_tabix_index
    // -----------------------------------------------------------------------

    #[test]
    fn test_has_tabix_index_exists() -> io::Result<()> {
        let dir = create_temp_dir("tabix-has-tabix")?;
        let (gz_path, _tbi_path) = write_gff3_with_tbi(&dir, None)?;

        assert!(has_tabix_index(&gz_path));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_has_tabix_index_missing_index() -> io::Result<()> {
        let dir = create_temp_dir("tabix-has-tabix-missing")?;
        let path = write_bgzf_file(
            &dir,
            vec!["chr1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1\n".to_string()],
            None,
        )?;

        assert!(!has_tabix_index(&path));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // associated_index_path
    // -----------------------------------------------------------------------

    #[test]
    fn test_associated_index_path_csi() {
        assert_eq!(
            associated_index_path(Path::new("annotations.gff3.gz"), "csi"),
            PathBuf::from("annotations.gff3.gz.csi")
        );
    }

    #[test]
    fn test_associated_index_path_tbi() {
        assert_eq!(
            associated_index_path(Path::new("sample.vcf.gz"), "tbi"),
            PathBuf::from("sample.vcf.gz.tbi")
        );
    }

    // -----------------------------------------------------------------------
    // TabixParser::reader
    // -----------------------------------------------------------------------
    #[test]
    fn test_tabix_parser_reader_initialization() -> io::Result<()> {
        let dir = create_temp_dir("tabix-test-reader-init")?;
        let (gz_path, _tbi_path) = write_gff3_with_tbi(&dir, None)?;

        let mut parser = TabixParser::new(
            Some(gz_path.to_string_lossy().to_string()),
            Some("tabix-nonexistent-binary-for-test".to_string()),
        )?;
        let reader = parser.reader()?;

        assert!(matches!(
            reader,
            TabixReader {
                backend: Backend::Tbi(_),
                buffered_line: None
            }
        ));
        Ok(())
    }

    #[test]
    fn test_tabix_parser_methods_without_open_reader() {
        Python::initialize();
        Python::attach(|py| {
            let mut parser = TabixParser::new(None, None).unwrap();

            assert_eq!(parser.format(), "tabix");
            assert!(parser.metadata(py).unwrap().is_empty());
            assert!(parser.header_lines().unwrap().is_empty());
            assert!(
                parser
                    .records(None)
                    .unwrap_err()
                    .to_string()
                    .contains("not open")
            );

            parser.close();
        });
    }

    #[test]
    fn test_tabix_parser_open_seek_records_and_close() -> io::Result<()> {
        let dir = create_temp_dir("tabix-test-parser-methods")?;
        let (gz_path, _tbi_path) = write_gff3_with_tbi(&dir, None)?;
        let file_path = gz_path.to_string_lossy().into_owned();

        let mut parser =
            TabixParser::new(None, Some("tabix-nonexistent-binary-for-test".to_string())).unwrap();
        parser.open(&file_path).unwrap();

        Python::initialize();
        Python::attach(|py| {
            let metadata = parser.metadata(py).unwrap();
            assert!(metadata.contains("header").unwrap());
        });

        assert!(parser.seek("chr1:1-300").unwrap());
        assert!(parser.next().unwrap().unwrap().contains("ID=gene1"));
        assert_eq!(parser.records(Some(2)).unwrap().len(), 2);
        assert_eq!(parser.records(None).unwrap(), Vec::<String>::new());
        assert_eq!(parser.records(Some(0)).unwrap(), Vec::<String>::new());

        parser.close();
        assert!(parser.next().unwrap_err().to_string().contains("not open"));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // QueryHeader::parse_end_position
    // -----------------------------------------------------------------------

    #[test]
    fn test_parse_end_position_valid() {
        let qh = gff_query_header();
        assert!(qh.parse_end_position("100").is_ok());

        let qh = bed_query_header();
        assert!(qh.parse_end_position("100").is_ok());
    }

    #[test]
    fn test_parse_end_position_zero() {
        // Position is 1-based so 0 is invalid
        let qh = gff_query_header();
        let err = qh.parse_end_position("0").unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidData);

        // BED end position can be 0 if start is also 0, inserion before 0 position
        // System tabix does not allow it but noodles backend does
        let qh = bed_query_header();
        assert_eq!(
            qh.parse_end_position("0").unwrap(),
            Position::try_from(1).unwrap()
        );
    }

    #[test]
    fn test_parse_end_position_non_numeric() {
        let qh = gff_query_header();
        let err = qh.parse_end_position("abc").unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidData);

        let qh = bed_query_header();
        let err = qh.parse_end_position("abc").unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidData);
    }

    #[test]
    fn test_parse_end_position_negative() {
        let qh = gff_query_header();
        let err = qh.parse_end_position("-1").unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidData);

        let qh = bed_query_header();
        let err = qh.parse_end_position("-1").unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidData);
    }

    // -----------------------------------------------------------------------
    // QueryHeader::parse_start_position
    // -----------------------------------------------------------------------

    #[test]
    fn test_parse_start_position_gff() {
        // GFF is 1-based — value passed through unchanged
        let qh = gff_query_header();
        let pos = qh.parse_start_position("10").unwrap();
        assert_eq!(usize::from(pos), 10);
    }

    #[test]
    fn test_parse_start_position_bed() {
        // BED is 0-based — value incremented to 1-based
        let qh = bed_query_header();
        let pos = qh.parse_start_position("9").unwrap();
        assert_eq!(usize::from(pos), 10);
    }

    #[test]
    fn test_parse_start_position_gff_zero_is_invalid() {
        // Position 0 is invalid in 1-based coordinates
        let qh = gff_query_header();
        let err = qh.parse_start_position("0").unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidData);
    }

    #[test]
    fn test_parse_start_position_non_numeric() {
        let qh = gff_query_header();
        let err = qh.parse_start_position("abc").unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidData);
    }

    // -----------------------------------------------------------------------
    // QueryHeader::line_intersects_region — GFF (1-based)
    // -----------------------------------------------------------------------

    #[test]
    fn test_line_intersects_region_gff_overlap() -> io::Result<()> {
        let qh = gff_query_header();
        let region = parse_region("chr1:15-25")?;
        assert!(qh.line_intersects_region("chr1\tsrc\tgene\t10\t20\t.\t+\t.\tID=gene1", &region)?);
        Ok(())
    }

    #[test]
    fn test_line_intersects_region_gff_no_overlap_before() -> io::Result<()> {
        let qh = gff_query_header();
        let region = parse_region("chr1:15-25")?;
        // feature ends at 14, region starts at 15 — no overlap
        assert!(!qh.line_intersects_region("chr1\tsrc\tgene\t5\t14\t.\t+\t.\tID=gene1", &region)?);
        Ok(())
    }

    #[test]
    fn test_line_intersects_region_gff_no_overlap_after() -> io::Result<()> {
        let qh = gff_query_header();
        let region = parse_region("chr1:15-25")?;
        // feature starts at 26, region ends at 25 — no overlap
        assert!(!qh.line_intersects_region("chr1\tsrc\tgene\t26\t35\t.\t+\t.\tID=gene1", &region)?);
        Ok(())
    }

    #[test]
    fn test_line_intersects_region_gff_exact_boundary() -> io::Result<()> {
        let qh = gff_query_header();
        let region = parse_region("chr1:15-25")?;
        // feature exactly matches region boundary
        assert!(qh.line_intersects_region("chr1\tsrc\tgene\t15\t25\t.\t+\t.\tID=gene1", &region)?);
        Ok(())
    }

    #[test]
    fn test_line_intersects_region_wrong_chrom() -> io::Result<()> {
        let qh = gff_query_header();
        let region = parse_region("chr1:15-25")?;
        assert!(!qh.line_intersects_region("chr2\tsrc\tgene\t10\t20\t.\t+\t.\tID=gene1", &region)?);
        Ok(())
    }

    #[test]
    fn test_line_intersects_region_comment_line() -> io::Result<()> {
        let qh = gff_query_header();
        let region = parse_region("chr1:15-25")?;
        assert!(!qh.line_intersects_region("#comment line", &region)?);
        Ok(())
    }

    #[test]
    fn test_line_intersects_region_missing_start_field() -> io::Result<()> {
        let qh = gff_query_header();
        let region = parse_region("chr1:15-25")?;
        // only 3 fields — start_position_index (3) is missing
        let err = qh
            .line_intersects_region("chr1\tsrc\tgene", &region)
            .unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidData);
        assert!(err.to_string().contains("start"), "unexpected error: {err}");
        Ok(())
    }

    #[test]
    fn test_line_intersects_region_missing_end_field() -> io::Result<()> {
        let qh = gff_query_header();
        let region = parse_region("chr1:15-25")?;
        // only 4 fields — end_position_index (4) is missing
        let err = qh
            .line_intersects_region("chr1\tsrc\tgene\t10", &region)
            .unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidData);
        assert!(err.to_string().contains("end"), "unexpected error: {err}");
        Ok(())
    }

    // -----------------------------------------------------------------------
    // QueryHeader::line_intersects_region — BED (0-based)
    // -----------------------------------------------------------------------

    #[test]
    fn test_line_intersects_region_bed_overlap() -> io::Result<()> {
        let qh = bed_query_header();
        // BED: chr1  9  20  (0-based half-open = 1-based 10..=20)
        let region = parse_region("chr1:15-25")?;
        assert!(qh.line_intersects_region("chr1\t9\t20\tfeature", &region)?);
        Ok(())
    }

    #[test]
    fn test_line_intersects_region_bed_no_overlap_before() -> io::Result<()> {
        let qh = bed_query_header();
        // BED: chr1  4  14  (0-based = 1-based 5..=14), region 15-25
        let region = parse_region("chr1:15-25")?;
        assert!(!qh.line_intersects_region("chr1\t4\t14\tfeature", &region)?);
        Ok(())
    }

    #[test]
    fn test_line_intersects_region_bed_no_overlap_after() -> io::Result<()> {
        let qh = bed_query_header();
        // BED: chr1  25  35  (0-based = 1-based 26..=35), region 15-25
        let region = parse_region("chr1:15-25")?;
        assert!(!qh.line_intersects_region("chr1\t25\t35\tfeature", &region)?);
        Ok(())
    }

    #[test]
    fn test_line_intersects_region_bed_exact_boundary() -> io::Result<()> {
        let qh = bed_query_header();
        // BED: chr1  24  35  (0-based = 1-based 25..=35), region 15-25
        let region = parse_region("chr1:15-25")?;
        assert!(qh.line_intersects_region("chr1\t24\t35\tfeature", &region)?);
        Ok(())
    }

    #[test]
    fn test_line_intersects_region_bed_wrong_chrom() -> io::Result<()> {
        let qh = bed_query_header();
        let region = parse_region("chr1:15-25")?;
        assert!(!qh.line_intersects_region("chr2\t9\t20\tfeature", &region)?);
        Ok(())
    }

    #[test]
    fn test_line_intersects_region_bed_comment_line() -> io::Result<()> {
        let qh = bed_query_header();
        let region = parse_region("chr1:15-25")?;
        assert!(!qh.line_intersects_region("#comment line", &region)?);
        Ok(())
    }

    #[test]
    fn test_line_intersects_region_bed_missing_start_field() -> io::Result<()> {
        let qh = bed_query_header();
        let region = parse_region("chr1:15-25")?;
        let err = qh.line_intersects_region("chr1", &region).unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidData);
        assert!(err.to_string().contains("start"), "unexpected error: {err}");
        Ok(())
    }

    #[test]
    fn test_line_intersects_region_bed_missing_end_field() -> io::Result<()> {
        let qh = bed_query_header();
        let region = parse_region("chr1:15-25")?;
        let err = qh.line_intersects_region("chr1\t9", &region).unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidData);
        assert!(err.to_string().contains("end"), "unexpected error: {err}");
        Ok(())
    }

    // -----------------------------------------------------------------------
    // QueryHeader::line_intersects_region — missing header configuration
    // -----------------------------------------------------------------------

    #[test]
    fn test_line_intersects_region_no_ref_seq_index() -> io::Result<()> {
        // invalid reference sequence name index (10) — line does not has 11 fields
        let header = Header::builder()
            .set_reference_sequence_name_index(10)
            .set_start_position_index(3)
            .set_end_position_index(Some(4))
            .build();
        let qh = QueryHeader::from(&header);
        let region = parse_region("chr1:15-25")?;
        let err = qh
            .line_intersects_region("chr1\tsrc\tgene\t15\t25\t.\t+\t.\tID=gene1", &region)
            .unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidData);
        assert!(
            err.to_string().contains("reference sequence"),
            "unexpected error: {err}"
        );
        Ok(())
    }

    #[test]
    fn test_line_intersects_region_no_end_position_index() -> io::Result<()> {
        // VCF does not have end position field
        let header = Header::builder()
            .set_reference_sequence_name_index(0)
            .set_start_position_index(2)
            .set_end_position_index(None)
            .build();
        let qh = QueryHeader::from(&header);
        let region = parse_region("chr1:15-25")?;
        assert!(qh.line_intersects_region("chr1\tvar1\t20\tA\tT", &region)?);
        Ok(())
    }

    // -----------------------------------------------------------------------
    // QueryHeader::classify_line — early-stop signal
    // -----------------------------------------------------------------------

    #[test]
    fn test_classify_line_past_end_on_matching_contig() -> io::Result<()> {
        let qh = gff_query_header();
        let region = parse_region("chr1:15-25")?;
        assert_eq!(
            qh.classify_line("chr1\tsrc\tgene\t30\t40\t.\t+\t.\tID=gene1", &region)?,
            LineRegionMatch::PastEnd
        );
        Ok(())
    }

    #[test]
    fn test_classify_line_not_past_end_when_region_unbounded() -> io::Result<()> {
        let qh = gff_query_header();
        // No end means "15 to end of contig" -- a line here can be excluded only for starting
        // *before* the region, never for starting "past" a nonexistent end.
        let region = parse_region("chr1:9000")?;
        assert_eq!(
            qh.classify_line("chr1\tsrc\tgene\t10\t20\t.\t+\t.\tID=gene1", &region)?,
            LineRegionMatch::NoMatch
        );
        Ok(())
    }

    #[test]
    fn test_classify_line_different_contig_is_no_match_not_past_end() -> io::Result<()> {
        let qh = gff_query_header();
        let region = parse_region("chr1:15-25")?;
        assert_eq!(
            qh.classify_line("chr2\tsrc\tgene\t9000\t9100\t.\t+\t.\tID=gene1", &region)?,
            LineRegionMatch::NoMatch
        );
        Ok(())
    }

    #[test]
    fn test_classify_line_before_start_is_no_match() -> io::Result<()> {
        let qh = gff_query_header();
        let region = parse_region("chr1:15-25")?;
        assert_eq!(
            qh.classify_line("chr1\tsrc\tgene\t1\t5\t.\t+\t.\tID=gene1", &region)?,
            LineRegionMatch::NoMatch
        );
        Ok(())
    }

    #[test]
    fn test_classify_line_intersects() -> io::Result<()> {
        let qh = gff_query_header();
        let region = parse_region("chr1:15-25")?;
        assert_eq!(
            qh.classify_line("chr1\tsrc\tgene\t10\t20\t.\t+\t.\tID=gene1", &region)?,
            LineRegionMatch::Intersects
        );
        Ok(())
    }

    // -----------------------------------------------------------------------
    // QueryHeader::classify_line — reference sequence column after start/end
    // -----------------------------------------------------------------------

    #[test]
    fn test_classify_line_finds_late_contig_column() -> io::Result<()> {
        // Regression: the single-pass scan used to stop as soon as it had read start and end,
        // even when the configured reference-name column (here, index 2) had not been reached
        // yet -- wrongly reporting the name as missing instead of reading three more columns.
        let qh = generic_late_contig_header();
        let region = parse_region("chr1:15-35")?;

        assert_eq!(
            qh.classify_line("10\t20\tchr1\tx", &region)?,
            LineRegionMatch::Intersects
        );
        assert_eq!(
            qh.classify_line("30\t40\tchr1\ty", &region)?,
            LineRegionMatch::Intersects
        );
        Ok(())
    }

    #[test]
    fn test_classify_line_late_contig_column_wrong_chrom_is_no_match() -> io::Result<()> {
        let qh = generic_late_contig_header();
        let region = parse_region("chr1:15-35")?;

        assert_eq!(
            qh.classify_line("10\t20\tchr2\tx", &region)?,
            LineRegionMatch::NoMatch
        );
        Ok(())
    }

    #[test]
    fn test_classify_line_late_contig_column_past_end() -> io::Result<()> {
        let qh = generic_late_contig_header();
        let region = parse_region("chr1:1-5")?;

        assert_eq!(
            qh.classify_line("10\t20\tchr1\tx", &region)?,
            LineRegionMatch::PastEnd
        );
        Ok(())
    }

    // -----------------------------------------------------------------------
    // NoodlesBackend::next — early stop on a past-end row
    // -----------------------------------------------------------------------

    /// Number of lines `NoodlesBackend::next` has physically read for the reader's active query.
    ///
    /// The record fixtures every other test shares happen to place their past-end row at (or
    /// very near) the end of the only chunk that query touches, so asserting on the *returned
    /// records* alone cannot tell early stop apart from a full scan that would have produced the
    /// same output anyway -- both read the same handful of lines regardless. This reads the
    /// test-only counter instead, to prove early stop actually reads fewer lines than a full
    /// scan of the query's rows would.
    fn scanned_lines(reader: &TabixReader) -> usize {
        match &reader.backend {
            Backend::Tbi(backend) => backend.scanned_lines,
            Backend::Csi(backend) => backend.scanned_lines,
            Backend::System(_) => 0,
        }
    }

    /// Six genes on chr1: two intersect "chr1:90-135", the third starts past its end, and three
    /// more sit further along the same contig -- read by a full scan but not by early stop.
    fn six_gene_single_contig_lines() -> Vec<String> {
        vec![
            "##gff-version 3".to_string(),
            "chr1\tsource\tgene\t100\t110\t.\t+\t.\tID=gene_a".to_string(),
            "chr1\tsource\tgene\t120\t130\t.\t+\t.\tID=gene_b".to_string(),
            "chr1\tsource\tgene\t140\t150\t.\t+\t.\tID=gene_c".to_string(),
            "chr1\tsource\tgene\t200\t210\t.\t+\t.\tID=gene_d".to_string(),
            "chr1\tsource\tgene\t300\t310\t.\t+\t.\tID=gene_e".to_string(),
            "chr1\tsource\tgene\t400\t410\t.\t+\t.\tID=gene_f".to_string(),
        ]
    }

    #[test]
    fn test_noodles_backend_next_early_stop_reads_fewer_lines_than_a_full_scan() -> io::Result<()> {
        let dir = create_temp_dir("tabix-test-noodles-early-stop-count")?;
        let (gz_path, _) = write_gff3_with_tbi(&dir, Some(six_gene_single_contig_lines()))?;
        let mut reader = TabixReader::open(&gz_path, "tabix-nonexistent-binary-for-test")?;

        assert!(reader.seek("chr1:90-135")?);
        let records: Vec<String> = std::iter::from_fn(|| reader.next().ok()?).collect();

        assert_eq!(records.len(), 2, "expected 2 records, got {records:?}");
        assert!(records[0].contains("ID=gene_a"));
        assert!(records[1].contains("ID=gene_b"));

        // Early stop reads exactly the two matches plus the one past-end row that triggers it --
        // never gene_d, gene_e or gene_f, which a full scan of the query's chunk(s) would still
        // have to read and discard one by one.
        assert_eq!(
            scanned_lines(&reader),
            3,
            "expected early stop to read only gene_a, gene_b and the triggering gene_c"
        );

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_noodles_backend_next_early_stop_ends_whole_query_not_just_chunk() -> io::Result<()> {
        // Same fixture and region as above, phrased as a query-level assertion: once the
        // past-end row is read, the query is fully spent -- a second `next()` call must not
        // resume into a later chunk and yield gene_d.
        let dir = create_temp_dir("tabix-test-noodles-early-stop-ends-query")?;
        let (gz_path, _) = write_gff3_with_tbi(&dir, Some(six_gene_single_contig_lines()))?;
        let mut reader = TabixReader::open(&gz_path, "tabix-nonexistent-binary-for-test")?;

        assert!(reader.seek("chr1:90-135")?);
        let _ = reader.next()?; // gene_a
        let _ = reader.next()?; // gene_b
        assert_eq!(reader.next()?, None);
        assert_eq!(reader.next()?, None, "a spent query must stay spent");

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_noodles_backend_next_stops_early_on_past_end_row_multi_contig_fixture() -> io::Result<()>
    {
        // The shared default fixture (chr1 genes 1-4, chr2 gene5) as a sanity check that the
        // early-stop behaviour above also holds against the fixture most other tests use, and
        // that it does not disturb the separate chr2 contig.
        let dir = create_temp_dir("tabix-test-noodles-early-stop-multi-contig")?;
        let (gz_path, _) = write_gff3_with_tbi(&dir, None)?;
        let mut reader = TabixReader::open(&gz_path, "tabix-nonexistent-binary-for-test")?;

        assert!(reader.seek("chr1:150-250")?);
        let records: Vec<String> = std::iter::from_fn(|| reader.next().ok()?).collect();

        assert_eq!(records.len(), 2, "expected 2 records, got {records:?}");
        assert!(records[0].contains("ID=gene2"));
        assert!(records[1].contains("ID=gene3"));

        assert!(reader.seek("chr2:100-200")?);
        let chr2_records: Vec<String> = std::iter::from_fn(|| reader.next().ok()?).collect();
        assert_eq!(chr2_records.len(), 1);
        assert!(chr2_records[0].contains("ID=gene5"));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_noodles_backend_next_long_feature_spanning_region_end_intersects() -> io::Result<()> {
        // A feature that starts inside the query and extends past its end is `Intersects`, not
        // `PastEnd`: only a feature that *starts* after the region end is past it.
        let dir = create_temp_dir("tabix-test-noodles-long-feature")?;
        let lines = vec![
            "##gff-version 3".to_string(),
            "chr1\tsource\tgene\t50\t1000\t.\t+\t.\tID=long_gene".to_string(),
        ];
        let (gz_path, _) = write_gff3_with_tbi(&dir, Some(lines))?;
        let mut reader = TabixReader::open(&gz_path, "tabix-nonexistent-binary-for-test")?;

        assert!(reader.seek("chr1:100-200")?);
        let records: Vec<String> = std::iter::from_fn(|| reader.next().ok()?).collect();

        assert_eq!(records.len(), 1);
        assert!(records[0].contains("ID=long_gene"));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_noodles_backend_next_open_ended_region_reads_every_downstream_match() -> io::Result<()>
    {
        // An unbounded region (no end) never classifies a row as `PastEnd` (there is no "end" to
        // be past), so every matching row downstream, however far along the contig, must still
        // be returned.
        let dir = create_temp_dir("tabix-test-noodles-open-ended")?;
        let (gz_path, _) = write_gff3_with_tbi(&dir, Some(six_gene_single_contig_lines()))?;
        let mut reader = TabixReader::open(&gz_path, "tabix-nonexistent-binary-for-test")?;

        assert!(reader.seek("chr1:90")?);
        let records: Vec<String> = std::iter::from_fn(|| reader.next().ok()?).collect();

        let ids: Vec<&str> = records
            .iter()
            .map(|r| r.rsplit("ID=").next().unwrap_or(""))
            .collect();
        assert_eq!(
            ids,
            vec!["gene_a", "gene_b", "gene_c", "gene_d", "gene_e", "gene_f"]
        );

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // NoodlesBackend::next — returned lines are independently owned
    // -----------------------------------------------------------------------

    /// A matched row is moved out of `query.line_buf` (not cloned): pin that each returned
    /// `String` is a fully independent value, and unaffected by whatever the reused buffer holds
    /// afterwards.
    #[test]
    fn test_noodles_backend_next_returned_lines_are_independently_owned() -> io::Result<()> {
        let dir = create_temp_dir("tabix-test-noodles-owned-lines")?;
        let (gz_path, _) = write_gff3_with_tbi(&dir, None)?;
        let mut reader = TabixReader::open(&gz_path, "tabix-nonexistent-binary-for-test")?;

        assert!(reader.seek("chr1")?);
        let first = reader.next()?.expect("expected first record");
        let second = reader.next()?.expect("expected second record");

        assert!(
            first.contains("ID=gene1"),
            "expected ID=gene1, got: {first}"
        );
        assert!(
            second.contains("ID=gene2"),
            "expected ID=gene2, got: {second}"
        );
        assert_ne!(first, second);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // NoodlesBackend::reference_sequence_id
    // -----------------------------------------------------------------------

    #[test]
    fn test_reference_sequence_id() -> io::Result<()> {
        let dir = create_temp_dir("tabix-test-reference-sequence-id")?;
        let (gz_file, tbi_file) = write_gff3_with_tbi(&dir, None)?;
        let index = tabix::fs::read(&tbi_file)?;
        let backend = NoodlesBackend::new(csi::io::IndexedReader::new(File::open(gz_file)?, index));

        assert_eq!(
            backend.reference_sequence_id(&parse_region("chr1")?),
            Some(0)
        );
        assert_eq!(
            backend.reference_sequence_id(&parse_region("chr2")?),
            Some(1)
        );
        assert_eq!(backend.reference_sequence_id(&parse_region("chr3")?), None);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // NoodlesBackend::reference_names
    // -----------------------------------------------------------------------

    #[test]
    fn test_noodles_backend_reference_names() -> io::Result<()> {
        let dir = create_temp_dir("tabix-test-reference-names")?;
        let (gz_file, tbi_file) = write_gff3_with_tbi(&dir, None)?;
        let index = tabix::fs::read(&tbi_file)?;
        let backend = NoodlesBackend::new(csi::io::IndexedReader::new(File::open(gz_file)?, index));

        assert_eq!(
            backend.reference_names(),
            vec!["chr1".to_string(), "chr2".to_string()]
        );

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_noodles_backend_reference_names_no_header() -> io::Result<()> {
        let dir = create_temp_dir("tabix-test-reference-names-no-header")?;
        let path = write_bgzf_file(
            &dir,
            vec!["##gff-version 3".to_string()],
            Some("test.gff3.gz"),
        )?;
        let index = csi::Index::default();
        let backend = NoodlesBackend::new(csi::io::IndexedReader::new(File::open(path)?, index));

        assert!(backend.reference_names().is_empty());

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // TabixReader::reference_names
    // -----------------------------------------------------------------------

    #[test]
    fn test_tabix_reader_reference_names_with_index() -> io::Result<()> {
        let dir = create_temp_dir("tabix-test-reader-reference-names")?;
        let (gz_path, _tbi_path) = write_gff3_with_tbi(&dir, None)?;
        let mut reader = TabixReader::open(&gz_path, "tabix-nonexistent-binary-for-test")?;

        assert_eq!(
            reader.reference_names()?,
            vec!["chr1".to_string(), "chr2".to_string()]
        );

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    #[cfg(unix)]
    fn test_tabix_reader_reference_names_system_backend() -> io::Result<()> {
        use std::os::unix::fs::PermissionsExt;

        let dir = create_temp_dir("tabix-test-reader-reference-names-system")?;

        let tabix_path = dir.join("tabix");
        fs::write(
            &tabix_path,
            "#!/bin/sh\n\
            if [ \"$1\" = '--version' ]; then echo 'tabix 1.0'; exit 0; fi\n\
            if [ \"$1\" = '-l' ]; then printf '%s\\n' chr1 chr2; exit 0; fi\n\
            exit 0\n",
        )?;
        let mut perms = fs::metadata(&tabix_path)?.permissions();
        perms.set_mode(0o755);
        fs::set_permissions(&tabix_path, perms)?;

        let data_path = dir.join("sample.gff3.gz");
        fs::write(&data_path, [])?;

        let mut reader = TabixReader::open(&data_path, tabix_path.to_str().unwrap())?;
        assert!(matches!(reader.backend, Backend::System(_)));
        assert_eq!(
            reader.reference_names()?,
            vec!["chr1".to_string(), "chr2".to_string()]
        );

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // NoodlesBackend::open — CSI branch
    // -----------------------------------------------------------------------

    #[test]
    fn test_backend_opens_csi_when_present() -> io::Result<()> {
        let dir = create_temp_dir("tabix-test-csi-open")?;

        let path = write_bgzf_file(
            &dir,
            vec!["##gff-version 3".to_string()],
            Some("test.gff3.gz"),
        )?;

        // Write a minimal CSI index
        let index = csi::Index::default();
        csi::fs::write(associated_index_path(&path, "csi"), &index)?;

        // Should open via CSI backend, not TBI
        let reader = TabixReader::open(&path, "tabix-nonexistent-binary-for-test")?;
        assert!(matches!(reader.backend, Backend::Csi(_)));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_backend_prefers_index_over_available_executable() -> io::Result<()> {
        let dir = create_temp_dir("tabix-test-index-over-executable")?;

        let path = write_bgzf_file(
            &dir,
            vec!["##gff-version 3".to_string()],
            Some("test.gff3.gz"),
        )?;
        csi::fs::write(associated_index_path(&path, "csi"), &csi::Index::default())?;

        // `sh` exists everywhere tests run, so the system backend would be available.
        let reader = TabixReader::open(&path, "sh")?;
        assert!(matches!(reader.backend, Backend::Csi(_)));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_backend_csi_index_missing() -> io::Result<()> {
        let dir = create_temp_dir("tabix-test-csi-index-missing")?;

        let path = write_bgzf_file(
            &dir,
            vec![
                "##gff-version 3".to_string(),
                "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=gene1".to_string(),
            ],
            Some("test.gff3.gz"),
        )?;

        let index = csi::Index::default();
        let file = File::open(&path)?;
        let mut backend = Backend::Csi(NoodlesBackend::new(csi::io::IndexedReader::new(
            file, index,
        )));

        assert_eq!(backend.next()?, None);
        assert_eq!(backend.header_lines()?, ["##gff-version 3".to_string()]);

        let err = backend.seek("chr1:1-100").unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidInput);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // NoodlesBackend::header_lines — reads from real BGZF data
    // -----------------------------------------------------------------------

    #[test]
    fn test_noodles_backend_header_lines() -> io::Result<()> {
        let dir = create_temp_dir("tabix-test-noodles-headers")?;

        let path = write_bgzf_file(
            &dir,
            vec![
                "##gff-version 3".to_string(),
                "##sequence-region chr1 1 1000".to_string(),
                "# a regular comment".to_string(),
                "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=gene1".to_string(),
            ],
            Some("test.gff3.gz"),
        )?;

        // Write a minimal CSI index
        let index = csi::Index::default();
        let file = File::open(&path)?;
        let mut backend = NoodlesBackend::new(csi::io::IndexedReader::new(file, index));

        // Override cached_headers to None (already None) and call directly.
        let headers = backend.header_lines()?;

        assert_eq!(headers.len(), 3);
        assert_eq!(headers[0], "##gff-version 3");
        assert_eq!(headers[1], "##sequence-region chr1 1 1000");
        assert_eq!(headers[2], "# a regular comment");

        // Second call should return cached result
        assert!(backend.cached_headers.is_some());
        let headers2 = backend.header_lines()?;
        assert_eq!(headers2.len(), 3);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_noodles_backend_header_lines_no_headers() -> io::Result<()> {
        let dir = create_temp_dir("tabix-test-noodles-no-headers")?;

        // File with no comment lines at all
        let path = write_bgzf_file(
            &dir,
            vec!["chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=gene1".to_string()],
            Some("test.gff3.gz"),
        )?;

        let index = csi::Index::default();
        let file = File::open(&path)?;
        let mut backend = NoodlesBackend::new(csi::io::IndexedReader::new(file, index));

        let headers = backend.header_lines()?;
        assert!(headers.is_empty());

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_noodles_backend_header_lines_only_headers() -> io::Result<()> {
        let dir = create_temp_dir("tabix-test-noodles-only-headers")?;

        // File that is entirely header/comment lines
        let path = write_bgzf_file(
            &dir,
            vec![
                "##gff-version 3".to_string(),
                "##sequence-region chr1 1 1000".to_string(),
            ],
            Some("test.gff3.gz"),
        )?;

        let index = csi::Index::default();
        let file = File::open(&path)?;
        let mut backend = NoodlesBackend::new(csi::io::IndexedReader::new(file, index));

        let headers = backend.header_lines()?;
        assert_eq!(headers.len(), 2);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    /// Verify that calling header_lines() while a query is active does not
    /// disturb the query — the saved/restored virtual position means the
    /// next call to next() still returns the expected record.
    #[test]
    fn test_noodles_backend_header_lines_preserves_bgzf_position() -> io::Result<()> {
        let dir = create_temp_dir("tabix-test-noodles-pos-preserve")?;

        let path = write_bgzf_file(
            &dir,
            vec![
                "##gff-version 3".to_string(),
                "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=gene1".to_string(),
            ],
            Some("test.gff3.gz"),
        )?;

        let index = csi::Index::default();
        let file = File::open(&path)?;
        let mut backend = NoodlesBackend::new(csi::io::IndexedReader::new(file, index));

        let pos_before = backend.reader.get_mut().virtual_position();
        backend.header_lines()?;
        let pos_after = backend.reader.get_mut().virtual_position();

        assert_eq!(pos_before, pos_after);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // NoodlesBackend::header_lines — with real tbi file
    // -----------------------------------------------------------------------

    #[test]
    fn test_noodles_backend_header_lines_with_tbi() -> io::Result<()> {
        let dir = create_temp_dir("tabix-test-noodles-headers-tbi")?;
        let (gz_path, _) = write_gff3_with_tbi(&dir, None)?;
        let mut reader = TabixReader::open(&gz_path, "tabix-nonexistent-binary-for-test")?;

        let headers = reader.header_lines()?;
        assert_eq!(headers.len(), 1);
        assert_eq!(headers[0], "##gff-version 3");

        // second call should return cached result
        let headers2 = reader.header_lines()?;
        assert_eq!(headers2.len(), 1);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // NoodlesBackend::seek — reference sequence absent from the index
    // -----------------------------------------------------------------------

    #[test]
    fn test_noodles_backend_seek_unknown_region_returns_no_rows() -> io::Result<()> {
        let dir = create_temp_dir("tabix-test-noodles-unknown-region")?;
        let (gz_path, _) = write_gff3_with_tbi(&dir, None)?;
        let mut reader = TabixReader::open(&gz_path, "tabix-nonexistent-binary-for-test")?;

        // "chrZ" is not in the index — no rows, not an error.
        assert!(!reader.seek("chrZ:1-100")?);
        assert_eq!(reader.next()?, None);

        // The reader stays usable for a subsequent, known contig.
        assert!(reader.seek("chr1:1-50")?);
        let rec = reader.next()?.expect("expected a record");
        assert!(rec.contains("ID=gene1"), "expected ID=gene1, got: {rec}");

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // NoodlesBackend::seek — invalid region string
    // -----------------------------------------------------------------------

    #[test]
    fn test_noodles_backend_seek_invalid_region_string() -> io::Result<()> {
        let dir = create_temp_dir("tabix-test-noodles-invalid-region")?;
        let (gz_path, _) = write_gff3_with_tbi(&dir, None)?;
        let mut reader = TabixReader::open(&gz_path, "tabix-nonexistent-binary-for-test")?;

        // "chr1:10-" fails region *syntax* parsing (trailing hyphen, no numeric end) -- distinct
        // from a syntactically valid but unindexed name, which now returns `false` instead (see
        // `test_noodles_backend_seek_unknown_region_returns_no_rows`).
        let err = reader.seek("chr1:10-").unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidInput);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // NoodlesBackend::next
    // -----------------------------------------------------------------------

    #[test]
    fn test_noodles_backend_seek_subregion() -> io::Result<()> {
        let dir = create_temp_dir("tabix-test-noodles-subregion")?;
        let (gz_path, _) = write_gff3_with_tbi(&dir, None)?;
        let mut reader = TabixReader::open(&gz_path, "tabix-nonexistent-binary-for-test")?;

        // overlaps records 2 and 3 only
        assert!(reader.seek("chr1:150-250")?);
        let rec1 = reader.next()?.expect("expected first record");
        assert!(rec1.contains("ID=gene2"), "expected ID=gene2, got: {rec1}");
        let rec2 = reader.next()?.expect("expected second record");
        assert!(rec2.contains("ID=gene3"), "expected ID=gene3, got: {rec2}");
        assert_eq!(reader.next()?, None);

        // re-seek: record 1 only
        assert!(reader.seek("chr1:1-50")?);
        let rec = reader.next()?.expect("expected record");
        assert!(rec.contains("ID=gene1"), "expected ID=gene1, got: {rec}");
        assert_eq!(reader.next()?, None);

        // different chromosome
        assert!(reader.seek("chr2:100-200")?);
        let rec = reader.next()?.expect("expected chr2 record");
        assert!(rec.contains("ID=gene5"), "expected ID=gene5, got: {rec}");
        assert_eq!(reader.next()?, None);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_noodles_backend_seek_full_chromosome() -> io::Result<()> {
        let dir = create_temp_dir("tabix-test-noodles-full-chrom")?;
        let (gz_path, _) = write_gff3_with_tbi(&dir, None)?;
        let mut reader = TabixReader::open(&gz_path, "tabix-nonexistent-binary-for-test")?;

        assert!(reader.seek("chr1")?);
        let records: Vec<String> = std::iter::from_fn(|| reader.next().ok()?).collect();

        assert_eq!(
            records.len(),
            4,
            "expected 4 chr1 records, got {}",
            records.len()
        );
        assert!(records[0].contains("ID=gene1"));
        assert!(records[1].contains("ID=gene2"));
        assert!(records[2].contains("ID=gene3"));
        assert!(records[3].contains("ID=gene4"));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_noodles_backend_seek_empty_region() -> io::Result<()> {
        let dir = create_temp_dir("tabix-test-noodles-empty-region")?;
        let (gz_path, _) = write_gff3_with_tbi(&dir, None)?;
        let mut reader = TabixReader::open(&gz_path, "tabix-nonexistent-binary-for-test")?;

        assert!(!reader.seek("chr1:600-700")?);
        assert_eq!(reader.next()?, None);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_noodles_backend_seek_without_end() -> io::Result<()> {
        let dir = create_temp_dir("tabix-test-noodles-seek-no-end")?;
        let (gz_path, _) = write_gff3_with_tbi(&dir, None)?;
        let mut reader = TabixReader::open(&gz_path, "tabix-nonexistent-binary-for-test")?;

        // seek with no end should return all records from start to end of chrom
        assert!(reader.seek("chr1:150")?);
        let records: Vec<String> = std::iter::from_fn(|| reader.next().ok()?).collect();

        assert_eq!(
            records.len(),
            3,
            "expected 3 chr1 records, got {}",
            records.len()
        );
        assert!(records[0].contains("ID=gene2"));
        assert!(records[1].contains("ID=gene3"));
        assert!(records[2].contains("ID=gene4"));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_noodles_backend_seek_start_greater_than_end() -> io::Result<()> {
        let dir = create_temp_dir("tabix-test-noodles-seek-start-gt-end")?;
        let (gz_path, _) = write_gff3_with_tbi(&dir, None)?;
        let mut reader = TabixReader::open(&gz_path, "tabix-nonexistent-binary-for-test")?;

        assert!(!reader.seek("chr1:400-200")?);
        assert_eq!(reader.next()?, None);

        assert!(!reader.seek("chr1:300-200")?);
        assert_eq!(reader.next()?, None);

        assert!(!reader.seek("chr1:300-100")?);
        assert_eq!(reader.next()?, None);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    /// Regression: an inverted-range seek used to return early without clearing a previous
    /// query, so leftover rows from a prior seek leaked through the next `next()` call.
    ///
    /// The first seek matches four rows; priming consumes only the first (`gene1`), leaving
    /// `gene2`..`gene4` in the query. Without the fix, the second (inverted-range) seek would
    /// leave that leftover query active, and the following `next()` would wrongly return
    /// `gene2` instead of `None`.
    #[test]
    fn test_noodles_backend_seek_start_greater_than_end_clears_previous_query() -> io::Result<()> {
        let dir = create_temp_dir("tabix-test-noodles-seek-clears-stale")?;
        let (gz_path, _) = write_gff3_with_tbi(&dir, None)?;
        let mut reader = TabixReader::open(&gz_path, "tabix-nonexistent-binary-for-test")?;

        assert!(reader.seek("chr1")?);
        // Re-seek with an inverted range before consuming the other three primed-query rows.
        assert!(!reader.seek("chr1:200-100")?);
        assert_eq!(
            reader.next()?,
            None,
            "stale rows from the earlier seek must not leak through"
        );

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_noodles_backend_next_clears_query_on_eof() -> io::Result<()> {
        let dir = create_temp_dir("tabix-test-noodles-eof")?;
        let path = write_bgzf_file(&dir, vec![], Some("empty.gff3.gz"))?;
        let file = File::open(&path)?;

        let mut backend =
            NoodlesBackend::new(csi::io::IndexedReader::new(file, csi::Index::default()));
        backend.query = Some(NoodlesQuery {
            region: parse_region("chr1:1-100")?,
            header: gff_query_header(),
            chunks: vec![Chunk::new(
                bgzf::VirtualPosition::default(),
                bgzf::VirtualPosition::from(1),
            )]
            .into_iter(),
            current_chunk_end: None,
            line_buf: String::new(),
        });

        assert_eq!(backend.next()?, None);
        assert!(backend.query.is_none());

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    #[cfg(unix)]
    fn test_noodles_backend_next_without_query_returns_none() -> io::Result<()> {
        let dir = create_temp_dir("tabix-test-noodles-no-query")?;
        let path = write_bgzf_file(&dir, vec![], Some("empty.gff3.gz"))?;
        let file = File::open(&path)?;
        let mut backend =
            NoodlesBackend::new(csi::io::IndexedReader::new(file, csi::Index::default()));

        assert_eq!(backend.next().unwrap(), None);
        Ok(())
    }

    // -----------------------------------------------------------------------
    // NoodlesBackend::next — edge cases
    // -----------------------------------------------------------------------

    #[test]
    fn test_noodles_backend_zero_based_positions() -> io::Result<()> {
        let dir = create_temp_dir("tabix-test-noodles-zero-based")?;
        let lines = vec![
            "chr1\t0\t0\tgene1".to_string(),
            "chr1\t99\t99\tgene2".to_string(),
            "chr1\t199\t299\tgene3".to_string(),
        ];
        let (gz_path, _) = write_bed_with_tbi(&dir, Some(lines))?;
        let mut reader = TabixReader::open(&gz_path, "tabix-nonexistent-binary-for-test")?;

        assert!(reader.seek("chr1:1-500")?);
        let records: Vec<String> = std::iter::from_fn(|| reader.next().ok()?).collect();

        assert_eq!(
            records.len(),
            3,
            "expected 3 chr1 records, got {}",
            records.len()
        );
        assert!(records[0].contains("gene1")); // Noodles backend return 0=..0 based records unlike system tabix
        assert!(records[1].contains("gene2"));
        assert!(records[2].contains("gene3"));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // SystemBackend::is_available — non-NotFound error propagates
    // -----------------------------------------------------------------------

    #[test]
    #[cfg(unix)]
    fn test_system_backend_is_available_not_found() {
        let result = SystemBackend::is_available("tabix-nonexistent-binary-for-test");
        assert!(result.is_ok());
        assert!(!result.unwrap());
    }

    #[test]
    #[cfg(unix)]
    fn test_system_backend_is_available_permission_error() -> io::Result<()> {
        use std::os::unix::fs::PermissionsExt;

        let dir = create_temp_dir("tabix-test-sys-perm")?;

        // Create a non-executable file — spawning it gives PermissionDenied,
        // not NotFound, so is_available should propagate the error rather than
        // returning Ok(false).
        let tabix_path = dir.join("tabix");
        fs::write(&tabix_path, b"#!/bin/sh\n")?;
        let mut perms = fs::metadata(&tabix_path)?.permissions();
        perms.set_mode(0o644); // not executable
        fs::set_permissions(&tabix_path, perms)?;

        let result = SystemBackend::is_available(tabix_path.to_str().unwrap());
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert_eq!(err.kind(), ErrorKind::PermissionDenied);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // SystemBackend::header_lines
    // -----------------------------------------------------------------------

    #[test]
    #[cfg(unix)]
    fn test_system_backend_header_lines() -> io::Result<()> {
        use std::os::unix::fs::PermissionsExt;

        let dir = create_temp_dir("tabix-test-sys-headers")?;

        let tabix_path = dir.join("tabix");
        fs::write(
            &tabix_path,
            "#!/bin/sh\n\
            if [ \"$1\" = '-H' ]; then printf '%s\\n' '##gff-version 3'; exit 0; fi\n\
            exit 0\n",
        )?;
        let mut perms = fs::metadata(&tabix_path)?.permissions();
        perms.set_mode(0o755);
        fs::set_permissions(&tabix_path, perms)?;

        let data_path = dir.join("sample.gff3.gz");
        fs::write(&data_path, [])?;

        let mut backend = SystemBackend::new(&data_path, tabix_path.to_str().unwrap());

        let headers = backend.header_lines()?;
        assert_eq!(headers.len(), 1);
        assert_eq!(headers[0], "##gff-version 3");

        // second call should return cached result
        assert!(backend.cached_headers.is_some());
        let headers2 = backend.header_lines()?;
        assert_eq!(headers2.len(), 1);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    #[cfg(unix)]
    fn test_system_backend_header_lines_empty() -> io::Result<()> {
        use std::os::unix::fs::PermissionsExt;

        let dir = create_temp_dir("tabix-test-sys-empty-headers")?;

        let tabix_path = dir.join("tabix");
        fs::write(
            &tabix_path,
            "#!/bin/sh\n\
            if [ \"$1\" = '-H' ]; then exit 0; fi\n\
            exit 0\n",
        )?;
        let mut perms = fs::metadata(&tabix_path)?.permissions();
        perms.set_mode(0o755);
        fs::set_permissions(&tabix_path, perms)?;

        let data_path = dir.join("sample.gff3.gz");
        fs::write(&data_path, [])?;

        let mut backend = SystemBackend::new(&data_path, tabix_path.to_str().unwrap());

        let headers = backend.header_lines()?;
        assert!(headers.is_empty());

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    #[cfg(unix)]
    fn test_system_backend_header_lines_command_fails() -> io::Result<()> {
        use std::os::unix::fs::PermissionsExt;

        let dir = create_temp_dir("tabix-test-sys-header-fail")?;

        let tabix_path = dir.join("tabix");
        fs::write(
            &tabix_path,
            "#!/bin/sh\n\
            if [ \"$1\" = '--version' ]; then echo 'tabix 1.0'; exit 0; fi\n\
            if [ \"$1\" = '-H' ]; then exit 2; fi\n",
        )?;
        let mut perms = fs::metadata(&tabix_path)?.permissions();
        perms.set_mode(0o755);
        fs::set_permissions(&tabix_path, perms)?;

        let data_path = dir.join("sample.gff3.gz");
        fs::write(&data_path, [])?;

        let mut backend = SystemBackend::new(&data_path, tabix_path.to_str().unwrap());
        let err = backend.header_lines().unwrap_err();
        assert!(
            err.to_string().contains("tabix -H failed"),
            "unexpected error: {err}"
        );

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // SystemBackend::reference_names
    // -----------------------------------------------------------------------

    #[test]
    #[cfg(unix)]
    fn test_system_backend_reference_names_are_cached() -> io::Result<()> {
        use std::os::unix::fs::PermissionsExt;

        let dir = create_temp_dir("tabix-test-sys-reference-names")?;

        // Each `-l` call appends to a log, so a second listing would show up as a second line.
        let calls_path = dir.join("calls");
        let tabix_path = dir.join("tabix");
        fs::write(
            &tabix_path,
            format!(
                "#!/bin/sh\n\
                if [ \"$1\" = '-l' ]; then echo l >> '{}'; printf '%s\\n' chr1 chrX; exit 0; fi\n\
                exit 0\n",
                calls_path.display()
            ),
        )?;
        let mut perms = fs::metadata(&tabix_path)?.permissions();
        perms.set_mode(0o755);
        fs::set_permissions(&tabix_path, perms)?;

        let data_path = dir.join("sample.gff3.gz");
        fs::write(&data_path, [])?;

        let mut backend = SystemBackend::new(&data_path, tabix_path.to_str().unwrap());
        assert_eq!(backend.reference_names()?, ["chr1", "chrX"]);
        assert_eq!(backend.reference_names()?, ["chr1", "chrX"]);
        assert_eq!(fs::read_to_string(&calls_path)?.lines().count(), 1);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    #[cfg(unix)]
    fn test_system_backend_reference_names_empty() -> io::Result<()> {
        use std::os::unix::fs::PermissionsExt;

        let dir = create_temp_dir("tabix-test-sys-reference-names-empty")?;

        let tabix_path = dir.join("tabix");
        fs::write(&tabix_path, "#!/bin/sh\nexit 0\n")?;
        let mut perms = fs::metadata(&tabix_path)?.permissions();
        perms.set_mode(0o755);
        fs::set_permissions(&tabix_path, perms)?;

        let data_path = dir.join("sample.gff3.gz");
        fs::write(&data_path, [])?;

        let mut backend = SystemBackend::new(&data_path, tabix_path.to_str().unwrap());
        assert!(backend.reference_names()?.is_empty());

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    #[cfg(unix)]
    fn test_system_backend_reference_names_command_fails() -> io::Result<()> {
        use std::os::unix::fs::PermissionsExt;

        let dir = create_temp_dir("tabix-test-sys-reference-names-fail")?;

        let tabix_path = dir.join("tabix");
        fs::write(
            &tabix_path,
            "#!/bin/sh\n\
            if [ \"$1\" = '-l' ]; then exit 2; fi\n\
            exit 0\n",
        )?;
        let mut perms = fs::metadata(&tabix_path)?.permissions();
        perms.set_mode(0o755);
        fs::set_permissions(&tabix_path, perms)?;

        let data_path = dir.join("sample.gff3.gz");
        fs::write(&data_path, [])?;

        let mut backend = SystemBackend::new(&data_path, tabix_path.to_str().unwrap());
        let err = backend.reference_names().unwrap_err();
        assert!(
            err.to_string().contains("tabix -l failed"),
            "unexpected error: {err}"
        );
        assert!(backend.cached_reference_names.is_none());

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    #[ignore = "skipped: requires tabix"]
    fn test_system_backend_reference_names_match_noodles_backend() -> io::Result<()> {
        let tabix_path = std::env::var("ENSEMBL_VEP_HTSLIB_BIN")
            .map(|bin| Path::new(&bin).join("tabix").to_string_lossy().into_owned())
            .unwrap_or_else(|_| "tabix".to_string());
        let dir = create_temp_dir("tabix-test-sys-reference-names-parity")?;
        let (gz_path, tbi_path) = write_gff3_with_tbi(&dir, None)?;
        let index = tabix::fs::read(&tbi_path)?;
        let noodles_backend =
            NoodlesBackend::new(csi::io::IndexedReader::new(File::open(&gz_path)?, index));
        let mut system_backend = SystemBackend::new(&gz_path, &tabix_path);

        assert_eq!(
            system_backend.reference_names()?,
            noodles_backend.reference_names()
        );

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // SystemBackend::seek — error propagates in next() after seek spawns child
    // -----------------------------------------------------------------------

    #[test]
    fn test_system_backend_next_without_query_returns_none() {
        let mut backend = SystemBackend::new(Path::new("sample.gff3.gz"), "tabix");
        assert_eq!(backend.next().unwrap(), None);
    }

    #[test]
    #[cfg(unix)]
    fn test_system_backend_next_failing_child() -> io::Result<()> {
        use std::os::unix::fs::PermissionsExt;

        let dir = create_temp_dir("tabix-test-sys-fail-child")?;

        // Mock tabix that succeeds for --version but exits with code 1 for queries
        let tabix_path = dir.join("tabix");
        fs::write(
            &tabix_path,
            "#!/bin/sh\n\
            if [ \"$1\" = '--version' ]; then echo 'tabix 1.0'; exit 0; fi\n\
            exit 1\n",
        )?;
        let mut perms = fs::metadata(&tabix_path)?.permissions();
        perms.set_mode(0o755);
        fs::set_permissions(&tabix_path, perms)?;

        let data_path = dir.join("sample.gff3.gz");
        fs::write(&data_path, [])?;

        let mut reader = TabixReader::open(&data_path, tabix_path.to_str().unwrap())?;
        reader.backend.seek("chr1:1-100")?; // seek succeeds (spawns child)

        // next() should propagate the non-zero exit as an error
        let err = reader.backend.next().unwrap_err();
        assert!(
            err.to_string().contains("tabix executable failed"),
            "unexpected error: {err}"
        );

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // SystemBackend: seek + next
    // -----------------------------------------------------------------------

    #[test]
    #[cfg(unix)]
    fn test_query_uses_system_tabix_without_readable_index() -> io::Result<()> {
        use std::os::unix::fs::PermissionsExt;

        let dir = create_temp_dir("tabix-test-sys")?;

        let tabix_path = dir.join("tabix");
        fs::write(
            &tabix_path,
            "#!/bin/sh\nprintf '%s\\n' 'chr1\ttabix\tgene\t8\t13\t.\t+\t.\tID=gene1'\n",
        )?;

        let mut permissions = fs::metadata(&tabix_path)?.permissions();
        permissions.set_mode(0o755);
        fs::set_permissions(&tabix_path, permissions)?;

        let data_path = dir.join("sample.gff3.gz");
        fs::write(&data_path, [])?;

        let mut reader = TabixReader::open(&data_path, tabix_path.to_str().unwrap())?;
        assert!(reader.seek("chr1:8-13")?);
        assert_eq!(
            reader.next()?.as_deref(),
            Some("chr1\ttabix\tgene\t8\t13\t.\t+\t.\tID=gene1")
        );
        assert_eq!(reader.next()?, None);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // SystemBackend::close_query + Drop — kill an active query
    // -----------------------------------------------------------------------

    #[test]
    #[cfg(unix)]
    fn test_system_backend_close_query_kills_child() -> io::Result<()> {
        use std::os::unix::fs::PermissionsExt;

        let dir = create_temp_dir("tabix-test-sys-close")?;

        // Mock tabix that execs a long-running process so the child PID can be
        // checked directly after close_query and Drop.
        let tabix_path = dir.join("tabix");
        fs::write(
            &tabix_path,
            "#!/bin/sh\n\
            if [ \"$1\" = '--version' ]; then echo 'tabix 1.0'; exit 0; fi\n\
            exec sleep 6000\n",
        )?;
        let mut perms = fs::metadata(&tabix_path)?.permissions();
        perms.set_mode(0o755);
        fs::set_permissions(&tabix_path, perms)?;

        let data_path = dir.join("sample.gff3.gz");
        fs::write(&data_path, [])?;

        let mut reader = TabixReader::open(&data_path, tabix_path.to_str().unwrap())?;
        reader.backend.seek("chr1:1-100")?;
        let first_child_pid = get_system_query_child_id(&reader);
        assert!(
            process_exists(first_child_pid)?,
            "expected process {first_child_pid} to exist"
        );

        reader.backend.seek("chr1:1-100")?;
        assert_process_is_gone(first_child_pid)?;

        let second_child_pid = get_system_query_child_id(&reader);
        assert!(
            process_exists(second_child_pid)?,
            "expected process {second_child_pid} to exist"
        );
        drop(reader);
        assert_process_is_gone(second_child_pid)?;

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // TabixReader::header_lines — system backend dispatch
    // -----------------------------------------------------------------------

    #[test]
    #[cfg(unix)]
    fn test_tabix_reader_system_header_lines_dispatch() -> io::Result<()> {
        use std::os::unix::fs::PermissionsExt;

        let dir = create_temp_dir("tabix-test-reader-sys-headers")?;

        let tabix_path = dir.join("tabix");
        fs::write(
            &tabix_path,
            "#!/bin/sh\n\
            if [ \"$1\" = '--version' ]; then echo 'tabix 1.0'; exit 0; fi\n\
            if [ \"$1\" = '-H' ]; then printf '%s\\n' '##gff-version 3'; exit 0; fi\n\
            exit 0\n",
        )?;
        let mut perms = fs::metadata(&tabix_path)?.permissions();
        perms.set_mode(0o755);
        fs::set_permissions(&tabix_path, perms)?;

        let data_path = dir.join("sample.gff3.gz");
        fs::write(&data_path, [])?;

        let mut reader = TabixReader::open(&data_path, tabix_path.to_str().unwrap())?;
        assert_eq!(reader.header_lines()?, ["##gff-version 3".to_string()]);

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // TabixReader::open — missing file / missing index
    // -----------------------------------------------------------------------

    #[test]
    fn test_tabix_reader_open_missing_file() {
        // No system tabix, no index file, and the gz itself doesn't exist —
        // the TBI builder should return NotFound.
        match TabixReader::open(
            Path::new("/nonexistent/path/file.gff3.gz"),
            "tabix-nonexistent-binary-for-test",
        ) {
            Err(err) => assert_eq!(err.kind(), ErrorKind::NotFound),
            Ok(_) => panic!("expected error, got Ok"),
        }
    }

    // -----------------------------------------------------------------------
    // TabixReader::seek — buffered_line cleared on re-seek
    // -----------------------------------------------------------------------

    #[test]
    fn test_tabix_reader_buffered_line_cleared_on_reseek() -> io::Result<()> {
        let dir = create_temp_dir("tabix-test-buf-clear")?;
        let (gz_path, _) = write_gff3_with_tbi(&dir, None)?;
        let mut reader = TabixReader::open(&gz_path, "tabix-nonexistent-binary-for-test")?;

        // First seek — primes buffered_line with ID=gene1
        assert!(reader.seek("chr1:1-50")?);

        // Second seek before consuming — buffered_line must be cleared,
        // so next() returns ID=gene2 not ID=gene1
        assert!(reader.seek("chr1:100-200")?);
        let rec = reader.next()?.expect("expected a record");
        assert!(
            rec.contains("ID=gene2"),
            "expected ID=gene2 after re-seek, got: {rec}"
        );

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // TabixReader::next — returns buffered_line first
    // -----------------------------------------------------------------------

    #[test]
    fn test_tabix_reader_next_returns_buffered_first() -> io::Result<()> {
        let dir = create_temp_dir("tabix-test-buf-next")?;
        let (gz_path, _) = write_gff3_with_tbi(&dir, None)?;
        let mut reader = TabixReader::open(&gz_path, "tabix-nonexistent-binary-for-test")?;

        assert!(reader.seek("chr1:100-200")?);

        // After seek the buffered_line holds the first record.
        // Calling next() twice should return record 2 then record 3.
        let first = reader.next()?.expect("first record from buffered_line");
        let second = reader.next()?.expect("second record from buffered_line");
        assert!(
            first.contains("ID=gene2"),
            "expected ID=gene2, got: {first}"
        );
        assert!(
            second.contains("ID=gene3"),
            "expected ID=gene3, got: {second}"
        );
        assert_eq!(reader.next()?, None);

        fs::remove_dir_all(dir)?;
        Ok(())
    }
}
