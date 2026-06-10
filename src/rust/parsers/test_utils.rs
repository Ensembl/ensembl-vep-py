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

use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_csi::binning_index::index::{
    Header,
    header::{Builder as HeaderBuilder, Format, format::CoordinateSystem},
    reference_sequence::bin::Chunk,
};
use noodles_tabix::{self as tabix, index::Indexer};
use std::{
    fs::File,
    io::{self, BufRead, Write},
    path::{Path, PathBuf},
};

struct IndexedRecord<'a> {
    chrom: &'a str,
    start: Position,
    end: Position,
}

fn parse_position(value: &str, context: &str) -> io::Result<Position> {
    let value = value
        .parse::<usize>()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("{context}: {e}")))?;
    Position::try_from(value).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn parse_gff3_record(line: &str) -> io::Result<Option<IndexedRecord<'_>>> {
    let fields: Vec<&str> = line.splitn(9, '\t').collect();
    if fields.len() < 5 {
        return Ok(None);
    }

    let start = parse_position(fields[3], "bad start")?;
    let end = parse_position(fields[4], "bad end")?;

    Ok(Some(IndexedRecord {
        chrom: fields[0],
        start,
        end,
    }))
}

fn parse_bed_record(line: &str) -> io::Result<Option<IndexedRecord<'_>>> {
    let fields: Vec<&str> = line.splitn(4, '\t').collect();
    if fields.len() < 3 {
        return Ok(None);
    }

    let start = Position::try_from(
        fields[1]
            .parse::<usize>()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("bad start: {e}")))?
            + 1,
    )
    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let end = Position::try_from(
        fields[2]
            .parse::<usize>()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("bad end: {e}")))?
            .max(1),
    )
    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok(Some(IndexedRecord {
        chrom: fields[0],
        start,
        end,
    }))
}

pub(crate) fn write_bgzf_file(
    dir: &Path,
    lines: Vec<String>,
    file_name: Option<&str>,
) -> io::Result<PathBuf> {
    let file_name = file_name.unwrap_or("test.gff3.gz");
    let path = dir.join(file_name);
    let file = File::create(&path)?;
    let mut writer = bgzf::io::Writer::new(file);
    for line in lines {
        writer.write_all(line.as_bytes())?;
        writer.write_all(b"\n")?;
    }
    writer.finish()?;
    Ok(path)
}

fn write_indexed_tabix_file<F>(
    dir: &Path,
    lines: Vec<String>,
    file_name: &str,
    header: Header,
    parse_record: F,
) -> io::Result<(PathBuf, PathBuf)>
where
    F: for<'a> Fn(&'a str) -> io::Result<Option<IndexedRecord<'a>>>,
{
    let gz_path = dir.join(file_name);
    let tbi_path = dir.join(format!("{file_name}.tbi"));

    write_bgzf_file(dir, lines, Some(file_name))?;

    let file = File::open(&gz_path)?;
    let mut reader = bgzf::io::Reader::new(file);

    let mut indexer = Indexer::default();
    indexer.set_header(header);

    loop {
        let chunk_start = reader.virtual_position();

        let mut buf = String::new();
        if reader.read_line(&mut buf)? == 0 {
            break;
        }

        let line = buf.trim_end_matches(['\n', '\r']);

        if line.starts_with('#') || line.is_empty() {
            continue;
        }

        let Some(record) = parse_record(line)? else {
            continue;
        };

        let chunk_end = reader.virtual_position();
        let chunk = Chunk::new(chunk_start, chunk_end);

        indexer.add_record(record.chrom, record.start, record.end, chunk)?;
    }

    let index = indexer.build();
    tabix::fs::write(&tbi_path, &index)?;

    Ok((gz_path, tbi_path))
}

pub(crate) fn write_gff3_with_tbi(
    dir: &Path,
    lines: Option<Vec<String>>,
) -> io::Result<(PathBuf, PathBuf)> {
    let default_lines = vec![
        "##gff-version 3".to_string(),
        "chr1\tsource\tgene\t1\t50\t.\t+\t.\tID=gene1".to_string(),
        "chr1\tsource\tgene\t100\t200\t.\t+\t.\tID=gene2".to_string(),
        "chr1\tsource\tgene\t200\t300\t.\t+\t.\tID=gene3".to_string(),
        "chr1\tsource\tgene\t400\t500\t.\t+\t.\tID=gene4".to_string(),
        "chr2\tsource\tgene\t100\t200\t.\t+\t.\tID=gene5".to_string(),
    ];
    let header = Header::builder()
        .set_format(Format::Generic(CoordinateSystem::Gff))
        .build();

    write_indexed_tabix_file(
        dir,
        lines.unwrap_or(default_lines),
        "test.gff3.gz",
        header,
        parse_gff3_record,
    )
}

pub(crate) fn write_bed_with_tbi(
    dir: &Path,
    lines: Option<Vec<String>>,
) -> io::Result<(PathBuf, PathBuf)> {
    let default_lines = vec![
        "chr1\t0\t50\tgene1".to_string(),
        "chr1\t99\t200\tgene2".to_string(),
        "chr1\t199\t300\tgene3".to_string(),
        "chr1\t399\t500\tgene4".to_string(),
        "chr2\t99\t200\tgene5".to_string(),
    ];
    let header = HeaderBuilder::bed().build();

    write_indexed_tabix_file(
        dir,
        lines.unwrap_or(default_lines),
        "test.bed.gz",
        header,
        parse_bed_record,
    )
}

#[cfg(test)]
mod tests {
    use std::{env, fs};

    use super::*;

    // ---------------------------------------------------------------------------
    // write_gff3_with_tbi
    // ---------------------------------------------------------------------------

    #[test]
    fn test_write_gff3_with_tbi_skips_non_record_lines() -> io::Result<()> {
        let dir = env::temp_dir().join(format!(
            "parser-test-utils-gff3-skip-{}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir)?;

        let lines = vec![
            "##gff-version 3".to_string(),
            "".to_string(),
            "chr1\tsource\tgene".to_string(),
            "chr1\tsource\tgene\t1\t50\t.\t+\t.\tID=gene1".to_string(),
        ];
        let (gz_path, tbi_path) = write_gff3_with_tbi(&dir, Some(lines))?;

        assert!(gz_path.exists());
        assert!(tbi_path.exists());

        // check first line from first chunk
        let mut tbi_reader = tabix::io::Reader::new(File::open(&tbi_path)?);
        let index = tbi_reader.read_index()?;

        // Get first chunk
        let reference_sequence = &index.reference_sequences()[0];
        let (_, bin) = reference_sequence.bins().iter().next().unwrap();
        let chunk = bin.chunks().first().unwrap();

        // Open BGZF file
        let mut reader = bgzf::io::Reader::new(File::open(&gz_path)?);

        // Seek to chunk start
        reader.seek(chunk.start())?;

        // Read first line from that position
        let mut line = String::new();
        reader.read_line(&mut line)?;
        assert!(line.starts_with("chr1"));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_write_gff3_with_tbi_rejects_bad_start() -> io::Result<()> {
        let dir = env::temp_dir().join(format!(
            "parser-test-utils-gff3-bad-start-{}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir)?;

        let lines = vec!["chr1\tsource\tgene\tbad\t50\t.\t+\t.\tID=gene1".to_string()];
        let err = write_gff3_with_tbi(&dir, Some(lines)).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("bad start"));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_write_gff3_with_tbi_rejects_bad_end() -> io::Result<()> {
        let dir = env::temp_dir().join(format!(
            "parser-test-utils-gff3-bad-end-{}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir)?;

        let lines = vec!["chr1\tsource\tgene\t1\tbad\t.\t+\t.\tID=gene1".to_string()];
        let err = write_gff3_with_tbi(&dir, Some(lines)).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("bad end"));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    // ---------------------------------------------------------------------------
    // write_bed_with_tbi
    // ---------------------------------------------------------------------------

    #[test]
    fn test_write_bed_with_tbi_skips_non_record_lines() -> io::Result<()> {
        let dir =
            env::temp_dir().join(format!("parser-test-utils-bed-skip-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir)?;

        let lines = vec![
            "# header".to_string(),
            "".to_string(),
            "chr1\t0".to_string(),
            "chr1\t0\t50\tgene1".to_string(),
        ];
        let (gz_path, tbi_path) = write_bed_with_tbi(&dir, Some(lines))?;

        assert!(gz_path.exists());
        assert!(tbi_path.exists());

        // read first line from first chunk
        let mut tbi_reader = tabix::io::Reader::new(File::open(&tbi_path)?);
        let index = tbi_reader.read_index()?;

        // Get first chunk
        let reference_sequence = &index.reference_sequences()[0];
        let (_, bin) = reference_sequence.bins().iter().next().unwrap();
        let chunk = bin.chunks().first().unwrap();

        // Open BGZF file
        let mut reader = bgzf::io::Reader::new(File::open(&gz_path)?);

        // Seek to chunk start
        reader.seek(chunk.start())?;

        // Read first line from that position
        let mut line = String::new();
        reader.read_line(&mut line)?;
        assert!(line.starts_with("chr1"));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_write_bed_with_tbi_rejects_bad_start() -> io::Result<()> {
        let dir = env::temp_dir().join(format!(
            "parser-test-utils-bed-bad-start-{}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir)?;

        let lines = vec!["chr1\tbad\t50\tgene1".to_string()];
        let err = write_bed_with_tbi(&dir, Some(lines)).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("bad start"));

        fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn test_write_bed_with_tbi_rejects_bad_end() -> io::Result<()> {
        let dir = env::temp_dir().join(format!(
            "parser-test-utils-bed-bad-end-{}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir)?;

        let lines = vec!["chr1\t0\tbad\tgene1".to_string()];
        let err = write_bed_with_tbi(&dir, Some(lines)).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("bad end"));

        fs::remove_dir_all(dir)?;
        Ok(())
    }
}
