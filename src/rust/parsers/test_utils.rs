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
    io::{self, BufRead, Write},
    path::{Path, PathBuf},
};

use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_csi::binning_index::index::{
    Header,
    header::{Builder as HeaderBuilder, Format, format::CoordinateSystem},
    reference_sequence::bin::Chunk,
};
use noodles_tabix::{self as tabix, index::Indexer};

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
    let lines = lines.unwrap_or(default_lines);

    let gz_path = dir.join("test.gff3.gz");
    let tbi_path = dir.join("test.gff3.gz.tbi");

    write_bgzf_file(dir, lines, Some("test.gff3.gz"))?;

    let file = File::open(&gz_path)?;
    let mut reader = bgzf::io::Reader::new(file);

    let header = Header::builder()
        .set_format(Format::Generic(CoordinateSystem::Gff))
        .build();

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

        let fields: Vec<&str> = line.splitn(9, '\t').collect();
        if fields.len() < 5 {
            continue;
        }

        let chunk_end = reader.virtual_position();
        let chrom = fields[0];

        let start: usize = fields[3]
            .parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("bad start: {e}")))?;
        let end: usize = fields[4]
            .parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("bad end: {e}")))?;

        let start =
            Position::try_from(start).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let end =
            Position::try_from(end).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        let chunk = Chunk::new(chunk_start, chunk_end);

        indexer.add_record(chrom, start, end, chunk)?;
    }

    let index = indexer.build();
    tabix::fs::write(&tbi_path, &index)?;

    Ok((gz_path, tbi_path))
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
    let lines = lines.unwrap_or(default_lines);

    let gz_path = dir.join("test.bed.gz");
    let tbi_path = dir.join("test.bed.gz.tbi");

    write_bgzf_file(dir, lines, Some("test.bed.gz"))?;

    let file = File::open(&gz_path)?;
    let mut reader = bgzf::io::Reader::new(file);

    let header = HeaderBuilder::bed().build();

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

        let fields: Vec<&str> = line.splitn(4, '\t').collect();
        if fields.len() < 3 {
            continue;
        }

        let chunk_end = reader.virtual_position();
        let chrom = fields[0];

        let start: usize = fields[1]
            .parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("bad start: {e}")))?;
        let end: usize = fields[2]
            .parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("bad end: {e}")))?;

        let start = Position::try_from(start + 1)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let end = Position::try_from(end.max(1))
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        let chunk = Chunk::new(chunk_start, chunk_end);

        indexer.add_record(chrom, start, end, chunk)?;
    }

    let index = indexer.build();
    tabix::fs::write(&tbi_path, &index)?;

    Ok((gz_path, tbi_path))
}

#[cfg(test)]
mod tests {
    use std::{env, fs};

    use super::*;

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
