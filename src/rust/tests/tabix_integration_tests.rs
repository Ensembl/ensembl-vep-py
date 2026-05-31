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

use ensembl_vep_core::parsers::tabix::TabixParser;
use pyo3::PyErr;
use std::{
    env, fs, io,
    path::{Path, PathBuf},
    process::{Command, Stdio},
};

const HEADER_LINES: [&str; 2] = [
    "##gff-version 3",
    "#seqid\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattributes",
];
const CHR1_GENE_ROW: &str = "chr1\ttabix\tgene\t10\t20\t.\t+\t.\tID=gene:1";
const CHR1_EXON_ROW: &str = "chr1\ttabix\texon\t12\t14\t.\t+\t.\tID=exon:1";
const CHR1_CDS_ROW: &str = "chr1\ttabix\tCDS\t15\t18\t.\t+\t0\tID=cds:1";
const CHR2_GENE_ROW: &str = "chr2\ttabix\tgene\t100\t120\t.\t-\t.\tID=gene:2";

struct HtslibToolchain {
    bgzip: PathBuf,
    tabix: PathBuf,
}

#[test]
#[ignore = "skipped: requires bgzip and tabix"]
fn test_tabix_parser_reads_records_with_noodles_backend() -> io::Result<()> {
    let toolchain = require_integration_toolchain();
    let dir = integration_temp_dir("tabix-parser-noodles")?;
    let path = write_tabix_gff3(&dir, &toolchain)?;
    let mut parser = TabixParser::new(
        Some(path.to_string_lossy().into_owned()),
        Some("tabix-nonexistent-binary-for-test".to_string()),
    )
    .map_err(py_err)?;

    assert_eq!(
        parser.header_lines().map_err(py_err)?,
        HEADER_LINES.map(str::to_owned)
    );
    assert!(parser.seek("chr1:10-20").map_err(py_err)?);
    assert_eq!(
        parser.next().map_err(py_err)?,
        Some(CHR1_GENE_ROW.to_string())
    );
    assert_eq!(
        parser.records(None).map_err(py_err)?,
        vec![CHR1_EXON_ROW.to_string(), CHR1_CDS_ROW.to_string()]
    );
    assert_eq!(parser.next().map_err(py_err)?, None);

    Ok(())
}

#[test]
#[ignore = "skipped: requires bgzip and tabix"]
fn test_tabix_parser_reads_records_with_system_tabix_backend() -> io::Result<()> {
    let toolchain = require_integration_toolchain();
    let dir = integration_temp_dir("tabix-parser-system")?;
    let path = write_tabix_gff3(&dir, &toolchain)?;
    let mut parser = TabixParser::new(
        Some(path.to_string_lossy().into_owned()),
        Some(toolchain.tabix.to_string_lossy().into_owned()),
    )
    .map_err(py_err)?;

    assert_eq!(
        parser.header_lines().map_err(py_err)?,
        HEADER_LINES.map(str::to_owned)
    );
    assert!(parser.seek("chr1:10-20").map_err(py_err)?);
    assert_eq!(
        parser.next().map_err(py_err)?,
        Some(CHR1_GENE_ROW.to_string())
    );
    assert_eq!(
        parser.records(None).map_err(py_err)?,
        vec![CHR1_EXON_ROW.to_string(), CHR1_CDS_ROW.to_string()]
    );
    assert_eq!(parser.next().map_err(py_err)?, None);

    Ok(())
}

fn require_integration_toolchain() -> HtslibToolchain {
    let bgzip =
        find_htslib_tool("bgzip").expect("bgzip executable is required for integration tests");
    let tabix =
        find_htslib_tool("tabix").expect("tabix executable is required for integration tests");

    HtslibToolchain { bgzip, tabix }
}

fn find_htslib_tool(name: &str) -> Option<PathBuf> {
    if let Ok(htslib_bin) = env::var("ENSEMBL_VEP_HTSLIB_BIN") {
        let candidate = PathBuf::from(htslib_bin).join(name);
        if candidate.is_file() {
            return Some(candidate);
        }
    }

    env::var_os("PATH").and_then(|paths| {
        env::split_paths(&paths)
            .map(|path| path.join(name))
            .find(|candidate| candidate.is_file())
    })
}

fn integration_temp_dir(prefix: &str) -> io::Result<PathBuf> {
    let dir = env::temp_dir().join(format!("{prefix}-{}", std::process::id()));
    if dir.exists() {
        fs::remove_dir_all(&dir)?;
    }
    fs::create_dir_all(&dir)?;
    Ok(dir)
}

fn write_tabix_gff3(dir: &Path, toolchain: &HtslibToolchain) -> io::Result<PathBuf> {
    let source_path = dir.join("sample.gff3");
    fs::write(
        &source_path,
        format!(
            "{}\n{}\n{}\n{}\n{}\n{}\n",
            HEADER_LINES[0],
            HEADER_LINES[1],
            CHR1_GENE_ROW,
            CHR1_EXON_ROW,
            CHR1_CDS_ROW,
            CHR2_GENE_ROW
        ),
    )?;

    let gz_path = dir.join("sample.gff3.gz");
    let output = fs::File::create(&gz_path)?;
    let status = Command::new(&toolchain.bgzip)
        .arg("-c")
        .arg(&source_path)
        .stdout(Stdio::from(output))
        .status()?;
    if !status.success() {
        return Err(io::Error::other(format!(
            "bgzip failed with status {status}"
        )));
    }

    let status = Command::new(&toolchain.tabix)
        .arg("-p")
        .arg("gff")
        .arg(&gz_path)
        .status()?;
    if !status.success() {
        return Err(io::Error::other(format!(
            "tabix failed with status {status}"
        )));
    }

    Ok(gz_path)
}

fn py_err(err: PyErr) -> io::Error {
    io::Error::other(err.to_string())
}
