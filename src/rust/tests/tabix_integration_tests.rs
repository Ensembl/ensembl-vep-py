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
use pyo3::{PyErr, Python};
use rstest::rstest;
use std::{
    env, fs, io,
    path::{Path, PathBuf},
    process::{Command, Stdio},
};

const GFF3_LINES: [&str; 6] = [
    "##gff-version 3",
    "#seqid\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattributes",
    "chr1\ttabix\tgene\t10\t20\t.\t+\t.\tID=gene:1",
    "chr1\ttabix\texon\t12\t14\t.\t+\t.\tID=exon:1",
    "chr1\ttabix\tCDS\t15\t18\t.\t+\t0\tID=cds:1",
    "chr2\ttabix\tgene\t100\t120\t.\t-\t.\tID=gene:2",
];
const BED_LINES: [&str; 5] = [
    "#chrom\tstart\tend\tname",
    "chr1\t9\t20\tfeature1",
    "chr1\t11\t14\tfeature2",
    "chr1\t14\t18\tfeature3",
    "chr2\t99\t120\tfeature4",
];
const BED_0_BASE_LINES: [&str; 1] = ["chr1\t0\t0\tfeature1"];
const BED_SAME_POS_LINES: [&str; 1] = ["chr1\t99\t99\tfeature1"];

struct HtslibToolchain {
    bgzip: PathBuf,
    tabix: PathBuf,
}

#[derive(Clone, Copy, Debug)]
enum Dataset {
    Gff3Default,
    BedDefault,
    Bed0Base,
    BedSamePos,
}

#[derive(Clone, Copy, Debug)]
enum Backend {
    Noodles,
    System,
}

#[rstest]
#[case::gff3_noodles(Dataset::Gff3Default, Backend::Noodles)]
#[case::gff3_system(Dataset::Gff3Default, Backend::System)]
#[case::bed_noodles(Dataset::BedDefault, Backend::Noodles)]
#[case::bed_system(Dataset::BedDefault, Backend::System)]
#[ignore = "skipped: requires bgzip and tabix"]
fn test_tabix_parser_header_lines(#[case] dataset: Dataset, #[case] backend: Backend) -> io::Result<()> {
    let toolchain = get_integration_toolchain();
    let dir = create_integration_temp_dir(&format!(
        "tabix-parser-header-lines-{:?}-{:?}",
        dataset, backend
    ))?;

    let parser = open_dataset_parser(&dir, &toolchain, dataset, backend)?;
    assert_header_lines(parser, expected_rows(dataset, 0))?;

    fs::remove_dir_all(dir)?;
    Ok(())
}

#[rstest]
#[case::gff3_noodles(Dataset::Gff3Default, Backend::Noodles)]
#[case::gff3_system(Dataset::Gff3Default, Backend::System)]
#[case::bed_noodles(Dataset::BedDefault, Backend::Noodles)]
#[case::bed_system(Dataset::BedDefault, Backend::System)]
#[ignore = "skipped: requires bgzip and tabix"]
fn test_tabix_parser_full_chromosome_query(
    #[case] dataset: Dataset,
    #[case] backend: Backend,
) -> io::Result<()> {
    let toolchain = get_integration_toolchain();
    let dir = create_integration_temp_dir(&format!(
        "tabix-parser-full-chromosome-{:?}-{:?}",
        dataset, backend
    ))?;

    let parser = open_dataset_parser(&dir, &toolchain, dataset, backend)?;
    assert_query_records_multistep(parser, "chr1", expected_rows(dataset, 1))?;

    fs::remove_dir_all(dir)?;
    Ok(())
}

#[rstest]
#[case::gff3_noodles(Dataset::Gff3Default, Backend::Noodles)]
#[case::gff3_system(Dataset::Gff3Default, Backend::System)]
#[case::bed_noodles(Dataset::BedDefault, Backend::Noodles)]
#[case::bed_system(Dataset::BedDefault, Backend::System)]
#[ignore = "skipped: requires bgzip and tabix"]
fn test_tabix_parser_sub_region_query(#[case] dataset: Dataset, #[case] backend: Backend) -> io::Result<()> {
    let toolchain = get_integration_toolchain();
    let dir = create_integration_temp_dir(&format!(
        "tabix-parser-sub-region-{:?}-{:?}",
        dataset, backend
    ))?;

    let parser = open_dataset_parser(&dir, &toolchain, dataset, backend)?;
    assert_query_records_multistep(parser, "chr1:10-20", expected_rows(dataset, 1))?;

    fs::remove_dir_all(dir)?;
    Ok(())
}

#[rstest]
#[case::gff3_noodles(Dataset::Gff3Default, Backend::Noodles)]
#[case::gff3_system(Dataset::Gff3Default, Backend::System)]
#[case::bed_noodles(Dataset::BedDefault, Backend::Noodles)]
#[case::bed_system(Dataset::BedDefault, Backend::System)]
#[ignore = "skipped: requires bgzip and tabix"]
fn test_tabix_parser_empty_region_query(#[case] dataset: Dataset, #[case] backend: Backend) -> io::Result<()> {
    let toolchain = get_integration_toolchain();
    let dir = create_integration_temp_dir(&format!(
        "tabix-parser-empty-region-{:?}-{:?}",
        dataset, backend
    ))?;

    let parser = open_dataset_parser(&dir, &toolchain, dataset, backend)?;
    assert_query_records(parser, "chr1:1000-2000", &[])?;

    fs::remove_dir_all(dir)?;
    Ok(())
}

#[rstest]
#[case::gff3_noodles(Dataset::Gff3Default, Backend::Noodles)]
#[case::gff3_system(Dataset::Gff3Default, Backend::System)]
#[case::bed_noodles(Dataset::BedDefault, Backend::Noodles)]
#[case::bed_system(Dataset::BedDefault, Backend::System)]
#[ignore = "skipped: requires bgzip and tabix"]
fn test_tabix_parser_query_just_before(#[case] dataset: Dataset, #[case] backend: Backend) -> io::Result<()> {
    let toolchain = get_integration_toolchain();
    let dir = create_integration_temp_dir(&format!(
        "tabix-parser-just-before-{:?}-{:?}",
        dataset, backend
    ))?;

    let parser = open_dataset_parser(&dir, &toolchain, dataset, backend)?;
    assert_query_records(parser, "chr2:99-99", &[])?;

    fs::remove_dir_all(dir)?;
    Ok(())
}

#[rstest]
#[case::gff3_noodles(Dataset::Gff3Default, Backend::Noodles)]
#[case::gff3_system(Dataset::Gff3Default, Backend::System)]
#[case::bed_noodles(Dataset::BedDefault, Backend::Noodles)]
#[case::bed_system(Dataset::BedDefault, Backend::System)]
#[ignore = "skipped: requires bgzip and tabix"]
fn test_tabix_parser_query_exactly_on_start(
    #[case] dataset: Dataset,
    #[case] backend: Backend,
) -> io::Result<()> {
    let toolchain = get_integration_toolchain();
    let dir = create_integration_temp_dir(&format!(
        "tabix-parser-exact-start-{:?}-{:?}",
        dataset, backend
    ))?;

    let parser = open_dataset_parser(&dir, &toolchain, dataset, backend)?;
    assert_query_records(parser, "chr2:100-100", expected_rows(dataset, 2))?;

    fs::remove_dir_all(dir)?;
    Ok(())
}

#[rstest]
#[case::gff3_noodles(Dataset::Gff3Default, Backend::Noodles)]
#[case::gff3_system(Dataset::Gff3Default, Backend::System)]
#[case::bed_noodles(Dataset::BedDefault, Backend::Noodles)]
#[case::bed_system(Dataset::BedDefault, Backend::System)]
#[ignore = "skipped: requires bgzip and tabix"]
fn test_tabix_parser_query_exactly_on_end(
    #[case] dataset: Dataset,
    #[case] backend: Backend,
) -> io::Result<()> {
    let toolchain = get_integration_toolchain();
    let dir = create_integration_temp_dir(&format!(
        "tabix-parser-exact-end-{:?}-{:?}",
        dataset, backend
    ))?;

    let parser = open_dataset_parser(&dir, &toolchain, dataset, backend)?;
    assert_query_records(parser, "chr2:120-120", expected_rows(dataset, 2))?;

    fs::remove_dir_all(dir)?;
    Ok(())
}

#[rstest]
#[case::gff3_noodles(Dataset::Gff3Default, Backend::Noodles)]
#[case::gff3_system(Dataset::Gff3Default, Backend::System)]
#[case::bed_noodles(Dataset::BedDefault, Backend::Noodles)]
#[case::bed_system(Dataset::BedDefault, Backend::System)]
#[ignore = "skipped: requires bgzip and tabix"]
fn test_tabix_parser_query_just_after(#[case] dataset: Dataset, #[case] backend: Backend) -> io::Result<()> {
    let toolchain = get_integration_toolchain();
    let dir = create_integration_temp_dir(&format!(
        "tabix-parser-just-after-{:?}-{:?}",
        dataset, backend
    ))?;

    let parser = open_dataset_parser(&dir, &toolchain, dataset, backend)?;
    assert_query_records(parser, "chr2:121-121", &[])?;

    fs::remove_dir_all(dir)?;
    Ok(())
}

#[rstest]
#[case::bed0base_noodles(Dataset::Bed0Base, Backend::Noodles)]
#[case::bed0base_system(Dataset::Bed0Base, Backend::System)]
#[ignore = "skipped: requires bgzip and tabix"]
fn test_tabix_parser_bed_0_position(#[case] dataset: Dataset, #[case] backend: Backend) -> io::Result<()> {
    let toolchain = get_integration_toolchain();
    let dir = create_integration_temp_dir(&format!(
        "tabix-parser-0-based-{:?}-{:?}",
        dataset, backend
    ))?;

    let parser = open_dataset_parser(&dir, &toolchain, dataset, backend)?;
    match backend {
        Backend::Noodles => assert_query_records(parser, "chr1:1-1", expected_rows(Dataset::Bed0Base, 0))?,
        Backend::System => assert_query_records(parser, "chr1:1-1", expected_rows(Dataset::Bed0Base, 1))?,
    }

    fs::remove_dir_all(dir)?;
    Ok(())
}

#[rstest]
#[case::bedsamepos_noodles(Dataset::BedSamePos, Backend::Noodles)]
#[case::bedsamepos_system(Dataset::BedSamePos, Backend::System)]
#[ignore = "skipped: requires bgzip and tabix"]
fn test_tabix_parser_bed_same_position(#[case] dataset: Dataset, #[case] backend: Backend) -> io::Result<()> {
    let toolchain = get_integration_toolchain();
    let dir = create_integration_temp_dir(&format!(
        "tabix-parser-same-position-{:?}-{:?}",
        dataset, backend
    ))?;
    
    let parser = open_dataset_parser(&dir, &toolchain, dataset, backend)?;
    assert_query_records(parser, "chr1:99-100", expected_rows(Dataset::BedSamePos, 0))?;

    fs::remove_dir_all(dir)?;
    Ok(())
}

#[rstest]
#[case::gff3_noodles(Dataset::Gff3Default, Backend::Noodles)]
#[case::gff3_system(Dataset::Gff3Default, Backend::System)]
#[case::bed_noodles(Dataset::BedDefault, Backend::Noodles)]
#[case::bed_system(Dataset::BedDefault, Backend::System)]
#[ignore = "skipped: requires bgzip and tabix"]
fn test_tabix_parser_close_and_reopen(#[case] dataset: Dataset, #[case] backend: Backend) -> io::Result<()> {
    let toolchain = get_integration_toolchain();
    let dir = create_integration_temp_dir(&format!("tabix-parser-reopen-{:?}-{:?}", dataset, backend))?;
    let path = write_dataset_using_system_tools(&dir, &toolchain, dataset)?;

    let parser = match backend {
        Backend::Noodles => open_noodles_parser(&path)?,
        Backend::System => open_system_parser(&path, &toolchain)?,
    };
    assert_close_and_reopen(parser, &path, expected_rows(dataset, 0), expected_rows(dataset, 1), expected_rows(dataset, 2))?;

    fs::remove_dir_all(dir)?;
    Ok(())
}

#[rstest]
#[case::gff3_noodles(Dataset::Gff3Default, Backend::Noodles)]
#[case::gff3_system(Dataset::Gff3Default, Backend::System)]
#[case::bed_noodles(Dataset::BedDefault, Backend::Noodles)]
#[case::bed_system(Dataset::BedDefault, Backend::System)]
#[ignore = "skipped: requires bgzip and tabix"]
fn test_tabix_parser_invalid_region(#[case] dataset: Dataset, #[case] backend: Backend) -> io::Result<()> {
    let toolchain = get_integration_toolchain();
    let dir = create_integration_temp_dir(&format!("tabix-parser-invalid-{:?}-{:?}", dataset, backend))?;
    let path = write_dataset_using_system_tools(&dir, &toolchain, dataset)?;
    let mut parser = match backend {
        Backend::Noodles => open_noodles_parser(&path)?,
        Backend::System => open_system_parser(&path, &toolchain)?,
    };

    match backend {
        Backend::Noodles => {
            let err = parser.seek(":::invalid:::").unwrap_err();
            assert!(
                err.to_string().contains("missing reference sequence name"),
                "unexpected error: {err}"
            );
        }
        Backend::System => match parser.seek(":::invalid:::") {
            Ok(found) => assert!(!found),
            Err(err) => assert!(
                err.to_string().contains("tabix executable failed"),
                "unexpected error: {err}"
            ),
        },
    }

    fs::remove_dir_all(dir)?;
    Ok(())
}

fn open_noodles_parser(path: &Path) -> io::Result<TabixParser> {
    open_tabix_parser(path, Path::new("tabix-nonexistent-binary-for-test"))
}

fn open_system_parser(path: &Path, toolchain: &HtslibToolchain) -> io::Result<TabixParser> {
    open_tabix_parser(path, &toolchain.tabix)
}

fn open_dataset_parser(
    dir: &Path,
    toolchain: &HtslibToolchain,
    dataset: Dataset,
    backend: Backend,
) -> io::Result<TabixParser> {
    let path = match dataset {
        Dataset::Gff3Default => write_gff3_using_system_tools(dir, toolchain, None)?,
        Dataset::BedDefault => write_bed_using_system_tools(dir, toolchain, None)?,
        Dataset::Bed0Base => write_bed_using_system_tools(dir, toolchain, Some(&BED_0_BASE_LINES))?,
        Dataset::BedSamePos => write_bed_using_system_tools(dir, toolchain, Some(&BED_SAME_POS_LINES))?,
    };

    match backend {
        Backend::Noodles => open_noodles_parser(&path),
        Backend::System => open_system_parser(&path, toolchain),
    }
}

fn write_dataset_using_system_tools(
    dir: &Path,
    toolchain: &HtslibToolchain,
    dataset: Dataset,
) -> io::Result<PathBuf> {
    match dataset {
        Dataset::Gff3Default => write_gff3_using_system_tools(dir, toolchain, None),
        Dataset::BedDefault => write_bed_using_system_tools(dir, toolchain, None),
        Dataset::Bed0Base => write_bed_using_system_tools(dir, toolchain, Some(&BED_0_BASE_LINES)),
        Dataset::BedSamePos => write_bed_using_system_tools(dir, toolchain, Some(&BED_SAME_POS_LINES)),
    }
}

fn expected_rows(dataset: Dataset, case_no: usize) -> &'static [&'static str] {
    match (dataset, case_no) {
        (Dataset::Gff3Default, 0) => &GFF3_LINES[0..2],
        (Dataset::BedDefault, 0) => &BED_LINES[0..1],
        (Dataset::Gff3Default, 1) => &GFF3_LINES[2..5],
        (Dataset::BedDefault, 1) => &BED_LINES[1..4],
        (Dataset::Gff3Default, 2) => &GFF3_LINES[5..6],
        (Dataset::BedDefault, 2) => &BED_LINES[4..5],
        (Dataset::Bed0Base, 0) => &BED_0_BASE_LINES,
        (Dataset::Bed0Base, 1) => &[],
        (Dataset::BedSamePos, 0) => &BED_SAME_POS_LINES,
        _ => &[],
    }
}

fn open_tabix_parser(path: &Path, tabix_path: &Path) -> io::Result<TabixParser> {
    Python::initialize();
    TabixParser::new(
        Some(path.to_string_lossy().into_owned()),
        Some(tabix_path.to_string_lossy().into_owned()),
    )
    .map_err(py_err)
}

fn assert_header_lines(mut parser: TabixParser, expected_headers: &[&str]) -> io::Result<()> {
    assert_eq!(
        parser.header_lines().map_err(py_err)?,
        expected_headers.iter().map(|line| (*line).to_owned()).collect::<Vec<_>>()
    );

    Ok(())
}

fn assert_query_records_multistep(
    mut parser: TabixParser,
    region: &str,
    expected_rows: &[&str],
) -> io::Result<()> {
    assert!(parser.seek(region).map_err(py_err)?);
    assert_eq!(parser.records(Some(0)).map_err(py_err)?, Vec::<String>::new());
    assert_eq!(
        parser.records(Some(2)).map_err(py_err)?,
        expected_rows[..2]
            .iter()
            .map(|row| (*row).to_owned())
            .collect::<Vec<_>>()
    );
    assert_eq!(
        parser.records(None).map_err(py_err)?,
        expected_rows[2..]
            .iter()
            .map(|row| (*row).to_owned())
            .collect::<Vec<_>>()
    );
    assert_eq!(parser.next().map_err(py_err)?, None);

    Ok(())
}

fn assert_query_records(
    mut parser: TabixParser,
    region: &str,
    expected_rows: &[&str],
) -> io::Result<()> {
    let found = parser.seek(region).map_err(py_err)?;
    assert_eq!(
        found, 
        !expected_rows.is_empty(), 
        "expected: {} records for region {region}, but got {}",
        expected_rows.len(),
        if found { "some" } else { "none" }
    );
    assert_eq!(
        parser.records(None).map_err(py_err)?,
        expected_rows
            .iter()
            .map(|row| (*row).to_owned())
            .collect::<Vec<_>>()
    );

    Ok(())
}

fn assert_close_and_reopen(
    mut parser: TabixParser,
    path: &Path,
    expected_headers: &[&str],
    expected_chr1_rows: &[&str],
    expected_chr2_rows: &[&str],
) -> io::Result<()> {
    assert!(parser.seek("chr2:100-120").map_err(py_err)?);
    assert_eq!(
        parser.next().map_err(py_err)?,
        expected_chr2_rows.first().map(|row| (*row).to_owned())
    );

    parser.close();
    let err = parser.next().unwrap_err();
    assert!(
        err.to_string().contains("not open"),
        "unexpected error: {err}"
    );

    parser
        .open(path.to_string_lossy().as_ref())
        .map_err(py_err)?;
    assert_eq!(
        parser.header_lines().map_err(py_err)?,
        expected_headers
            .iter()
            .map(|row| (*row).to_owned())
            .collect::<Vec<_>>()
    );
    assert!(parser.seek("chr1:10-20").map_err(py_err)?);
    assert_eq!(
        parser.records(None).map_err(py_err)?,
        expected_chr1_rows
            .iter()
            .map(|row| (*row).to_owned())
            .collect::<Vec<_>>()
    );

    Ok(())
}

fn get_integration_toolchain() -> HtslibToolchain {
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

fn create_integration_temp_dir(prefix: &str) -> io::Result<PathBuf> {
    let dir = env::temp_dir().join(format!("{prefix}-{}", std::process::id()));
    if dir.exists() {
        fs::remove_dir_all(&dir)?;
    }
    fs::create_dir_all(&dir)?;
    Ok(dir)
}

fn write_gff3_using_system_tools(
    dir: &Path,
    toolchain: &HtslibToolchain,
    lines: Option<&[&str]>,
) -> io::Result<PathBuf> {
    write_tabix_fixture_using_system_tools(
        dir,
        toolchain,
        "sample.gff3",
        "gff",
        lines.unwrap_or(&GFF3_LINES),
    )
}

fn write_bed_using_system_tools(
    dir: &Path,
    toolchain: &HtslibToolchain,
    lines: Option<&[&str]>,
) -> io::Result<PathBuf> {
    write_tabix_fixture_using_system_tools(
        dir,
        toolchain,
        "sample.bed",
        "bed",
        lines.unwrap_or(&BED_LINES),
    )
}

fn write_tabix_fixture_using_system_tools(
    dir: &Path,
    toolchain: &HtslibToolchain,
    file_name: &str,
    tabix_preset: &str,
    lines: &[&str],
) -> io::Result<PathBuf> {
    let source_path = write_tabix_source_file(dir, file_name, lines)?;
    let gz_path = dir.join(format!("{file_name}.gz"));

    compress_with_bgzip(&source_path, &gz_path, &toolchain.bgzip)?;
    index_with_tabix(&gz_path, tabix_preset, &toolchain.tabix)?;

    Ok(gz_path)
}

fn write_tabix_source_file(dir: &Path, file_name: &str, lines: &[&str]) -> io::Result<PathBuf> {
    let source_path = dir.join(file_name);
    let mut content = lines.join("\n");
    content.push('\n');
    fs::write(&source_path, content)?;
    Ok(source_path)
}

fn compress_with_bgzip(source_path: &Path, gz_path: &Path, bgzip_path: &Path) -> io::Result<()> {
    let output = fs::File::create(gz_path)?;
    let status = Command::new(bgzip_path)
        .arg("-c")
        .arg(source_path)
        .stdout(Stdio::from(output))
        .status()?;
    if !status.success() {
        return Err(io::Error::other(format!(
            "bgzip failed with status {status}"
        )));
    }

    Ok(())
}

fn index_with_tabix(gz_path: &Path, tabix_preset: &str, tabix_path: &Path) -> io::Result<()> {
    let status = Command::new(tabix_path)
        .arg("-p")
        .arg(tabix_preset)
        .arg(gz_path)
        .status()?;
    if !status.success() {
        return Err(io::Error::other(format!(
            "tabix failed with status {status}"
        )));
    }

    Ok(())
}

fn py_err(err: PyErr) -> io::Error {
    io::Error::other(err.to_string())
}
