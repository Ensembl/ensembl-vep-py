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
    env, fs, io,
    path::{Path, PathBuf},
};

#[cfg(unix)]
pub fn process_exists(pid: u32) -> std::io::Result<bool> {
    let status = std::process::Command::new("kill")
        .arg("-0")
        .arg(pid.to_string())
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status()?;

    Ok(status.success())
}

#[cfg(unix)]
pub fn assert_process_is_gone(pid: u32) -> std::io::Result<()> {
    for _ in 0..10 {
        if !process_exists(pid)? {
            return Ok(());
        }

        std::thread::sleep(std::time::Duration::from_millis(10));
    }

    panic!("expected process {pid} to be gone");
}

#[cfg(unix)]
pub fn create_temp_dir(prefix: &str) -> io::Result<PathBuf> {
    let dir = env::temp_dir().join(format!("{prefix}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir)?;
    Ok(dir)
}

pub fn write_sample_file(dir: &Path, file_name: &str, lines: &[u8]) -> io::Result<PathBuf> {
    let path = dir.join(file_name);
    fs::write(&path, lines)?;
    Ok(path)
}

#[cfg(all(test, unix))]
mod tests {
    use std::{io, process::Command};

    use super::*;

    // -----------------------------------------------------------------------
    // write_sample_file
    // -----------------------------------------------------------------------

    #[test]
    fn test_write_sample_file() -> io::Result<()> {
        let dir = create_temp_dir("test-write-sample-file")?;
        let path = write_sample_file(&dir, "sample.txt", b"sample lines\n")?;

        assert_eq!(fs::read(&path)?, b"sample lines\n");

        fs::remove_dir_all(&dir)?;
        Ok(())
    }

    // ---------------------------------------------------------------------------
    // create_temp_dir
    // ---------------------------------------------------------------------------

    #[test]
    fn test_create_temp_dir() -> io::Result<()> {
        let dir = create_temp_dir("test")?;
        assert!(dir.exists());
        fs::remove_dir_all(&dir)?;
        assert!(!dir.exists());
        Ok(())
    }

    // ---------------------------------------------------------------------------
    // process_exists
    // ---------------------------------------------------------------------------

    #[test]
    fn test_process_exists_with_running_child() -> io::Result<()> {
        let mut child = Command::new("sleep").arg("60").spawn()?;
        let pid = child.id();

        assert!(process_exists(pid)?);

        let _ = child.kill();
        let _ = child.wait();
        assert_process_is_gone(pid)
    }

    // ---------------------------------------------------------------------------
    // assert_process_is_gone
    // ---------------------------------------------------------------------------

    #[test]
    #[should_panic(expected = "expected process")]
    fn test_assert_process_is_gone_panics_for_current_process() {
        assert_process_is_gone(std::process::id()).unwrap();
    }
}
