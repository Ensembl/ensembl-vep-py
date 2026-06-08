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
}
