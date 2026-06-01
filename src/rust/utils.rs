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

pub(crate) fn trim_newline(line: &mut String) {
    if line.ends_with('\n') {
        line.pop();
        if line.ends_with('\r') {
            line.pop();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trim_newline_lf() {
        let mut line = "hello\n".to_string();
        trim_newline(&mut line);
        assert_eq!(line, "hello");
    }

    #[test]
    fn test_trim_newline_crlf() {
        let mut line = "hello\r\n".to_string();
        trim_newline(&mut line);
        assert_eq!(line, "hello");
    }

    #[test]
    fn test_trim_newline_no_newline() {
        let mut line = "hello".to_string();
        trim_newline(&mut line);
        assert_eq!(line, "hello");
    }

    #[test]
    fn test_trim_newline_empty() {
        let mut line = String::new();
        trim_newline(&mut line);
        assert_eq!(line, "");
    }
}
