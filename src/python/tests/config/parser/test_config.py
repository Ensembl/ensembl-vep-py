# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Tests for reading configuration from a file."""

from __future__ import annotations

from pathlib import Path

import pytest

from ensembl.vep.config import Config
from ensembl.vep.config.parser.config import read_config_file, read_legacy_config_file


@pytest.fixture(name="toml_file")
def fixture_toml_file(tmp_path: Path) -> Path:
    """TOML file setting two corpus options."""
    path = tmp_path / "vep.toml"
    path.write_text("[corpus]\nfuzz = 111\nstep = 222\n", encoding="utf-8")
    return path


class TestConfigFiles:
    """File loading."""

    def test_toml_tables_map_onto_sections(self, toml_file: Path) -> None:
        """TOML tables map to sections."""
        layer = read_config_file(toml_file, Config.section_types_for(("corpus",)))
        assert layer == {"corpus": {"fuzz": 111, "step": 222}}

    def test_a_perl_style_flat_file_is_read(self, tmp_path: Path) -> None:
        """A flat legacy file is read unchanged."""
        path = tmp_path / "vep.conf"
        path.write_text("# comment\nfuzz 777\n\n", encoding="utf-8")
        assert read_config_file(path, Config.section_types_for(("corpus",))) == {"corpus": {"fuzz": 777}}

    def test_a_valueless_perl_key_reads_as_a_switch(self, tmp_path: Path) -> None:
        """Bare legacy flags read as switches."""
        path = tmp_path / "vep.conf"
        path.write_text("no_stats\n", encoding="utf-8")
        assert read_config_file(path, Config.section_types_for(("run",))) == {"run": {"no_stats": True}}

    def test_unknown_perl_keys_are_ignored(self, tmp_path: Path) -> None:
        """Unknown legacy keys are ignored."""
        path = tmp_path / "vep.conf"
        path.write_text("sift b\nfuzz 5\n", encoding="utf-8")
        assert read_config_file(path, Config.section_types_for(("corpus",))) == {"corpus": {"fuzz": 5}}

    def test_unknown_toml_sections_are_ignored(self, tmp_path: Path) -> None:
        """Unknown TOML sections are ignored."""
        path = tmp_path / "vep.toml"
        path.write_text("[nope]\nx = 1\n[corpus]\nfuzz = 2\n", encoding="utf-8")
        assert read_config_file(path, Config.section_types_for(("corpus",))) == {"corpus": {"fuzz": 2}}

    def test_a_flat_file_may_use_an_equals_sign(self, tmp_path: Path) -> None:
        """A flat file may use an equals sign."""
        path = tmp_path / "vep.conf"
        path.write_text("fuzz = 3\n", encoding="utf-8")
        assert read_config_file(path, Config.section_types_for(("corpus",))) == {"corpus": {"fuzz": 3}}

    def test_the_legacy_reader_can_be_called_directly(self, tmp_path: Path) -> None:
        """The legacy reader is public."""
        path = tmp_path / "vep.conf"
        path.write_text("fuzz 12\n", encoding="utf-8")
        types = Config.section_types_for(("corpus",))
        assert read_legacy_config_file(path, types) == {"corpus": {"fuzz": 12}}

    def test_a_missing_file_raises(self, tmp_path: Path) -> None:
        """Missing files fail."""
        with pytest.raises(FileNotFoundError):
            read_config_file(tmp_path / "absent.toml", Config.section_types_for(("corpus",)))
