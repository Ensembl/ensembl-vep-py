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

import logging
from pathlib import Path

import pytest

from ensembl.vep.config import Config, load
from ensembl.vep.config.parser.config import read_config_file, read_legacy_config_file


@pytest.fixture(name="toml_file")
def fixture_toml_file(tmp_path: Path) -> Path:
    """TOML file setting two corpus options."""
    path = tmp_path / "vep.toml"
    path.write_text("[corpus]\nfuzz = 111\nstep = 222\n", encoding="utf-8")
    return path


PERL_CONFIG = """\
###############################################################################
# Example VEP configuration, in the layout of an Ensembl VEP installation's own
###############################################################################

# --- input ---
input_file\tvariants.vcf
--fasta "/data/ref genome/GRCh38.fa"
gff /data/GRCh38.gff3.gz
format vcf
synonyms "/data/chr synonyms/chr_synonyms.txt

# --- output ---
output_file /scratch/run#2/annotated.vcf
output_file /scratch/ignored.vcf
no_stats
log_level debug
log_level error

# --- annotation ---
distance 5000,2000
custom /data/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG
custom "/data/my annotations.bed.gz" MyBed bed overlap 0

# --- bench end-to-end ---
batch_size 512

# options this implementation has not grown yet
no_clean
sift b
symbol
cache_version 110
"""


@pytest.fixture(name="perl_config")
def fixture_perl_config(tmp_path: Path) -> Path:
    """A Perl VEP config file with tabs, dashes, quotes, repeats and unknown options."""
    path = tmp_path / "vep.ini"
    path.write_text(PERL_CONFIG, encoding="utf-8")
    return path


@pytest.fixture(name="perl_layer")
def fixture_perl_layer(perl_config: Path) -> dict[str, dict[str, object]]:
    """The layer read from the Perl config file."""
    types = Config.section_types_for(("source", "run", "logging"))
    return read_legacy_config_file(perl_config, types)


class TestLegacyConfigFiles:
    """Reading a Perl VEP config file as Perl VEP does."""

    def test_comments_blanks_and_unknown_options_are_skipped(self, perl_layer: dict) -> None:
        """Only declared options reach the layer."""
        assert set(perl_layer) == {"source", "run", "logging"}

    def test_a_tab_separates_key_and_value(self, perl_layer: dict) -> None:
        """Tab-separated lines are read."""
        assert perl_layer["source"]["input_file"] == Path("variants.vcf")

    def test_leading_dashes_on_a_key_are_dropped(self, perl_layer: dict) -> None:
        """A key may be written as a flag."""
        assert perl_layer["source"]["fasta"] == Path("/data/ref genome/GRCh38.fa")

    def test_a_quoted_value_keeps_its_spaces_and_loses_its_quotes(self, perl_layer: dict) -> None:
        """Quotes group words and are removed."""
        assert "ref genome" in str(perl_layer["source"]["fasta"])
        assert '"' not in str(perl_layer["source"]["fasta"])

    def test_a_missing_end_quote_is_not_an_error(self, perl_layer: dict) -> None:
        """Opening quote without closing quote still returns value."""
        assert "chr synonyms" in str(perl_layer["source"]["synonyms"])
        assert '"' not in str(perl_layer["source"]["synonyms"])

    def test_a_hash_inside_a_value_is_part_of_the_value(self, perl_layer: dict) -> None:
        """Only a line that starts with a hash is a comment."""
        assert perl_layer["run"]["output_file"] == Path("/scratch/run#2/annotated.vcf")

    def test_the_first_occurrence_of_a_repeated_key_wins(self, perl_layer: dict) -> None:
        """A repeated key does not override an earlier one in the same file."""
        assert perl_layer["logging"]["log_level"] == "debug"

    def test_a_repeated_text_tuple_option_accumulates(self, perl_layer: dict) -> None:
        """Each custom line adds one item, kept whole."""
        assert perl_layer["source"]["custom"] == (
            "/data/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG",
            "/data/my annotations.bed.gz MyBed bed overlap 0",
        )

    def test_a_numeric_tuple_option_still_splits_on_commas(self, perl_layer: dict) -> None:
        """A comma-separated number list is one value, not an accumulating option."""
        assert perl_layer["source"]["distance"] == (5000, 2000)

    def test_a_bare_key_is_a_switch(self, perl_layer: dict) -> None:
        """A key with no value reads as true."""
        assert perl_layer["run"]["no_stats"] is True

    def test_no_option_not_recognised(self, perl_config: Path) -> None:
        """Config files do not currently support the automatic creation of an ``no-`` option
        with opposite effect for a boolean option."""
        types = Config.section_types_for(("bench_endtoend",))
        perl_layer = read_legacy_config_file(perl_config, types)
        assert "clean" not in perl_layer["bench_endtoend"]

    def test_a_hash_line_with_leading_space_is_not_a_comment_but_names_no_option(
        self, tmp_path: Path
    ) -> None:
        """Only a hash in the first column is a comment; an indented one is just an unknown key."""
        path = tmp_path / "vep.ini"
        path.write_text("  # fuzz 1\nfuzz 2\n", encoding="utf-8")
        types = Config.section_types_for(("corpus",))
        assert read_legacy_config_file(path, types) == {"corpus": {"fuzz": 2}}

    def test_a_later_layer_still_overrides_an_earlier_one(self, perl_config: Path) -> None:
        """First-wins is per file: the command line still beats the file."""
        config = load_with_file(["--log_level", "info"], perl_config)
        assert config.logging.log_level == "info"

    def test_a_bad_choice_in_a_legacy_file_is_rejected(self, tmp_path: Path) -> None:
        """Choices are checked for the legacy reader too."""
        path = tmp_path / "vep.ini"
        path.write_text("log_level verbose\n", encoding="utf-8")
        with pytest.raises(ValueError, match="logging.log_level: 'verbose' is not one of"):
            read_legacy_config_file(path, Config.section_types_for(("logging",)))


class TestTomlFiles:
    """TOML detection, errors and unknown keys."""

    def test_a_syntax_error_in_a_toml_file_is_an_error(self, tmp_path: Path) -> None:
        """A broken TOML file is not reread as legacy config."""
        path = tmp_path / "vep.toml"
        path.write_text('[stats]\noutput_file = "r.html"\nprofile = yes\n', encoding="utf-8")
        with pytest.raises(ValueError, match=r"vep\.toml: invalid TOML.*line 3"):
            read_config_file(path, Config.section_types_for(("run",)))

    def test_dotted_keys_make_a_toml_file(self, tmp_path: Path) -> None:
        """A document with a table but no header is still TOML."""
        path = tmp_path / "vep.toml"
        path.write_text("corpus.fuzz = 4\n", encoding="utf-8")
        assert read_config_file(path, Config.section_types_for(("corpus",))) == {"corpus": {"fuzz": 4}}

    def test_a_value_outside_its_choices_is_an_error_naming_the_option(self, tmp_path: Path) -> None:
        """Choices are enforced for file values."""
        path = tmp_path / "vep.toml"
        path.write_text('[run]\nstats_format = "pdf"\n', encoding="utf-8")
        with pytest.raises(ValueError, match="run.stats_format: 'pdf' is not one of html, text, json"):
            read_config_file(path, Config.section_types_for(("run",)))

    def test_an_unregistered_annotator_is_an_error(self, tmp_path: Path) -> None:
        """Callable choices are resolved for file values."""
        path = tmp_path / "vep.toml"
        path.write_text('[annotation]\nannotator = "consequence.typo"\n', encoding="utf-8")
        with pytest.raises(ValueError, match="annotation.annotator: 'consequence.typo' is not one of"):
            read_config_file(path, Config.section_types_for(("annotation",)))

    def test_an_unknown_key_in_a_known_table_is_warned_about(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """A misspelt key is reported rather than dropped silently."""
        path = tmp_path / "vep.toml"
        path.write_text("[corpus]\nfuzz = 2\nfuzzy = 3\n", encoding="utf-8")
        with caplog.at_level(logging.WARNING, logger="ensembl.vep"):
            layer = read_config_file(path, Config.section_types_for(("corpus",)))
        assert layer == {"corpus": {"fuzz": 2}}
        assert "unknown option 'fuzzy' in table [corpus]" in caplog.text

    def test_a_key_in_the_wrong_table_is_warned_about(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """A key that belongs to another section is not applied silently."""
        path = tmp_path / "vep.toml"
        path.write_text("[run]\nfuzz = 2\n", encoding="utf-8")
        with caplog.at_level(logging.WARNING, logger="ensembl.vep"):
            layer = read_config_file(path, Config.section_types_for(("run", "corpus")))
        assert not layer
        assert "unknown option 'fuzz' in table [run]" in caplog.text


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

    def test_a_perl_cache_line_reads_as_a_switch_beside_dir_cache(self, tmp_path: Path) -> None:
        """Perl's ``cache 1`` is a switch; the directory comes from ``dir_cache``."""
        path = tmp_path / "vep.conf"
        path.write_text("cache 1\ndir_cache /data/vep\n", encoding="utf-8")

        layer = read_config_file(path, Config.section_types_for(("source",)))

        assert layer == {"source": {"cache": True, "dir_cache": Path("/data/vep")}}

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

    def test_unknown_toml_no_options_are_ignored(self, tmp_path: Path) -> None:
        """Unknown ``no-`` prefixed options in TOML files are ignored."""
        path = tmp_path / "vep.toml"
        path.write_text(
            "[bench]\nno_gate = true\n[bench_endtoend]\nbatch_size = 512\nno_clean = true\n", encoding="utf-8"
        )
        assert read_config_file(path, Config.section_types_for(("bench", "bench_endtoend"))) == {
            "bench": {"no_gate": True},
            "bench_endtoend": {"batch_size": 512},
        }

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


def load_with_file(argv: list[str], config_file: Path) -> Config:
    """Load a logging-only configuration with a file layer and a command line."""
    return load(argv, section_names=("logging",), config_file=config_file, environ={})
