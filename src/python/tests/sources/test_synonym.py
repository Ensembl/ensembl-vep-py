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
"""Tests for sequence-region synonyms and resolving names against one source."""

from __future__ import annotations
import logging

from pathlib import Path

import pytest

from ensembl.vep.sources.synonym import ChromosomeSynonyms, SequenceRegionResolver

CHR1_ALL_SYNONYMS_TEXT = """\
1\tchr1
1\tCM000663.2
1\tNC_000001.11
CM000663.2\tchr1
CM000663.2\tNC_000001.11
NC_000001.11\tchr1
"""

CHR1_NOT_ALL_COMINATIONS_TEXT = """\
1\tCM000663.2
CM000663.2\tNC_000001.11
"""


def _synonyms(reference_name: str, *synonyms: str) -> ChromosomeSynonyms:
    """Return a map of one region's names, as reading a one-line synonyms file would give."""

    return ChromosomeSynonyms(
        by_name={
            reference_name: frozenset(synonyms),
            **{synonym: frozenset({reference_name}) for synonym in synonyms},
        }
    )


class TestChromosomeSynonyms:
    """Reading Perl VEP's synonyms file."""

    def test_reads_names_in_both_directions(self, tmp_path: Path) -> None:
        """Each synonym maps to the reference name and the reference name to all of them."""

        path = tmp_path / "synonyms.txt"
        path.write_text("1\tchr1\n1\tNC_000001.11\n2\tchr2\n", encoding="utf-8")

        synonyms = ChromosomeSynonyms.from_file(path)

        assert synonyms.synonyms_of("1") == frozenset({"chr1", "NC_000001.11"})
        assert synonyms.synonyms_of("NC_000001.11") == frozenset({"1"})
        assert synonyms.synonyms_of("chr2") == frozenset({"2"})

    def test_does_not_pair_two_synonyms_of_one_name(self) -> None:
        """Perl records reference-to-synonym pairs only, and so do we."""

        synonyms = _synonyms("1", "chr1", "NC_000001.11")

        assert synonyms.synonyms_of("chr1") == frozenset({"1"})

    def test_ignores_a_line_naming_nothing_to_map(self, tmp_path: Path) -> None:
        """A blank line, or one lone name, describes no pair."""

        path = tmp_path / "synonyms.txt"
        path.write_text("\n1\n  \n1\tchr1\n", encoding="utf-8")

        synonyms = ChromosomeSynonyms.from_file(path)

        assert synonyms.by_name == {"1": frozenset({"chr1"}), "chr1": frozenset({"1"})}

    def test_an_empty_map_has_no_synonyms(self) -> None:
        """The default is no synonyms at all, which every resolver can take."""

        assert not ChromosomeSynonyms().synonyms_of("1")

    def test_a_missing_file_is_reported_with_its_path(self, tmp_path: Path) -> None:
        """A synonyms file the caller asked for and that is not there stops the run."""

        with pytest.raises(OSError, match="absent.txt"):
            ChromosomeSynonyms.from_file(tmp_path / "absent.txt")

    def test_more_than_two_columns_gives_a_warning(
        self,
        tmp_path: Path,
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """A synonyms file with more than two columns triggers a warning."""
        path = tmp_path / "synonyms.txt"
        path.write_text("1\tchr1\t1\tNC_000001.11\n2\tchr2\n", encoding="utf-8")

        with caplog.at_level(logging.WARNING):
            _ = ChromosomeSynonyms.from_file(path)
            assert (
                "Only two columns are expected in the synonyms file" in caplog.text
            ), "Incorrect warning in ChromosomeSynonyms.from_file"


class TestSequenceRegionResolver:
    """Uploaded-to-source name resolution."""

    def test_keeps_an_exact_source_name(self) -> None:
        """Prefer an exact source name to an explicit synonym."""

        resolver = SequenceRegionResolver(
            source_names=frozenset({"1", "chr1"}),
            synonyms=_synonyms("1", "chr1"),
        )

        assert resolver.resolve("1") == "1"
        assert resolver.resolve("1") == "1"

    def test_adds_a_chr_prefix_case_insensitively(self) -> None:
        """Resolve a bare uploaded name to a source chr name."""

        resolver = SequenceRegionResolver(source_names=frozenset({"Chr1"}))

        assert resolver.resolve("1") == "Chr1"

    def test_removes_a_chr_prefix_case_insensitively(self) -> None:
        """Resolve a prefixed uploaded name to a bare source name."""

        resolver = SequenceRegionResolver(source_names=frozenset({"1"}))

        assert resolver.resolve("cHr1") == "1"

    def test_uses_an_explicit_synonym_before_chr_adjustment(self) -> None:
        """Resolve a configured RefSeq-style alias exactly as mapped."""

        resolver = SequenceRegionResolver(
            source_names=frozenset({"chr1"}),
            synonyms=_synonyms("NC_000001.11", "chr1"),
        )

        assert resolver.resolve("NC_000001.11") == "chr1"

    def test_resolves_a_synonym_in_either_direction(self) -> None:
        """One map serves both: the FASTA wanting "chr1" and a transcript file wanting the accession."""

        synonyms = _synonyms("NC_000001.11", "chr1")
        to_accession = SequenceRegionResolver(source_names=frozenset({"NC_000001.11"}), synonyms=synonyms)

        assert to_accession.resolve("chr1") == "NC_000001.11"

    def test_picks_the_lowest_declared_synonym_when_several_fit(self) -> None:
        """Perl takes whichever hash key comes first; sorting makes the answer repeatable."""

        resolver = SequenceRegionResolver(
            source_names=frozenset({"chr1", "NC_000001.11"}),
            synonyms=ChromosomeSynonyms(by_name={"1": frozenset({"chr1", "NC_000001.11"})}),
        )

        assert resolver.resolve("1") == "NC_000001.11"

    def test_ignores_a_synonym_the_source_does_not_declare(self) -> None:
        """A synonym for a contig this source has no row for leaves the name as uploaded.

        The map is the run's, not the source's: it legitimately names contigs a given FASTA or
        transcript file does not carry, and that is not an error in the map.
        """

        resolver = SequenceRegionResolver(
            source_names=frozenset({"chr1"}),
            synonyms=_synonyms("MT", "chrM"),
        )

        assert resolver.resolve("MT") == "MT"

    def test_keeps_an_unresolved_uploaded_name(self) -> None:
        """Fall back deterministically to the original spelling."""

        resolver = SequenceRegionResolver(source_names=frozenset({"chr1"}))

        assert resolver.resolve("unplaced") == "unplaced"

    def test_keeps_resolution_caches_per_instance(self) -> None:
        """Avoid sharing one source's aliases with another source."""

        first = SequenceRegionResolver(source_names=frozenset({"chr1"}))
        second = SequenceRegionResolver(source_names=frozenset({"1"}))

        assert first.resolve("1") == "chr1"
        assert second.resolve("1") == "1"

    def test_rejects_case_folded_source_name_collisions(self) -> None:
        """Require case-insensitive matching to have one possible result."""

        with pytest.raises(ValueError, match="differ only by case"):
            SequenceRegionResolver(source_names=frozenset({"chrM", "CHRM"}))

    def test_chr_synonym_from_contig(self, tmp_path: Path) -> None:
        """Get contig name for chromosome name and vice-versa."""

        path = tmp_path / "synonyms.txt"
        path.write_text("1\tchr1\n1\tNC_000001.11\nchr1\tNC_000001.11\n2\tchr2\n", encoding="utf-8")
        synonyms = ChromosomeSynonyms.from_file(path)
        contig_resolver = SequenceRegionResolver(source_names=frozenset({"NC_000001.11"}), synonyms=synonyms)
        assert contig_resolver.resolve("chr1") == "NC_000001.11"

        chr_resolver = SequenceRegionResolver(source_names=frozenset({"chr1"}), synonyms=synonyms)
        assert chr_resolver.resolve("NC_000001.11") == "chr1"

    def test_all_chr1_synonyms(self, tmp_path: Path) -> None:
        """Test that all possible combinations exist in the resolver."""

        path = tmp_path / "synonyms.txt"
        path.write_text(CHR1_ALL_SYNONYMS_TEXT, encoding="utf-8")
        synonyms = ChromosomeSynonyms.from_file(path)

        chr_num_resolver = SequenceRegionResolver(source_names=frozenset({"1"}), synonyms=synonyms)
        assert chr_num_resolver.resolve("NC_000001.11") == "1"
        assert chr_num_resolver.resolve("chr1") == "1"
        assert chr_num_resolver.resolve("CM000663.2") == "1"

        chr_resolver = SequenceRegionResolver(source_names=frozenset({"chr1"}), synonyms=synonyms)
        assert chr_resolver.resolve("NC_000001.11") == "chr1"
        assert chr_resolver.resolve("1") == "chr1"
        assert chr_resolver.resolve("CM000663.2") == "chr1"

        contig_resolver = SequenceRegionResolver(source_names=frozenset({"NC_000001.11"}), synonyms=synonyms)
        assert contig_resolver.resolve("1") == "NC_000001.11"
        assert contig_resolver.resolve("chr1") == "NC_000001.11"
        assert contig_resolver.resolve("CM000663.2") == "NC_000001.11"

        accession_resolver = SequenceRegionResolver(source_names=frozenset({"CM000663.2"}), synonyms=synonyms)
        assert accession_resolver.resolve("chr1") == "CM000663.2"
        assert accession_resolver.resolve("1") == "CM000663.2"
        assert accession_resolver.resolve("NC_000001.11") == "CM000663.2"

    def test_not_all_synonym_links_in_file(self, tmp_path: Path) -> None:
        """Test whether all combinations need to be present in synonyms file"""

        path = tmp_path / "synonyms.txt"
        path.write_text(CHR1_NOT_ALL_COMINATIONS_TEXT, encoding="utf-8")
        synonyms = ChromosomeSynonyms.from_file(path)

        # 1 is connected to CM000663.2
        accession_resolver = SequenceRegionResolver(source_names=frozenset({"CM000663.2"}), synonyms=synonyms)
        assert accession_resolver.resolve("1") == "CM000663.2"

        contig_resolver = SequenceRegionResolver(source_names=frozenset({"NC_000001.11"}), synonyms=synonyms)
        # CM000663.2 is connected to NC_000001.11
        assert contig_resolver.resolve("CM000663.2") == "NC_000001.11"
        # but 1 is not connected to NC_000001.11
        assert contig_resolver.resolve("1") == "1"
