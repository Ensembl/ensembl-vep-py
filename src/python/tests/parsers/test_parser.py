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
"""Unit tests for the parser interface."""

from attrs import define
import pytest

from ensembl.vep.parsers.parser import Parser


_RECORDS: list[str] = []


@define(kw_only=True, auto_attribs=True)
class ConcreteParser(Parser[str]):
    """Concrete parser used to test Parser interface helpers."""

    def open(self, file_path: str) -> "ConcreteParser":
        """Open a file for parsing."""

        return self

    def close(self) -> None:
        """Close the parser."""

        return None

    def seek(
        self,
        region: str | None = None,
        *,
        seqid: str | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> bool:
        """Seek to a region."""

        return bool(region or seqid or start or end)

    def next(self) -> str | None:
        """Return the next row."""

        if not _RECORDS:
            return None

        return _RECORDS.pop(0)


def _parser_with_records(records: list[str]) -> ConcreteParser:
    """Create a parser with records for Parser.records tests."""

    global _RECORDS  # pylint: disable=global-statement

    _RECORDS = records.copy()
    return ConcreteParser(format="dummy")


class TestParserInitialization:
    """Test Parser initialization."""

    def test_abstract_parser_cannot_be_instantiated(self) -> None:
        """Test Parser cannot be instantiated directly."""

        with pytest.raises(TypeError, match="abstract class"):
            _ = Parser(format="tabix")  # type: ignore[abstract], pylint: disable=abstract-class-instantiated

    def test_concrete_parser_sets_format(self) -> None:
        """Test concrete parser stores the parser format."""

        parser = ConcreteParser(format="dummy")

        assert parser.format == "dummy"


class TestParserMetadata:
    """Tests for Parser.metadata."""

    def test_default_metadata_is_empty(self) -> None:
        """Test default metadata is an empty dictionary."""

        parser = ConcreteParser(format="dummy")

        assert not parser.metadata


class TestParserRecords:
    """Tests for Parser.records."""

    def test_records_returns_all_records(self) -> None:
        """Test records consumes rows until next returns None."""

        parser = _parser_with_records(["one", "two", "three"])

        assert parser.records() == ["one", "two", "three"]
        assert not parser.records()

    def test_records_respects_limit(self) -> None:
        """Test records stops once the limit is reached."""

        parser = _parser_with_records(["one", "two", "three"])

        assert parser.records(limit=2) == ["one", "two"]
        assert parser.records() == ["three"]

    def test_records_zero_limit_does_not_call_next(self) -> None:
        """Test zero limit returns no records without consuming rows."""

        parser = _parser_with_records(["one"])

        assert not parser.records(limit=0)
        assert _RECORDS == ["one"]

    def test_records_stops_when_next_returns_none_before_limit(self) -> None:
        """Test records stops early when the parser is exhausted."""

        parser = _parser_with_records(["one"])

        assert parser.records(limit=3) == ["one"]
        assert not _RECORDS
