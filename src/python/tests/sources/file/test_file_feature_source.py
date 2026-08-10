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

# pylint: disable=redefined-outer-name
"""Unit tests for FileSource."""

from __future__ import annotations

from typing import Iterator
from attrs import define, field
import pytest

from ensembl.vep.common import SourceType, Strand, format_region_string
from ensembl.vep.models import Feature
from ensembl.vep.parsers import Parser
from ensembl.vep.sources.sequence_source import SequenceSource
from ensembl.vep.sources.file.file_feature_source import FileFeatureSource

_RECORDS: list[int] = []
_SEEKABLE_REGION: list[str] = []


@define(kw_only=True, auto_attribs=True)
class ConcreteParser(Parser[int]):
    """Concrete parser used to test FileSource."""

    opened: bool = False
    closed: bool = False

    def open(self, file_path: str, mode: str | None = None) -> ConcreteParser:
        """Open a file for parsing."""

        self.opened = True
        return self

    def close(self) -> None:
        """Close the parser."""

        self.closed = True

    def seek(
        self,
        region: str | None = None,
        *,
        seq_region_name: str | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> bool:
        """Seek to a region."""

        parsed_region = format_region_string(region, seq_region_name=seq_region_name, start=start, end=end)
        return parsed_region in _SEEKABLE_REGION

    def next(self) -> int | None:
        """Return the next row."""

        if not _RECORDS:
            return None

        return _RECORDS.pop(0)

    @staticmethod
    def parse_file(file_path: str, limit: int | None = None) -> list[int]:
        """Parse the entire file and return all records."""

        return _RECORDS.copy()


def _get_seekable_parser(records: list[int], seekable_region: list[str] | None = None) -> ConcreteParser:
    """Create a parser with records."""

    global _RECORDS  # pylint: disable=global-statement
    global _SEEKABLE_REGION  # pylint: disable=global-statement

    _RECORDS = records.copy()
    _SEEKABLE_REGION = seekable_region.copy() if seekable_region is not None else []

    return ConcreteParser(format="dummy")


@define(kw_only=True, auto_attribs=True)
class DummyFileFeatureSource(FileFeatureSource[ConcreteParser, Feature]):
    """Concrete file source used for testing FileFeatureSource interface helpers."""

    closed: bool = False
    sequence_source: SequenceSource | None = field(default=None)  # ignore[assignment]

    def _records_to_features(self, records: list[int]) -> list[Feature]:
        return [Feature(start=value, end=value, strand=Strand.FORWARD, slice=None) for value in records]

    def close(self) -> None:
        """Close the source."""

        self.closed = True


class TestFileFeatureSourceInitialization:
    """Test FileFeatureSource initialization."""

    def test_abstract_file_feature_source_cannot_be_instantiated(self) -> None:
        """Test FileFeatureSource cannot be instantiated directly."""

        with pytest.raises(TypeError, match="abstract class"):
            FileFeatureSource(  # type: ignore[abstract]  # pylint: disable=abstract-class-instantiated
                file_path="dummy.txt", adaptor=ConcreteParser(format="dummy"), sequence_source=None
            )

    def test_concrete_file_feature_source_sets_source_type(self) -> None:
        """Test concrete FileFeatureSource sets the source_type default to SourceType.FILE."""

        source = DummyFileFeatureSource(file_path="dummy.txt", adaptor=ConcreteParser(format="dummy"))
        assert source.source_type == SourceType.FILE


class TestFileFeatureSourceGetAllFeatures:
    """Test get_all_features."""

    def test_parse_file_returns_features(self) -> None:
        """Test features are built from adaptor.parse_file()."""

        parser = _get_seekable_parser([1, 2, 3], seekable_region=["1:1-100"])
        source = DummyFileFeatureSource(file_path="dummy.txt", adaptor=parser)

        features = source.get_all_features()

        assert [f.start for f in features] == [1, 2, 3]

    def test_adaptor_is_not_opened(self) -> None:
        """Test parse_file is a one-shot call that bypasses open()/close()."""

        parser = _get_seekable_parser([1, 2, 3], seekable_region=["1:1-100"])
        source = DummyFileFeatureSource(file_path="dummy.txt", adaptor=parser)

        source.get_all_features()

        assert parser.opened is False


class TestFileFeatureSourceGetFeaturesByRegion:
    """Test get_features_by_region."""

    def test_missing_region_args_returns_empty_list(self, caplog: pytest.LogCaptureFixture) -> None:
        """Test an empty list is returned when no region info is given."""
        parser = _get_seekable_parser([1, 2, 3], seekable_region=["1:1-100"])
        source = DummyFileFeatureSource(file_path="dummy.txt", adaptor=parser)

        with caplog.at_level("DEBUG"):
            assert source.get_features_by_region() == []
            assert "Cannot parse region" in caplog.text

    def test_conflicting_region_and_seq_region_name_returns_empty_list(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        """Test an empty list is returned when region and seq_region_name conflict."""

        parser = _get_seekable_parser([1, 2, 3], seekable_region=["1:1-100"])
        source = DummyFileFeatureSource(file_path="dummy.txt", adaptor=parser)

        with caplog.at_level("DEBUG"):
            assert source.get_features_by_region("1:1-100", seq_region_name="1") == []
            assert "Cannot parse region" in caplog.text

    def test_seek_failure_returns_empty_list(self) -> None:
        """Test an empty list is returned when the adaptor can't seek to the region."""

        parser = _get_seekable_parser([1, 2, 3], seekable_region=["1:1-100"])
        source = DummyFileFeatureSource(file_path="dummy.txt", adaptor=parser)

        assert source.get_features_by_region("1:100-200") == []

    def test_seek_success_returns_features(self) -> None:
        """Test features are built from adaptor.records() after a successful seek."""

        parser = _get_seekable_parser([10, 20], seekable_region=["1:1-100"])
        source = DummyFileFeatureSource(file_path="dummy.txt", adaptor=parser)

        features = source.get_features_by_region("1:1-100")

        assert [f.start for f in features] == [10, 20]

    def test_seek_using_seq_region_name_start_end(self) -> None:
        """Test seek using seq_region_name, start, and end."""

        parser = _get_seekable_parser([10], seekable_region=["1:1-100"])
        source = DummyFileFeatureSource(file_path="dummy.txt", adaptor=parser)

        features = source.get_features_by_region(seq_region_name="1", start=1, end=100)

        assert [f.start for f in features] == [10]

    def test_opens_and_closes_adaptor(self) -> None:
        """Test the adaptor is opened and closed around the region query."""

        parser = _get_seekable_parser([10, 20], seekable_region=["1:1-100"])
        source = DummyFileFeatureSource(file_path="dummy.txt", adaptor=parser)

        source.get_features_by_region("1:1-100")

        assert parser.opened is True
        assert parser.closed is True

    def test_closes_adaptor_even_when_seek_fails(self) -> None:
        """Test the adaptor is still closed when seek() returns False."""

        parser = _get_seekable_parser([1, 2, 3], seekable_region=["1:1-100"])
        source = DummyFileFeatureSource(file_path="dummy.txt", adaptor=parser)

        source.get_features_by_region("1:15-200")

        assert parser.closed is True


class TestFileFeatureSourceIterAllFeatures:
    """Test iter_all_features."""

    def test_returns_an_iterator(self) -> None:
        """Test iter_all_features returns a lazy iterator, not a list."""

        parser = _get_seekable_parser([1, 2, 3], seekable_region=["1:1-100"])
        source = DummyFileFeatureSource(file_path="dummy.txt", adaptor=parser)

        result = source.iter_all_features()

        assert isinstance(result, Iterator)

    def test_yields_features_one_at_a_time(self) -> None:
        """Test all features are yielded in order."""

        records = [1, 2, 3]
        parser = _get_seekable_parser(records, seekable_region=["1:1-100"])
        source = DummyFileFeatureSource(file_path="dummy.txt", adaptor=parser)

        for i, feature in enumerate(source.iter_all_features()):
            assert feature.start == records[i]

    def test_opens_adaptor_before_first_yield(self) -> None:
        """Test the adaptor is opened once, lazily, on first consumption."""

        parser = _get_seekable_parser([1, 2, 3], seekable_region=["1:1-100"])
        source = DummyFileFeatureSource(file_path="dummy.txt", adaptor=parser)

        iterator = source.iter_all_features()
        assert parser.opened is False

        first = next(iterator)
        assert first.start == 1
        assert parser.opened is True

    def test_keeps_adaptor_open_across_multiple_next_calls(self) -> None:
        """Test the adaptor session spans the whole generator, not just one record."""

        records = [1, 2, 3]
        parser = _get_seekable_parser(records, seekable_region=["1:1-100"])
        source = DummyFileFeatureSource(file_path="dummy.txt", adaptor=parser)

        for i, feature in enumerate(source.iter_all_features()):
            assert feature.start == records[i]
            assert parser.closed is False

    def test_closes_adaptor_once_generator_is_exhausted(self) -> None:
        """Test the adaptor is closed once all records have been consumed."""

        parser = _get_seekable_parser([1, 2, 3], seekable_region=["1:1-100"])
        source = DummyFileFeatureSource(file_path="dummy.txt", adaptor=parser)

        list(source.iter_all_features())

        assert parser.closed is True

    def test_adaptor_remains_open_outside_of_context_manager(self) -> None:
        """Test the adaptor remains open if the iterator is not used within a context-manager."""

        parser = _get_seekable_parser([1, 2, 3], seekable_region=["1:1-100"])
        source = DummyFileFeatureSource(file_path="dummy.txt", adaptor=parser)

        iterator = source.iter_all_features()
        first = next(iterator)
        assert first.start == 1
        assert parser.opened is True
        assert parser.closed is False


class TestFileFeatureSourceIterFeaturesByRegion:
    """Test iter_features_by_region."""

    def test_missing_region_arg_returns_nothing(self) -> None:
        """Test nothing is yielded when no region info is given."""

        parser = _get_seekable_parser([1, 2, 3], seekable_region=["1:1-100"])
        source = DummyFileFeatureSource(file_path="dummy.txt", adaptor=parser)

        assert not list(source.iter_features_by_region())

    def test_seek_failure_yields_nothing(self) -> None:
        """Test nothing is yielded when the adaptor can't seek to the region."""

        parser = _get_seekable_parser([1, 2, 3], seekable_region=["1:1-100"])
        source = DummyFileFeatureSource(file_path="dummy.txt", adaptor=parser)

        assert not list(source.iter_features_by_region("1:100-200"))

    def test_seek_success_yields_features(self) -> None:
        """Test features are yielded after a successful seek."""

        records = [10, 20]
        parser = _get_seekable_parser(records, seekable_region=["1:1-100"])
        source = DummyFileFeatureSource(file_path="dummy.txt", adaptor=parser)

        for i, feature in enumerate(source.iter_features_by_region("1:1-100")):
            assert feature.start == records[i]


class TestFileFeatureSourceClose:
    """Test close."""

    def test_close_is_a_no_op(self) -> None:
        """Test close does not raise and does not touch the adaptor.

        FileFeatureSource opens and closes a fresh adaptor session per call, so
        there is no persistent state for close() to release.
        """

        parser = _get_seekable_parser([1], seekable_region=["1:1-100"])
        source = DummyFileFeatureSource(file_path="dummy.txt", adaptor=parser)

        source.close()

        assert parser.closed is False


class TestFileFeatureSourceContextManager:
    """Test the inherited FeatureSource context manager support."""

    def test_enter_returns_self(self) -> None:
        """Test __enter__ returns the source instance itself."""

        parser = _get_seekable_parser([1], seekable_region=["1:1-100"])
        source = DummyFileFeatureSource(file_path="dummy.txt", adaptor=parser)

        with source as s:
            assert s is source

    def test_exit_calls_close(self) -> None:
        """Test __exit__ calls close on normal exit."""

        parser = _get_seekable_parser([1], seekable_region=["1:1-100"])
        source = DummyFileFeatureSource(file_path="dummy.txt", adaptor=parser)

        with source:
            pass

        assert source.closed is True
