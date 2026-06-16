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

from typing import Iterator

import pytest
from attrs import define

from ensembl.vep.common import SourceType, Strand
from ensembl.vep.models import Feature
from ensembl.vep.sources.file.file_source import FileSource


class _FakeParser:
    """Minimal test double standing in for a real Parser/adaptor.

    Mimics the open()-as-context-manager / seek() / next() / records() /
    parse_file() contract that FileSource relies on, without needing a real
    parser implementation.
    """

    def __init__(self, records: list[int], seekable_region: str | None = None) -> None:
        self.records_source = records
        self.seekable_region = seekable_region
        self.opened = False
        self.closed = False
        self._cursor = 0

    def parse_file(self, file_path: str) -> list[int]:
        """Return all records in one shot, independent of open()/seek()."""
        return list(self.records_source)

    def open(self, file_path: str, mode: str = "auto") -> "_FakeParser":
        """Open the parser and reset the read cursor."""
        self.opened = True
        self._cursor = 0
        return self

    def __enter__(self) -> "_FakeParser":
        return self

    def __exit__(self, *exc_info: object) -> None:
        self.closed = True

    def seek(self, region: str) -> bool:
        """Succeed only when asked for the configured seekable region."""
        return region == self.seekable_region

    def next(self) -> int | None:
        """Return the next record, or None when exhausted."""
        if self._cursor >= len(self.records_source):
            return None
        value = self.records_source[self._cursor]
        self._cursor += 1
        return value

    def records(self, limit: int | None = None) -> list[int]:
        """Return all records from the current cursor position onward."""
        remaining = self.records_source[self._cursor :]
        self._cursor = len(self.records_source)
        return remaining


@define(kw_only=True, auto_attribs=True)
class _DummyFileSource(FileSource[_FakeParser, Feature]):
    """Concrete FileSource for exercising the abstract base in isolation."""

    def _records_to_features(self, records: list[int]) -> list[Feature]:
        return [Feature(start=value, end=value, strand=Strand.FORWARD, slice=None) for value in records]


class TestFileSourceType:
    """Test default attributes."""

    def test_source_type_is_file(self) -> None:
        """Test source_type defaults to SourceType.FILE."""
        source = _DummyFileSource(file_path="dummy.txt", adaptor=_FakeParser([1]))
        assert source.source_type == SourceType.FILE


class TestFileSourceGetAllFeatures:
    """Test get_all_features."""

    def test_returns_features_from_parse_file(self) -> None:
        """Test features are built from adaptor.parse_file()."""
        parser = _FakeParser([1, 2, 3])
        source = _DummyFileSource(file_path="dummy.txt", adaptor=parser)

        features = source.get_all_features()

        assert [f.start for f in features] == [1, 2, 3]

    def test_does_not_open_adaptor(self) -> None:
        """Test parse_file is a one-shot call that bypasses open()/seek()."""
        parser = _FakeParser([1, 2, 3])
        source = _DummyFileSource(file_path="dummy.txt", adaptor=parser)

        source.get_all_features()

        assert parser.opened is False


class TestFileSourceGetFeaturesByRegion:
    """Test get_features_by_region."""

    def test_returns_empty_list_when_region_args_missing(self) -> None:
        """Test an empty list is returned when no region info is given."""
        parser = _FakeParser([1, 2, 3])
        source = _DummyFileSource(file_path="dummy.txt", adaptor=parser)

        assert source.get_features_by_region() == []

    def test_returns_empty_list_when_region_and_seq_region_name_both_given(self) -> None:
        """Test an empty list is returned when region and seq_region_name conflict."""
        parser = _FakeParser([1, 2, 3])
        source = _DummyFileSource(file_path="dummy.txt", adaptor=parser)

        assert source.get_features_by_region("1:1-100", seq_region_name="1") == []

    def test_returns_empty_list_when_seek_fails(self) -> None:
        """Test an empty list is returned when the adaptor can't seek to the region."""
        parser = _FakeParser([1, 2, 3], seekable_region=None)
        source = _DummyFileSource(file_path="dummy.txt", adaptor=parser)

        assert source.get_features_by_region("1:1-100") == []

    def test_returns_features_when_seek_succeeds(self) -> None:
        """Test features are built from adaptor.records() after a successful seek."""
        parser = _FakeParser([10, 20], seekable_region="1:1-100")
        source = _DummyFileSource(file_path="dummy.txt", adaptor=parser)

        features = source.get_features_by_region("1:1-100")

        assert [f.start for f in features] == [10, 20]

    def test_builds_region_from_seq_region_name_start_end(self) -> None:
        """Test seq_region_name/start/end are combined into a region string."""
        parser = _FakeParser([10], seekable_region="1:5-15")
        source = _DummyFileSource(file_path="dummy.txt", adaptor=parser)

        features = source.get_features_by_region(seq_region_name="1", start=5, end=15)

        assert [f.start for f in features] == [10]

    def test_opens_and_closes_adaptor(self) -> None:
        """Test the adaptor is opened and closed around the region query."""
        parser = _FakeParser([10, 20], seekable_region="1:1-100")
        source = _DummyFileSource(file_path="dummy.txt", adaptor=parser)

        source.get_features_by_region("1:1-100")

        assert parser.opened is True
        assert parser.closed is True

    def test_closes_adaptor_even_when_seek_fails(self) -> None:
        """Test the adaptor is still closed when seek() returns False."""
        parser = _FakeParser([1, 2, 3], seekable_region=None)
        source = _DummyFileSource(file_path="dummy.txt", adaptor=parser)

        source.get_features_by_region("1:1-100")

        assert parser.closed is True


class TestFileSourceIterAllFeatures:
    """Test iter_all_features."""

    def test_returns_an_iterator(self) -> None:
        """Test iter_all_features returns a lazy iterator, not a list."""
        parser = _FakeParser([1, 2, 3])
        source = _DummyFileSource(file_path="dummy.txt", adaptor=parser)

        result = source.iter_all_features()

        assert isinstance(result, Iterator)
        assert parser.opened is False  # nothing happens until consumed

    def test_yields_features_one_at_a_time(self) -> None:
        """Test all features are yielded in order."""
        parser = _FakeParser([1, 2, 3])
        source = _DummyFileSource(file_path="dummy.txt", adaptor=parser)

        features = list(source.iter_all_features())

        assert [f.start for f in features] == [1, 2, 3]

    def test_opens_adaptor_before_first_yield(self) -> None:
        """Test the adaptor is opened once, lazily, on first consumption."""
        parser = _FakeParser([1, 2, 3])
        source = _DummyFileSource(file_path="dummy.txt", adaptor=parser)

        iterator = source.iter_all_features()
        first = next(iterator)

        assert first.start == 1
        assert parser.opened is True

    def test_keeps_adaptor_open_across_multiple_next_calls(self) -> None:
        """Test the adaptor session spans the whole generator, not just one record."""
        parser = _FakeParser([1, 2, 3])
        source = _DummyFileSource(file_path="dummy.txt", adaptor=parser)

        iterator = source.iter_all_features()
        next(iterator)

        assert parser.closed is False

    def test_closes_adaptor_once_generator_is_exhausted(self) -> None:
        """Test the adaptor is closed once all records have been consumed."""
        parser = _FakeParser([1, 2, 3])
        source = _DummyFileSource(file_path="dummy.txt", adaptor=parser)

        list(source.iter_all_features())

        assert parser.closed is True


class TestFileSourceIterFeaturesByRegion:
    """Test iter_features_by_region."""

    def test_yields_nothing_when_region_args_missing(self) -> None:
        """Test nothing is yielded when no region info is given."""
        parser = _FakeParser([1, 2, 3])
        source = _DummyFileSource(file_path="dummy.txt", adaptor=parser)

        assert list(source.iter_features_by_region()) == []

    def test_yields_nothing_when_seek_fails(self) -> None:
        """Test nothing is yielded when the adaptor can't seek to the region."""
        parser = _FakeParser([1, 2, 3], seekable_region=None)
        source = _DummyFileSource(file_path="dummy.txt", adaptor=parser)

        assert list(source.iter_features_by_region("1:1-100")) == []

    def test_yields_features_when_seek_succeeds(self) -> None:
        """Test features are yielded after a successful seek."""
        parser = _FakeParser([10, 20], seekable_region="1:1-100")
        source = _DummyFileSource(file_path="dummy.txt", adaptor=parser)

        features = list(source.iter_features_by_region("1:1-100"))

        assert [f.start for f in features] == [10, 20]


class TestFileSourceClose:
    """Test close."""

    def test_close_is_a_no_op(self) -> None:
        """Test close does not raise and does not touch the adaptor.

        FileSource opens/closes a fresh adaptor session per call, so there is
        no persistent state for close() to release.
        """
        parser = _FakeParser([1])
        source = _DummyFileSource(file_path="dummy.txt", adaptor=parser)

        source.close()

        assert parser.closed is False


class TestFileSourceContextManager:
    """Test the inherited FeatureSource context manager support."""

    def test_enter_returns_self(self) -> None:
        """Test __enter__ returns the source instance itself."""
        source = _DummyFileSource(file_path="dummy.txt", adaptor=_FakeParser([1]))

        with source as entered:
            assert entered is source

    def test_exit_calls_close(self) -> None:
        """Test __exit__ calls close on normal exit."""
        source = _DummyFileSource(file_path="dummy.txt", adaptor=_FakeParser([1]))
        calls: list[bool] = []
        source.close = lambda: calls.append(True)  # type: ignore[method-assign]

        with source:
            pass

        assert calls == [True]
