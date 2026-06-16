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
"""Unit tests for DBSource."""

import pytest
from attrs import define

from ensembl.vep.common import SourceType, Strand
from ensembl.vep.models import Feature
from ensembl.vep.models.sequence.slice import Slice
from ensembl.vep.sources.db.db_source import DBSource


class _FakeDBAdaptor:
    """Minimal test double standing in for a real DB adaptor."""

    def __init__(self, rows: list[int]) -> None:
        self.rows = rows
        self.fetch_by_region_calls: list[dict] = []

    def fetch_all(self) -> list[int]:
        """Return all rows, independent of region filtering."""
        return list(self.rows)

    def fetch_by_region(
        self, *, seq_region_name: str | None, start: int | None, end: int | None
    ) -> list[int]:
        """Record the call and return all configured rows."""
        self.fetch_by_region_calls.append(
            {"seq_region_name": seq_region_name, "start": start, "end": end}
        )
        return list(self.rows)


@define(kw_only=True, auto_attribs=True)
class _DummyDBSource(DBSource[_FakeDBAdaptor, Feature]):
    """Concrete DBSource for exercising the abstract base in isolation."""

    def _rows_to_features(self, rows: list[int]) -> list[Feature]:
        return [Feature(start=row, end=row, strand=Strand.FORWARD, slice=None) for row in rows]


class TestDBSourceType:
    """Test default attributes."""

    def test_source_type_is_db(self) -> None:
        """Test source_type defaults to SourceType.DB."""
        source = _DummyDBSource(adaptor=_FakeDBAdaptor([1]))
        assert source.source_type == SourceType.DB


class TestDBSourceGetAllFeatures:
    """Test get_all_features."""

    def test_returns_features_from_fetch_all(self) -> None:
        """Test features are built from adaptor.fetch_all()."""
        adaptor = _FakeDBAdaptor([1, 2, 3])
        source = _DummyDBSource(adaptor=adaptor)

        features = source.get_all_features()

        assert [f.start for f in features] == [1, 2, 3]


class TestDBSourceGetFeaturesByRegion:
    """Test get_features_by_region."""

    def test_returns_features_from_fetch_by_region(self) -> None:
        """Test features are built from adaptor.fetch_by_region()."""
        adaptor = _FakeDBAdaptor([10, 20])
        source = _DummyDBSource(adaptor=adaptor)

        features = source.get_features_by_region("1", start=100, end=200)

        assert [f.start for f in features] == [10, 20]

    def test_passes_seq_region_name_start_end_to_adaptor(self) -> None:
        """Test region arguments are forwarded to the adaptor unchanged."""
        adaptor = _FakeDBAdaptor([10])
        source = _DummyDBSource(adaptor=adaptor)

        source.get_features_by_region("1", start=100, end=200)

        assert adaptor.fetch_by_region_calls == [{"seq_region_name": "1", "start": 100, "end": 200}]

    def test_accepts_slice_as_seq_region_name(self) -> None:
        """Test a Slice's seq_region_name is extracted before calling the adaptor."""
        adaptor = _FakeDBAdaptor([10])
        source = _DummyDBSource(adaptor=adaptor)
        genomic_slice = Slice(
            coord_system="chromosome",
            seq_region_name="1",
            start=1,
            end=1000,
            strand=Strand.FORWARD,
        )

        source.get_features_by_region(genomic_slice, start=100, end=200)

        assert adaptor.fetch_by_region_calls == [{"seq_region_name": "1", "start": 100, "end": 200}]


class TestDBSourceIteration:
    """Test that streaming iteration is not yet supported for DB sources."""

    def test_iter_all_features_raises_not_implemented(self) -> None:
        """Test iter_all_features raises NotImplementedError."""
        source = _DummyDBSource(adaptor=_FakeDBAdaptor([1]))
        with pytest.raises(NotImplementedError):
            source.iter_all_features()

    def test_iter_features_by_region_raises_not_implemented(self) -> None:
        """Test iter_features_by_region raises NotImplementedError."""
        source = _DummyDBSource(adaptor=_FakeDBAdaptor([1]))
        with pytest.raises(NotImplementedError):
            source.iter_features_by_region("1", start=1, end=10)


class TestDBSourceClose:
    """Test close."""

    def test_close_is_a_no_op(self) -> None:
        """Test close does not raise."""
        source = _DummyDBSource(adaptor=_FakeDBAdaptor([1]))
        source.close()


class TestDBSourceContextManager:
    """Test the inherited FeatureSource context manager support."""

    def test_enter_returns_self(self) -> None:
        """Test __enter__ returns the source instance itself."""
        source = _DummyDBSource(adaptor=_FakeDBAdaptor([1]))

        with source as entered:
            assert entered is source

    def test_exit_calls_close(self) -> None:
        """Test __exit__ calls close on normal exit."""
        source = _DummyDBSource(adaptor=_FakeDBAdaptor([1]))
        calls: list[bool] = []
        source.close = lambda: calls.append(True)  # type: ignore[method-assign]

        with source:
            pass

        assert calls == [True]
