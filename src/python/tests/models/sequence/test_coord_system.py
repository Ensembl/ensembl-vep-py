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
"""Unit tests for CoordSystem class."""

import pytest

from ensembl.vep.models.sequence.coord_system import CoordSystem
from ensembl.vep.common import CoordSystemLevel


class TestCoordSystemInitialization:
    """Test CoordSystem initialization."""

    def test_initialization_without_keywords(self) -> None:
        """Test CoordSystem initialization without keywords raises expected error."""
        with pytest.raises(TypeError, match=r"missing \d+ required keyword-only arguments?"):
            _ = CoordSystem()  # type: ignore[call-arg]


class TestCoordSystemProperties:
    """Test property accessors."""

    @pytest.fixture
    def coord(self) -> CoordSystem:
        """Create a CoordSystem instance for testing."""
        return CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME)

    def test_core_properties(self, coord: CoordSystem) -> None:
        """Test core properties are set correctly."""
        assert coord.name == "chromosome"
        assert coord.level == CoordSystemLevel.CHROMOSOME
        assert coord.default is False
        assert coord.alias_to is None
        assert coord.rank is None

    def test_is_some_level(self, coord: CoordSystem) -> None:
        """Test is_chromosome_level property."""
        assert coord.is_chromosome_level is True
        assert coord.is_sequence_level is False
        assert coord.is_toplevel is False

    def test_is_some_level_immutable(self, coord: CoordSystem) -> None:
        """Test is_chromosome_level property is cached and immutable."""
        with pytest.raises(
            AttributeError, match="property 'is_chromosome_level' of 'CoordSystem' object has no setter"
        ):
            coord.is_chromosome_level = False  # type: ignore[misc]

        with pytest.raises(
            AttributeError, match="property 'is_sequence_level' of 'CoordSystem' object has no setter"
        ):
            coord.is_sequence_level = True  # type: ignore[misc]

        with pytest.raises(
            AttributeError, match="property 'is_toplevel' of 'CoordSystem' object has no setter"
        ):
            coord.is_toplevel = True  # type: ignore[misc]

    def test_is_some_level_update(self, coord: CoordSystem) -> None:
        """Test is_chromosome_level property update."""
        assert coord.is_chromosome_level is True
        assert coord.is_sequence_level is False
        assert coord.is_toplevel is False

        coord.level = CoordSystemLevel.CONTIG
        assert coord.is_chromosome_level is False
        assert coord.is_sequence_level is True
        assert coord.is_toplevel is False

        coord.level = CoordSystemLevel.TOPLEVEL
        assert coord.is_chromosome_level is False
        assert coord.is_sequence_level is False
        assert coord.is_toplevel is True

        coord.level = CoordSystemLevel.CHROMOSOME  # switch back to chromosome
        assert coord.is_chromosome_level is True
        assert coord.is_sequence_level is False
        assert coord.is_toplevel is False


class TestCoordSystemEquality:
    """Test equality comparison."""

    @pytest.fixture
    def coord(self) -> CoordSystem:
        """Create a CoordSystem instance for testing."""
        return CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME)

    def test_equal_coord_systems(self, coord: CoordSystem) -> None:
        """Test equal CoordSystem instances are considered equal."""
        coord2 = CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME)
        assert coord == coord2

    def test_different_name(self, coord: CoordSystem) -> None:
        """Test two CoordSystems with different names are not equal."""
        coord2 = CoordSystem(name="chr1", level=CoordSystemLevel.CHROMOSOME)
        assert coord != coord2

    def test_different_level(self, coord: CoordSystem) -> None:
        """Test two CoordSystems with different levels are not equal."""
        coord2 = CoordSystem(name="chromosome", level=CoordSystemLevel.CONTIG)
        assert coord != coord2

    def test_different_default(self, coord: CoordSystem) -> None:
        """Test two CoordSystems with different default values are still equal."""
        coord2 = CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME, default=True)
        assert coord == coord2

    def test_different_alias_to(self, coord: CoordSystem) -> None:
        """Test two CoordSystems with different alias_to values are still equal."""
        coord2 = CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME, alias_to="NC_000001.11")
        assert coord == coord2

    def test_different_rank(self, coord: CoordSystem) -> None:
        """Test two CoordSystems with different rank values are not equal."""
        coord2 = CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME, rank=5)
        assert coord != coord2

    def test_different_type(self, coord: CoordSystem) -> None:
        """Test comparing CoordSystem to different type returns False."""

        class AnotherCoord(CoordSystem):
            """CoordSystem subclass for testing."""

        coord2 = AnotherCoord(name="chromosome", level=CoordSystemLevel.CHROMOSOME)
        assert coord != coord2
