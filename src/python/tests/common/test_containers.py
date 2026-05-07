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
"""Unit tests for common containers."""

import pytest

from ensembl.vep.common.containers import Locus, Attribute, ExternalDB, ExternalReference

# Define test cases: (class, same_kwargs, differing_kwargs)
CONTAINER_TEST_CASES = [
    (
        Locus,
        {"start": 100, "end": 200, "allele_string": "A/T"},
        {"start": 100, "end": 200, "allele_string": "A/G"},  # one field differs
    ),
    (
        Attribute,
        {"name": "gene_symbol", "value": "BRCA1", "code": "symbol", "description": "Gene symbol"},
        {"name": "gene_symbol", "value": "TP53", "code": "symbol", "description": "Gene symbol"},
    ),
    (
        ExternalDB,
        {
            "id": "NCBI",
            "name": "NCBI",
            "description": "National Center for Biotechnology Information",
            "url": "https://www.ncbi.nlm.nih.gov/",
            "release": "2024-06",
        },
        {
            "id": "NCBI",
            "name": "NCBI",
            "description": "National Center for Biotechnology Information",
            "url": "https://www.ncbi.nlm.nih.gov/",
            "release": "2023-12",
        },
    ),
    (
        ExternalReference,
        {
            "accession_id": "HGNC:1100",
            "name": "HGNC",
            "description": "HUGO Gene Nomenclature Committee",
            "assignment_method": "DIRECT",
            "url": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:1100",
        },
        {
            "accession_id": "HGNC:1101",
            "name": "HGNC",
            "description": "HUGO Gene Nomenclature Committee",
            "assignment_method": "DIRECT",
            "url": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:1101",
        },
    ),
]


@pytest.mark.parametrize("cls, kwargs, differing_kwargs", CONTAINER_TEST_CASES)
class TestContainers:  # type: ignore[pytest-classes]
    """General tests for all frozen container classes in ensembl.vep.common.containers."""

    def test_initialization(self, cls, kwargs, differing_kwargs) -> None:  # type: ignore[no-untyped-def]  # pylint: disable=unused-argument
        """Test that the class can be initialized with the correct attributes."""
        obj = cls(**kwargs)
        for key, value in kwargs.items():
            assert getattr(obj, key) == value

    def test_immutability(self, cls, kwargs, differing_kwargs) -> None:  # type: ignore[no-untyped-def]  # pylint: disable=unused-argument
        """Test that attributes cannot be modified after initialization."""
        obj = cls(**kwargs)
        with pytest.raises((AttributeError, TypeError)):
            setattr(obj, next(iter(kwargs)), "mutated")

    def test_equality(self, cls, kwargs, differing_kwargs) -> None:  # type: ignore[no-untyped-def]  # pylint: disable=unused-argument
        """Test that two instances with the same attributes are equal."""
        assert cls(**kwargs) == cls(**kwargs)

    def test_inequality(self, cls, kwargs, differing_kwargs) -> None:  # type: ignore[no-untyped-def]
        """Test that two instances with different attributes are not equal."""
        assert cls(**kwargs) != cls(**differing_kwargs)
