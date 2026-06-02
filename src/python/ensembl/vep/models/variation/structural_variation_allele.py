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
"""StructuralVariationAllele class."""

from __future__ import annotations

import logging
from attrs import define

from ensembl.vep.models.variation import VariationAllele

logger = logging.getLogger(__name__)


@define(kw_only=True, auto_attribs=True)
class StructuralVariationAllele(VariationAllele):
    """Represents a specific allele of a structural variation."""

    def minimise(self, force: bool = False) -> None:
        """Minimises the variation allele by computing the minimal representation of the variant.

        Structural variants are not minimised, so this method simply returns the current representation
        of the structural variant.
        """
        msg = "Structural Variation alleles are currently not minimised. "
        if force:
            msg += "force=True doesn't change this behaviour. "
        msg += "Returning the original representation"
        logger.warning(msg)

        original_locus = self._original_locus or self._locus

        setattr(self, "_original_locus", self._original_locus or original_locus)
        setattr(self, "_minimised_locus", self._locus)
        setattr(self, "_minimised", True)

    def shift(self, force: bool = False) -> None:
        """Shifts the variation allele to the leftmost position in the sequence.

        Currently, this method simply returns the current representation of the structural variant.
        """
        msg = "Structural Variation alleles are currently not shifted. "
        if force:
            msg += "force=True doesn't change this behaviour. "
        msg += "Returning the original representation"
        logger.warning(msg)

        original_locus = self._original_locus or self._locus

        setattr(self, "_original_locus", self._original_locus or original_locus)
        setattr(self, "_shifted_locus", self._locus)
        setattr(self, "_shifted", True)
