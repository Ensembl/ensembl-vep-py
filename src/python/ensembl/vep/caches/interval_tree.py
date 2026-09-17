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
"""IntervalTree class."""

import logging
from threading import RLock

from typing import Generic, TypeVar

from ailist import AIList
from attr import define, field

from ensembl.vep.models import Feature

logger = logging.getLogger(__name__)

F = TypeVar("F", bound=Feature)


@define(kw_only=True, auto_attribs=True)
class IntervalTree(Generic[F]):
    """
    Provides methods to hold intervals in a data structure for calculating overlaps with other intervals.

    Attributes:
        features: A list of Feature objects to add
    """

    features: list[F] = field(factory=list)
    _intervals: AIList = field(factory=AIList, init=False)
    _lock: RLock = field(factory=RLock, init=False, repr=False)

    def __attrs_post_init__(self) -> None:
        """Populate interval index for features supplied at initialization."""
        self._add_features_to_intervals(self.features, start_idx=0)

    def add_features(self, features: list[F]) -> None:
        """Add the features to the IntervalTree object.

        Args:
            features: A list of Feature objects to add.

        Returns
            None

        Raises:
            ValueError: If the coordinates are invalid.
        """

        self._add_features_to_intervals(features, start_idx=len(self.features))
        self.features.extend(features)

    def _add_features_to_intervals(self, features: list[F], start_idx: int) -> None:
        """For each feature, check the coordinates and then add to the ``features`` and
        ``_intervals`` attributes.

        Args:
            features: A list of Feature objects to add.
            start_idx: The next index to add the first interval under

        Returns
            None

        Raises:
            Warning: If the coordinates are invalid.
        """
        for idx, feature in enumerate(features):
            try:
                start, end = self._check_interval(feature.start, feature.end)
                self._intervals.add(start, end + 1, start_idx + idx)
            except ValueError:
                logger.warning("Feature has invalid coordinates")

    def get_overlapping_features(self, start: int, end: int) -> list[F]:
        """Query the IntervalTree object with an interval and return overlapping Features.

        Args:
            start: start coordinate of the interval
            end: end coordinate of the interval

        Returns
            A list of Feature objects that overlap the interval

        Raises:
            ValueError: If the coordinates are invalid.
        """
        with self._lock:
            start, end = self._check_interval(start, end)
            return [self.features[idx] for idx in self._intervals.intersect_ids(start, end + 1)]

    def get_overlapping_features_by_feature(self, feature: Feature) -> list[F]:
        """Query the IntervalTree object with a Feature object and return overlapping Features.

        Args:
            feature: Feature object to check for overlaps

        Returns
            A list of Feature objects that overlap the supplied Feature

        Raises:
            ValueError: If the coordinates are invalid.
        """

        return self.get_overlapping_features(start=feature.start, end=feature.end)

    def get_nearest_features(self, start: int, end: int, k: int = 10) -> list[F]:
        """Query the IntervalTree object with an interval and return the ``k`` nearest Features.

        Args:
            start: start coordinate of the interval
            end: end coordinate of the interval

        Returns
            A list of the ``k`` nearest Feature objects to the supplied interval

        Raises:
            ValueError: If the coordinates are invalid.
        """
        with self._lock:
            start, end = self._check_interval(start, end)
            closest = self._intervals.closest(start, end + 1, k=k)
            print(f"Closest ids: {closest.ids}")
            return [self.features[idx] for idx in closest.ids]

    def get_nearest_features_by_feature(self, feature: Feature, k: int = 10) -> list[F]:
        """Query the IntervalTree object with a Feature object and return the nearest ``k`` Features.

        Args:
            feature: Feature object to check for overlaps

        Returns
            A list of the ``k`` nearest Feature objects to the supplied Feature

        Raises:
            ValueError: If the coordinates are invalid.
        """

        return self.get_nearest_features(start=feature.start, end=feature.end, k=k)

    def _check_interval(self, start: int, end: int) -> tuple[int, int]:
        """Helper method to check start and end values of an interval.

        Args:
            start: start coordinate of the interval
            end: end coordinate of the interval

        Returns
            None

        Raises:
            ValueError: If the coordinates are invalid.
        """

        if start < 1:
            raise ValueError("start must be 1-based: cannot be less than 1")
        if start > end:
            start, end = end, start
        return start, end
