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
"""Interactive variant progress reporting."""

from __future__ import annotations

import sys
from time import monotonic
from typing import TextIO

from attrs import define, field

__all__ = ["VariantProgress"]


def _format_duration(seconds: int) -> str:
    """Format a duration as hours, minutes, and seconds."""
    hours, remainder = divmod(seconds, 3600)
    minutes, remaining_seconds = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{remaining_seconds:02d}"


@define(kw_only=True, auto_attribs=True)
class VariantProgress:
    """Progress renderer for completed input variants.

    Attributes:
        total (int | None): Number of input variant records, or None when the input declares no
            count. Without one there is nothing to measure against, so the line shows what has
            been annotated so far and leaves out the bar, the percentage and the ETA.
        enabled (bool): Whether updates should be rendered.
        stream (TextIO): Interactive output stream.
        width (int): Character width of the progress bar.
    """

    total: int | None
    enabled: bool = True
    stream: TextIO = field(factory=lambda: sys.stderr)
    width: int = 24
    _completed: int = field(default=0, init=False)
    _started_at: float = field(init=False)
    _rendered_width: int = field(default=0, init=False)

    def __attrs_post_init__(self) -> None:
        if self.total is not None and self.total < 0:
            raise ValueError("progress total must not be negative")
        if self.width < 1:
            raise ValueError("progress width must be positive")
        self._started_at = monotonic()

    def __enter__(self) -> VariantProgress:
        if self.enabled:
            self._render()
        return self

    def __exit__(self, *_error: object) -> None:
        if self.enabled:
            self.stream.write("\n")
            self.stream.flush()

    def advance(self, variants: int) -> None:
        """Advance by a non-negative number of completed variants."""
        if variants < 0:
            raise ValueError("progress increment must not be negative")
        completed = self._completed + variants
        self._completed = completed if self.total is None else min(self.total, completed)
        if self.enabled:
            self._render()

    def _render(self) -> None:
        """Redraw the current variant count."""
        elapsed_seconds = max(0.0, monotonic() - self._started_at)
        rate = self._completed / elapsed_seconds if elapsed_seconds > 0 else 0.0

        if self.total is None:
            line = self._unsized_line(rate, elapsed_seconds)
        else:
            line = self._sized_line(self.total, rate, elapsed_seconds)

        padding = " " * max(0, self._rendered_width - len(line))
        self.stream.write(f"\r{line}{padding}")
        self.stream.flush()
        self._rendered_width = len(line)

    def _sized_line(self, total: int, rate: float, elapsed_seconds: float) -> str:
        """Return the progress line for an input whose record count is known."""
        percentage = 100 if total == 0 else int(100 * self._completed / total)
        filled = self.width if total == 0 else int(self.width * self._completed / total + 0.5)
        progress_bar = "█" * filled + "░" * (self.width - filled)

        if total in (0, self._completed):
            eta = "00:00:00"
            elapsed = "00:00:00"
        elif rate == 0:
            eta = "--:--:--"
            elapsed = "--:--:--"
        else:
            remaining_seconds = int((total - self._completed) / rate + 0.5)
            eta = _format_duration(remaining_seconds)
            elapsed = _format_duration(int(elapsed_seconds + 0.5))

        return (
            f"Annotating variants: [{progress_bar}] {percentage}% | "
            f"{self._completed:,}/{total:,} | {int(rate + 0.5):,} variants/s | {elapsed}>{eta}"
        )

    def _unsized_line(self, rate: float, elapsed_seconds: float) -> str:
        """Return the progress line for an input that declares no record count.

        No total means no fraction of it, so the bar, the percentage and the ETA are all left
        out rather than drawn against a guess.
        """
        elapsed = _format_duration(int(elapsed_seconds + 0.5))

        return f"Annotating variants: {self._completed:,} | {int(rate + 0.5):,} variants/s | {elapsed}"
