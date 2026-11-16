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
    """Render completed variants against the input record count.

    Attributes:
        total (int): Number of input variant records.
        enabled (bool): Whether updates should be rendered.
        stream (TextIO): Interactive output stream.
        width (int): Character width of the progress bar.
    """

    total: int
    enabled: bool = True
    stream: TextIO = field(factory=lambda: sys.stderr)
    width: int = 24
    _completed: int = field(default=0, init=False)
    _started_at: float = field(init=False)
    _rendered_width: int = field(default=0, init=False)

    def __attrs_post_init__(self) -> None:
        if self.total < 0:
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
        self._completed = min(self.total, self._completed + variants)
        if self.enabled:
            self._render()

    def _render(self) -> None:
        """Redraw the current variant count."""
        percentage = 100 if self.total == 0 else int(100 * self._completed / self.total)
        filled = self.width if self.total == 0 else int(self.width * self._completed / self.total + 0.5)
        progress_bar = "█" * filled + "░" * (self.width - filled)

        elapsed_seconds = max(0.0, monotonic() - self._started_at)
        rate = self._completed / elapsed_seconds if elapsed_seconds > 0 else 0.0
        rounded_rate = int(rate + 0.5)
        if self.total in (0, self._completed):
            eta = "00:00:00"
            elapsed = "00:00:00"
        elif rate == 0:
            eta = "--:--:--"
            elapsed = "--:--:--"
        else:
            remaining_seconds = int((self.total - self._completed) / rate + 0.5)
            eta = _format_duration(remaining_seconds)
            elapsed = _format_duration(int(elapsed_seconds + 0.5))

        line = (
            f"Annotating variants: [{progress_bar}] {percentage}% | "
            f"{self._completed:,}/{self.total:,} | {rounded_rate:,} variants/s | {elapsed}>{eta}"
        )
        padding = " " * max(0, self._rendered_width - len(line))
        self.stream.write(f"\r{line}{padding}")
        self.stream.flush()
        self._rendered_width = len(line)
