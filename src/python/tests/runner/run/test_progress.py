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
"""Tests for variant progress reporting."""

from __future__ import annotations

from io import StringIO
from unittest.mock import Mock

import pytest

from ensembl.vep.runner.run import progress as progress_module
from ensembl.vep.runner.run.progress import VariantProgress


class TestVariantProgress:
    """Interactive progress rendering."""

    def test_renders_counts_rate_and_eta(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Show exact counts, average rate, and estimated completion time."""
        stream = StringIO()
        clock = Mock(side_effect=[0.0, 0.0, 3240 / 1847])
        monkeypatch.setattr(progress_module, "monotonic", clock)

        with VariantProgress(total=5220, stream=stream) as progress:
            progress.advance(3240)

        rendered = stream.getvalue()
        assert (
            "Annotating variants: [███████████████░░░░░░░░░] "
            "62% | 3,240/5,220 | 1,847 variants/s | 00:00:02>00:00:01"
        ) in rendered
        assert rendered.endswith("\n")

    def test_completed_progress_has_no_remaining_time(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Finish at one hundred percent with a zero ETA."""
        stream = StringIO()
        clock = Mock(side_effect=[0.0, 0.0, 2.0])
        monkeypatch.setattr(progress_module, "monotonic", clock)

        with VariantProgress(total=4, stream=stream, width=4) as progress:
            progress.advance(4)

        assert "[████] 100% | 4/4 | 2 variants/s | 00:00:00>00:00:00" in stream.getvalue()

    def test_disabled_progress_writes_nothing(self) -> None:
        """Keep redirected, quiet, and explicitly disabled runs silent."""
        stream = StringIO()

        with VariantProgress(total=4, enabled=False, stream=stream) as progress:
            progress.advance(4)

        assert stream.getvalue() == ""

    def test_empty_input_is_complete(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Represent an empty input without dividing by zero."""
        stream = StringIO()
        monkeypatch.setattr(progress_module, "monotonic", Mock(return_value=0.0))

        with VariantProgress(total=0, stream=stream):
            pass

        assert "100% | 0/0 | 0 variants/s | 00:00:00>00:00:00" in stream.getvalue()

    def test_pending_progress_has_no_eta(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Leave the ETA unknown until a measurable rate is available."""
        stream = StringIO()
        monkeypatch.setattr(progress_module, "monotonic", Mock(return_value=0.0))

        with VariantProgress(total=4, stream=stream):
            pass

        assert "0% | 0/4 | 0 variants/s | --:--:-->--:--:--" in stream.getvalue()

    def test_an_unknown_total_renders_counts_without_a_bar(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Report what has been annotated when the input declares no record count."""
        stream = StringIO()
        clock = Mock(side_effect=[0.0, 0.0, 2.0])
        monkeypatch.setattr(progress_module, "monotonic", clock)

        with VariantProgress(total=None, stream=stream) as progress:
            progress.advance(3000)

        rendered = stream.getvalue()
        assert "Annotating variants: 3,000 | 1,500 variants/s | 00:00:02" in rendered
        assert "%" not in rendered
        assert "[" not in rendered

    def test_an_unknown_total_does_not_clamp_the_count(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Keep counting past any bound, there being none to stop at."""
        stream = StringIO()
        monkeypatch.setattr(progress_module, "monotonic", Mock(return_value=0.0))

        with VariantProgress(total=None, stream=stream) as progress:
            progress.advance(2)
            progress.advance(3)

        assert "Annotating variants: 5 | " in stream.getvalue()

    @pytest.mark.parametrize(
        ("total", "width", "message"),
        [
            pytest.param(-1, 24, "total", id="negative-total"),
            pytest.param(1, 0, "width", id="zero-width"),
        ],
    )
    def test_rejects_invalid_dimensions(self, total: int, width: int, message: str) -> None:
        """Reject progress dimensions that cannot be rendered."""
        with pytest.raises(ValueError, match=message):
            VariantProgress(total=total, width=width)

    def test_rejects_a_negative_increment(self) -> None:
        """Prevent completed progress from moving backwards."""
        progress = VariantProgress(total=1, enabled=False)

        with pytest.raises(ValueError, match="increment"):
            progress.advance(-1)
