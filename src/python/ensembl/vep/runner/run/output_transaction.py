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
"""Atomic output transactions for the default annotation runner."""

from __future__ import annotations

import logging
import os
import secrets
import stat
from pathlib import Path

from attrs import define, field

__all__ = ["OutputTransaction"]

logger = logging.getLogger(__name__)


def _temporary_sibling(destination: Path, suffix: str) -> Path:
    """Create a unique sibling with normal output permissions."""
    try:
        destination_mode = stat.S_IMODE(destination.stat().st_mode)
    except FileNotFoundError:
        destination_mode = None

    for _attempt in range(100):
        temporary = destination.with_name(f".{destination.name}.{secrets.token_hex(8)}.{suffix}")
        try:
            descriptor = os.open(
                temporary,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL,
                0o666,
            )
        except FileExistsError:
            continue
        os.close(descriptor)
        if destination_mode is not None:
            try:
                temporary.chmod(destination_mode)
            except BaseException:
                temporary.unlink(missing_ok=True)
                raise
        return temporary
    raise FileExistsError(f"could not create a temporary file beside {destination}")


@define(kw_only=True, auto_attribs=True)
class OutputTransaction:
    """Atomic transaction for files produced by one run.

    Attributes:
        replace_existing (bool): Whether commit may replace existing destinations.
    """

    replace_existing: bool
    _staged_files: list[tuple[Path, Path]] = field(factory=list, init=False)

    def stage(self, destination: Path) -> Path:
        """Create and track a temporary file beside ``destination``."""
        temporary = _temporary_sibling(destination, "tmp")
        self._staged_files.append((temporary, destination))
        return temporary

    def commit(self) -> None:
        """Commit every staged file or restore the previous generation."""
        if not self.replace_existing:
            self._commit_new_files()
            return
        self._replace_files()

    def cleanup(self) -> None:
        """Remove temporary files that were not consumed by commit."""
        for temporary, _destination in self._staged_files:
            temporary.unlink(missing_ok=True)

    def _commit_new_files(self) -> None:
        """Commit files without replacing destinations created during the run."""
        committed: list[Path] = []
        try:
            for temporary, destination in self._staged_files:
                os.link(temporary, destination)
                committed.append(destination)
        except BaseException:
            for destination in reversed(committed):
                destination.unlink(missing_ok=True)
            raise

    def _replace_files(self) -> None:
        """Replace destinations and restore all predecessors if any replacement fails."""
        backups: dict[Path, Path] = {}
        committed: set[Path] = set()
        try:
            for temporary, destination in self._staged_files:
                if destination.exists():
                    backup = _temporary_sibling(destination, "backup")
                    try:
                        destination.replace(backup)
                    except BaseException:
                        backup.unlink(missing_ok=True)
                        raise
                    backups[destination] = backup
                temporary.replace(destination)
                committed.add(destination)
        except BaseException as commit_failure:
            try:
                self._rollback(backups, committed)
            except OSError as rollback_failure:
                raise rollback_failure from commit_failure
            raise

        for backup in backups.values():
            try:
                backup.unlink(missing_ok=True)
            except OSError as failure:
                logger.warning("could not remove output transaction backup %s: %s", backup, failure)

    def _rollback(self, backups: dict[Path, Path], committed: set[Path]) -> None:
        """Restore destinations changed by an incomplete commit."""
        failures: list[OSError] = []
        for _temporary, destination in reversed(self._staged_files):
            if destination in committed:
                try:
                    destination.unlink(missing_ok=True)
                except OSError as failure:
                    failures.append(failure)
            backup = backups.get(destination)
            if backup is not None:
                try:
                    backup.replace(destination)
                except OSError as failure:
                    failures.append(failure)
        if failures:
            raise OSError(f"could not restore output transaction: {failures[0]}")
