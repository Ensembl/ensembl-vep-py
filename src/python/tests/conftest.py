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
"""Shared pytest configuration."""

from dataclasses import dataclass
import os
from pathlib import Path
import re
import shutil
import subprocess
import tarfile
import urllib.request

import pytest

RUN_INTEGRATION_ENV = "ENSEMBL_VEP_RUN_INTEGRATION_TESTS"
NO_INSTALL_HTSLIB_ENV = "ENSEMBL_VEP_NO_INSTALL_HTSLIB"
HTSLIB_BIN_ENV = "ENSEMBL_VEP_HTSLIB_BIN"
HTSLIB_VERSION_ENV = "ENSEMBL_VEP_HTSLIB_VERSION"
DEFAULT_HTSLIB_VERSION = "1.23"


@dataclass(frozen=True)
class HtslibToolchain:
    """Paths to HTSlib command line tools used by integration tests."""

    bgzip: str
    tabix: str


def pytest_addoption(parser: pytest.Parser) -> None:
    """Add integration test command line options."""

    parser.addoption(
        "--run-integration",
        action="store_true",
        default=_truthy(os.environ.get(RUN_INTEGRATION_ENV)),
        help=f"run integration tests marked with @pytest.mark.integration; also set {RUN_INTEGRATION_ENV}=1",
    )
    parser.addoption(
        "--htslib-bin",
        action="store",
        default=os.environ.get(HTSLIB_BIN_ENV),
        help=f"directory containing bgzip and tabix; also set {HTSLIB_BIN_ENV}",
    )
    parser.addoption(
        "--htslib-version",
        action="store",
        default=os.environ.get(HTSLIB_VERSION_ENV),
        help=f"HTSlib version to use when installing tools; also set {HTSLIB_VERSION_ENV}",
    )
    parser.addoption(
        "--no-install-htslib",
        action="store_true",
        default=_truthy(os.environ.get(NO_INSTALL_HTSLIB_ENV)),
        help=f"skip automatic HTSlib source installation; also set {NO_INSTALL_HTSLIB_ENV}=1",
    )


def pytest_configure(config: pytest.Config) -> None:
    """Register custom pytest markers."""

    config.addinivalue_line(
        "markers",
        "integration: tests that require external tools such as bgzip and tabix",
    )


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    """Skip integration tests unless explicitly requested."""

    if config.getoption("--run-integration"):
        return

    skip_integration = pytest.mark.skip(reason="skipped: integration test requires --run-integration")
    for item in items:
        if "integration" in item.keywords:
            item.add_marker(skip_integration)


@pytest.fixture(scope="session")
def htslib_toolchain(
    pytestconfig: pytest.Config, tmp_path_factory: pytest.TempPathFactory
) -> HtslibToolchain:
    """Return bgzip and tabix paths for integration tests."""

    if not pytestconfig.getoption("--run-integration"):
        pytest.skip("skipped: integration test requires --run-integration")

    requested_version = pytestconfig.getoption("--htslib-version")
    htslib_bin = pytestconfig.getoption("--htslib-bin")
    toolchain = _find_htslib_toolchain(htslib_bin)

    if toolchain is not None and _tabix_version_matches(toolchain.tabix, requested_version):
        return toolchain

    if pytestconfig.getoption("--no-install-htslib"):
        pytest.skip("skipped: bgzip and tabix were not found and automatic HTSlib installation is disabled")

    version = requested_version or DEFAULT_HTSLIB_VERSION
    try:
        return _install_htslib(version, tmp_path_factory)
    except (OSError, subprocess.CalledProcessError, tarfile.TarError) as exc:
        pytest.skip(f"skipped: could not install HTSlib {version}: {exc}")

    raise RuntimeError("unreachable")


def _truthy(value: str | None) -> bool:
    """Return whether an environment value means true."""

    return value is not None and value.lower() in {"1", "true", "yes", "on"}


def _find_htslib_toolchain(htslib_bin: str | None) -> HtslibToolchain | None:
    """Find bgzip and tabix from an explicit directory or PATH."""

    bgzip = _find_executable("bgzip", htslib_bin)
    tabix = _find_executable("tabix", htslib_bin)
    if bgzip is None or tabix is None:
        return None

    return HtslibToolchain(bgzip=bgzip, tabix=tabix)


def _find_executable(name: str, htslib_bin: str | None) -> str | None:
    """Find an executable in htslib_bin or PATH."""

    if htslib_bin is not None:
        candidate = Path(htslib_bin) / name
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)

    return shutil.which(name)


def _tabix_version_matches(tabix: str, requested_version: str | None) -> bool:
    """Return whether tabix satisfies the requested HTSlib version."""

    if requested_version is None:
        return True

    actual_version = _tabix_version(tabix)
    if actual_version is None:
        return False

    return actual_version == requested_version or actual_version.startswith(f"{requested_version}.")


def _tabix_version(tabix: str) -> str | None:
    """Return the version reported by tabix."""

    output = subprocess.run([tabix, "--version"], check=False, capture_output=True, text=True)
    version_match = re.search(r"\b(\d+\.\d+(?:\.\d+)?)\b", output.stdout + output.stderr)
    if version_match is None:
        return None

    return version_match.group(1)


def _install_htslib(version: str, tmp_path_factory: pytest.TempPathFactory) -> HtslibToolchain:
    """Install HTSlib from a release tarball and return bgzip/tabix paths."""

    install_root = tmp_path_factory.mktemp(f"htslib-{version.replace('.', '-')}")
    prefix = install_root / "install"
    archive_path = install_root / f"htslib-{version}.tar.bz2"
    source_url = f"https://github.com/samtools/htslib/releases/download/{version}/htslib-{version}.tar.bz2"

    urllib.request.urlretrieve(source_url, archive_path)  # noqa: S310
    with tarfile.open(archive_path) as archive:
        archive.extractall(install_root, filter="data")

    source_dir = install_root / f"htslib-{version}"
    jobs = str(os.cpu_count() or 2)

    subprocess.run(["sh", "configure", f"--prefix={prefix}"], cwd=source_dir, check=True)
    subprocess.run(["make", "-j", jobs, "bgzip", "tabix"], cwd=source_dir, check=True)
    subprocess.run(["make", "install"], cwd=source_dir, check=True)

    toolchain = _find_htslib_toolchain(str(prefix / "bin"))
    if toolchain is None:
        raise OSError(f"HTSlib {version} installed without bgzip and tabix")

    return toolchain
