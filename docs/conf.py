# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Sphinx configuration for the ensembl-vep documentation.

All shared defaults live in :mod:`ensembl.utils.docs`; this file only sets project-specific values and any
local overrides. It doubles as the first consumer of the shared module, so every docs build exercises it.
"""

from pathlib import Path

import ensembl.utils
from ensembl.utils.docs import configure
import ensembl.vep

coverage_root = Path(__file__).parent / "reports"
configure(
    globals(),
    project="ensembl-vep",
    repo_url="https://github.com/Ensembl/ensembl-vep-py",
    release=ensembl.vep.__version__,
    docs_base_url="https://ensembl.github.io/ensembl-vep-py",
    coverage_root=coverage_root if coverage_root.exists() else None,
    include_entrypoints=True,
    add_pypi_icon=True,
    napoleon_attr_annotations=True,
    napoleon_use_rtype=True,
    autodoc_typehints="description",
    typehints_use_rtype=True,
    autodoc_default_options={"member-order": "alphabetical"},
)
