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

# Use the VEP tool icon from https://www.ensembl.org/tools/vep rather than the shared Ensembl mark.
# The local CSS adapts the dark artwork for the dark navigation bar.
html_logo = "_static/vep-icon.svg"
html_favicon = "_static/vep-icon.svg"
html_theme_options["logo"]["text"] = ""
html_static_path = ["_static", "auto/static"]
html_css_files.append("vep.css")

# The entry-points table lives on the CLI reference landing page, which `make clidoc` generates with the
# table markers in place.
entrypoints_target_file = "auto/reference/cli/index.md"

# The coverage report lives in the Developer Guide. The shared configuration hides both sidebars for it by
# its old top-level name, so give the embedded reports the full page width under the new name too.
COVERAGE_PAGE = "developer_guide/coverage_report"
html_sidebars[COVERAGE_PAGE] = []
html_theme_options["secondary_sidebar_items"][COVERAGE_PAGE] = []
