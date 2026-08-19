# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

.PHONY: install-* apidoc docs docs-clean

install-uv:
	curl -LsSf https://astral.sh/uv/install.sh | sh

install-cargo:
	curl https://sh.rustup.rs -sSf | sh

docs-clean:
	rm -rf docs/_build docs/reference docs/reports

apidoc:
	uv run sphinx-apidoc -o docs/reference/ src/python/ensembl --force --module-first \
	--no-toc --implicit-namespaces --templatedir docs/_templates/apidoc
	rm -f docs/reference/ensembl.rst

docs: docs-clean apidoc
	uv run sphinx-build -b html docs/ docs/_build/html

coverage:
	uv run coverage run -m pytest
	uv run coverage html -d docs/reports/htmlcov
	uv run coverage xml -o docs/reports/coverage.xml
	uv run genbadge coverage -i docs/reports/coverage.xml -o docs/reports/htmlcov/coverage-badge.svg -n "Python"
	rm docs/reports/coverage.xml
