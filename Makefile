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

coverage:
	mkdir -p docs/reports/coverage/{python,rust}
	uv run coverage run -m pytest
	uv run coverage html -d docs/reports/coverage/python/htmlcov
	uv run coverage xml -o docs/reports/coverage/python/coverage.xml
	uv run genbadge coverage -i docs/reports/coverage/python/coverage.xml -o docs/reports/coverage/python/coverage-badge.svg -n "Python"
	cp docs/reports/coverage/python/coverage-badge.svg docs/img/python-coverage-badge.svg
	cargo llvm-cov --ignore-filename-regex '(^|/)src/rust/lib\.rs$$' --html --output-dir docs/reports/coverage/rust/
	cargo llvm-cov report --cobertura --output-path docs/reports/coverage/rust/coverage.xml
	uv run genbadge coverage -i docs/reports/coverage/rust/coverage.xml -o docs/reports/coverage/rust/coverage-badge.svg -n "Rust"
	cp docs/reports/coverage/rust/coverage-badge.svg docs/img/rust-coverage-badge.svg
	rm docs/reports/coverage/python/coverage.xml docs/reports/coverage/rust/coverage.xml

apidoc:
	uv run sphinx-apidoc -o docs/reference/ src/python/ensembl --force --module-first \
	--no-toc --implicit-namespaces --templatedir docs/_templates/apidoc src/python/ensembl/vep/ensembl_vep_core.cpython-312-darwin.so
	rm -f docs/reference/ensembl.rst

docs: apidoc
	uv run sphinx-build -b html docs/ docs/_build/html

docs-clean:
	rm -rf docs/_build docs/reference docs/reports
