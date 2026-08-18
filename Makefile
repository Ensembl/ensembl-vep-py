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

build:
	uv build

.PHONY: install-*

install-uv:
	curl -LsSf https://astral.sh/uv/install.sh | sh

install-cargo:
	curl https://sh.rustup.rs -sSf | sh

install-package:
	uv pip install -e .

.PHONY: format* lint* type-check test* checks*

format-py-check:
	uv run black --config pyproject.toml --target-version py312 --check .

format-py-diff:
	uv run black --config pyproject.toml --target-version py312 --diff .

format-py:
	uv run black --config pyproject.toml --target-version py312 .

format-rust:
	cargo fmt --all -- --check

format: format-py format-rust

lint-py:
	uv run pylint --rcfile pyproject.toml src/python/ensembl
	uv run pylint --rcfile pyproject.toml src/python/tests

lint-rust:
	cargo clippy --all --all-features --tests -- -D warnings

type-check:
	uv run mypy --config-file pyproject.toml src/python/ensembl
	uv run mypy --config-file pyproject.toml src/python/tests

test-py: install-package
	uv run coverage run -m pytest

test-py-cov: test-py
	uv run coverage report -m

test-py-int: install-package
	uv run coverage run -m pytest --run-integration src/python
	uv run coverage report -m

test-rust:
	cargo test --all-features

test: test-py test-rust

checks-py: format-py-check lint-py type-check test-py

checks-rust: format-rust lint-rust test-rust

checks-all: checks-rust format-py-check lint-py type-check test-py-int
