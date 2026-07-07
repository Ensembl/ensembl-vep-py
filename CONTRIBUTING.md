> **_NOTE:_**  This project is in an initial development phase. We are not accepting any contribution at this moment. The following content is provided as an example for future guideline purposes.

# Contributing to Ensembl VEP

Thank you for your interest in contributing to **Ensembl Variant Effect Predictor (VEP)**! We welcome contributions that improve functionality, documentation, performance, or usability.

## 🔬 Project Overview

Ensembl VEP (Variant Effect Predictor) predicts the functional effects of genomic variants. Contributions should maintain high standards for scientific accuracy, performance, and reproducibility.

## Pre-requisites
While most code in this repository is developer-platform agnostic, ensembl-vep is assumed to be run mainly on *NIX systems (linux or MacOS).
If you are developing on Windows, you are advised to use a linux distribution (such as Ubuntu) through WSL as your development environment,
to ensure full compatibility with all tools used by this repository (such as `make`).

The following packages are pre-requisites to be able to contribute to the Ensembl VEP code:
 * [Cargo](https://doc.rust-lang.org/cargo/): Rust package manager
   For convenience, Cargo can be installed using the `Makefile` in this repository by running `make install-cargo`.
 * [UV](https://docs.astral.sh/uv/): python package manager
   For convenience, uv can be installed using the `Makefile` in this repository by running `make install-uv`.
 * Python 3.12 dev library. On Ubuntu, this can be installed with `apt install libpython-3.12-dev`.

## 🚀 Getting Started

1. **Clone the repository**:
   ```bash
   git clone https://github.com/Ensembl/ensembl-vep-py.git
   cd ensembl-vep-py
   ```
2. **Set up uv venv**:
    ```bash
    uv venv
    source .venv/bin/activate
    ```

3. **Install dependencies**:
   ```bash
   # Rust components
   cargo install
   
   # Python environment
   uv sync --extra cicd,build
   ```

4. **Building the Project**:
    Create a `MANIFEST.in` file with following content -
   ```bash
   include Cargo.toml
   recursive-include src/rust *.rs
   ```

   and, run -
   ```bash
   uv pip install -e .
   ```

5. **Run tests**:
   ```bash
   ## Rust tests with coverage (install first via: cargo install cargo-llvm-cov)
   cargo llvm-cov --ignore-filename-regex '(^|/)src/rust/lib\.rs$' --show-missing-lines
   ## Python tests with coverage
   uv run coverage run -m pytest
   uv run coverage report -m

   ## Integration test for tabix (requires tabix to be installed and available in PATH)
   cargo test --test tabix_integration -- --ignored
   uv run coverage run -m pytest --run-integration src/python
   ```

## 📋 Development Workflow

### 1. Branching Strategy
```
main          # Stable releases
release/xxx   # Release branches  
feature/xxx   # New features (e.g., feature/faster-vcf-parsing)
bugfix/xxx    # Bug fixes
chore/xxx     # Maintenance tasks (e.g., chore/update-dependencies)
task/xxx      # Miscellaneous tasks (e.g., task/add-new-tests)
```

#### Workflow for Contributors

1. **Clone the repository** to your local machine

2. **Create your development branch** from the remote main branch `origin/main`

3. **Make your changes** and commit them with clear messages

4. **Pre commit checks** do pre-commit checks (see [Pre-commit Checks](#3-pre-commit-checks)) locally and fix any failing tests

5. **Create a Pull Request** from your development branch to `origin/main` and ensure all CI checks pass

6. **After review and approval**, maintainers will merge your PR

7. **Releases** are created from `main` by tagging specific commits as `release/xxx`

#### Branch Flow Diagram
```
Upstream Repo                         
───────────── 

                main ◄────┐
                .         │
                .         │
                .         │ 
   release/1.1 ─┤         │ PR (your contribution)
                │         │           
                ├─ feature/xxx 
                │
                │
                │  ◄──────┐
                │         │
                │         │ PR (your contribution)
                ├─ bugfix/xxx
                │
                │
                │
                │
   release/1.0 ─┤
                │
                .
                .
                .
```

#### Key Points

- ✅ **DO**: Ensure all CI checks pass before creating PR
- ✅ **DO**: Keep your branch up-to-date with upstream
- ❌ **DON'T**: Create PR with failing tests

#### Keeping Your Branch Updated
```bash
# Add upstream remote (one time only)
git remote add upstream https://github.com/original-org/repo.git

# Update your branch with latest changes
git fetch origin
git pull --rebase origin/main

# If you already have a PR open, try doing a merge instead
git fetch origin/main
git merge origin/main
```

### 2. Commit Messages (TBD)


### 3. Pre-commit Checks

1. Python checks
```bash
# Unit and integration tests
uv run coverage run -m pytest
uv run coverage run -m pytest --run-integration src/python

# Static analysis
uv run pylint --rcfile pyproject.toml src/python
uv run mypy --config-file pyproject.toml src/python

# Formatting
uv run black --config pyproject.toml .
```

2. Rust checks
```bash
# Unit and integration tests
cargo llvm-cov --ignore-filename-regex '(^|/)src/rust/lib\.rs$' --show-missing-lines
cargo test --test tabix_integration -- --ignored

# Linting
cargo clippy

# Formatting
cargo fmt
```

## 🧪 Testing Requirements
All PRs must:

- ✅ Pass all existing tests
- ✅ Add tests for new functionality (90%+ coverage)
- ✅ Not decrease test coverage
- ✅ Include type annotations (Python)
- ✅ Handle edge cases (multi-allelic variants, structural variants)

## 📊 Code Coverage

Coverage reports are automatically generated and posted to PRs. Maintain:
- **Rust**: ≥90% line coverage
- **Python**: ≥90% line coverage

Rust coverage is generated with stable Rust through `cargo llvm-cov`. The crate
root `src/rust/lib.rs` is intentionally limited to PyO3 module registration glue
and is excluded with `--ignore-filename-regex`; add new testable Rust logic in
submodules instead.

## 🔧 Code Style

- For Python see `pyproject.toml`
- For Rust see `Cargo.toml`

## 🐛 Bug Reports (TBD)
1. Use the [GitHub Issues](https://github.com/ensembl-vep-py/issues) template
2. Include:
   - VEP version
   - Input data sample (VCF/JSON)
   - Expected vs actual output
   - Steps to reproduce
   - Environment (OS, Python/Rust versions)

## ✨ Feature Requests (TBD)
1. Open an issue with:
   - Use case (scientific/biological motivation)
   - Proposed API changes
   - Performance requirements
   - Compatibility constraints

2. Discuss before implementing major features

## 📚 Documentation (TBD)

- Update `README.md` for user-facing changes
- Add examples to `docs/examples/`
- Update API docs in code (Rust `///` comments, Python docstrings)
- Plugin documentation in `docs/plugins/`

## 🔌 Plugins (TBD)

New plugins should:
1. Follow `plugins/template/` structure
2. Include comprehensive tests
3. Document input/output format
4. Provide example usage


## 🚫 What we don't accept

- Untested code
- Breaking changes without migration path
- Deprecated bioinformatics formats without justification
- Performance regressions

## 🌍 Community Guidelines

- Be respectful and constructive
- Focus on code, not people
- Scientific accuracy > performance

## 💬 Support (TBD)

- Bugs: [Issues](https://github.com/ensembl-vep-py/issues)
- Features: [Discussions → Ideas](https://github.com/ensembl-vep-py/discussions/categories/ideas)

## 🎉 Thank You!

Your contributions help researchers worldwide analyze genomic variants faster and more accurately. Welcome aboard! 🚀

---

*Last updated: July 2026*
