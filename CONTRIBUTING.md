> **_NOTE:_**  This project is in an initial development phase. We are not accepting any contribution at this moment. The following content is provided as an example for future guideline purposes.

# Contributing to Ensembl VEP

Thank you for your interest in contributing to **Ensembl Variant Effect Predictor (VEP)**! We welcome contributions that improve functionality, documentation, performance, or usability.

## 🔬 Project Overview

Ensembl VEP (Variant Effect Predictor) predicts the functional effects of genomic variants. Contributions should maintain high standards for scientific accuracy, performance, and reproducibility.

## Pre-requisites
The following packages are pre-requisites to be able to contribute to the Ensembl VEP code:
 * [Cargo](https://doc.rust-lang.org/cargo/): Rust package manager
   For convenience, Cargo can be installed using the `Makefile` in this repository by running `make install-cargo`.
 * [UV](https://docs.astral.sh/uv/): python package manager
   For convenience, uv can be installed using the `Makefile` in this repository by running `make install-uv`.
 * Python 3.12 dev library. On Ubuntu, this can be installed with `apt install libpython-3.12-dev`.

## 🚀 Getting Started

1. **Fork the repository**
2. **Clone your fork**:
   ```bash
   git clone https://github.com/yourusername/ensembl-vep-py.git
   cd ensembl-vep-py
   ```
3. **Install dependencies**:
   ```bash
   # Rust components
   cargo install
   
   # Python environment
   uv sync --extra dev,test
   ```

3. **Building the Project**:
    Create a `MANIFEST.in` file with following content -
   ```bash
   include Cargo.toml
   recursive-include src/rust *.rs
   ```

   and, run -
   ```bash
   uv build
   ```

## 📋 Development Workflow

### 1. Branching Strategy
```
main          # Stable releases
release/xxx   # Release branches  
feature/xxx   # New features (e.g., feature/faster-vcf-parsing)
bugfix/xxx    # Bug fixes
```

#### Workflow for Contributors

1. **Fork the repository** to your GitHub account

2. **Create your development branch** from the relevant upstream branch:
```bash
   # For new features
   git checkout -b my-feature-name origin/feature/xxx
   
   # For bug fixes
   git checkout -b my-bugfix-name origin/bugfix/xxx
```

3. **Make your changes** and commit them with clear messages

4. **Ensure CI passes** locally and in your fork:
   - Push to your fork and verify all GitHub Actions pass
   - Fix any failing tests before creating a PR

5. **Create a Pull Request**:
   - **From:** your fork's branch (`my-feature-name`)
   - **To:** upstream's `feature/xxx` or `bugfix/xxx` branch
   - NOT directly to `main`

6. **After review and approval**, maintainers will merge your PR to the feature/bugfix branch

7. **Maintainers** will create a separate PR from `feature/xxx` or `bugfix/xxx` to `main` when ready

8. **Releases** are created from `main` by tagging specific commits as `release/xxx`

#### Branch Flow Diagram
```
Upstream Repo                          Your Fork
─────────────                          ─────────

                main ◄────┐
                .         │
                .         │
                .         │ 
   release/1.1 ─┤         │ PR (reviewed by maintainers)
                │         │           
                ├─ feature/xxx ◄─────┐
                │                    │
                │                    │ PR (your contribution)
                │                    │
                │                    └─── my-feature-name (your branch)
                │
                ├─ bugfix/xxx ◄──────┐
                │                    │
                │                    │ PR (your contribution)
                │                    │
                │                    └─── my-bugfix-name (your branch)
   release/1.0 ─┤
                │
                .
                .
                .


Timeline Flow:
──────────────

1. Fork repo
2. Create branch from feature/xxx or bugfix/xxx
3. Make changes → Run tests → Push to fork
4. ✓ CI passes in your fork
5. Create PR to upstream feature/xxx or bugfix/xxx
6. Code review & approval
7. Merge to feature/xxx or bugfix/xxx
8. Maintainers merge to main
9. Release created from main
```

#### Key Points

- ✅ **DO**: Create PR to `feature/xxx` or `bugfix/xxx` branches
- ✅ **DO**: Ensure all CI checks pass before creating PR
- ✅ **DO**: Keep your branch up-to-date with upstream
- ❌ **DON'T**: Create PR directly to `main`
- ❌ **DON'T**: Create PR with failing tests

#### Keeping Your Branch Updated
```bash
# Add upstream remote (one time only)
git remote add upstream https://github.com/original-org/repo.git

# Update your branch with latest changes
git fetch upstream
git rebase upstream/feature/xxx
```

### 2. Commit Messages
(TBD)
Follow [Conventional Commits](https://www.conventionalcommits.org/):
```
feat: add support for GRCh38.p14
fix: resolve VCF parsing edge case for multi-allelic variants
docs: update plugin installation instructions
perf: optimize transcript lookup by 25%
refactor: extract variant parser into separate module
```

### 3. Pre-commit Checks
(TBD)
```bash
# Format code
cargo fmt
ruff format src/python/

# Lint
cargo clippy
ruff check src/python/
mypy src/python/

# Test
cargo test --all-features
pytest src/python/ --cov=src/python
```

## 🧪 Testing Requirements
(TBD)
All PRs must:

- ✅ Pass all existing tests
- ✅ Add tests for new functionality (90%+ coverage)
- ✅ Not decrease test coverage
- ✅ Include type annotations (Python)
- ✅ Handle edge cases (multi-allelic variants, structural variants)

**Rust tests**: `cargo test --lib --tests`
**Python tests**: `pytest src/python/ -v`

## 📊 Code Coverage

Coverage reports are automatically generated and posted to PRs. Maintain:
- **Rust**: ≥90% line coverage
- **Python**: ≥90% line coverage

## 🔧 Code Style

- For Python see `pyproject.toml`
- For Rust see `Cargo.toml`

## 🐛 Bug Reports
(TBD)
1. Use the [GitHub Issues](https://github.com/ensembl-vep-py/issues) template
2. Include:
   - VEP version
   - Input data sample (VCF/JSON)
   - Expected vs actual output
   - Steps to reproduce
   - Environment (OS, Python/Rust versions)

## ✨ Feature Requests
(TBD)
1. Open an issue with:
   - Use case (scientific/biological motivation)
   - Proposed API changes
   - Performance requirements
   - Compatibility constraints

2. Discuss before implementing major features

## 🔬 Pull Request Checklist

(TBD)
```markdown
- [ ] Passes `cargo test --all-features`
- [ ] Passes `pytest --cov-fail-under=90`
- [ ] Passes `cargo fmt` and `ruff format`
- [ ] Passes `cargo clippy` and `ruff check`
- [ ] Passes `mypy src/python/`
- [ ] Updated documentation
- [ ] Added/changed tests
- [ ] Maintains 90% Rust / ≥90% Python coverage
- [ ] Squash merged commits
```

## 📚 Documentation

- Update `README.md` for user-facing changes
- Add examples to `docs/examples/`
- Update API docs in code (Rust `///` comments, Python docstrings)
- Plugin documentation in `docs/plugins/`

## 🔌 Plugins

New plugins should:
1. Follow `plugins/template/` structure
2. Include comprehensive tests
3. Document input/output format
4. Provide example usage

## 🤝 Code Owners

```
# CODEOWNERS
# Core VEP logic
/src/rust/ @ensembl-vep-core

# Python bindings  
/src/python/ @ensembl-vep-python

# Documentation
/docs/ @ensembl-vep-docs
```

## 🚫 What we don't accept

- Untested code
- Breaking changes without migration path
- Deprecated bioinformatics formats without justification
- Performance regressions

## 🌍 Community Guidelines

- Be respectful and constructive
- Focus on code, not people
- Scientific accuracy > cleverness
- Performance > micro-optimizations

## 💬 Support

- Bugs: [Issues](https://github.com/ensembl-vep-py/issues)
- Features: [Discussions → Ideas](https://github.com/ensembl-vep-py/discussions/categories/ideas)

## 🎉 Thank You!

Your contributions help researchers worldwide analyze genomic variants faster and more accurately. Welcome aboard! 🚀

---

*Last updated: January 2026*