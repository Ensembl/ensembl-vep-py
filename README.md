# Ensembl VEP

[![License](https://img.shields.io/badge/license-Apache_2.0-blue.svg)](https://github.com/Ensembl/ensembl-vep-py/blob/main/LICENSE)
[![Python-Coverage](https://img.shields.io/endpoint?url=https://ensembl.github.io/ensembl-vep-py/assets/badges/python_coverage.json)](https://ensembl.github.io/ensembl-vep-py/coverage)
[![Rust-Coverage](https://img.shields.io/endpoint?url=https://ensembl.github.io/ensembl-vep-py/assets/badges/rust_coverage.json)](https://ensembl.github.io/ensembl-vep-py/coverage)
[![CI](https://github.com/Ensembl/ensembl-vep-py/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/Ensembl/ensembl-utilstemplate-python/actions/workflows/ci.yml)
[![Release](https://img.shields.io/pypi/v/ensembl-vep)](https://pypi.org/project/ensembl-vep)

Ensembl VEP (Variant Effect Predictor) predicts the functional effects of genomic variants.

> **_NOTE:_**  This project is in an initial development phase. Many features are incomplete, APIs may change without notice, and documentation is partial. Use at your own risk and feel free to open issues or contribute; your feedback is welcome.

## Documentation

Full documentation is available at: https://ensembl.github.io/ensembl-vep-py/

## Installation

### Requirements

- Python 3.12 or higher
- Dependencies listed in `pyproject.toml` and `Cargo.toml`

### Using pip

```bash
pip install ensembl-vep
```

### Using uv (recommended)

```bash
uv pip install ensembl-vep
```

### Development Installation

(Until the repo is private)
1. fork the repository and clone it
    ```bash
    git clone https://github.com/yourusername/ensembl-vep-py.git
    git remote add upstream https://github.com/Ensembl/ensembl-vep-py.git
    ```

2. Create a branch with upstream to a feature/** or bugfix/** branch
    ```bash
    cd ensembl-vep-py
    git checkout -b some_developement upstream/feature/some_feature
    ```

3. Develop
    ```bash
    git clone https://github.com/Ensembl/ensembl-vep-py.git
    cd ensembl-vep-py
    cargo install
    uv sync --extra dev,test
    ```

See more details in the [Contributing Guide](CONTRIBUTING.md) guide.

## Quick Start

### Basic Usage

```python
from ensembl_vep import hello

print(hello("world!"))
```

### Command-Line Interface

```bash
vep world!
```

## Configuration

### Cache Files
TODO

## Testing

Run the test suite:

```bash
# Run pytest with coverage
uv run coverage run -m pytest

# Run Rust llvm-cov
cargo llvm-cov --ignore-filename-regex '(^|/)src/rust/lib\.rs$'
```

The Rust crate root `src/rust/lib.rs` is reserved for PyO3 module registration
glue. Keep testable Rust logic in submodules; coverage runs ignore `lib.rs` so
the project can use stable Rust without nightly-only coverage attributes.

## Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details.

## Citation

If you use this software, please cite:

McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GR, Thormann A, Flicek P, Cunningham F.  
The Ensembl Variant Effect Predictor.  
*Genome Biology* Jun 6;17(1):122. (2016)  
doi:[10.1186/s13059-016-0974-4](https://doi.org/10.1186/s13059-016-0974-4)

## Support

- Bugs: [Issues](https://github.com/ensembl-vep-py/issues)
- Features: [Discussions → Ideas](https://github.com/ensembl-vep-py/discussions/categories/ideas)

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for version history and updates.

## Acknowledgments

TODO
