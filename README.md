# Ensembl VEP

[![License](https://img.shields.io/badge/license-Apache_2.0-blue.svg)](https://github.com/Ensembl/ensembl-vep-py/blob/main/LICENSE)
[![Python-Coverage](https://ensembl.github.io/ensembl-vep-py/latest/coverage/python/coverage-badge.svg)](https://ensembl.github.io/ensembl-vep-py/latest/coverage/python/htmlcov/index.html)
[![Rust-Coverage](https://ensembl.github.io/ensembl-vep-py/latest/coverage/rust/coverage-badge.svg)](https://ensembl.github.io/ensembl-vep-py/latest/coverage/rust/html/index.html)
[![CI](https://github.com/Ensembl/ensembl-vep-py/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/Ensembl/ensembl-utilstemplate-python/actions/workflows/ci.yml)
[![Release](https://img.shields.io/pypi/v/ensembl-vep)](https://pypi.org/project/ensembl-vep)

Ensembl VEP (Variant Effect Predictor) predicts the functional effects of genomic variants.

> **_NOTE:_**  This project is in an initial development phase. Many features are incomplete, APIs may change without notice, and documentation is partial. Use at your own risk and feel free to open issues or contribute; your feedback is welcome.

## Documentation

Full documentation is available at: https://ensembl.github.io/ensembl-vep-py/

## Installation

### Requirements

- Python 3.12 or higher
- Rust 1.95.0 or higher (TBD - not strictly tested)
- Dependencies listed in `pyproject.toml` and `Cargo.toml`

### Using pip (TBD - after PyPI release)

```bash
pip install ensembl-vep
```

### Using uv (recommended)  (TBD - after PyPI release)

```bash
uv pip install ensembl-vep
```

### Development Installation

1. clone the repository
    ```bash
    git clone https://github.com/Ensembl/ensembl-vep-py.git
    cd ensembl-vep-py
    ```

2. Create a branch against remote main branch
    ```bash
    git checkout -b feature/some_feature origin/main
    ```

3. Install the package in editable mode
    ```bash
    uv pip install -e .
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

### Cache Files (TBD)

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

## Citation (TBD)

If you use this software, please cite:

McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GR, Thormann A, Flicek P, Cunningham F.  
The Ensembl Variant Effect Predictor.  
*Genome Biology* Jun 6;17(1):122. (2016)  
doi:[10.1186/s13059-016-0974-4](https://doi.org/10.1186/s13059-016-0974-4)

## Support (TBD)

- Bugs: [Issues](https://github.com/ensembl-vep-py/issues)
- Features: [Discussions → Ideas](https://github.com/ensembl-vep-py/discussions/categories/ideas)

## Changelog (TBD)

See [CHANGELOG.md](CHANGELOG.md) for version history and updates.

## Acknowledgments (TBD)
