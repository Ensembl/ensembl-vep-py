# Ensembl VEP

[![License](https://img.shields.io/badge/license-Apache_2.0-blue.svg)](https://github.com/Ensembl/ensembl-vep-py/blob/main/LICENSE)
[![Python-Coverage](https://ensembl.github.io/ensembl-vep-py/latest/_images/python-coverage-badge.svg)](https://ensembl.github.io/ensembl-vep-py/latest/coverage/python/htmlcov/index.html)
[![Rust-Coverage](https://ensembl.github.io/ensembl-vep-py/latest/_images/rust-coverage-badge.svg)](https://ensembl.github.io/ensembl-vep-py/latest/coverage/rust/html/index.html)
[![CI](https://github.com/Ensembl/ensembl-vep-py/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/Ensembl/ensembl-utilstemplate-python/actions/workflows/ci.yml)
[![Release](https://img.shields.io/pypi/v/ensembl-vep)](https://pypi.org/project/ensembl-vep)

Ensembl VEP (Variant Effect Predictor) predicts the functional effects of genomic variants.

> **_NOTE:_**  This project is in an initial development phase. Many features are incomplete, APIs may change without notice, and documentation is partial. Use at your own risk and feel free to open issues or contribute; your feedback is welcome.

## Documentation

Full documentation is available at: https://ensembl.github.io/ensembl-vep-py/

## Installation

### From PyPI

```bash
pip install ensembl-vep
```

Or alternatively, using uv (recommended)

```bash
uv pip install ensembl-vep
```

### Using source 

#### Requirements

- Git
- uv (https://github.com/Ensembl/uv)    
- a Rust toolchain, including Cargo.

uv installs Python 3.12 when it is not already available. Rust and Cargo are required because the source distribution includes a native Rust extension.

#### Download and install

1. clone the repository
    ```bash
    git clone https://github.com/Ensembl/ensembl-vep-py.git
    cd ensembl-vep-py
    ```

2. Create virtual environment
    ```bash
    uv venv --python 3.12
    source .venv/bin/activate
    ```

3. Install the package
    ```bash
    uv pip install .
    ```

For more details, see the [Download and install](https://ensembl.github.io/ensembl-vep-py/latest/user_guide/install.html) page of the documentation.

## Quick Start

### Basic Usage

Run Ensembl VEP with the example VCF, GTF and FASTA files:

```bash
vep run \
  --input_file examples/homo_sapiens/homo_sapiens_GRCh38.vcf \
  --gtf examples/homo_sapiens/genes.gtf.bgz \
  --fasta examples/homo_sapiens/unmasked.fa.bgz \
  --output_file quickstart_output.vcf \
  --stats_file quickstart_summary.html
```

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

- Report bugs or request features: [Issues](https://github.com/ensembl-vep-py/issues)

## Changelog (TBD)

See [CHANGELOG.md](CHANGELOG.md) for version history and updates.
