# Setup

This page describes how to set up a local checkout for development. If you only want to use Ensembl VEP, follow the [Download and install](../user_guide/install.md) instructions instead.

## Pre-requisites

Development is supported on Linux and macOS. On Windows, use a Linux distribution through [WSL](https://learn.microsoft.com/en-us/windows/wsl/) so the same tools and Make targets are available.

The development environment requires:

- Git and Make;
- Python 3.12 or later, including the Python development headers;
- [uv](https://docs.astral.sh/uv/);
- a Rust toolchain, including Cargo.

The integration tests also require `bgzip` and `tabix` from [HTSlib](https://www.htslib.org/). They are not required for the unit tests.

The Makefile provides convenience targets for installing uv and Rust when they are not already available:

```bash
make install-uv
make install-cargo
```

Restart the shell after installing either tool so the updated `PATH` is available.

## Clone the repository

```bash
git clone https://github.com/Ensembl/ensembl-vep-py.git
cd ensembl-vep-py
```

Create a branch from the remote `main` branch:

```bash
git fetch origin
git switch --create feature/my-change origin/main
```

## Create the development environment

Create and activate a virtual environment:

```bash
uv venv --python 3.12
source .venv/bin/activate
```

Install the development, build and documentation dependencies:

```bash
uv sync --extra cicd --extra build --extra docs
```

Install Ensembl VEP in editable mode so changes in the checkout are immediately available to Python:

```bash
make install-package
```

## Check the setup

Run the Python and Rust unit tests:

```bash
make test
```

Run all formatting, linting, type-checking and integration checks before opening a pull request:

```bash
make checks-all
```

The integration checks require `bgzip` and `tabix` to be available on `PATH`. See the [Contributing Guide](https://github.com/Ensembl/ensembl-vep-py/blob/main/CONTRIBUTING.md) for the branching and review workflow.
