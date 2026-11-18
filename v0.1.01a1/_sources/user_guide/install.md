# Download and install

## Pre-requisites

To install Ensembl VEP from source you need:

- Git;
- [uv](https://docs.astral.sh/uv/);
- a Rust toolchain, including Cargo.

uv installs Python 3.12 when it is not already available. Rust and Cargo are required because the source distribution includes a native Rust extension.

## Docker

> **_NOTE:_** A Docker image is not available yet. Download and installation instructions will be added here when it is published.

## PyPI

Ensembl VEP is published on PyPI as [`ensembl-vep`](https://pypi.org/project/ensembl-vep) and requires Python 3.12 or later. Install it into a virtual environment with pip:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install ensembl-vep
```

Or, using uv (recommended):

```bash
uv venv
source .venv/bin/activate
uv pip install ensembl-vep
```

If no pre-built wheel is available for your platform, the package is built from its source distribution, which needs a Rust toolchain as described in [Pre-requisites](#pre-requisites).

Check that the command is available:

```bash
vep --help
```

## Source

### Install uv

On macOS and Linux:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

On Windows, run the following command in PowerShell:

```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

See the [uv installation guide](https://docs.astral.sh/uv/getting-started/installation/) for package-manager and other installation options. Restart the terminal after installation, then check that uv is available:

```bash
uv --version
```
(You may need to restart your terminal for the changes to take effect or add uv to your PATH manually.)

### Install Rust and Cargo

On macOS and Linux, install Rust with `rustup`:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

On Windows, download and run [`rustup-init.exe`](https://rustup.rs/). The installer may also prompt you to install the Visual Studio C++ Build Tools.

Restart the terminal after installation, then check that both tools are available:

```bash
rustc --version
cargo --version
```
(You may need to restart your terminal for the changes to take effect or add the binaries to your PATH manually.)

### Build and install Ensembl VEP

Clone the repository and change into its directory:

```bash
git clone https://github.com/Ensembl/ensembl-vep-py.git
cd ensembl-vep-py
```

Create a virtual environment using Python 3.12:

```bash
uv venv --python 3.12
source .venv/bin/activate
```

Build and install Ensembl VEP into `.venv`:

```bash
uv pip install .
```

This is a non-editable installation. uv builds a wheel before installing it, so the Rust extension is compiled with Cargo's release profile rather than its debug profile.

Check that the command is available:

```bash
vep --help
```

Continue with the [Quickstart guide](quickstart.md) to annotate the example variants included in the repository.
