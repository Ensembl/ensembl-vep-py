# How to install this library

This Python library requires Python 3.12+ to work.

## Basic installation

To install it, just run the following command:

```bash
pip3 install git+https://github.com/Ensembl/ensembl-vep-py.git
```

## Development-oriented installation

To perform a development installation it is advised to use a [Python virtual environment](https://docs.python.org/3/library/venv.html)
(and activate it):

```bash
python3 -m venv path/to/virtual_env
source path/to/virtual_env/bin/activate
```

And then clone the GitHub repository and install it in edit mode:

```bash
git clone https://github.com/Ensembl/ensembl-vep-py.git
pip3 install -e ensembl-vep-py/.[cicd,doc]
```
