# Ensembl Python project template

[![License](https://img.shields.io/badge/license-Apache_2.0-blue.svg)](https://github.com/Ensembl/ensembl-vep-py/blob/main/LICENSE)
[![Rust-Coverage](https://ensembl.github.io/template-python/coverage/coverage-badge.svg)](https://ensembl.github.io/ensembl-vep-py/coverage)
[![Py-Coverage](https://ensembl.github.io/template-python/coverage/coverage-badge.svg)](https://ensembl.github.io/ensembl-vep-py/coverage)
[![CI](https://github.com/Ensembl/template-python/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/Ensembl/ensembl-utilstemplate-python/actions/workflows/ci.yml)
[![Release](https://img.shields.io/pypi/v/ensembl-vep)](https://pypi.org/project/ensembl-vep)

This repo structure can be forked and used as the base template for developing a Python package (app or library). It contains boilerplate code and configuration for unit testing, linting, type checking, automatic documentation and CI/CD.

For more information, please consult this repository's [GitHub pages](https://ensembl.github.io/template-python/).

## Getting Started

This template is a simple but fully functioning Python project.

### Template Project Layout

```plain
.
├── .github               <- GitHub Actions workflows
│   ├── actions
│   │   └── python_build
│   │       └── action.yml
│   ├── workflows
│   │   ├── ci_dev.yml
│   │   ├── ci.yml
│   │   ├── docs.yml
│   │   ├── notice.yml
│   │   ├── publish.yml
│   │   └── python_ci.yml
│   └── dependabot.yml
├── docs                  <- Folder with boilerplate for auto-generated mkdocs documentation
│   ├── img
│   │   └── ebang.png
│   ├── scripts
│   │   └── gen_ref_pages.sh  <- Script to generate reference pages and build the documentation
│   ├── stylesheets
│   │   └── extra.css
│   ├── code_of_conduct.md
│   ├── index.md
│   ├── install.md
│   └── usage.md
├── src                   <- Root folder for this project source code
│   └── ensembl
│       └── template
│           ├── __init__.py
│           └── hello_world.py
├── tests                 <- Root folder for tests code
│   └── test_hello_world.py
├── .gitignore            <- File patterns to be ignored by Git
├── .licenserc.yml        <- License header template expected in every file in this repository
├── LICENSE               <- License text
├── mkdocs.yml            <- mkdocs configuration file
├── NOTICE                <- Notice file (required by License)
├── pyproject.toml        <- This project and Python tools packaging configuration
├── README.md             <- Information about this project, e.g. description, typical installation and usage
└── setup.py              <- Setuptools stub required to install the project in edit mode (`pip install -e`)
```

### Creating a new project

The best and easiest way to use this template is to go to [GitHub](https://github.com/Ensembl/template-python) and on the top right click on the button "Use this template -> Create a new repository" and follow the steps to create your own repository based on the code found here.

Alternatively, you can also clone this repository:

```bash
git clone --depth 1 --branch main https://github.com/Ensembl/template-python.git
```

Next, remove the `git` repository data and rename the main folder:

```bash
rm -rf template-python/.git
mv template-python <NEW_PROJECT_NAME>
cd <NEW_PROJECT_NAME>
git init
```

Later, the following files (or folder names) will need to be modified:

- `NOTICE`: Substitute "Ensembl Python Template" with your project's name
- `pyproject.toml`: Change package name and other settings
- `mkdocs.yml`: Change accordingly to project's name
- `docs/*.md`: Change documentation accordingly to new project folder structure and instructions
- `README.md`: Write a meaningful README for your project

Finally, once done with the basic customisation, create the initial commit:

```bash
cd <NEW_PROJECT_NAME>
git add .
git commit -m 'Initial commit'
```

### Installing the development environment with `venv`

```bash
python -m venv <VIRTUAL_ENVIRONMENT_NAME>
source <VIRTUAL_ENVIRONMENT_NAME>/bin/activate
pip install -e .[cicd,doc]
```

### Testing with `pytest`

Run test suite:

```bash
cd <NEW_PROJECT_NAME>
pytest
```

To run tests, calculate and display testing coverage stats:

```bash
cd <NEW_PROJECT_NAME>
coverage run -m pytest
coverage report -m
```

### Generate documentation via `mkdocs`

```bash
cd <NEW_PROJECT_NAME>
mkdocs build
```

Open automatically generated documentation page at `site/index.html`.

### Automatic formatting (PEP8 compliance)

```bash
cd <NEW_PROJECT_NAME>
black --check src tests
```

Use `--diff` to print a diff of what Black would change, without actually changing the files.

To actually reformat all files contained in `src` and `test`:

```bash
cd <NEW_PROJECT_NAME>
black src tests
```

### Linting and type checking

```bash
cd <NEW_PROJECT_NAME>
pylint src tests
mypy src tests
```

`pylint` will check the code for syntax, name errors and formatting style. `mypy` will use type hints to statically type check the code.

It should be relatively easy (and definitely useful) to integrate both `pylint` and `mypy` in your IDE/Text editor.

## Useful resources

### Python Documentation
- [Official Python Docs](https://docs.python.org/3/)

### Python virtual environments management
- [Python venv](https://docs.python.org/3/library/venv.html)

### Auto-generating documentation
- [mkdocs](https://www.mkdocs.org)

### Linting, type checking and formatting
- [Pylint](https://www.pylint.org/)
- [Mypy](https://mypy.readthedocs.io/en/stable/)
- [Black](https://black.readthedocs.io/en/stable/)
- [Google Python Style Guide](https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings)

### Testing
- [pytest](https://docs.pytest.org/en/6.2.x/)
- [Coverage](https://coverage.readthedocs.io/)

### CI/CD
- [GitLab CI/CD](https://docs.gitlab.com/ee/ci/quick_start/)

### Development tools
- [IPython](https://ipython.org/)

### Distributing
- [Packaging Python](https://packaging.python.org/tutorials/packaging-projects/)
- [Apache 2.0 License](https://www.apache.org/licenses/LICENSE-2.0#apply)
- [Semantic Versioning](https://semver.org/)
