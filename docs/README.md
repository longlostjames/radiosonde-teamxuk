# Documentation

This directory contains the Sphinx documentation for radiosonde-teamxuk.

## Building the Documentation

### Prerequisites

Install the documentation dependencies:

```bash
pip install -e ".[docs]"
```

Or install them separately:

```bash
pip install sphinx sphinx-rtd-theme sphinx-autodoc-typehints myst-parser
```

### Build Locally

To build the HTML documentation:

```bash
cd docs
make html
```

The built documentation will be in `docs/_build/html/`. Open `index.html` in a browser to view it.

To clean the build directory:

```bash
make clean
```

### Live Preview

For live reloading during development, install and use `sphinx-autobuild`:

```bash
pip install sphinx-autobuild
sphinx-autobuild docs docs/_build/html
```

Then open http://127.0.0.1:8000 in your browser.

## Documentation Structure

- `index.rst` - Main documentation index
- `installation.md` - Installation instructions
- `quickstart.md` - Quick start guide
- `usage.md` - Detailed usage guide
- `api.rst` - API reference (auto-generated from docstrings)
- `contributing.md` - Contributing guidelines
- `conf.py` - Sphinx configuration
- `_static/` - Static files (CSS, images, etc.)
- `_templates/` - Custom templates

## ReadTheDocs

The documentation is automatically built and hosted on ReadTheDocs when changes are pushed to the repository.

Configuration: `.readthedocs.yaml` in the project root.

Documentation URL: https://radiosonde-teamxuk.readthedocs.io

## Writing Documentation

### Markdown Files

Most content uses Markdown (`.md`) via MyST Parser. See the [MyST documentation](https://myst-parser.readthedocs.io/) for syntax.

### reStructuredText Files

Some files use reStructuredText (`.rst`) for advanced features like API documentation. See the [Sphinx documentation](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html) for syntax.

### API Documentation

API documentation is automatically generated from docstrings in the Python code using Sphinx autodoc. Make sure to write comprehensive docstrings in your code.

## Troubleshooting

If you encounter build errors:

1. Check that all dependencies are installed
2. Ensure Python code can be imported (package is installed)
3. Review docstring formatting in Python files
4. Check for syntax errors in Markdown/RST files

Some docstring formatting warnings are expected and don't prevent the documentation from building successfully.
