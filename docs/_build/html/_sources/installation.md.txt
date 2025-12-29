# Installation

## Requirements

- Python 3.9 or higher
- pip (Python package installer)

## Installation Methods

### From GitHub (Recommended for Users)

Install the latest version directly from GitHub:

```bash
pip install git+https://github.com/longlostjames/radiosonde-teamxuk.git
```

### From Source (For Development)

Clone the repository and install in editable mode:

```bash
git clone https://github.com/longlostjames/radiosonde-teamxuk.git
cd radiosonde-teamxuk
pip install -e .
```

The `-e` flag installs the package in "editable" mode, allowing you to modify the code and see changes immediately without reinstalling.

### From Local Directory

If you have downloaded the source code:

```bash
cd radiosonde-teamxuk
pip install .
```

## Dependencies

The package will automatically install the following dependencies:

- numpy
- matplotlib
- metpy
- cartopy
- netCDF4
- ncas-amof-netcdf-template
- python-eccodes

## Verifying Installation

After installation, verify that the command-line tools are available:

```bash
process-teamxuk-radiosondes --version
generate-teamxuk-quicklooks --version
```

You should see version `0.1.0` displayed.

## Troubleshooting

### Command Not Found

If you get a "command not found" error, the scripts may not be in your PATH. You can:

1. **Add to PATH**: Add `~/.local/bin` to your PATH:
   ```bash
   export PATH="$HOME/.local/bin:$PATH"
   ```

2. **Use Python module**: Run the scripts as Python modules:
   ```bash
   python -m radiosonde_teamxuk.process_radiosondes
   python -m radiosonde_teamxuk.generate_quicklooks
   ```

3. **Reinstall**: Try reinstalling with `--force-reinstall`:
   ```bash
   pip install --force-reinstall -e .
   ```

### Missing ECCODES

If you encounter issues with ECCODES, you may need to install it separately:

```bash
# On Ubuntu/Debian
sudo apt-get install libeccodes-dev

# On macOS with Homebrew
brew install eccodes
```

Then install the Python bindings:

```bash
pip install python-eccodes
```
