# Installation Guide

## Install from source

### Development Installation (Editable)

To install the package in development mode (changes to the code will be reflected immediately):

```bash
cd /path/to/radiosonde-teamxuk
pip install -e .
```

### Regular Installation

To install the package normally:

```bash
cd /path/to/radiosonde-teamxuk
pip install .
```

### Install from Git Repository

You can also install directly from the GitHub repository:

```bash
pip install git+https://github.com/longlostjames/radiosonde-teamxuk.git
```

## After Installation

Once installed, you can use the command-line tools:

```bash
# Process EDT files to NetCDF
process-teamxuk-radiosondes <input_directory> <output_directory> [metadata_directory]

# Generate quicklook plots
generate-teamxuk-radiosonde-quicklooks <input_directory> <output_directory>
```

Or import the package in Python:

```python
from radiosonde_teamxuk import read_edt_file, save_to_ncas_netcdf
from radiosonde_teamxuk import create_quicklook

# Your code here
```

## Dependencies

All required dependencies will be automatically installed with the package:
- numpy
- matplotlib
- metpy
- cartopy
- netCDF4
- ncas-amof-netcdf-template
- eccodes

## Uninstalling

To remove the package:

```bash
pip uninstall radiosonde-teamxuk
```
