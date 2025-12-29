# Quick Start Guide

This guide will help you get started with processing radiosonde data using the radiosonde-teamxuk package.

## Basic Usage

### Processing Radiosonde Data

Convert EDT format files to NCAS-AMOF NetCDF format:

```bash
process-teamxuk-radiosondes input_directory/ output_directory/
```

**Example:**
```bash
process-teamxuk-radiosondes test_edt_output/ processed_data/
```

This will:
- Read all EDT files from `test_edt_output/`
- Convert them to NetCDF format
- Save the output files in `processed_data/`
- Use default metadata from the package

### Using Custom Metadata

Specify a custom metadata directory:

```bash
process-teamxuk-radiosondes input_dir/ output_dir/ metadata_dir/
```

**Example:**
```bash
process-teamxuk-radiosondes test_edt_output/ processed_data/ my_metadata/
```

The metadata directory should contain:
- `metadata_ncas-radiosonde-*.json` - Instrument metadata files

### Generating Quicklook Plots

Create visualization plots from NetCDF files:

```bash
generate-teamxuk-quicklooks netcdf_directory/
```

**Example:**
```bash
generate-teamxuk-quicklooks processed_data/
```

This generates plots showing:
- Temperature profile
- Relative humidity profile
- Wind speed and direction
- Trajectory map

### Adding Stability Analysis

Include atmospheric stability plots:

```bash
generate-teamxuk-quicklooks netcdf_directory/ --stability
```

**Example:**
```bash
generate-teamxuk-quicklooks processed_data/ --stability
```

This adds additional plots for:
- Skew-T log-P diagram
- Stability indices

## Python API Usage

You can also use the package directly in Python:

### Processing Files

```python
from radiosonde_teamxuk import process_radiosondes

# Process a single file
process_radiosondes.process_file(
    'edt_files/edt1sdataforv217_20250224_1101.txt',
    'processed_data/',
    metadata_dir='metadata/'
)
```

### Generating Plots

```python
from radiosonde_teamxuk import generate_quicklooks

# Generate quicklook for a NetCDF file
generate_quicklooks.generate_quicklook(
    'processed_data/ncas-radiosonde-2_sterzing_20250224-110100_sonde_teamxuk_v1.0.nc',
    stability=True
)
```

## Example Workflow

Here's a complete workflow from raw data to visualizations:

```bash
# 1. Process all EDT files
process-teamxuk-radiosondes edt_files/ processed_data/

# 2. Generate quicklook plots with stability analysis
generate-teamxuk-quicklooks processed_data/ --stability

# 3. Review the output
ls processed_data/
ls processed_data/quicklooks/
```

## Output Files

### Processed NetCDF Files

Files follow the NCAS-AMOF naming convention:
```
ncas-radiosonde-2_sterzing_YYYYMMDD-HHMMSS_sonde_teamxuk_v1.0.nc
```

### Quicklook Plots

Plots are saved in the `quicklooks/` subdirectory:
```
processed_data/
├── ncas-radiosonde-2_sterzing_20250224-110100_sonde_teamxuk_v1.0.nc
└── quicklooks/
    └── ncas-radiosonde-2_sterzing_20250224-110100_sonde_teamxuk_v1.0.png
```

## Next Steps

- See [Usage](usage.md) for detailed information on all features
- Check the [API Reference](api.rst) for programmatic usage
- Read [Contributing](contributing.md) to help improve the package
