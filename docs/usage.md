# Usage Guide

Complete guide to using the radiosonde-teamxuk package.

## Command-Line Tools

### process-teamxuk-radiosondes

Process EDT format radiosonde data files and convert them to NCAS-AMOF NetCDF format.

#### Syntax

```bash
process-teamxuk-radiosondes <input_directory> <output_directory> [metadata_directory]
```

#### Arguments

- `input_directory` (required): Directory containing EDT format files (`.nc`)
- `output_directory` (required): Directory where NetCDF files will be saved
- `metadata_directory` (optional): Directory containing custom metadata JSON files

#### Examples

Basic processing:
```bash
process-teamxuk-radiosondes test_edt_output/ processed_data/
```

With custom metadata:
```bash
process-teamxuk-radiosondes test_edt_output/ processed_data/ my_metadata/
```

Process specific campaign data:
```bash
process-teamxuk-radiosondes /data/teamxuk/raw/ /data/teamxuk/processed/
```

#### Output

The tool creates NetCDF files following the NCAS-AMOF naming convention:
```
ncas-radiosonde-2_teamx-uk_sterzing_YYYYMMDD-HHMMSS_v1.0.nc
```

Each file contains:
- Atmospheric profiles (pressure, temperature, humidity)
- Wind data (speed and direction)
- GPS position data
- Derived variables (potential temperature, equivalent potential temperature, mixing ratio)
- Comprehensive metadata

### generate-teamxuk-quicklooks

Generate quicklook visualization plots from processed NetCDF files.

#### Syntax

```bash
generate-teamxuk-quicklooks <netcdf_directory> [--stability]
```

#### Arguments

- `netcdf_directory` (required): Directory containing NetCDF files
- `--stability` (optional): Include atmospheric stability analysis plots

#### Examples

Basic quicklook generation:
```bash
generate-teamxuk-quicklooks processed_data/
```

With stability analysis:
```bash
generate-teamxuk-quicklooks processed_data/ --stability
```

#### Output

Plots are saved in `<netcdf_directory>/quicklooks/` with naming:
```
<original_filename>_quicklook.png
```

**Standard Quicklook** includes:
- Temperature profile vs altitude
- Relative humidity profile
- Wind speed and direction
- Trajectory map with launch location

**With --stability flag** adds:
- Skew-T log-P diagram
- CAPE (Convective Available Potential Energy)
- CIN (Convective Inhibition)
- Lifted Index
- Other stability indices

## Python API

### Processing Module

```python
from radiosonde_teamxuk import process_radiosondes
```

#### Process a Single File

```python
process_radiosondes.process_file(
    edt_file='test_edt_output/edt1sdata_20250224_1101.nc',
    output_dir='processed_data/',
    metadata_dir=None  # Use default metadata
)
```

#### Process Multiple Files

```python
import glob

edt_files = glob.glob('test_edt_output/*.nc')
for edt_file in edt_files:
    process_radiosondes.process_file(
        edt_file=edt_file,
        output_dir='processed_data/',
        metadata_dir='metadata/'
    )
```

### Quicklooks Module

```python
from radiosonde_teamxuk import generate_quicklooks
```

#### Generate Single Quicklook

```python
generate_quicklooks.generate_quicklook(
    netcdf_file='processed_data/ncas-radiosonde-2_teamx-uk_sterzing_20250224-110100_v1.0.nc',
    stability=True
)
```

#### Batch Processing

```python
import glob

nc_files = glob.glob('processed_data/*.nc')
for nc_file in nc_files:
    try:
        generate_quicklooks.generate_quicklook(nc_file, stability=True)
        print(f"Generated quicklook for {nc_file}")
    except Exception as e:
        print(f"Error processing {nc_file}: {e}")
```

## Data Format

### Input: EDT Files

EDT (Enhanced Data Transfer) files are NetCDF format files produced by Vaisala ground stations. They contain:
- Raw radiosonde measurements
- GPS data
- Telemetry information

### Output: NCAS-AMOF NetCDF Files

The processed files follow the NCAS-AMOF standard and include:

**Coordinates:**
- `altitude` - Height above mean sea level (m)
- `time` - Time since launch (seconds)

**Variables:**
- `air_pressure` - Atmospheric pressure (Pa)
- `air_temperature` - Temperature (K)
- `relative_humidity` - Relative humidity (%)
- `wind_speed` - Wind speed (m/s)
- `wind_from_direction` - Wind direction (degrees)
- `eastward_wind` - U component (m/s)
- `northward_wind` - V component (m/s)
- `latitude` - GPS latitude (degrees)
- `longitude` - GPS longitude (degrees)

**Derived Variables:**
- `air_potential_temperature` - Potential temperature (K)
- `equivalent_potential_temperature` - Equivalent potential temperature (K)
- `humidity_mixing_ratio` - Water vapor mixing ratio (kg/kg)

## Metadata Configuration

### Default Metadata

The package includes default metadata files:
- `metadata_ncas-radiosonde-1.json`
- `metadata_ncas-radiosonde-2.json`

### Custom Metadata

To use custom metadata:

1. Create a directory (e.g., `my_metadata/`)
2. Copy and modify the default JSON files
3. Specify the directory when processing:

```bash
process-teamxuk-radiosondes input/ output/ my_metadata/
```

### Metadata Fields

Key metadata fields include:
- `platform` - Launch site information
- `instrument` - Radiosonde model and serial number
- `campaign` - Campaign name and description
- `contacts` - Principal investigators and data managers

## Best Practices

### File Organization

Organize your data with a clear directory structure:

```
project/
├── raw_edt/          # Original EDT files
├── processed/        # Converted NetCDF files
│   └── quicklooks/   # Generated plots
└── metadata/         # Custom metadata (if needed)
```

### Batch Processing

For large datasets, process in batches and handle errors:

```python
import os
import glob
from radiosonde_teamxuk import process_radiosondes

input_dir = 'raw_edt/'
output_dir = 'processed/'

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Process all files
edt_files = glob.glob(f'{input_dir}/*.nc')
for edt_file in edt_files:
    try:
        process_radiosondes.process_file(edt_file, output_dir)
        print(f"✓ Processed {os.path.basename(edt_file)}")
    except Exception as e:
        print(f"✗ Error processing {os.path.basename(edt_file)}: {e}")
```

### Quality Control

Always review quicklook plots for data quality:

```bash
# Process data
process-teamxuk-radiosondes raw_edt/ processed/

# Generate quicklooks with stability analysis
generate-teamxuk-quicklooks processed/ --stability

# Review plots
ls processed/quicklooks/
```

## Troubleshooting

### Missing Data

If processed files have missing data:
- Check that EDT files are valid
- Verify radiosonde successfully completed its ascent
- Review raw EDT file for data gaps

### Plot Generation Failures

If quicklook generation fails:
- Ensure NetCDF file contains required variables
- Check that coordinate dimensions are valid
- Verify matplotlib backend is properly configured

### Memory Issues

For large batches:
- Process files individually rather than all at once
- Close NetCDF datasets explicitly after reading
- Use batch processing scripts with error handling
