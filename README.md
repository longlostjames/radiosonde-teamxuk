# Radiosonde TEAMxUK Processing

Python tools for processing Vaisala RS41-SGP radiosonde data from the TEAMxUK campaign, converting EDT files to NCAS-AMOF compliant NetCDF format and generating quicklook plots.

## Overview

This repository contains two main scripts:
- **`process_radiosondes.py`**: Converts EDT files to NCAS-AMOF NetCDF format
- **`generate_quicklooks.py`**: Creates quicklook plots from NetCDF files

## Requirements

- Python 3.9+
- Required packages:
  - `numpy`
  - `matplotlib`
  - `metpy`
  - `cartopy`
  - `netCDF4`
  - `ncas-amof-netcdf-template`

Install dependencies:
```bash
pip install -r requirements.txt
```

## Usage

### 1. Processing EDT Files to NetCDF

Convert Vaisala EDT (Extended Data) files to NCAS-AMOF compliant NetCDF format:

```bash
python process_radiosondes.py <input_directory> <output_directory> [metadata_directory]
```

**Example:**
```bash
python process_radiosondes.py ./edt_files ./output
```

**Features:**
- Automatically detects instrument type (ncas-radiosonde-1 or ncas-radiosonde-2) from hardware version in EDT files
- Selects appropriate metadata file based on instrument
- Extracts MW41 software version dynamically from EDT files
- Creates directory structure: `output/instrument/project/YYYY/MM/`
- Generates NCAS-compliant filenames: `instrument_platform_YYYYMMDDHHmmss_sonde_project_v1.0.nc`

**Metadata Files:**
The script looks for instrument-specific metadata files:
- `metadata_ncas-radiosonde-1.json` for ncas-radiosonde-1
- `metadata_ncas-radiosonde-2.json` for ncas-radiosonde-2

These files contain project information, PI details, instrument specifications, etc.

**EDT File Format:**
Only processes EDT files starting with `edt1sdataforv217*`. The script reads:
- Launch location (lat/lon)
- Release date/time
- Sonde serial number
- Ground check device hardware version (0x70 → radiosonde-1, 0x60 → radiosonde-2)
- MW41 software version
- Station name
- Profile data (pressure, temperature, humidity, winds, GPS trajectory)

### 2. Generating Quicklook Plots

Create quicklook plots from NetCDF files:

```bash
python generate_quicklooks.py <input_directory> [--stability]
```

**Example:**
```bash
python generate_quicklooks.py /path/to/netcdf/files
python generate_quicklooks.py /path/to/netcdf/files --stability
```

**Features:**
- Searches recursively for all `.nc` files in input directory
- Creates plots with:
  - Skew-T log-P diagram
  - Hodograph
  - Trajectory map with topography
- Optional CAPE/CIN shading with `--stability` flag
- Consistent map bounds and altitude scaling across all flights
- Saves plots to `input_directory/quicklooks/` in flat structure

**Output:**
- PNG files saved as: `input_directory/quicklooks/filename.png`
- Resolution: 150 DPI
- All plots from subdirectories combined in single quicklooks folder

## File Naming Convention

NetCDF files follow NCAS naming convention:
```
{instrument}_{platform}_{YYYYMMDDHHmmss}_sonde_{project}_v1.0.nc
```

Example:
```
ncas-radiosonde-2_sterzing_20250225-080054_sonde_teamxuk_v1.0.nc
```

## Directory Structure

### Input (EDT files):
```
input_dir/
  ├── edt1sdataforv217_20250225_1100.nc
  ├── edt1sdataforv217_20250225_2300.nc
  └── ...
```

### Output (NetCDF files):
```
output_dir/
  ├── ncas-radiosonde-1/
  │   └── teamxuk/
  │       └── 2025/
  │           ├── 02/
  │           │   ├── ncas-radiosonde-1_sterzing_20250211-230026_sonde_teamxuk_v1.0.nc
  │           │   └── ...
  │           └── quicklooks/
  │               ├── ncas-radiosonde-1_sterzing_20250211-230026_sonde_teamxuk_v1.0.png
  │               └── ...
  └── ncas-radiosonde-2/
      └── teamxuk/
          └── 2025/
              ├── 02/
              │   ├── ncas-radiosonde-2_sterzing_20250225-080054_sonde_teamxuk_v1.0.nc
              │   └── ...
              └── quicklooks/
                  ├── ncas-radiosonde-2_sterzing_20250225-080054_sonde_teamxuk_v1.0.png
                  └── ...
```

## Data Variables

NetCDF files include:
- **Coordinates**: time, latitude, longitude, altitude
- **Atmospheric**: air_pressure, air_temperature, relative_humidity
- **Wind**: wind_speed, wind_from_direction
- **Balloon**: upward_balloon_velocity, elapsed_time
- **Metadata**: instrument info, project details, geospatial bounds, serial numbers

All data follows CF-1.8 conventions and NCAS-AMOF standards.

## Instrument Mapping

Hardware version detection:
- `0x70` → ncas-radiosonde-1 (Vaisala Sounding Station unit 1)
- `0x60` → ncas-radiosonde-2 (Vaisala Sounding Station unit 2)
- Unknown → defaults to ncas-radiosonde-2

## Authors

- Chris Walden (chris.walden@ncas.ac.uk)
- Processing software: https://github.com/longlostjames/radiosonde-teamx

## Project

TEAMxUK: Quantifying atmospheric processes in mountainous regions
- PI: Charles Chemel (charles.chemel@ncas.ac.uk)
- Instrument Scientists: Hugo Ricketts, Philip Rosenberg

## License

This software is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

### Data License

The radiosonde data produced by this software is made available under the [UK Open Government Licence (OGL)](http://www.nationalarchives.gov.uk/doc/open-government-licence/).

Acknowledgement of NCAS as the data provider is required whenever and wherever these data are used.
