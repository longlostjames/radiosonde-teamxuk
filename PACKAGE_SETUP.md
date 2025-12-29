# Package Installation Summary

Your radiosonde-teamxuk project is now pip installable! 🎉

## What Was Done

### 1. Created Package Structure
- Moved Python scripts into `radiosonde_teamxuk/` package directory
- Moved JSON metadata files into the package
- Created `__init__.py` with proper imports and exports

### 2. Created Installation Files

#### `pyproject.toml` (Modern Python packaging)
- Package metadata (name, version, description, authors)
- Dependencies automatically installed
- Command-line entry points:
  - `process-teamxuk-radiosondes` → runs process_radiosondes.py
  - `generate-teamxuk-quicklooks` → runs generate_quicklooks.py

#### `requirements.txt`
- Lists all dependencies with version constraints
- Can be used separately if needed

#### `MANIFEST.in`
- Ensures JSON metadata files are included in distribution

#### `INSTALL.md`
- Detailed installation instructions
- Usage examples for both CLI and Python import

### 3. Updated Code
- Added importlib.resources support to find metadata files in the package
- Updated README with installation instructions
- Updated test_trajectory.py to import from the package

## How to Install

### For Users

**From source directory:**
```bash
cd /path/to/radiosonde-teamxuk
pip install .
```

**From GitHub (when pushed):**
```bash
pip install git+https://github.com/longlostjames/radiosonde-teamxuk.git
```

### For Development

```bash
cd /path/to/radiosonde-teamxuk
pip install -e .
```

The `-e` (editable) flag means changes to the code are immediately reflected without reinstalling.

## How to Use

### Command-Line Tools

After installation, you can run:

```bash
# Process EDT files to NetCDF
process-teamxuk-radiosondes <input_directory> <output_directory> [metadata_directory]

# Generate quicklook plots
generate-teamxuk-quicklooks <input_directory> [--stability]
```

**Note:** If you see warnings about commands not being on PATH, add `~/.local/bin` to your PATH:
```bash
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

### Python Library

```python
import radiosonde_teamxuk

# Read and process EDT files
profile = radiosonde_teamxuk.read_edt_file('path/to/file.edt')
output_file = radiosonde_teamxuk.save_to_ncas_netcdf(profile, './output')

# Create quicklook plot
radiosonde_teamxuk.create_quicklook(profile, stability=True)

# Process entire directories
radiosonde_teamxuk.process_edt_files('./input', './output')
radiosonde_teamxuk.process_netcdf_files('./netcdf', show_stability=True)
```

## Package Contents

```
radiosonde-teamxuk/
├── radiosonde_teamxuk/           # Main package
│   ├── __init__.py               # Package initialization
│   ├── process_radiosondes.py    # EDT to NetCDF conversion
│   ├── generate_quicklooks.py    # Quicklook plot generation
│   ├── metadata_ncas-radiosonde-1.json
│   ├── metadata_ncas-radiosonde-2.json
│   └── AMF_product_sonde_variable.json
├── pyproject.toml                # Package configuration
├── requirements.txt              # Dependencies
├── MANIFEST.in                   # Include data files
├── README.md                     # Documentation
├── INSTALL.md                    # Installation guide
└── LICENSE                       # MIT License
```

## Verification

Run the test script to verify installation:
```bash
python test_install.py
```

Expected output:
```
✓ Package imported successfully
  Version: 0.1.0
✓ process_radiosondes functions imported
✓ generate_quicklooks functions imported
✓ Found metadata file: metadata_ncas-radiosonde-1.json
✓ Found metadata file: metadata_ncas-radiosonde-2.json
✓ Found metadata file: AMF_product_sonde_variable.json

✓ All imports successful!
```

## Next Steps

1. **Test the installation** with your actual data
2. **Push to GitHub** to enable installation via `pip install git+...`
3. **Optional:** Publish to PyPI for `pip install radiosonde-teamxuk`
4. **Optional:** Add GitHub Actions for automated testing

## Publishing to PyPI (Optional)

To make it installable via `pip install radiosonde-teamxuk`:

1. Create accounts on PyPI and TestPyPI
2. Build the package:
   ```bash
   pip install build twine
   python -m build
   ```
3. Upload to TestPyPI first (test):
   ```bash
   python -m twine upload --repository testpypi dist/*
   ```
4. Upload to PyPI (production):
   ```bash
   python -m twine upload dist/*
   ```

## Troubleshooting

- **Import errors:** Reinstall with `pip install -e .`
- **Command not found:** Add `~/.local/bin` to PATH
- **Missing dependencies:** Run `pip install -r requirements.txt`
- **Metadata files not found:** Check that MANIFEST.in is included in build
