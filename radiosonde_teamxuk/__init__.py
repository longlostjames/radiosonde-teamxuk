"""
Radiosonde TEAMxUK Processing Package

Python tools for processing Vaisala RS41-SGP radiosonde data from the TEAMxUK campaign,
converting EDT files to NCAS-AMOF compliant NetCDF format and generating quicklook plots.
"""

__version__ = "0.2.0"
__author__ = "Chris Walden"

from .process_radiosondes import (
    read_edt_file,
    save_to_ncas_netcdf,
    process_edt_files
)
from .generate_quicklooks import (
    plot_combined_analysis as create_quicklook,
    process_netcdf_files
)

__all__ = [
    'read_edt_file',
    'save_to_ncas_netcdf',
    'process_edt_files',
    'create_quicklook',
    'process_netcdf_files',
]
