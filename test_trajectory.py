#!/usr/bin/env python
"""Test script to verify trajectory-based lat/lon"""
import sys
from radiosonde_teamxuk.process_radiosondes import read_edt_file, save_to_ncas_netcdf
from netCDF4 import Dataset

# Read EDT file
edt_file = '/gws/pw/j07/ncas_obs_vol1/amf/raw_data/ncas-radiosonde-2/incoming/teamx/edt_files/edt1sdata_20250211_2300.txt'
profiles = read_edt_file(edt_file)

if profiles:
    profile = profiles[0]
    print(f"Profile has {len(profile['latitude_profile'])} GPS points")
    print(f"Latitude range: {profile['latitude_profile'].min():.4f} to {profile['latitude_profile'].max():.4f}")
    print(f"Longitude range: {profile['longitude_profile'].min():.4f} to {profile['longitude_profile'].max():.4f}")
    
    # Save to NetCDF
    output_file = 'test_edt_output/test_trajectory.nc'
    save_to_ncas_netcdf(profile, output_file, metadata_file='metadata.json')
    print(f"\nSaved to {output_file}")
    
    # Check the NetCDF file
    nc = Dataset(output_file)
    print(f"\nNetCDF structure:")
    print(f"  latitude dimensions: {nc.variables['latitude'].dimensions}")
    print(f"  latitude shape: {nc.variables['latitude'].shape}")
    print(f"  longitude dimensions: {nc.variables['longitude'].dimensions}")
    print(f"  longitude shape: {nc.variables['longitude'].shape}")
    print(f"  Latitude range in file: {nc.variables['latitude'][:].min():.4f} to {nc.variables['latitude'][:].max():.4f}")
    print(f"  Longitude range in file: {nc.variables['longitude'][:].min():.4f} to {nc.variables['longitude'][:].max():.4f}")
    print(f"  geospatial_bounds: {nc.geospatial_bounds}")
    nc.close()
else:
    print("Failed to read EDT file")
