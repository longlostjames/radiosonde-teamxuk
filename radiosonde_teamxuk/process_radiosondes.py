#!/usr/bin/env python3
"""
Process radiosonde files from TEAMxUK campaign
Usage: python process_radiosondes.py <input_directory> <output_directory>
"""

import os
import sys
import glob
from pathlib import Path
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from metpy.plots import SkewT, Hodograph
from metpy.units import units
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from eccodes import codes_bufr_new_from_file, codes_get, codes_set, codes_release
from netCDF4 import Dataset
import ncas_amof_netcdf_template as nant
try:
    from importlib.resources import files
except ImportError:
    from importlib_resources import files

# Import package version
try:
    from importlib.metadata import version
    PACKAGE_VERSION = version("radiosonde-teamxuk")
except Exception:
    try:
        from radiosonde_teamxuk import __version__ as PACKAGE_VERSION
    except ImportError:
        PACKAGE_VERSION = "0.2.0"  # Fallback to hardcoded version


def wind_components(speed, direction):
    """
    Calculate u and v wind components from speed and direction
    
    Parameters:
    -----------
    speed : array-like
        Wind speed in m/s
    direction : array-like
        Wind direction in degrees (meteorological convention)
        
    Returns:
    --------
    u, v : tuple of arrays
        u (east-west) and v (north-south) components
    """
    direction_rad = np.deg2rad(direction)
    u = -speed * np.sin(direction_rad)
    v = -speed * np.cos(direction_rad)
    return u, v


def read_rs41_bufr(filename):
    """
    Read RS41-SGP radiosonde data from BUFR file
    
    Parameters:
    -----------
    filename : str
        Path to BUFR file
        
    Returns:
    --------
    list of dict
        List containing profile dictionary with atmospheric data
    """
    profiles = []
    
    with open(filename, 'rb') as f:
        # Read the first message
        bufr = codes_bufr_new_from_file(f)
        
        if bufr is None:
            print(f"No BUFR messages found in {filename}")
            return profiles
        
        # Unpack the data
        codes_set(bufr, 'unpack', 1)
        
        # Get the number of levels
        try:
            num_levels = codes_get(bufr, 'numberOfSubsets')
        except:
            num_levels = 1
        
        # Initialize arrays
        pressure = []
        temperature = []
        humidity = []
        dewpoint = []
        height = []
        wind_speed = []
        wind_direction = []
        latitude_profile = []  # GPS latitude at each level (if available)
        longitude_profile = []  # GPS longitude at each level (if available)
        latitude = None
        longitude = None
        timestamp = None
        
        # Try to get station info (launch location)
        try:
            latitude = codes_get(bufr, '#1#latitude')
            longitude = codes_get(bufr, '#1#longitude')
        except:
            try:
                latitude = codes_get(bufr, 'latitude')
                longitude = codes_get(bufr, 'longitude')
            except:
                latitude = 0.0
                longitude = 0.0
        
        # Try to get timestamp
        try:
            year = codes_get(bufr, '#1#year')
            month = codes_get(bufr, '#1#month')
            day = codes_get(bufr, '#1#day')
            hour = codes_get(bufr, '#1#hour')
            minute = codes_get(bufr, '#1#minute')
            timestamp = f"{year:04d}-{month:02d}-{day:02d} {hour:02d}:{minute:02d}"
        except:
            timestamp = "Unknown"
        
        # Read level-by-level data using numbered keys
        level = 1
        while True:
            try:
                # Try to read pressure for this level
                p = codes_get(bufr, f'#{level}#pressure')
                
                # If pressure is missing, break
                if p == codes_get(bufr, 'missingValue') or p < -1e10:
                    break
                
                pressure.append(p / 100.0)  # Convert Pa to hPa
                
                # Temperature
                try:
                    t = codes_get(bufr, f'#{level}#airTemperature')
                    if t > -1e10:
                        temperature.append(t - 273.15)  # Convert K to C
                    else:
                        temperature.append(np.nan)
                except:
                    temperature.append(np.nan)
                
                # Humidity
                try:
                    rh = codes_get(bufr, f'#{level}#relativeHumidity')
                    if rh > -1e10:
                        humidity.append(rh)
                    else:
                        humidity.append(np.nan)
                except:
                    humidity.append(np.nan)
                
                # Dewpoint
                try:
                    td = codes_get(bufr, f'#{level}#dewpointTemperature')
                    if td > -1e10:
                        dewpoint.append(td - 273.15)  # Convert K to C
                    else:
                        dewpoint.append(np.nan)
                except:
                    dewpoint.append(np.nan)
                
                # Height
                try:
                    h = codes_get(bufr, f'#{level}#nonCoordinateGeopotentialHeight')
                    if h > -1e10:
                        height.append(h)
                    else:
                        height.append(np.nan)
                except:
                    height.append(np.nan)
                
                # Wind speed
                try:
                    ws = codes_get(bufr, f'#{level}#windSpeed')
                    if ws > -1e10:
                        wind_speed.append(ws)
                    else:
                        wind_speed.append(np.nan)
                except:
                    wind_speed.append(np.nan)
                
                # Wind direction
                try:
                    wd = codes_get(bufr, f'#{level}#windDirection')
                    if wd > -1e10:
                        wind_direction.append(wd)
                    else:
                        wind_direction.append(np.nan)
                except:
                    wind_direction.append(np.nan)
                
                # Try to get GPS position at this level (if available)
                try:
                    lat_level = codes_get(bufr, f'#{level}#latitude')
                    lon_level = codes_get(bufr, f'#{level}#longitude')
                    if lat_level > -1e10 and lon_level > -1e10:
                        latitude_profile.append(lat_level)
                        longitude_profile.append(lon_level)
                    else:
                        latitude_profile.append(np.nan)
                        longitude_profile.append(np.nan)
                except:
                    latitude_profile.append(np.nan)
                    longitude_profile.append(np.nan)
                
                level += 1
                
            except Exception as e:
                # End of data
                break
        
        # Release the message
        codes_release(bufr)
        
        # Create profile dictionary
        if len(pressure) > 0:
            # Check if we have GPS trajectory data
            has_gps_trajectory = len(latitude_profile) > 0 and not all(np.isnan(latitude_profile))
            
            profile = {
                'pressure': np.array(pressure),
                'temperature': np.array(temperature),
                'humidity': np.array(humidity),
                'dewpoint': np.array(dewpoint),
                'height': np.array(height),
                'wind_speed': np.array(wind_speed),
                'wind_direction': np.array(wind_direction),
                'latitude': latitude,  # Launch location
                'longitude': longitude,  # Launch location
                'latitude_profile': np.array(latitude_profile) if has_gps_trajectory else None,  # GPS track
                'longitude_profile': np.array(longitude_profile) if has_gps_trajectory else None,  # GPS track
                'timestamp': timestamp
            }
            profiles.append(profile)
    
    return profiles


def read_edt_file(filename):
    """
    Read Vaisala EDT (Extended Data) file
    
    Parameters:
    -----------
    filename : str
        Path to EDT file
        
    Returns:
    --------
    list of dict
        List containing profile dictionary with atmospheric data
    """
    profiles = []
    
    try:
        with open(filename, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()
        
        # Parse header for metadata
        latitude = None
        longitude = None
        timestamp = None
        sonde_serial = None
        hardware_version = None
        station_name = None
        mw41_software_version = None
        
        # Debug: Print lines that might contain hardware version
        debug_hw_search = False
        if debug_hw_search:
            print(f"  Searching for hardware version in first 80 lines:")
            for i, line in enumerate(lines[:80]):
                if 'hardware' in line.lower() or 'ground' in line.lower() or 'device' in line.lower():
                    print(f"    Line {i}: {line.strip()}")
        
        for line in lines[:80]:  # Check first 80 lines for metadata (hardware version is around line 67, MW41 version nearby)
            if 'Release point latitude' in line:
                # Extract latitude - handle various degree symbols
                parts = line.split()
                lat_found = False
                for part in reversed(parts):
                    # Remove any degree symbols and direction indicators
                    clean_part = part
                    for char in ['�', '°', '\xb0', 'N', 'S', 'n', 's']:
                        clean_part = clean_part.replace(char, '')
                    try:
                        latitude = float(clean_part)
                        if any(c in part.upper() for c in ['S']):
                            latitude = -latitude
                        lat_found = True
                        break
                    except:
                        continue
            elif 'Release point longitude' in line:
                # Extract longitude - handle various degree symbols  
                parts = line.split()
                lon_found = False
                for part in reversed(parts):
                    # Remove any degree symbols and direction indicators
                    clean_part = part
                    for char in ['�', '°', '\xb0', 'E', 'W', 'e', 'w']:
                        clean_part = clean_part.replace(char, '')
                    try:
                        longitude = float(clean_part)
                        if any(c in part.upper() for c in ['W']):
                            longitude = -longitude
                        lon_found = True
                        break
                    except:
                        continue
            elif 'Balloon release date' in line:
                date_str = line.split()[-1]  # e.g., "18/02/25"
            elif 'Balloon release time' in line:
                time_str = line.split()[-1]  # e.g., "23:07:28"
                try:
                    # Combine date and time with seconds
                    day, month, year = date_str.split('/')
                    hour, minute, second = time_str.split(':')
                    year = f"20{year}" if int(year) < 50 else f"19{year}"
                    timestamp = f"{year}-{month}-{day} {hour}:{minute}:{second}"
                except:
                    timestamp = None
            elif 'Sonde serial number' in line:
                sonde_serial = line.split()[-1]
            elif 'Ground check device hardware version' in line:
                # Extract hex value like 0x70 or 0x60
                hw_str = line.split()[-1]
                hardware_version = hw_str  # Store the hex string
            elif 'Software version' in line and 'MW41' in line:
                # Extract MW41 software version from line like "Software version                                    MW41 2.17.0"
                # The format is "MW41 X.Y.Z" so we need to get the version number after MW41
                parts = line.split()
                if 'MW41' in parts:
                    mw41_idx = parts.index('MW41')
                    if mw41_idx + 1 < len(parts):
                        mw41_software_version = parts[mw41_idx + 1]
                        print(f"  Debug: Found MW41 version line: {line.strip()}")
                        print(f"  Debug: Extracted MW41 version: {mw41_software_version}")
            elif 'Station name' in line:
                # Extract station name (e.g., "Sterzing")
                station_name = line.split()[-1].lower()  # Store lowercase
        
        # Find the data section (starts after header row with column names)
        data_start_idx = None
        for i, line in enumerate(lines):
            if line.strip().startswith('Elapsed time') and 'TimeUTC' in line:
                data_start_idx = i + 2  # Skip header and units line
                break
        
        if data_start_idx is None:
            print(f"No data section found in {filename}")
            return profiles
        
        # Read data rows
        pressure = []
        temperature = []
        humidity = []
        dewpoint = []
        height = []
        wind_speed = []
        wind_direction = []
        latitude_profile = []
        longitude_profile = []
        elapsed_time = []
        ascent_rate = []
        
        for line in lines[data_start_idx:]:
            line = line.strip()
            if not line or line.startswith('//'):
                continue
            
            parts = line.split()
            if len(parts) < 20:  # Need at least the main variables
                continue
            
            try:
                # Parse each column based on EDT format
                # 0: Elapsed time (s)
                # 1: TimeUTC
                # 2: P (hPa)
                # 3: Temp (°C)
                # 4: RH (%)
                # 5: Dewp (°C)
                # 6: Speed (m/s)
                # 7: Dir (°)
                # 8: Ecomp (m/s)
                # 9: Ncomp (m/s)
                # 10: Lat (°)
                # 11: Lon (°)
                # 12: AscRate (m/s)
                # 13: HeightMSL (m)
                
                elapsed_time.append(float(parts[0]))
                pressure.append(float(parts[2]) if parts[2] != '/////' else np.nan)
                temperature.append(float(parts[3]) if parts[3] != '/////' else np.nan)
                humidity.append(float(parts[4]) if parts[4] != '/////' else np.nan)
                dewpoint.append(float(parts[5]) if parts[5] != '/////' else np.nan)
                wind_speed.append(float(parts[6]) if parts[6] != '/////' else np.nan)
                wind_direction.append(float(parts[7]) if parts[7] != '/////' else np.nan)
                latitude_profile.append(float(parts[10]) if parts[10] != '/////' else np.nan)
                longitude_profile.append(float(parts[11]) if parts[11] != '/////' else np.nan)
                ascent_rate.append(float(parts[12]) if parts[12] != '/////' else np.nan)
                height.append(float(parts[13]) if parts[13] != '/////' else np.nan)
                
            except (ValueError, IndexError):
                continue
        
        # Create profile dictionary
        if len(pressure) > 0:
            profile = {
                'pressure': np.array(pressure),
                'temperature': np.array(temperature),
                'humidity': np.array(humidity),
                'dewpoint': np.array(dewpoint),
                'height': np.array(height),
                'wind_speed': np.array(wind_speed),
                'wind_direction': np.array(wind_direction),
                'latitude': latitude,  # Launch location
                'longitude': longitude,  # Launch location
                'latitude_profile': np.array(latitude_profile),  # GPS trajectory
                'longitude_profile': np.array(longitude_profile),  # GPS trajectory
                'ascent_rate': np.array(ascent_rate),  # Upward balloon velocity
                'elapsed_time': np.array(elapsed_time),  # Actual elapsed time from file
                'timestamp': timestamp,
                'sonde_serial': sonde_serial,
                'hardware_version': hardware_version,  # Ground check device hardware version
                'mw41_software_version': mw41_software_version,  # MW41 software version from EDT
                'station_name': station_name,  # Ground station name
                'source_file': os.path.basename(filename)  # Original EDT filename
            }
            profiles.append(profile)
    
    except Exception as e:
        print(f"Error reading EDT file {filename}: {str(e)}")
        import traceback
        traceback.print_exc()
    
    return profiles


def save_to_ncas_netcdf(profile, output_dir, instrument=None, metadata_dir=None):
    """
    Save radiosonde profile data to NCAS-AMOF NetCDF format
    
    Parameters:
    -----------
    profile : dict
        Dictionary containing radiosonde data
    output_dir : str
        Directory where output NetCDF file will be saved (uses NCAS naming convention)
    instrument : str, optional
        Instrument identifier. If None, determined from hardware_version in profile
        (0x70 -> ncas-radiosonde-1, 0x60 -> ncas-radiosonde-2, default: ncas-radiosonde-2)
    metadata_dir : str, optional
        Directory containing instrument-specific metadata files
        (metadata_ncas-radiosonde-1.json and metadata_ncas-radiosonde-2.json)
    """
    import os
    from netCDF4 import Dataset
    
    try:
        # Determine instrument from hardware version if not specified
        if instrument is None:
            hardware_version = profile.get('hardware_version')
            print(f"  Determining instrument from EDT metadata...")
            if hardware_version == '0x70':
                instrument = 'ncas-radiosonde-1'
                print(f"  ✓ Ground check device hardware version: {hardware_version}")
                print(f"  ✓ Mapped to instrument: {instrument}")
            elif hardware_version == '0x60':
                instrument = 'ncas-radiosonde-2'
                print(f"  ✓ Ground check device hardware version: {hardware_version}")
                print(f"  ✓ Mapped to instrument: {instrument}")
            else:
                # Default to ncas-radiosonde-2 if version not recognized
                instrument = 'ncas-radiosonde-2'
                if hardware_version:
                    print(f"  ⚠️  Ground check device hardware version: {hardware_version} (unknown)")
                    print(f"  ⚠️  Using default instrument: {instrument}")
                else:
                    print(f"  ⚠️  No hardware version found in EDT file")
                    print(f"  ⚠️  Using default instrument: {instrument}")
        else:
            print(f"  Using explicitly specified instrument: {instrument}")
        
        # Select metadata file based on instrument
        metadata_file = None
        if metadata_dir:
            metadata_file = os.path.join(metadata_dir, f'metadata_{instrument}.json')
            if os.path.exists(metadata_file):
                print(f"  ✓ Using metadata file: {os.path.basename(metadata_file)}")
            else:
                # Try to load from package resources
                try:
                    pkg_files = files('radiosonde_teamxuk')
                    pkg_metadata = pkg_files / f'metadata_{instrument}.json'
                    if pkg_metadata.is_file():
                        metadata_file = str(pkg_metadata)
                        print(f"  ✓ Using metadata file from package: {os.path.basename(metadata_file)}")
                    else:
                        print(f"  ⚠️  Metadata file not found: metadata_{instrument}.json")
                        metadata_file = None
                except:
                    print(f"  ⚠️  Metadata file not found: {os.path.basename(metadata_file)}")
                    metadata_file = None
        
        # Parse timestamp to get date and time with seconds
        timestamp_str = profile.get('timestamp', 'Unknown')
        if timestamp_str != "Unknown":
            try:
                # Try parsing with seconds first
                dt = datetime.strptime(timestamp_str, "%Y-%m-%d %H:%M:%S")
                date_str = dt.strftime("%Y%m%d")
                datetime_str = dt.strftime("%Y%m%d-%H%M%S")
            except:
                try:
                    # Fall back to without seconds
                    dt = datetime.strptime(timestamp_str, "%Y-%m-%d %H:%M")
                    date_str = dt.strftime("%Y%m%d")
                    datetime_str = dt.strftime("%Y%m%d-%H%M%S")
                except:
                    dt = datetime.now()
                    date_str = dt.strftime("%Y%m%d")
                    datetime_str = dt.strftime("%Y%m%d-%H%M%S")
        else:
            dt = datetime.now()
            date_str = dt.strftime("%Y%m%d")
            datetime_str = dt.strftime("%Y%m%d-%H%M%S")
        
        # Get platform/location name from profile or metadata for filename
        platform = profile.get('station_name')  # Try EDT station name first
        
        # Get project identifier from metadata (e.g., "teamxuk" from "TEAMxUK: ...")
        project_id = None
        if metadata_file and os.path.exists(metadata_file):
            try:
                import json
                with open(metadata_file, 'r') as f:
                    metadata = json.load(f)
                    if 'project' in metadata:
                        # Extract project ID from project string
                        # e.g., "TEAMxUK: Quantifying..." -> "teamxuk"
                        project_str = metadata['project']
                        if ':' in project_str:
                            project_id = project_str.split(':')[0].strip().lower()
                        else:
                            project_id = project_str.strip().lower()
                        print(f"  ✓ Extracted project ID: {project_id}")
                    if 'platform' in metadata and not platform:
                        platform = metadata['platform'].lower().replace(' ', '-')
            except Exception as e:
                print(f"  ⚠️  Error reading metadata: {e}")
                pass
        
        # If no platform in metadata, use the NCAS template default for sonde product
        if platform is None:
            platform = 'mobile'  # NCAS default for sonde product
        
        # Create output directory structure: output_dir/instrument/project-id/YYYY/MM/
        instrument_output_dir = os.path.join(output_dir, instrument)
        if project_id:
            # Extract year and month from datetime string (format: YYYYMMDD-HHMMSS)
            year = datetime_str[:4]
            month = datetime_str[4:6]
            instrument_output_dir = os.path.join(instrument_output_dir, project_id, year, month)
        os.makedirs(instrument_output_dir, exist_ok=True)
        
        # Number of measurements along the balloon trajectory
        n_measurements = len(profile['pressure'])
        
        # Create NCAS-AMOF NetCDF file using sonde product
        # The sonde product treats data as a time series along the balloon trajectory
        # This creates a file with NCAS naming convention automatically
        nc = nant.create_netcdf.main(
            instrument,
            date=date_str,
            dimension_lengths={'time': n_measurements},
            products='sonde',
            loc='land',
            file_location=instrument_output_dir,
            verbose=0
        )
        
        # Close the template file so we can modify its structure
        nc_path = nc.filepath()
        nc.close()
        
        # Rename to custom NCAS convention with location, project, and time
        # Format: ncas-radiosonde-2_sterzing_20250216-110000_sonde_teamx-uk_v1.0.nc
        nc_path_obj = Path(nc_path)
        print(f"  Creating filename with: instrument={instrument}, platform={platform}, datetime={datetime_str}, project_id={project_id}")
        if project_id:
            new_filename = f"{instrument}_{platform}_{datetime_str}_sonde_{project_id}_v1.0.nc"
        else:
            new_filename = f"{instrument}_{platform}_{datetime_str}_v1.0.nc"
        print(f"  Generated filename: {new_filename}")
        new_path = nc_path_obj.parent / new_filename
        
        # If a file with this name already exists, add a suffix
        if new_path.exists():
            counter = 1
            while new_path.exists():
                if project_id:
                    new_filename = f"{instrument}_{platform}_{datetime_str}_sonde_{project_id}_v1.0_{counter}.nc"
                else:
                    new_filename = f"{instrument}_{platform}_{datetime_str}_v1.0_{counter}.nc"
                new_path = nc_path_obj.parent / new_filename
                counter += 1
        
        os.rename(nc_path, new_path)
        nc_path = str(new_path)
        
        # Reopen with netCDF4 in write mode to modify lat/lon dimensions
        nc = Dataset(nc_path, 'r+')
        
        # Switch to define mode to allow structural changes
        nc.set_fill_off()
        
        # Remove template's latitude and longitude variables (they use wrong dimensions)
        if 'latitude' in nc.variables:
            # Save attributes before deletion
            lat_attrs = {attr: nc.variables['latitude'].getncattr(attr) 
                        for attr in nc.variables['latitude'].ncattrs()}
            # Must be in define mode to delete
            nc.renameVariable('latitude', 'latitude_old_temp')
        else:
            lat_attrs = {}
            
        if 'longitude' in nc.variables:
            lon_attrs = {attr: nc.variables['longitude'].getncattr(attr) 
                        for attr in nc.variables['longitude'].ncattrs()}
            nc.renameVariable('longitude', 'longitude_old_temp')
        else:
            lon_attrs = {}
        
        # Create new latitude and longitude variables with time dimension
        lat_var = nc.createVariable('latitude', 'f4', ('time',), fill_value=-999.0)
        for attr, val in lat_attrs.items():
            if attr not in ['_FillValue']:
                try:
                    lat_var.setncattr(attr, val)
                except:
                    pass
        # Set required attributes
        lat_var.standard_name = 'latitude'
        lat_var.long_name = 'Latitude of radiosonde along trajectory'
        lat_var.units = 'degree_north'
        lat_var.axis = 'Y'
        lat_var.valid_min = np.float32(-90.0)
        lat_var.valid_max = np.float32(90.0)
        lat_var.coordinates = 'time latitude longitude altitude'
        
        lon_var = nc.createVariable('longitude', 'f4', ('time',), fill_value=-999.0)
        for attr, val in lon_attrs.items():
            if attr not in ['_FillValue']:
                try:
                    lon_var.setncattr(attr, val)
                except:
                    pass
        lon_var.standard_name = 'longitude'
        lon_var.long_name = 'Longitude of radiosonde along trajectory'
        lon_var.units = 'degree_east'
        lon_var.axis = 'X'
        lon_var.valid_min = np.float32(-180.0)
        lon_var.valid_max = np.float32(180.0)
        lon_var.coordinates = 'time latitude longitude altitude'
        
        # Create time array - estimate time for each measurement
        # Use elapsed time from file if available (EDT files), otherwise estimate from height
        if profile.get('elapsed_time') is not None:
            elapsed_seconds = profile['elapsed_time']
        else:
            # Estimate from ascent rate
            ascent_rate = 5.0  # m/s
            elapsed_seconds = profile['height'] / ascent_rate
            elapsed_seconds[np.isnan(elapsed_seconds)] = 0
        
        # Ensure elapsed_seconds has the same length as n_measurements
        if len(elapsed_seconds) != n_measurements:
            print(f"  ⚠️  Warning: elapsed_seconds length ({len(elapsed_seconds)}) doesn't match n_measurements ({n_measurements})")
            print(f"  Truncating to shortest length")
            min_len = min(len(elapsed_seconds), n_measurements)
            elapsed_seconds = elapsed_seconds[:min_len]
            # Also truncate other arrays
            for key in ['pressure', 'temperature', 'humidity', 'height', 'wind_speed', 
                       'wind_direction', 'latitude_profile', 'longitude_profile']:
                if key in profile and profile[key] is not None:
                    profile[key] = profile[key][:min_len]
            n_measurements = min_len
        
        # Create datetime objects for each measurement point
        times = [dt + np.timedelta64(int(s), 's') for s in elapsed_seconds]
        
        # Get times in required formats
        unix_times, day_of_year, years, months, days, hours, minutes, seconds, \
            time_coverage_start_unix, time_coverage_end_unix, file_date = nant.util.get_times(times)
        
        # Update coordinate variables
        altitude_masked = np.ma.masked_invalid(profile['height'])
        nant.util.update_variable(nc, 'altitude', altitude_masked)
        nant.util.update_variable(nc, 'time', unix_times)
        nant.util.update_variable(nc, 'day_of_year', day_of_year)
        nant.util.update_variable(nc, 'year', years)
        nant.util.update_variable(nc, 'month', months)
        nant.util.update_variable(nc, 'day', days)
        nant.util.update_variable(nc, 'hour', hours)
        nant.util.update_variable(nc, 'minute', minutes)
        nant.util.update_variable(nc, 'second', seconds)
        
        # Update location variables for trajectory data
        # Use GPS trajectory if available (EDT files), otherwise use launch location for all points
        if 'latitude_profile' in profile and 'longitude_profile' in profile:
            # Use GPS trajectory data
            lat_trajectory = np.ma.masked_invalid(profile['latitude_profile'])
            lon_trajectory = np.ma.masked_invalid(profile['longitude_profile'])
        else:
            # Use fixed launch location for all points (BUFR files)
            lat_trajectory = np.full(n_measurements, profile['latitude'])
            lon_trajectory = np.full(n_measurements, profile['longitude'])
        
        # Update latitude and longitude with trajectory data
        nc.variables['latitude'][:] = lat_trajectory
        nc.variables['longitude'][:] = lon_trajectory
        
        # Update atmospheric variables with proper masking for NaN values
        # Air pressure (convert hPa to Pa)
        pressure_pa = profile['pressure'] * 100
        pressure_pa = np.ma.masked_invalid(pressure_pa)
        nant.util.update_variable(nc, 'air_pressure', pressure_pa)
        
        # Air temperature (convert C to K)
        temp_k = profile['temperature'] + 273.15
        temp_k = np.ma.masked_invalid(temp_k)
        nant.util.update_variable(nc, 'air_temperature', temp_k)
        
        # Relative humidity
        rh = np.ma.masked_invalid(profile['humidity'])
        nant.util.update_variable(nc, 'relative_humidity', rh)
        
        # Wind speed
        wspd = np.ma.masked_invalid(profile['wind_speed'])
        nant.util.update_variable(nc, 'wind_speed', wspd)
        
        # Wind direction
        wdir = np.ma.masked_invalid(profile['wind_direction'])
        nant.util.update_variable(nc, 'wind_from_direction', wdir)
        
        # Elapsed time since launch (in seconds)
        elapsed_time = np.ma.masked_invalid(elapsed_seconds)
        nant.util.update_variable(nc, 'elapsed_time', elapsed_time)
        
        # Upward balloon velocity (ascent rate)
        if 'ascent_rate' in profile:
            ascent_rate = np.ma.masked_invalid(profile['ascent_rate'])
            nant.util.update_variable(nc, 'upward_balloon_velocity', ascent_rate)
        
        # Add metadata from file if provided (but skip variables that might cause issues)
        if metadata_file and os.path.exists(metadata_file):
            try:
                nant.util.add_metadata_to_netcdf(nc, metadata_file)
            except Exception as e:
                # Silently handle metadata update issues (usually from renamed old lat/lon variables)
                # Manually add the important metadata as global attributes
                import json
                with open(metadata_file, 'r') as f:
                    metadata = json.load(f)
                for key, value in metadata.items():
                    if key not in nc.variables:  # Only set global attributes
                        try:
                            nc.setncattr(key, value)
                        except:
                            pass
        
        # Override instrument_software_version with value from EDT file if available
        if 'mw41_software_version' in profile and profile['mw41_software_version']:
            nc.setncattr('instrument_software', 'MW41')
            nc.setncattr('instrument_software_version', profile['mw41_software_version'])
            print(f"  ✓ Set MW41 software version from EDT: {profile['mw41_software_version']}")
        else:
            print(f"  Debug: MW41 version not found in profile or is None/empty")
            print(f"  Debug: profile keys: {list(profile.keys())}")
            if 'mw41_software_version' in profile:
                print(f"  Debug: mw41_software_version value: {repr(profile['mw41_software_version'])}")
        
        # Set processing software version from package
        nc.setncattr('processing_software', 'radiosonde-teamxuk')
        nc.setncattr('processing_software_version', PACKAGE_VERSION)
        nc.setncattr('processing_software_url', 'https://github.com/longlostjames/radiosonde-teamxuk')
        
        # Append source file information to existing comment
        source_file = profile.get('source_file', 'unknown')
        existing_comment = nc.getncattr('comment') if 'comment' in nc.ncattrs() else ''
        if existing_comment:
            nc.setncattr('comment', f'{existing_comment}. Radiosonde profile data converted from EDT file: {source_file}')
        else:
            nc.setncattr('comment', f'Radiosonde profile data converted from EDT file: {source_file}')
        
        # Set time coverage
        nc.setncattr('time_coverage_start',
                     datetime.fromtimestamp(time_coverage_start_unix, datetime.now().astimezone().tzinfo).strftime("%Y-%m-%dT%H:%M:%S"))
        nc.setncattr('time_coverage_end',
                     datetime.fromtimestamp(time_coverage_end_unix, datetime.now().astimezone().tzinfo).strftime("%Y-%m-%dT%H:%M:%S"))
        
        # Set geospatial bounds based on trajectory (SW corner, NE corner format)
        if 'latitude_profile' in profile and 'longitude_profile' in profile:
            valid_lats = profile['latitude_profile'][~np.isnan(profile['latitude_profile'])]
            valid_lons = profile['longitude_profile'][~np.isnan(profile['longitude_profile'])]
            if len(valid_lats) > 0 and len(valid_lons) > 0:
                lat_min, lat_max = valid_lats.min(), valid_lats.max()
                lon_min, lon_max = valid_lons.min(), valid_lons.max()
                # Format: "SW_lat SW_lon, NE_lat NE_lon"
                geobounds = f"{lat_min:.6f}N {lon_min:.6f}E, {lat_max:.6f}N {lon_max:.6f}E"
            else:
                geobounds = f"{profile['latitude']:.6f}N {profile['longitude']:.6f}E"
        else:
            geobounds = f"{profile['latitude']:.6f}N {profile['longitude']:.6f}E"
        nc.setncattr('geospatial_bounds', geobounds)
        
        # Set instrument serial number from EDT file if available
        if 'sonde_serial' in profile and profile['sonde_serial']:
            nc.setncattr('instrument_serial_number', profile['sonde_serial'])
        
        # Close the file
        nc.close()
        
        # Clean up: remove temporary variables and unused dimensions
        # We'll copy to a clean file excluding the unwanted variables/dimensions
        temp_path = nc_path + '.tmp'
        
        try:
            with Dataset(nc_path, 'r') as src:
                with Dataset(temp_path, 'w', format='NETCDF4') as dst:
                    # Copy global attributes
                    dst.setncatts({attr: src.getncattr(attr) for attr in src.ncattrs()})
                    
                    # Copy dimensions (exclude the singleton lat/lon dimensions)
                    for name, dimension in src.dimensions.items():
                        if name not in ['latitude', 'longitude']:
                            dst.createDimension(
                                name, (len(dimension) if not dimension.isunlimited() else None))
                    
                    # Copy variables (exclude the temporary renamed ones)
                    exclude_vars = ['latitude_old_temp', 'longitude_old_temp']
                    for name, variable in src.variables.items():
                        if name not in exclude_vars:
                            # Create variable
                            dst_var = dst.createVariable(
                                name, variable.datatype, variable.dimensions,
                                fill_value=variable._FillValue if hasattr(variable, '_FillValue') else None
                            )
                            # Copy attributes, excluding problematic ones
                            for attr in variable.ncattrs():
                                if attr not in ['_FillValue', 'valid_min', 'valid_max']:
                                    try:
                                        dst_var.setncattr(attr, variable.getncattr(attr))
                                    except:
                                        pass
                            # Copy valid_min and valid_max if they're compatible with data type
                            for attr in ['valid_min', 'valid_max']:
                                if hasattr(variable, attr):
                                    try:
                                        val = variable.getncattr(attr)
                                        # Only set if compatible with variable dtype
                                        if variable.dtype.kind in ['f', 'i', 'u']:  # numeric types
                                            dst_var.setncattr(attr, val)
                                    except:
                                        pass
                            # Copy data (suppress warnings about valid_min/max)
                            import warnings
                            with warnings.catch_warnings():
                                warnings.filterwarnings('ignore', category=UserWarning, message='.*valid_min.*')
                                warnings.filterwarnings('ignore', category=UserWarning, message='.*valid_max.*')
                                dst_var[:] = variable[:]
            
            # Replace original file with cleaned version
            import os
            os.replace(temp_path, nc_path)
        except Exception as e:
            print(f"Warning: Could not clean up temporary variables: {e}")
            # Remove temp file if it exists
            if os.path.exists(temp_path):
                os.remove(temp_path)
        
        # Find the created file and use NCAS naming convention
        # The template creates files like: ncas-radiosonde-2_cao_20250211_sonde_v1.0.nc
        pattern = f"{instrument}_*_{date_str}_*.nc"
        created_files = list(Path(output_dir).glob(pattern))
        
        if created_files:
            # Return the path to the created NetCDF file
            return str(Path(nc_path).name)
        
        return str(Path(nc_path).name)
        
    except Exception as e:
        print(f"Error creating NetCDF file: {str(e)}")
        import traceback
        traceback.print_exc()
        return None


def plot_combined_analysis(profiles, filename=None):
    """
    Combined plot: Skew-T on left, hodograph and map stacked on right
    
    Parameters:
    -----------
    profiles : list
        List of radiosonde profiles
    filename : str, optional
        Name of the source file to display in plot
    """
    if not profiles or len(profiles) == 0:
        print("No profiles to plot")
        return None
    
    profile = profiles[0]
    lat = profile['latitude']
    lon = profile['longitude']
    
    # Create figure
    fig = plt.figure(figsize=(20, 12))
    
    # ============ SKEW-T (Left side) ============
    # Use rect parameter for proper skew projection - (left, bottom, width, height)
    skew = SkewT(fig, rotation=45, rect=(0.06, 0.06, 0.42, 0.87))
    
    # Plot Skew-T data
    for i, prof in enumerate(profiles):
        mask = ~(np.isnan(prof['temperature']) | np.isnan(prof['pressure']) | np.isnan(prof['dewpoint']))
        
        if np.sum(mask) > 0:
            p = prof['pressure'][mask] * units.hPa
            T = prof['temperature'][mask] * units.degC
            Td = prof['dewpoint'][mask] * units.degC
            
            skew.plot(p, T, color='#DE8F05', linewidth=2.5, label="Temperature")
            skew.plot(p, Td, color='#0173B2', linewidth=2.5, label="Dewpoint")
            
            # Wind barbs
            wind_mask = ~(np.isnan(prof['wind_speed']) | np.isnan(prof['wind_direction']) | np.isnan(prof['pressure']))
            if np.sum(wind_mask) > 10:
                wind_subsample = wind_mask.copy()
                indices = np.where(wind_mask)[0]
                keep_indices = indices[::40]
                wind_subsample[:] = False
                wind_subsample[keep_indices] = True
                
                p_wind = prof['pressure'][wind_subsample] * units.hPa
                u, v = wind_components(prof['wind_speed'][wind_subsample], 
                                      prof['wind_direction'][wind_subsample])
                skew.plot_barbs(p_wind, u * units('m/s'), v * units('m/s'))
    
    skew.plot_dry_adiabats(alpha=0.25, color='#CC78BC')
    skew.plot_moist_adiabats(alpha=0.25, color='#029E73')
    skew.plot_mixing_lines(alpha=0.25, color='#949494')
    
    skew.ax.set_ylim(1000, 100)
    skew.ax.set_xlim(-40, 40)
    skew.ax.set_xlabel('Temperature (°C)', fontsize=12, fontweight='bold')
    skew.ax.set_ylabel('Pressure (hPa)', fontsize=12, fontweight='bold')
    skew.ax.set_title('Skew-T Log-P Diagram', fontsize=13, fontweight='bold', pad=10)
    skew.ax.legend(loc='upper left', fontsize=11)
    
    # ============ HODOGRAPH (Top right) ============
    hodo_ax = fig.add_axes([0.55, 0.52, 0.40, 0.40])
    h = Hodograph(hodo_ax, component_range=60)
    h.add_grid(increment=10)
    
    wind_mask = ~(np.isnan(profile['wind_speed']) | np.isnan(profile['wind_direction']) | np.isnan(profile['height']))
    
    if np.sum(wind_mask) > 10:
        u, v = wind_components(profile['wind_speed'][wind_mask], 
                              profile['wind_direction'][wind_mask])
        heights = profile['height'][wind_mask]
        
        u = u * units('m/s')
        v = v * units('m/s')
        
        # Color gradient based on height (consistent with map)
        max_height = max(heights)
        colors = plt.cm.plasma(heights / max_height)
        for j in range(len(u) - 1):
            h.ax.plot([u[j].magnitude, u[j+1].magnitude], 
                     [v[j].magnitude, v[j+1].magnitude], 
                     color=colors[j], linewidth=2.5, alpha=0.8)
        
        # Height markers
        height_markers = [0, 1000, 3000, 6000, 9000, 12000]
        for hm in height_markers:
            idx = np.argmin(np.abs(heights - hm))
            if np.abs(heights[idx] - hm) < 500:
                h.ax.plot(u[idx].magnitude, v[idx].magnitude, 'ko', markersize=7)
                h.ax.text(u[idx].magnitude + 1.5, v[idx].magnitude + 1.5, 
                        f'{int(heights[idx]/1000)}km', 
                        fontsize=9, fontweight='bold')
        
        # Colorbar
        sm = plt.cm.ScalarMappable(cmap=plt.cm.plasma, 
                                   norm=plt.Normalize(vmin=0, vmax=max(heights)/1000))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=hodo_ax, pad=0.05, fraction=0.046, aspect=15)
        cbar.set_label('Height (km)', fontsize=11, fontweight='bold')
    
    hodo_ax.set_title('Hodograph', fontsize=13, fontweight='bold', pad=10)
    
    # ============ MAP (Bottom right) ============
    if np.sum(wind_mask) < 10:
        print("Not enough wind data for map")
        return fig
    
    u_vals, v_vals = wind_components(profile['wind_speed'][wind_mask], 
                                     profile['wind_direction'][wind_mask])
    heights = profile['height'][wind_mask]
    
    # Calculate drift trajectory
    ascent_rate = 5.0
    dt = np.diff(heights) / ascent_rate
    dt = np.concatenate([[0], dt])
    
    cumulative_u = np.cumsum(u_vals * dt)
    cumulative_v = np.cumsum(v_vals * dt)
    
    meters_per_deg_lat = 111000
    meters_per_deg_lon = 111000 * np.cos(np.radians(lat))
    
    drift_lons = lon + cumulative_u / meters_per_deg_lon
    drift_lats = lat + cumulative_v / meters_per_deg_lat
    
    # Map extent
    lon_padding = max(0.3, (drift_lons.max() - drift_lons.min()) * 0.2)
    lat_padding = max(0.3, (drift_lats.max() - drift_lats.min()) * 0.2)
    lon_min = min(drift_lons.min(), lon) - lon_padding
    lon_max = max(drift_lons.max(), lon) + lon_padding
    lat_min = min(drift_lats.min(), lat) - lat_padding
    lat_max = max(drift_lats.max(), lat) + lat_padding
    map_extent = [lon_min, lon_max, lat_min, lat_max]
    
    map_ax = fig.add_axes([0.55, 0.06, 0.40, 0.40], projection=ccrs.PlateCarree())
    map_ax.set_extent(map_extent, crs=ccrs.PlateCarree())
    map_ax.set_facecolor('#f0f0f0')
    
    # Map features
    map_ax.add_feature(cfeature.LAND, facecolor='#e8e8e8', edgecolor='none', zorder=1)
    map_ax.add_feature(cfeature.OCEAN, facecolor='#d4e7f5', edgecolor='none', zorder=1)
    map_ax.add_feature(cfeature.COASTLINE, linewidth=1.2, edgecolor='#666666', zorder=2)
    map_ax.add_feature(cfeature.BORDERS, linewidth=0.8, linestyle='--', edgecolor='#999999', zorder=2)
    
    # Gridlines
    gl = map_ax.gridlines(draw_labels=True, linewidth=0.8, color='gray', 
                          alpha=0.4, linestyle=':', zorder=2)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 10}
    gl.ylabel_style = {'size': 10}
    
    # Trajectory
    map_ax.plot(drift_lons, drift_lats, 
               color='#555555', linewidth=4, alpha=0.3,
               transform=ccrs.PlateCarree(), zorder=3,
               solid_capstyle='round')
    
    max_height = max(heights)
    colors = plt.cm.plasma(heights / max_height)
    for j in range(len(drift_lons) - 1):
        map_ax.plot([drift_lons[j], drift_lons[j+1]], 
                   [drift_lats[j], drift_lats[j+1]], 
                   color=colors[j], linewidth=3.5, alpha=0.9,
                   transform=ccrs.PlateCarree(), zorder=4,
                   solid_capstyle='round')
    
    # Launch site
    map_ax.plot(lon, lat, marker='*', color='#ff3333', markersize=25, 
                transform=ccrs.PlateCarree(), zorder=10, 
                markeredgecolor='white', markeredgewidth=2,
                label='Launch Site')
    
    # Altitude markers
    height_markers = [0, 3000, 6000, 9000, 12000, 15000]
    for hm in height_markers:
        idx = np.argmin(np.abs(heights - hm))
        if idx < len(drift_lons) and np.abs(heights[idx] - hm) < 800:
            map_ax.plot(drift_lons[idx], drift_lats[idx], 'o', 
                       color=colors[idx], markersize=10,
                       markeredgecolor='white', markeredgewidth=2,
                       transform=ccrs.PlateCarree(), zorder=8)
            
            map_ax.text(drift_lons[idx], drift_lats[idx], 
                       f'  {int(heights[idx]/1000)}km  ', 
                       fontsize=10, fontweight='bold',
                       ha='left', va='bottom',
                       bbox=dict(boxstyle='round,pad=0.3', 
                                facecolor='white', 
                                edgecolor=colors[idx],
                                linewidth=1.5,
                                alpha=0.95),
                       transform=ccrs.PlateCarree(), zorder=9)
    
    # Colorbar for map
    sm_map = plt.cm.ScalarMappable(cmap=plt.cm.plasma, 
                                   norm=plt.Normalize(vmin=0, vmax=max(heights)/1000))
    sm_map.set_array([])
    cbar_map = plt.colorbar(sm_map, ax=map_ax, pad=0.02, fraction=0.046, aspect=15)
    cbar_map.set_label('Altitude (km)', fontsize=11, fontweight='bold')
    
    # Calculate drift
    total_drift = np.sqrt((drift_lons[-1] - lon)**2 + (drift_lats[-1] - lat)**2) * 111
    
    map_ax.set_title(f'Wind Trajectory (Total Drift: {total_drift:.1f} km)', 
                     fontsize=13, fontweight='bold', pad=10)
    map_ax.legend(loc='upper left', fontsize=10, framealpha=0.95, 
                 edgecolor='gray', fancybox=True)
    
    # Overall title
    timestamp = profile['timestamp'] if profile['timestamp'] else ''
    title = f'RS41-SGP Radiosonde Analysis: {timestamp}  (Launch: {lat:.3f}°N, {lon:.3f}°E)'
    if filename:
        title += f'\n{filename}'
    fig.suptitle(title, fontsize=16, fontweight='bold')
    
    return fig


def process_edt_files(input_dir, output_dir, metadata_dir=None):
    """
    Process all EDT files in input directory and save NetCDF files to output directory
    
    Parameters:
    -----------
    input_dir : str
        Directory containing EDT files
    output_dir : str
        Directory to save output NetCDF files
    metadata_dir : str, optional
        Directory containing instrument-specific metadata files for NetCDF generation
    """
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Metadata files will be selected per instrument based on hardware version
    # Look for instrument-specific metadata files in common locations
    if metadata_dir is None:
        # Try package directory first
        try:
            pkg_files = files('radiosonde_teamxuk')
            if (pkg_files / 'metadata_ncas-radiosonde-1.json').is_file() or \
               (pkg_files / 'metadata_ncas-radiosonde-2.json').is_file():
                metadata_dir = str(pkg_files)
                print(f"Found metadata files in package: {metadata_dir}")
        except:
            pass
        
        # Fall back to searching in current directory and input directory
        if metadata_dir is None:
            for search_dir in ['.', input_dir]:
                if os.path.exists(os.path.join(search_dir, 'metadata_ncas-radiosonde-1.json')) or \
                   os.path.exists(os.path.join(search_dir, 'metadata_ncas-radiosonde-2.json')):
                    metadata_dir = os.path.abspath(search_dir)
                    print(f"Found metadata files in: {metadata_dir}")
                    break
    
    # Find all EDT files
    edt_files = glob.glob(os.path.join(input_dir, '*.txt'))
    edt_files += glob.glob(os.path.join(input_dir, '*.edt'))
    
    # Filter EDT files to only include those starting with "edt1sdataforv217"
    edt_files = [f for f in edt_files if os.path.basename(f).startswith('edt1sdataforv217')]
    
    all_files = [(f, 'edt') for f in edt_files]
    all_files.sort()
    
    if not all_files:
        print(f"No EDT files (starting with 'edt1sdataforv217') found in {input_dir}")
        return
    
    print(f"Found {len(edt_files)} EDT files (edt1sdataforv217*)")
    
    # Process each file
    for i, (data_file, file_type) in enumerate(all_files, 1):
        filename = os.path.basename(data_file)
        print(f"\n[{i}/{len(all_files)}] Processing {filename}...")
        
        try:
            # Read EDT file
            profiles = read_edt_file(data_file)
            
            if not profiles:
                print(f"  ⚠️  No data found in {filename}")
                continue
            
            # Save NetCDF file for each profile (uses NCAS naming convention from template)
            nc_filenames = []
            for j, profile in enumerate(profiles):
                # Template generates NCAS-compliant filenames like:
                # ncas-radiosonde-2_sterzing_20250218-230728_v1.0.nc
                nc_filename = save_to_ncas_netcdf(profile, output_dir, metadata_dir=metadata_dir)
                if nc_filename:
                    nc_filenames.append(nc_filename)
                    print(f"  ✓ Saved NetCDF: {nc_filename}")
                else:
                    print(f"  ⚠️  Failed to save NetCDF")
            
        except Exception as e:
            print(f"  ✗ Error processing {filename}: {str(e)}")
            continue
    
    print(f"\n✓ Processing complete! Output saved to {output_dir}")


def main():
    """Main entry point for script"""
    # Check for --version flag
    if len(sys.argv) == 2 and sys.argv[1] in ['--version', '-v']:
        print(f"radiosonde-teamxuk version {PACKAGE_VERSION}")
        sys.exit(0)
    
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print("Usage: process-teamxuk-radiosondes <input_directory> <output_directory> [metadata_directory]")
        print("\nExample:")
        print("  process-teamxuk-radiosondes ./edt_files ./")
        print("\nThe script will automatically select metadata files based on instrument:")
        print("  - metadata_ncas-radiosonde-1.json for ncas-radiosonde-1")
        print("  - metadata_ncas-radiosonde-2.json for ncas-radiosonde-2")
        print("\nIf metadata_directory is not provided, the script will search in")
        print("the current directory and input directory.")
        print("\nOptions:")
        print("  --version, -v  Show version number")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    metadata_dir = sys.argv[3] if len(sys.argv) == 4 else None
    
    if not os.path.isdir(input_dir):
        print(f"Error: Input directory '{input_dir}' does not exist")
        sys.exit(1)
    
    if metadata_dir and not os.path.isdir(metadata_dir):
        print(f"Warning: Metadata directory '{metadata_dir}' does not exist")
        metadata_dir = None
    
    process_edt_files(input_dir, output_dir, metadata_dir)


if __name__ == "__main__":
    main()
