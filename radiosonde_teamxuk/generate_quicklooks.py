#!/usr/bin/env python3
"""
Generate quicklook plots from NetCDF radiosonde files
Usage: python generate_quicklooks.py <input_directory> <output_directory>
"""

import os
import sys
import glob
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from metpy.plots import SkewT, Hodograph
from metpy.units import units
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
from netCDF4 import Dataset


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


def read_netcdf_profile(nc_file):
    """
    Read radiosonde profile from NetCDF file
    
    Parameters:
    -----------
    nc_file : str
        Path to NetCDF file
        
    Returns:
    --------
    dict
        Profile dictionary with atmospheric data
    """
    try:
        with Dataset(nc_file, 'r') as nc:
            # Read variables
            pressure = nc.variables['air_pressure'][:] / 100  # Pa to hPa
            temperature = nc.variables['air_temperature'][:] - 273.15  # K to C
            humidity = nc.variables['relative_humidity'][:]
            height = nc.variables['altitude'][:]
            wind_speed = nc.variables['wind_speed'][:]
            wind_direction = nc.variables['wind_from_direction'][:]
            
            # Get launch location (first valid lat/lon or from attributes)
            if 'latitude' in nc.variables:
                lats = nc.variables['latitude'][:]
                lons = nc.variables['longitude'][:]
                # Get first valid position
                valid_mask = ~(np.ma.getmaskarray(lats) | np.ma.getmaskarray(lons))
                if np.any(valid_mask):
                    idx = np.where(valid_mask)[0][0]
                    latitude = float(lats[idx])
                    longitude = float(lons[idx])
                else:
                    latitude = 0.0
                    longitude = 0.0
            else:
                latitude = 0.0
                longitude = 0.0
            
            # Get timestamp
            timestamp = nc.getncattr('time_coverage_start') if 'time_coverage_start' in nc.ncattrs() else 'Unknown'
            
            # Convert masked arrays to regular arrays with NaN
            profile = {
                'pressure': np.ma.filled(pressure, np.nan),
                'temperature': np.ma.filled(temperature, np.nan),
                'humidity': np.ma.filled(humidity, np.nan),
                'height': np.ma.filled(height, np.nan),
                'wind_speed': np.ma.filled(wind_speed, np.nan),
                'wind_direction': np.ma.filled(wind_direction, np.nan),
                'latitude': latitude,
                'longitude': longitude,
                'timestamp': timestamp
            }
            
            return profile
            
    except Exception as e:
        print(f"Error reading NetCDF file {nc_file}: {str(e)}")
        return None


def plot_combined_analysis(profile, filename=None, map_extent=None, max_altitude=None, show_stability=False):
    """
    Combined plot: Skew-T on left, hodograph and map stacked on right
    
    Parameters:
    -----------
    profile : dict
        Radiosonde profile dictionary
    filename : str, optional
        Name of the source file to display in plot
    map_extent : list, optional
        Map extent [lon_min, lon_max, lat_min, lat_max] for consistent bounds
    max_altitude : float, optional
        Maximum altitude across all flights for consistent colormap scaling
    show_stability : bool, optional
        If True, display CAPE/CIN shading and LCL/LFC markers (default: False)
    """
    if not profile:
        print("No profile to plot")
        return None
    
    lat = profile['latitude']
    lon = profile['longitude']
    
    # Create figure
    fig = plt.figure(figsize=(20, 12))
    
    # ============ SKEW-T (Left side) ============
    # SkewT expects fig and rect parameters, not an axes object
    skew = SkewT(fig, rotation=45, rect=(0.05, 0.06, 0.42, 0.86))
    
    # Filter valid data for plotting
    valid_mask = ~(np.isnan(profile['pressure']) | 
                   np.isnan(profile['temperature']))
    
    if np.sum(valid_mask) > 0:
        p = profile['pressure'][valid_mask] * units.hPa
        T = profile['temperature'][valid_mask] * units.degC
        
        # Dewpoint
        rh = profile['humidity'][valid_mask]
        rh_valid = ~np.isnan(rh)
        Td = None
        if np.sum(rh_valid) > 0:
            from metpy.calc import dewpoint_from_relative_humidity
            Td = dewpoint_from_relative_humidity(T[rh_valid], rh[rh_valid] * units.percent)
            skew.plot(p[rh_valid], Td, color='#0173B2', linewidth=2.5, label='Dewpoint', alpha=0.8)
        
        skew.plot(p, T, color='#DE8F05', linewidth=2.5, label='Temperature', alpha=0.8)
        
        # Calculate and plot CAPE/CIN (if enabled)
        if show_stability and Td is not None and len(p) > 10:
            try:
                from metpy.calc import parcel_profile, surface_based_cape_cin, lcl, lfc
                
                # Remove duplicate pressure levels to avoid warnings
                # Need to work with magnitudes for np.unique, then reattach units
                p_magnitude = p.magnitude
                _, unique_indices = np.unique(p_magnitude, return_index=True)
                unique_indices = np.sort(unique_indices)
                p_unique = p[unique_indices]
                T_unique = T[unique_indices]
                Td_unique = Td[unique_indices]
                
                # Calculate LCL and LFC first
                lcl_pressure, lcl_temperature = lcl(p_unique[0], T_unique[0], Td_unique[0])
                lfc_pressure, lfc_temperature = lfc(p_unique, T_unique, Td_unique)
                
                # Create parcel profile using LCL-aware method
                # The parcel follows dry adiabat to LCL, then moist adiabat above
                prof = parcel_profile(p_unique, T_unique[0], Td_unique[0])
                
                # Plot the parcel temperature trace to visualize the parcel path
                skew.plot(p_unique, prof, 'k--', linewidth=1.5, label='Parcel Path', alpha=0.6)
                
                # Create pressure masks for proper CAPE/CIN shading
                # CIN: between LCL and LFC (where saturated parcel is colder than environment)
                # CAPE: above LFC (where parcel is warmer than environment)
                
                # Shade CIN (only between LCL and LFC)
                if not np.isnan(lfc_pressure.magnitude):
                    # CIN exists between LCL and LFC
                    p_between = (p_unique <= lcl_pressure) & (p_unique >= lfc_pressure)
                    if np.any(p_between):
                        p_cin = p_unique[p_between]
                        T_cin = T_unique[p_between]
                        prof_cin = prof[p_between]
                        Td_cin = Td_unique[p_between]
                        skew.shade_cin(p_cin, T_cin, prof_cin, Td_cin, alpha=0.3, label='CIN')
                    
                    # Shade CAPE (only above LFC)
                    p_above_lfc = p_unique <= lfc_pressure
                    if np.any(p_above_lfc):
                        p_cape = p_unique[p_above_lfc]
                        T_cape = T_unique[p_above_lfc]
                        prof_cape = prof[p_above_lfc]
                        skew.shade_cape(p_cape, T_cape, prof_cape, alpha=0.3, label='CAPE')
                
                # Calculate CAPE/CIN values
                cape, cin = surface_based_cape_cin(p_unique, T_unique, Td_unique)
                
                # Plot LCL (Lifting Condensation Level) - only if different from LFC
                lcl_label = f'LCL ({lcl_pressure.magnitude:.0f} hPa)'
                lfc_label = f'LFC ({lfc_pressure.magnitude:.0f} hPa)' if not np.isnan(lfc_pressure.magnitude) else None
                
                # Check if LCL and LFC are at same level (within 1 hPa)
                if lfc_label and abs(lcl_pressure.magnitude - lfc_pressure.magnitude) < 1.0:
                    # Same level - plot once with combined label
                    skew.plot(lcl_pressure, lcl_temperature, 'ko', markersize=8, markerfacecolor='purple',
                             label=f'LCL/LFC ({lcl_pressure.magnitude:.0f} hPa)')
                else:
                    # Different levels - plot separately
                    skew.plot(lcl_pressure, lcl_temperature, 'ko', markersize=8, markerfacecolor='black', 
                             label=lcl_label)
                    if not np.isnan(lfc_pressure.magnitude):
                        skew.plot(lfc_pressure, lfc_temperature, 'ko', markersize=8, markerfacecolor='red',
                                 label=lfc_label)
                
                # Add text annotation with values (lower left to avoid legend)
                skew.ax.text(0.02, 0.25, f'CAPE: {cape.magnitude:.0f} J/kg\nCIN: {cin.magnitude:.0f} J/kg',
                           transform=skew.ax.transAxes, fontsize=10, verticalalignment='top',
                           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            except Exception as e:
                print(f"  Warning: Could not calculate CAPE/CIN: {e}")
        
        # Wind barbs (subsample)
        wind_mask = ~(np.isnan(profile['wind_speed']) | np.isnan(profile['wind_direction']))
        combined_mask = valid_mask & wind_mask
        
        if np.sum(combined_mask) > 0:
            indices = np.where(combined_mask)[0]
            if len(indices) > 40:
                keep_indices = indices[::len(indices)//40]
            else:
                keep_indices = indices
            
            p_wind = profile['pressure'][keep_indices] * units.hPa
            u, v = wind_components(profile['wind_speed'][keep_indices], 
                                  profile['wind_direction'][keep_indices])
            skew.plot_barbs(p_wind, u * units('m/s'), v * units('m/s'))
    
    skew.plot_dry_adiabats(alpha=0.25, color='#CC78BC')
    skew.plot_moist_adiabats(alpha=0.25, color='#029E73')
    skew.plot_mixing_lines(alpha=0.25, color='#949494')
    
    skew.ax.set_ylim(1000, 100)
    skew.ax.set_xlim(-40, 40)
    skew.ax.set_xlabel('Temperature (°C)', fontsize=12, fontweight='bold')
    skew.ax.set_ylabel('Pressure (hPa)', fontsize=12, fontweight='bold')
    
    title = 'Skew-T Log-P Diagram'
    if filename:
        title += f'\n{filename}'
    if profile.get('timestamp'):
        title += f'\n{profile["timestamp"]}'
    skew.ax.set_title(title, fontsize=13, fontweight='bold', pad=10)
    
    # Legend in upper left
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
        color_max_height = max_altitude if max_altitude else max(heights)
        colors = plt.cm.plasma(heights / color_max_height)
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
        color_max_height = max_altitude if max_altitude else max(heights)
        sm = plt.cm.ScalarMappable(cmap=plt.cm.plasma, 
                                   norm=plt.Normalize(vmin=0, vmax=color_max_height/1000))
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
    
    # Map extent - use provided extent or calculate from this flight
    if map_extent is None:
        lon_padding = max(0.3, (drift_lons.max() - drift_lons.min()) * 0.2)
        lat_padding = max(0.3, (drift_lats.max() - drift_lats.min()) * 0.2)
        lon_min = min(drift_lons.min(), lon) - lon_padding
        lon_max = max(drift_lons.max(), lon) + lon_padding
        lat_min = min(drift_lats.min(), lat) - lat_padding
        lat_max = max(drift_lats.max(), lat) + lat_padding
        map_extent = [lon_min, lon_max, lat_min, lat_max]
    
    map_ax = fig.add_axes([0.55, 0.06, 0.40, 0.40], projection=ccrs.PlateCarree())
    map_ax.set_extent(map_extent, crs=ccrs.PlateCarree())
    
    # Add OpenTopoMap tiles with faded appearance and caching
    tiles_loaded = False
    try:
        # Create custom OpenTopoMap tile source with cache
        class OpenTopoMap(cimgt.GoogleTiles):
            def _image_url(self, tile):
                x, y, z = tile
                url = f'https://tile.opentopomap.org/{z}/{x}/{y}.png'
                return url
        
        # Enable caching in user's home directory
        cache_dir = os.path.expanduser('~/.cache/cartopy_tiles')
        os.makedirs(cache_dir, exist_ok=True)
        opentopo = OpenTopoMap(cache=cache_dir)
        
        # Set timeout and cache settings to handle server issues
        import warnings
        import socket
        # Set a socket timeout for tile downloads
        socket.setdefaulttimeout(10)  # 10 second timeout per tile
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # Use zoom level 9 (fewer tiles = faster, more reliable)
            map_ax.add_image(opentopo, 9, alpha=0.4, interpolation='bilinear')
        tiles_loaded = True
    except Exception:
        # Silently fall back to basic features if tiles fail (server down, no internet, etc.)
        pass
    finally:
        # Reset socket timeout
        import socket
        socket.setdefaulttimeout(None)
    
    if not tiles_loaded:
        # Use basic cartographic features as fallback
        map_ax.set_facecolor('#f0f0f0')
        map_ax.add_feature(cfeature.LAND, facecolor='#e8e8e8', edgecolor='none', zorder=1)
        map_ax.add_feature(cfeature.OCEAN, facecolor='#d4e7f5', edgecolor='none', zorder=1)
        map_ax.add_feature(cfeature.COASTLINE, linewidth=1.2, edgecolor='#666666', zorder=2)
        map_ax.add_feature(cfeature.BORDERS, linewidth=1.0, edgecolor='#888888', linestyle='--', zorder=2)
    
    gl = map_ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.6, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    # Plot trajectory with height-based colors
    color_max_height = max_altitude if max_altitude else max(heights)
    colors_traj = plt.cm.plasma(heights / color_max_height)
    for j in range(len(drift_lons) - 1):
        map_ax.plot([drift_lons[j], drift_lons[j+1]], 
                   [drift_lats[j], drift_lats[j+1]], 
                   color=colors_traj[j], linewidth=2.5, 
                   transform=ccrs.PlateCarree(), zorder=3, alpha=0.8)
    
    # Launch location
    map_ax.plot(lon, lat, marker='*', color='#ff3333', markersize=20, 
               transform=ccrs.PlateCarree(), label='Launch', zorder=4, 
               markeredgecolor='black', markeredgewidth=1.5)
    
    map_ax.set_title('Estimated Trajectory', fontsize=13, fontweight='bold', pad=10)
    map_ax.legend(loc='upper left', fontsize=11, framealpha=0.95)
    
    # Add OpenTopoMap attribution
    map_ax.text(0.99, 0.01, '© OpenTopoMap (CC-BY-SA)', 
               transform=map_ax.transAxes, fontsize=7, 
               verticalalignment='bottom', horizontalalignment='right',
               bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7, edgecolor='none'))
    
    return fig


def compute_overall_bounds(nc_files):
    """
    Compute overall geospatial bounds and maximum altitude from all NetCDF files
    
    Parameters:
    -----------
    nc_files : list
        List of NetCDF file paths
        
    Returns:
    --------
    tuple
        (map_extent, max_altitude) where map_extent is [lon_min, lon_max, lat_min, lat_max]
        or (None, None) if no valid data
    """
    all_lats = []
    all_lons = []
    all_altitudes = []
    
    print("Computing overall geospatial bounds and maximum altitude...")
    
    for nc_file in nc_files:
        try:
            with Dataset(nc_file, 'r') as nc:
                # Check if geospatial_bounds attribute starts with "CHANGE"
                skip_for_bounds = False
                if 'geospatial_bounds' in nc.ncattrs():
                    geospatial_bounds = nc.getncattr('geospatial_bounds')
                    if isinstance(geospatial_bounds, str) and geospatial_bounds.startswith('CHANGE'):
                        skip_for_bounds = True
                        print(f"  Skipping {os.path.basename(nc_file)} for bounds calculation (geospatial_bounds starts with 'CHANGE')")
                
                if not skip_for_bounds and 'latitude' in nc.variables and 'longitude' in nc.variables:
                    lats = nc.variables['latitude'][:]
                    lons = nc.variables['longitude'][:]
                    
                    # Extract valid (non-masked) values
                    valid_lats = lats[~np.ma.getmaskarray(lats)]
                    valid_lons = lons[~np.ma.getmaskarray(lons)]
                    
                    # Validate that lat/lon values are within valid ranges
                    if len(valid_lats) > 0 and len(valid_lons) > 0:
                        # Check for reasonable latitude values (-90 to 90)
                        lat_valid = np.all((valid_lats >= -90) & (valid_lats <= 90))
                        # Check for reasonable longitude values (-180 to 180)
                        lon_valid = np.all((valid_lons >= -180) & (valid_lons <= 180))
                        
                        if lat_valid and lon_valid:
                            all_lats.extend(valid_lats)
                            all_lons.extend(valid_lons)
                        else:
                            print(f"  Skipping {os.path.basename(nc_file)}: invalid lat/lon values")
                            continue
                
                # Get altitude/height data (always collect, even if skipped for bounds)
                if 'altitude' in nc.variables:
                    alts = nc.variables['altitude'][:]
                    valid_alts = alts[~np.ma.getmaskarray(alts)]
                    if len(valid_alts) > 0:
                        all_altitudes.extend(valid_alts)
        except Exception as e:
            print(f"  Warning: Could not read data from {os.path.basename(nc_file)}: {e}")
            continue
    
    if not all_lats or not all_lons:
        print("  No valid geospatial data found")
        return None, None
    
    # Calculate overall bounds with padding
    lat_min, lat_max = min(all_lats), max(all_lats)
    lon_min, lon_max = min(all_lons), max(all_lons)
    
    lon_padding = max(0.3, (lon_max - lon_min) * 0.15)
    lat_padding = max(0.3, (lat_max - lat_min) * 0.15)
    
    map_extent = [
        lon_min - lon_padding,
        lon_max + lon_padding,
        lat_min - lat_padding,
        lat_max + lat_padding
    ]
    
    # Get maximum altitude
    max_altitude = max(all_altitudes) if all_altitudes else None
    
    print(f"  Overall bounds: Lon [{map_extent[0]:.3f}, {map_extent[1]:.3f}], Lat [{map_extent[2]:.3f}, {map_extent[3]:.3f}]")
    if max_altitude:
        print(f"  Maximum altitude: {max_altitude:.1f} m ({max_altitude/1000:.1f} km)")
    
    return map_extent, max_altitude


def process_netcdf_files(input_dir, output_dir=None, show_stability=False):
    """
    Process all NetCDF files in input directory tree and generate quicklook plots
    
    Searches recursively for NetCDF files in subdirectories (e.g., instrument/project/YYYY/MM/)
    and saves quicklooks to quicklooks directories at the instrument/project level.
    
    Parameters:
    -----------
    input_dir : str
        Root directory containing NetCDF files (searches recursively)
    output_dir : str, optional
        Not used - kept for backwards compatibility
    show_stability : bool, optional
        If True, display CAPE/CIN shading and LCL/LFC markers (default: False)
    """
    # Find all NetCDF files recursively
    nc_pattern = os.path.join(input_dir, "**", "*.nc")
    nc_files = glob.glob(nc_pattern, recursive=True)
    
    if not nc_files:
        print(f"No NetCDF files found in {input_dir}")
        return
    
    print(f"Found {len(nc_files)} NetCDF file(s)")
    
    # Compute overall geospatial bounds and max altitude from all files
    map_extent, max_altitude = compute_overall_bounds(nc_files)
    
    for nc_file in sorted(nc_files):
        print(f"\nProcessing: {os.path.basename(nc_file)}")
        
        # Read profile from NetCDF
        profile = read_netcdf_profile(nc_file)
        if profile is None:
            print(f"  Skipping {os.path.basename(nc_file)} - could not read profile")
            continue
        
        # Generate plot with consistent map bounds and colormap scaling
        try:
            fig = plot_combined_analysis(profile, filename=os.path.basename(nc_file), 
                                        map_extent=map_extent, max_altitude=max_altitude,
                                        show_stability=show_stability)
            
            if fig is not None:
                # Place quicklooks directory directly under input_dir
                # All plots go in a flat structure: input_dir/quicklooks/file.png
                nc_path = Path(nc_file)
                quicklooks_dir = Path(input_dir) / 'quicklooks'
                
                os.makedirs(quicklooks_dir, exist_ok=True)
                plot_file = quicklooks_dir / (nc_path.stem + '.png')
                
                fig.savefig(str(plot_file), dpi=150, bbox_inches='tight')
                plt.close(fig)
                print(f"  Saved plot: {plot_file}")
            else:
                print(f"  Could not generate plot for {os.path.basename(nc_file)}")
        except Exception as e:
            print(f"  Error generating plot: {str(e)}")
            import traceback
            traceback.print_exc()
    
    print(f"\nProcessing complete.")


def main():
    """Main entry point for script"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Generate quicklook plots from radiosonde NetCDF files')
    parser.add_argument('input_dir', help='Directory containing NetCDF files')
    parser.add_argument('output_dir', nargs='?', default=None, help='(Not used - quicklooks are created alongside NetCDF files)')
    parser.add_argument('--stability', action='store_true', 
                       help='Show CAPE/CIN shading and LCL/LFC markers')
    
    args = parser.parse_args()
    
    if not os.path.isdir(args.input_dir):
        print(f"Error: Input directory '{args.input_dir}' does not exist")
        sys.exit(1)
    
    process_netcdf_files(args.input_dir, args.output_dir, show_stability=args.stability)


if __name__ == "__main__":
    main()
