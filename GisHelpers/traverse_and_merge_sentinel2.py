#!/usr/bin/env python3
"""
Traverse Sentinel-2 data structure and extract band information for 20m and 60m resolutions
"""

import os
import sys
import re
from pathlib import Path
import rasterio
from rasterio.merge import merge


def extract_band_from_filename(filename):
    """
    Extract band name from Sentinel-2 .jp2 filename.
    
    Sentinel-2 filenames typically follow the pattern:
    T[UTM_ZONE][LATITUDE_BAND][GRID_SQUARE]_[ACQUISITION_DATE]_[BAND]_[RESOLUTION].jp2
    Example: T32TMS_20220615T104031_B02_20m.jp2 -> Band B02
    Example: T29SQB_20230829T110621_AOT_20m.jp2 -> AOT
    Example: T29SQB_20230829T110621_B8A_20m.jp2 -> B8A
    
    Args:
        filename (str): The .jp2 filename
        
    Returns:
        str: Band name (e.g., 'B02', 'B03', 'AOT', 'TCI', 'SCL', 'WVP') or None if not found
    """
    # Check for auxiliary data layers first (AOT, TCI, SCL, WVP)
    aux_match = re.search(r'_(AOT|TCI|SCL|WVP)_', filename)
    if aux_match:
        return aux_match.group(1)
    
    # Check for B8A band (special case)
    if '_B8A_' in filename:
        return 'B8A'
    
    # Try to find standard band pattern like B01, B02, B03, etc.
    band_match = re.search(r'_B(\d{2})_', filename)
    if band_match:
        return f"B{band_match.group(1)}"
    
    # Alternative pattern for some Sentinel-2 files
    band_match = re.search(r'_B(\d{2})\d*_', filename)
    if band_match:
        return f"B{band_match.group(1)}"
    
    # Another common pattern
    band_match = re.search(r'B(\d{2})', filename)
    if band_match:
        return f"B{band_match.group(1)}"
    
    return None


def extract_resolution_from_path_or_filename(file_path, filename):
    """
    Extract resolution from directory name (R20m, R60m) or filename.
    
    Args:
        file_path (Path): Full path to the file
        filename (str): The .jp2 filename
        
    Returns:
        str: Resolution ('20m', '60m') or None if not found or not relevant
    """
    # First try to get resolution from parent directory name
    parent_dir = file_path.parent.name
    if parent_dir.startswith('R'):
        resolution = parent_dir[1:]  # Remove 'R' prefix
        if resolution in ['20m', '60m']:
            return resolution
    
    # If not in directory name, try to extract from filename
    res_match = re.search(r'_(\d{2})m\.jp2$', filename)
    if res_match:
        resolution = f"{res_match.group(1)}m"
        if resolution in ['20m', '60m']:
            return resolution
    
    return None


def traverse_sentinel2_data(root_path):
    """
    Traverse the Sentinel-2 data directory structure and organize files by resolution and band.
    
    Args:
        root_path (str): Path to the root directory containing .SAFE folders
        
    Returns:
        tuple: (bands_20m_dict, bands_60m_dict) where keys are band names and values are filenames
    """
    root_path = Path(root_path)
    
    if not root_path.exists():
        print(f"Error: Directory {root_path} does not exist")
        return {}, {}
    
    print(f"Traversing Sentinel-2 data in: {root_path}")
    print("=" * 80)
    
    # Find all .SAFE directories
    safe_dirs = list(root_path.glob("*.SAFE"))
    
    if not safe_dirs:
        print("No .SAFE directories found")
        return {}, {}
    
    print(f"Found {len(safe_dirs)} .SAFE directories:")
    print("-" * 40)
    
    # Dictionaries to store band:list_of_file_paths mappings
    bands_20m = {}  # band_name: list of full file paths
    bands_60m = {}  # band_name: list of full file paths
    
    for safe_dir in sorted(safe_dirs):
        print(f"\nProcessing: {safe_dir.name}")
        
        # Look for GRANULE directories
        granule_dir = safe_dir / "GRANULE"
        
        if not granule_dir.exists():
            print(f"  Warning: No GRANULE directory found in {safe_dir.name}")
            continue
        
        # Find all subdirectories in GRANULE
        granule_subdirs = [d for d in granule_dir.iterdir() if d.is_dir()]
        
        for granule_subdir in granule_subdirs:
            print(f"  GRANULE: {granule_subdir.name}")
            
            # Look for IMG_DATA directory
            img_data_dir = granule_subdir / "IMG_DATA"
            
            if not img_data_dir.exists():
                print(f"    Warning: No IMG_DATA directory found in {granule_subdir.name}")
                continue
            
            # Check for resolution directories (R10m, R20m, R60m)
            resolution_dirs = [d for d in img_data_dir.iterdir() if d.is_dir() and d.name.startswith('R')]
            
            if resolution_dirs:
                # Process resolution directories
                for res_dir in sorted(resolution_dirs):
                    # Only process 20m and 60m resolutions
                    if res_dir.name not in ['R20m', 'R60m']:
                        print(f"    Skipping resolution: {res_dir.name}")
                        continue
                    
                    print(f"    Processing resolution: {res_dir.name}")
                    
                    # Find all .jp2 files
                    jp2_files_in_dir = list(res_dir.glob("*.jp2"))
                    
                    for jp2_file in sorted(jp2_files_in_dir):
                        print(f"      {jp2_file.name}")
                        
                        # Extract band and resolution information
                        band_name = extract_band_from_filename(jp2_file.name)
                        resolution = extract_resolution_from_path_or_filename(jp2_file, jp2_file.name)
                        
                        if band_name and resolution:
                            if resolution == '20m':
                                if band_name not in bands_20m:
                                    bands_20m[band_name] = []
                                bands_20m[band_name].append(str(jp2_file))
                                print(f"        -> Added to 20m dict: {band_name} -> {jp2_file.name}")
                            elif resolution == '60m':
                                if band_name not in bands_60m:
                                    bands_60m[band_name] = []
                                bands_60m[band_name].append(str(jp2_file))
                                print(f"        -> Added to 60m dict: {band_name} -> {jp2_file.name}")
                        else:
                            print(f"        -> Could not extract band/resolution info from {jp2_file.name}")
            else:
                # Sometimes .jp2 files are directly in IMG_DATA
                print(f"    Checking IMG_DATA directory directly")
                jp2_files_in_dir = list(img_data_dir.glob("*.jp2"))
                
                for jp2_file in sorted(jp2_files_in_dir):
                    print(f"      {jp2_file.name}")
                    
                    # Extract band and resolution information
                    band_name = extract_band_from_filename(jp2_file.name)
                    resolution = extract_resolution_from_path_or_filename(jp2_file, jp2_file.name)
                    
                    if band_name and resolution:
                        if resolution == '20m':
                            if band_name not in bands_20m:
                                bands_20m[band_name] = []
                            bands_20m[band_name].append(str(jp2_file))
                            print(f"        -> Added to 20m dict: {band_name} -> {jp2_file.name}")
                        elif resolution == '60m':
                            if band_name not in bands_60m:
                                bands_60m[band_name] = []
                            bands_60m[band_name].append(str(jp2_file))
                            print(f"        -> Added to 60m dict: {band_name} -> {jp2_file.name}")
                    else:
                        print(f"        -> Could not extract band/resolution info from {jp2_file.name}")
    
    print("\n" + "=" * 80)
    print("RESULTS:")
    print(f"20m resolution bands found: {len(bands_20m)}")
    print(f"60m resolution bands found: {len(bands_60m)}")
    
    if bands_20m:
        print("\n20m Resolution Band Dictionary:")
        print("-" * 40)
        for band, file_list in sorted(bands_20m.items()):
            print(f"  {band}: {len(file_list)} files")
            for file_path in file_list:
                filename = Path(file_path).name
                print(f"    - {filename}")
    
    if bands_60m:
        print("\n60m Resolution Band Dictionary:")
        print("-" * 40)
        for band, file_list in sorted(bands_60m.items()):
            print(f"  {band}: {len(file_list)} files")
            for file_path in file_list:
                filename = Path(file_path).name
                print(f"    - {filename}")
    
    return bands_20m, bands_60m


def merge_sentinel_bands(root_path, bands_20m, bands_60m, output_dir="combined"):
    """
    Simple merging of Sentinel-2 bands by resolution, saving as JP2 files.
    
    Args:
        root_path (Path): Path to the root directory containing .SAFE folders
        bands_20m (dict): Dictionary of 20m resolution bands
        bands_60m (dict): Dictionary of 60m resolution bands
        output_dir (str): Output directory name
    """
    root_path = Path(root_path)
    output_path = root_path.parent / output_dir
    output_path.mkdir(exist_ok=True)
    
    print(f"\n{'='*80}")
    print("MERGING BANDS")
    print(f"Output directory: {output_path}")
    print(f"{'='*80}")
    
    # Process each resolution
    for resolution, bands_dict in [("20m", bands_20m), ("60m", bands_60m)]:
        if not bands_dict:
            print(f"No {resolution} bands to process")
            continue
            
        print(f"\nProcessing {resolution} resolution bands...")
        res_output_dir = output_path / resolution
        res_output_dir.mkdir(exist_ok=True)
        
        # Process each band
        for band_name, file_paths_list in bands_dict.items():
            print(f"  Processing band {band_name}...")
            
            # Use the file paths directly from our collected data
            band_files = [Path(file_path) for file_path in file_paths_list]
            
            if not band_files:
                print(f"    Warning: No files found for band {band_name}")
                continue
                
            print(f"    Found {len(band_files)} files for band {band_name}")
            
            # Simple merging - if multiple files, merge them; if single file, copy
            try:
                if len(band_files) == 1:
                    # Single file - read and write as JP2
                    with rasterio.open(band_files[0]) as src:
                        data = src.read()
                        profile = src.profile
                        file_crs = src.crs
                    
                    # Report CRS information for single file
                    print(f"    CRS for merged {band_name}: {file_crs}")
                    if file_crs is None:
                        print(f"    Warning: No CRS information found in {band_files[0].name}")
                        
                    # Try different JP2 drivers in order of preference
                    jp2_drivers = ['j2k', 'JP2ECW', 'JP2KAK']
                    output_file = res_output_dir / f"merged_{band_name}_{resolution}.jp2"
                    
                    saved = False
                    for driver in jp2_drivers:
                        try:
                            profile_copy = profile.copy()
                            profile_copy.update(driver=driver)
                            with rasterio.open(output_file, 'w', **profile_copy) as dst:
                                dst.write(data)
                            saved = True
                            print(f"    Saved using {driver} driver: {output_file.name}")
                            break
                        except Exception as driver_error:
                            if "No such driver registered" not in str(driver_error):
                                print(f"    Failed with {driver}: {driver_error}")
                            continue
                    
                    # Fallback: save as GeoTIFF with JP2 extension for compatibility
                    if not saved:
                        try:
                            profile_copy = profile.copy()
                            profile_copy.update(driver='GTiff')
                            with rasterio.open(output_file, 'w', **profile_copy) as dst:
                                dst.write(data)
                            print(f"    Saved as GeoTIFF with .jp2 extension (JP2 drivers unavailable): {output_file.name}")
                            saved = True
                        except Exception as e:
                            print(f"    Error: Could not save {band_name}: {e}")
                    
                    if not saved:
                        print(f"    Warning: Failed to save {band_name}")
                    
                else:
                    # Multiple files - merge them
                    datasets = []
                    try:
                        for file_path in band_files:
                            ds = rasterio.open(file_path)
                            datasets.append(ds)
                        
                        # Group datasets by CRS to avoid merge errors
                        crs_groups = {}
                        for ds in datasets:
                            crs_key = str(ds.crs) if ds.crs else 'None'
                            if crs_key not in crs_groups:
                                crs_groups[crs_key] = []
                            crs_groups[crs_key].append(ds)
                        
                        # Report CRS information
                        print(f"    CRS analysis for {len(datasets)} files:")
                        for crs_key, group in crs_groups.items():
                            if crs_key == 'None':
                                print(f"      No CRS: {len(group)} files")
                            else:
                                print(f"      {crs_key}: {len(group)} files")
                        
                        # Process the largest CRS group
                        largest_group = max(crs_groups.values(), key=len)
                        main_crs = str(largest_group[0].crs) if largest_group[0].crs else 'None'
                        print(f"    Using CRS '{main_crs}' for merged {band_name} ({len(largest_group)} files)")
                        
                        # Report CRS mismatches with specific file names
                        if len(largest_group) < len(datasets):
                            print(f"    Warning: Using {len(largest_group)} of {len(datasets)} files due to CRS differences")
                            print(f"    CRS mismatches detected:")
                            for crs_key, group in crs_groups.items():
                                if group != largest_group:
                                    for ds in group:
                                        file_name = ds.name.split('/')[-1] if '/' in ds.name else ds.name
                                        if crs_key == 'None':
                                            print(f"      {file_name}: No CRS")
                                        else:
                                            print(f"      {file_name}: {crs_key}")
                        
                        # Merge the datasets from the largest CRS group
                        merged_data, merged_transform = merge(largest_group)
                        
                        # Get profile from first dataset and update for JP2
                        profile = largest_group[0].profile
                        profile.update(
                            height=merged_data.shape[1],
                            width=merged_data.shape[2],
                            transform=merged_transform
                        )
                        
                        output_file = res_output_dir / f"merged_{band_name}_{resolution}.jp2"
                        
                        # Try different JP2 drivers
                        jp2_drivers = ['j2k', 'JP2ECW', 'JP2KAK']
                        saved = False
                        for driver in jp2_drivers:
                            try:
                                profile_copy = profile.copy()
                                profile_copy.update(driver=driver)
                                with rasterio.open(output_file, 'w', **profile_copy) as dst:
                                    dst.write(merged_data)
                                saved = True
                                print(f"    Merged {len(largest_group)} files using {driver} -> {output_file.name}")
                                break
                            except Exception as driver_error:
                                if "No such driver registered" not in str(driver_error):
                                    print(f"    Failed with {driver}: {driver_error}")
                                continue
                        
                        # Fallback: save as GeoTIFF with JP2 extension for compatibility
                        if not saved:
                            try:
                                profile_copy = profile.copy()
                                profile_copy.update(driver='GTiff')
                                with rasterio.open(output_file, 'w', **profile_copy) as dst:
                                    dst.write(merged_data)
                                print(f"    Merged {len(largest_group)} files as GeoTIFF with .jp2 extension (JP2 drivers unavailable) -> {output_file.name}")
                                saved = True
                            except Exception as e:
                                print(f"    Error: Could not save merged {band_name}: {e}")
                        
                        if not saved:
                            print(f"    Warning: Failed to save merged {band_name}")
                        
                    finally:
                        # Close all datasets
                        for ds in datasets:
                            try:
                                ds.close()
                            except:
                                pass
                                
            except Exception as e:
                print(f"    Error processing band {band_name}: {e}")
                continue
    
    print(f"\n{'='*80}")
    print(f"MERGING COMPLETED - Results saved to: {output_path}")
    print(f"{'='*80}")


def main():
    """Main function to run the traversal and merging"""
    # Use the path based on the project structure discovered
    # The path has a space in it: "Data/Tharsis_all_Sentinel2 /raw"
    base_path = Path(__file__).parent.parent  # Go up from GisHelpers to project root
    sentinel_path = base_path / "Data" / "Tharsis_all_Sentinel2 " / "raw"
    
    # Alternative path without space in case the space was a display artifact
    if not sentinel_path.exists():
        sentinel_path = base_path / "Data" / "Tharsis_all_Sentinel2" / "raw"
    
    # If still not found, try relative path as specified in issue description
    if not sentinel_path.exists():
        sentinel_path = Path("../Data/Tharsis_all_Sentinel2/raw")
    
    print(f"Attempting to traverse: {sentinel_path}")
    
    # Run traversal to get band dictionaries
    bands_20m, bands_60m = traverse_sentinel2_data(sentinel_path)
    
    # Run merging if any bands were found
    if bands_20m or bands_60m:
        merge_sentinel_bands(sentinel_path, bands_20m, bands_60m)
    else:
        print("No bands found to merge")


if __name__ == "__main__":
    main()