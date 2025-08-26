#!/usr/bin/env python3
"""
Traverse Sentinel-2 data structure and extract band information for 20m and 60m resolutions
"""

import os
import sys
import re
from pathlib import Path


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
    
    # Dictionaries to store band:filename mappings
    bands_20m = {}  # band_name: filename
    bands_60m = {}  # band_name: filename
    
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
                                bands_20m[band_name] = jp2_file.name
                                print(f"        -> Added to 20m dict: {band_name} -> {jp2_file.name}")
                            elif resolution == '60m':
                                bands_60m[band_name] = jp2_file.name
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
                            bands_20m[band_name] = jp2_file.name
                            print(f"        -> Added to 20m dict: {band_name} -> {jp2_file.name}")
                        elif resolution == '60m':
                            bands_60m[band_name] = jp2_file.name
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
        for band, filename in sorted(bands_20m.items()):
            print(f"  {band}: {filename}")
    
    if bands_60m:
        print("\n60m Resolution Band Dictionary:")
        print("-" * 40)
        for band, filename in sorted(bands_60m.items()):
            print(f"  {band}: {filename}")
    
    return bands_20m, bands_60m


def main():
    """Main function to run the traversal"""
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
    traverse_sentinel2_data(sentinel_path)


if __name__ == "__main__":
    main()