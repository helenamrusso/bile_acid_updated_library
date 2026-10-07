# Root of the local project folder; edit to run elsewhere.
PROJECT_DIR = '/media/dorrestein/Helena2/Bile_acids_update'

#!/usr/bin/env python3
"""
Script to process MGF files:
- Add SCANS= line after CLUSTER
- Increment scan numbers for each spectrum
- Save with _scan suffix
"""

import os
from pathlib import Path
import re

def process_mgf_file(input_path, output_path):
    """
    Process a single MGF file to add SCANS= line after CLUSTER.
    
    Args:
        input_path: Path to input MGF file
        output_path: Path to output MGF file
    
    Returns:
        Number of spectra processed
    """
    scan_number = 1
    spectra_count = 0
    
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        for line in infile:
            # Write the current line
            outfile.write(line)
            
            # If we find a CLUSTER line, add SCANS line after it
            if line.startswith('CLUSTER='):
                outfile.write(f'SCANS={scan_number}\n')
                scan_number += 1
            
            # Count spectra
            if line.startswith('BEGIN IONS'):
                spectra_count += 1
    
    return spectra_count

def process_all_mgf_files(folder_path):
    """
    Process all MGF files in a folder.
    
    Args:
        folder_path: Path to folder containing MGF files
    """
    folder = Path(folder_path)
    
    if not folder.exists():
        print(f"❌ Error: Folder does not exist: {folder_path}")
        return
    
    # Find all MGF files
    mgf_files = list(folder.glob('*.mgf'))
    
    if not mgf_files:
        print(f"⚠️  No MGF files found in: {folder_path}")
        return
    
    print(f"Found {len(mgf_files)} MGF file(s) to process\n")
    
    total_spectra = 0
    processed_files = 0
    errors = []
    
    for mgf_file in mgf_files:
        # Create output filename with _scan suffix
        output_name = mgf_file.stem + '_scan.mgf'
        output_path = mgf_file.parent / output_name
        
        try:
            print(f"Processing: {mgf_file.name}")
            spectra_count = process_mgf_file(mgf_file, output_path)
            total_spectra += spectra_count
            processed_files += 1
            print(f"  ✓ Processed {spectra_count} spectra")
            print(f"  ✓ Saved to: {output_name}\n")
            
        except Exception as e:
            error_msg = f"Error processing {mgf_file.name}: {str(e)}"
            print(f"  ❌ {error_msg}\n")
            errors.append(error_msg)
    
    # Summary
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Total MGF files processed: {processed_files}/{len(mgf_files)}")
    print(f"Total spectra processed: {total_spectra}")
    print(f"Errors: {len(errors)}")
    
    if errors:
        print("\nErrors encountered:")
        for error in errors:
            print(f"  - {error}")
    
    print(f"\nProcessed files saved in: {folder_path}")

if __name__ == "__main__":
    # Folder containing MGF files
    mgf_folder = f"{PROJECT_DIR}/Falcon/all_falcon_outputs"
    
    print("Starting MGF processing...\n")
    print(f"Input folder: {mgf_folder}\n")
    
    process_all_mgf_files(mgf_folder)
    
    print("\nOperation complete!")