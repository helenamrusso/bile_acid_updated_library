#!/usr/bin/env python3
"""
Script to copy all files from multiple falcon_output folders to a single destination.
"""

import os
import shutil
from pathlib import Path

# Root of the local project folder; edit to run elsewhere.
PROJECT_DIR = '/media/dorrestein/Helena2/Bile_acids_update'

# Define source folders
source_folders = [
    f"{PROJECT_DIR}/Queries_final_mgf/C24_final/Nonhydroxy/falcon_output",
    f"{PROJECT_DIR}/Queries_final_mgf/C24_final/Monohydroxy/falcon_output",
    f"{PROJECT_DIR}/Queries_final_mgf/C24_final/Dihydroxy/falcon_output",
    f"{PROJECT_DIR}/Queries_final_mgf/C24_final/Trihydroxy/falcon_output",
    f"{PROJECT_DIR}/Queries_final_mgf/C24_final/Tetrahydroxy/falcon_output",
    f"{PROJECT_DIR}/Queries_final_mgf/C24_final/Pentahydroxy/falcon_output",
    f"{PROJECT_DIR}/Queries_final_mgf/C23/C23_monohydroxy/falcon_output",
    f"{PROJECT_DIR}/Queries_final_mgf/C23/C23_dihydroxy/falcon_output",
    f"{PROJECT_DIR}/Queries_final_mgf/C23/C23_trihydroxy/falcon_output",
    f"{PROJECT_DIR}/Queries_final_mgf/C23/C23_tetrahydroxy/falcon_output",
    f"{PROJECT_DIR}/Queries_final_mgf/C23/C23_pentahydroxy/falcon_output",
    f"{PROJECT_DIR}/Queries_final_mgf/C27/C27_monohydroxy/falcon_output",
    f"{PROJECT_DIR}/Queries_final_mgf/C27/C27_dihydroxy/falcon_output",
    f"{PROJECT_DIR}/Queries_final_mgf/C27/C27_trihydroxy/falcon_output",
    f"{PROJECT_DIR}/Queries_final_mgf/C27/C27_tetrahydroxy/falcon_output",
    f"{PROJECT_DIR}/Queries_final_mgf/C27/C27_pentahydroxy/falcon_output",
    f"{PROJECT_DIR}/Queries_final_mgf/Aceto/Dihydroxy/falcon_output",
    f"{PROJECT_DIR}/Queries_final_mgf/C24_alkylamine/C24_alkylamine_monohydroxy/falcon_output",
    f"{PROJECT_DIR}/Queries_final_mgf/C24_alkylamine/C24_alkylamine_dihydroxy/falcon_output",
    f"{PROJECT_DIR}/Queries_final_mgf/C24_alkylamine/C24_alkylamine_trihydroxy/falcon_output",
    f"{PROJECT_DIR}/Queries_final_mgf/C24_alkylamine/C24_alkylamine_tetrahydroxy/falcon_output",
    f"{PROJECT_DIR}/Queries_final_mgf/C24_alkylamine/C24_alkylamine_pentahydroxy/falcon_output",
    f"{PROJECT_DIR}/Queries_final_mgf/C28_alcohol/One_OH/falcon_output",
    f"{PROJECT_DIR}/Queries_final_mgf/C28_alcohol/Two_OH/falcon_output",
    f"{PROJECT_DIR}/Queries_final_mgf/C28_alcohol/Three_OH/falcon_output",
]

# Define destination folder
destination = f"{PROJECT_DIR}/Falcon/all_falcon_outputs"

def copy_files(source_folders, destination):
    """
    Copy all files from source folders to destination folder.
    
    Args:
        source_folders: List of source folder paths
        destination: Destination folder path
    """
    # Create destination folder if it doesn't exist
    dest_path = Path(destination)
    dest_path.mkdir(parents=True, exist_ok=True)
    
    total_files = 0
    total_copied = 0
    skipped_files = 0
    errors = []
    
    print(f"Destination folder: {destination}\n")
    
    for source_folder in source_folders:
        source_path = Path(source_folder)
        
        # Check if source folder exists
        if not source_path.exists():
            print(f"⚠️  Source folder does not exist: {source_folder}")
            errors.append(f"Folder not found: {source_folder}")
            continue
        
        if not source_path.is_dir():
            print(f"⚠️  Not a directory: {source_folder}")
            errors.append(f"Not a directory: {source_folder}")
            continue
        
        # Get all files in the source folder
        files = [f for f in source_path.iterdir() if f.is_file()]
        
        if not files:
            print(f"📁 Empty folder: {source_folder}")
            continue
        
        print(f"📁 Processing: {source_folder}")
        print(f"   Found {len(files)} file(s)")
        
        for file_path in files:
            total_files += 1
            dest_file = dest_path / file_path.name
            
            # Handle filename conflicts
            if dest_file.exists():
                # Check if files are identical
                if file_path.stat().st_size == dest_file.stat().st_size:
                    print(f"   ⏭️  Skipping (already exists): {file_path.name}")
                    skipped_files += 1
                    continue
                else:
                    # Add suffix to filename to avoid overwriting
                    counter = 1
                    stem = file_path.stem
                    suffix = file_path.suffix
                    while dest_file.exists():
                        new_name = f"{stem}_{counter}{suffix}"
                        dest_file = dest_path / new_name
                        counter += 1
                    print(f"   📝 Renamed to avoid conflict: {dest_file.name}")
            
            try:
                shutil.copy2(file_path, dest_file)
                total_copied += 1
                print(f"   ✓ Copied: {file_path.name}")
            except Exception as e:
                error_msg = f"Error copying {file_path}: {str(e)}"
                print(f"   ❌ {error_msg}")
                errors.append(error_msg)
        
        print()  # Empty line between folders
    
    # Summary
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Total files found: {total_files}")
    print(f"Files copied: {total_copied}")
    print(f"Files skipped (already exist): {skipped_files}")
    print(f"Errors: {len(errors)}")
    
    if errors:
        print("\nErrors encountered:")
        for error in errors:
            print(f"  - {error}")
    
    print(f"\nAll files copied to: {destination}")

if __name__ == "__main__":
    print("Starting file copy operation...\n")
    copy_files(source_folders, destination)
    print("\nOperation complete!")