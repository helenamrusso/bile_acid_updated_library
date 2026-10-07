#!/usr/bin/env python3
"""
Sort and split MGF files by PEPMASS.

Usage:
    python sort_and_split_mgf.py input.mgf --output-dir ./output --num-parts 3
    python sort_and_split_mgf.py input.mgf --output sorted.mgf --no-split
"""

import argparse
import os
from pathlib import Path
from tqdm import tqdm


class MGFBlock:
    """Represents a single MGF spectrum block."""
    
    def __init__(self, lines):
        self.lines = lines
        self._pepmass = None
        self._scan = None
    
    @property
    def pepmass(self):
        """Extract MASSQL_PEPMASS value from block."""
        if self._pepmass is None:
            for line in self.lines:
                if line.startswith("MASSQL_PEPMASS="):
                    self._pepmass = float(line.split("=", 1)[1].strip())
                    break
            if self._pepmass is None:
                raise ValueError("MASSQL_PEPMASS missing in block")
        return self._pepmass
    
    def to_string(self, scan_number=None):
        """Convert block to string, optionally renumbering scan."""
        if scan_number is None:
            return "".join(self.lines)
        
        result = []
        for line in self.lines:
            if line.startswith("SCANS="):
                result.append(f"SCANS={scan_number}\n")
            else:
                result.append(line)
        return "".join(result)


def read_mgf_blocks(filepath):
    """
    Read MGF file and return list of MGFBlock objects.
    
    Args:
        filepath: Path to input MGF file
        
    Returns:
        List of MGFBlock objects
    """
    blocks = []
    current_block = []
    
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("BEGIN IONS"):
                current_block = [line]
            elif line.startswith("END IONS"):
                current_block.append(line)
                blocks.append(MGFBlock(current_block))
                current_block = []
            elif current_block:  # Only append if we're inside a block
                current_block.append(line)
    
    return blocks


def sort_mgf_by_pepmass(input_file, output_file):
    """
    Sort MGF file by MASSQL_PEPMASS and renumber scans.
    
    Args:
        input_file: Path to input MGF file
        output_file: Path to output sorted MGF file
        
    Returns:
        Number of blocks processed
    """
    print(f"Reading and parsing {input_file}...")
    blocks = read_mgf_blocks(input_file)
    
    print(f"Found {len(blocks)} spectra")
    print("Sorting by PEPMASS...")
    
    # Sort by pepmass
    blocks.sort(key=lambda b: b.pepmass)
    
    print(f"Writing sorted MGF to {output_file}...")
    with open(output_file, "w") as f:
        for i, block in enumerate(tqdm(blocks, desc="Writing"), 1):
            f.write(block.to_string(scan_number=i))
    
    print(f"Successfully sorted {len(blocks)} spectra")
    return len(blocks)


def split_mgf(input_file, output_dir, num_parts, prefix="part"):
    """
    Split MGF file into equal parts with renumbered scans.
    
    Args:
        input_file: Path to input MGF file
        output_dir: Directory for output files
        num_parts: Number of parts to split into
        prefix: Prefix for output filenames
        
    Returns:
        List of output file paths
    """
    if num_parts < 1:
        raise ValueError("num_parts must be at least 1")
    
    print(f"Reading {input_file}...")
    blocks = read_mgf_blocks(input_file)
    total_scans = len(blocks)
    
    print(f"Total scans: {total_scans}")
    print(f"Splitting into {num_parts} part(s)...")
    
    # Calculate split sizes
    base_size = total_scans // num_parts
    remainder = total_scans % num_parts
    
    # Create output directory if needed
    os.makedirs(output_dir, exist_ok=True)
    
    # Calculate sizes for each part (distribute remainder across first parts)
    part_sizes = [base_size + (1 if i < remainder else 0) for i in range(num_parts)]
    
    # Split blocks
    output_files = []
    start_idx = 0
    
    for part_num, size in enumerate(part_sizes, 1):
        end_idx = start_idx + size
        part_blocks = blocks[start_idx:end_idx]
        
        # Generate output filename
        input_name = Path(input_file).stem
        output_file = os.path.join(output_dir, f"{input_name}_{prefix}{part_num}.mgf")
        
        print(f"Part {part_num}: {size} scans → {output_file}")
        
        # Write part with renumbered scans
        with open(output_file, "w") as f:
            for i, block in enumerate(part_blocks, 1):
                f.write(block.to_string(scan_number=i))
        
        output_files.append(output_file)
        start_idx = end_idx
    
    print(f"\nSuccessfully split into {num_parts} part(s)")
    return output_files


def main():
    parser = argparse.ArgumentParser(
        description="Sort MGF file by PEPMASS and optionally split into parts",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Sort only (no splitting)
  python sort_and_split_mgf.py input.mgf --output sorted.mgf --no-split
  
  # Sort and split into 3 parts
  python sort_and_split_mgf.py input.mgf --output-dir ./output --num-parts 3
  
  # Sort and split into 5 parts with custom prefix
  python sort_and_split_mgf.py input.mgf --output-dir ./output --num-parts 5 --prefix section
        """
    )
    
    parser.add_argument("input", help="Input MGF file")
    parser.add_argument("--output", help="Output file for sorted MGF (used with --no-split)")
    parser.add_argument("--output-dir", help="Output directory for split files")
    parser.add_argument("--num-parts", type=int, default=3, 
                       help="Number of parts to split into (default: 3)")
    parser.add_argument("--prefix", default="part", 
                       help="Prefix for split filenames (default: 'part')")
    parser.add_argument("--no-split", action="store_true",
                       help="Only sort, don't split")
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.no_split:
        if not args.output:
            parser.error("--output is required when using --no-split")
        
        # Just sort
        sort_mgf_by_pepmass(args.input, args.output)
    else:
        if not args.output_dir:
            parser.error("--output-dir is required when splitting")
        
        # Create temporary sorted file
        temp_sorted = os.path.join(args.output_dir, "_temp_sorted.mgf")
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Sort first
        sort_mgf_by_pepmass(args.input, temp_sorted)
        
        # Then split
        split_mgf(temp_sorted, args.output_dir, args.num_parts, args.prefix)
        
        # Clean up temp file
        os.remove(temp_sorted)
        print(f"\nCleaned up temporary file: {temp_sorted}")


if __name__ == "__main__":
    main()

 # Example of use:   
 # Example: python split_mgf.py C24_alkylamine_trihydroxy_input_falcon_charge.mgf --output-dir script_split_mgf --num-parts 2 --prefix section