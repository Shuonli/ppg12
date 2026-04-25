#!/usr/bin/env python3
"""
Script to compare segment-number combinations between dst_fittinglist and dstLists directories.
Reads the content of each .list file and extracts segment-number patterns like 00047289-00001.
"""

import os
import re
from pathlib import Path

def extract_segment_numbers_from_files(directory):
    """
    Extract segment-number combinations from all .list files in a directory.
    Reads file contents and extracts patterns like 00047289-00001.

    Args:
        directory: Path to directory containing .list files

    Returns:
        set: Set of segment-number combinations found in all files
    """
    segments = set()

    if not os.path.exists(directory):
        print(f"Warning: Directory {directory} does not exist!")
        return segments

    # Pattern to match segment-number in file content (e.g., -00047289-00001.root)
    pattern = re.compile(r'-(\d{8}-\d{5})\.root')

    list_files = [f for f in os.listdir(directory) if f.endswith('.list')]

    for filename in list_files:
        filepath = os.path.join(directory, filename)
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line:
                        match = pattern.search(line)
                        if match:
                            segment_num = match.group(1)
                            segments.add(segment_num)
        except Exception as e:
            print(f"Error reading {filepath}: {e}")

    return segments

def main():
    # Base directory
    base_dir = "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/data/ana521"

    # Subdirectories
    dst_fitting_dir = os.path.join(base_dir, "dst_fittinglist")
    dst_lists_dir = os.path.join(base_dir, "dstLists")

    print("=" * 80)
    print("Segment-Number Comparison Tool")
    print("=" * 80)
    print()

    # Extract segments from both directories
    print("Extracting segment-numbers from dst_fittinglist files...")
    fitting_segments = extract_segment_numbers_from_files(dst_fitting_dir)
    print(f"Found {len(fitting_segments)} unique segment-number combinations in dst_fittinglist")
    print()

    print("Extracting segment-numbers from dstLists files...")
    dst_segments = extract_segment_numbers_from_files(dst_lists_dir)
    print(f"Found {len(dst_segments)} unique segment-number combinations in dstLists")
    print()

    # Find differences
    print("=" * 80)
    print("Analysis Results")
    print("=" * 80)
    print()

    # Segments in dst_fittinglist but not in dstLists
    only_in_fitting = fitting_segments - dst_segments
    print(f"Segment-numbers in dst_fittinglist but NOT in dstLists: {len(only_in_fitting)}")
    if only_in_fitting:
        sorted_fitting = sorted(only_in_fitting)
        print("First 50 entries:" if len(only_in_fitting) > 50 else "All entries:")
        for seg in sorted_fitting[:50]:
            print(f"  {seg}")
        if len(only_in_fitting) > 50:
            print(f"  ... and {len(only_in_fitting) - 50} more")
    print()

    # Segments in dstLists but not in dst_fittinglist
    only_in_dst = dst_segments - fitting_segments
    print(f"Segment-numbers in dstLists but NOT in dst_fittinglist: {len(only_in_dst)}")
    if only_in_dst:
        sorted_dst = sorted(only_in_dst)
        print("First 50 entries:" if len(only_in_dst) > 50 else "All entries:")
        for seg in sorted_dst[:50]:
            print(f"  {seg}")
        if len(only_in_dst) > 50:
            print(f"  ... and {len(only_in_dst) - 50} more")
    print()

    # Segments in both
    common_segments = fitting_segments & dst_segments
    print(f"Segment-numbers in BOTH directories: {len(common_segments)}")
    print()

    # Summary
    print("=" * 80)
    print("Summary")
    print("=" * 80)
    print(f"Total unique segment-numbers in dst_fittinglist: {len(fitting_segments)}")
    print(f"Total unique segment-numbers in dstLists:        {len(dst_segments)}")
    print(f"Common segment-numbers:                          {len(common_segments)}")
    print(f"Only in dst_fittinglist:                         {len(only_in_fitting)}")
    print(f"Only in dstLists:                                {len(only_in_dst)}")
    print()

    # Option to save results to file
    if only_in_fitting or only_in_dst:
        output_file = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/segment_differences.txt"
        with open(output_file, 'w') as f:
            f.write("Segment-Number Differences\n")
            f.write("=" * 80 + "\n\n")

            if only_in_fitting:
                f.write(f"In dst_fittinglist but NOT in dstLists ({len(only_in_fitting)} entries):\n")
                for seg in sorted(only_in_fitting):
                    f.write(f"{seg}\n")
                f.write("\n")

            if only_in_dst:
                f.write(f"In dstLists but NOT in dst_fittinglist ({len(only_in_dst)} entries):\n")
                for seg in sorted(only_in_dst):
                    f.write(f"{seg}\n")
                f.write("\n")

        print(f"Full difference list saved to: {output_file}")
        print()

if __name__ == "__main__":
    main()
