#!/usr/bin/env python3
# Peter L. Morrell - St. Paul, MN 17 October 2016
# Improved version with error handling and directory flexibility

import os
import sys

"""A script to rename files, especially FASTQ files downloaded 
from the NCBI Sequence Read Archive."""

Usage = """
file_rename.py - version 2
Reads in a file containing a list of old and new file names 
in two columns. It checks if the file exists in the specified 
(or current) directory and renames it if permissions allow.

Usage:
    file_rename.py [file_list.txt] [optional_directory]
"""

# Check for correct arguments
if len(sys.argv) < 2:
    print(Usage)
    sys.exit(1)

file_list = sys.argv[1]
directory = sys.argv[2] if len(sys.argv) > 2 else os.getcwd()

# Ensure the directory exists
if not os.path.exists(directory):
    print(f"Error: Directory '{directory}' does not exist.")
    sys.exit(1)

# Create a dictionary to hold old and new file names
file_data = {}
try:
    with open(file_list) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            names = line.strip().split()
            if len(names) != 2:
                print(f"Skipping invalid line: {line.strip()}")
                continue
            file_data[names[0]] = names[1]
except IOError as ex:
    print("File open issue:", ex.strerror)
    sys.exit(1)

# Get list of files in the directory
files = os.listdir(directory)

# Filter files that need renaming
re_names = {k: file_data[k] for k in files if k in file_data}

# Rename files with error handling
for old_name, new_name in re_names.items():
    try:
        os.rename(os.path.join(directory, old_name), os.path.join(directory, new_name))
        print(f"Renamed: {old_name} -> {new_name}")
    except OSError as ex:
        print(f"Error renaming '{old_name}': {ex.strerror}")

print("File renaming complete.")
