#!/usr/bin/env python3

# Peter L. Morrell - A tool for creating tarred and bzipped archives

import os
import sys
import tarfile

USAGE = """
Usage:
    python tarchiver2.py directory_name
    
Creates a compressed tar archive (.tar.bz2) of the specified directory.
"""

def make_tarfile(output_filename, source_dir):
    """Create a compressed tar archive of the specified directory."""
    with tarfile.open(output_filename, "w:bz2") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))
    print(f"Archive created: {output_filename}")

def main():
    if len(sys.argv) < 2:
        print(USAGE)
        sys.exit(1)
    
    try:
        directory_name = sys.argv[1]
        
        # Verify directory exists
        if not os.path.isdir(directory_name):
            print(f"Error: '{directory_name}' is not a directory or doesn't exist")
            sys.exit(1)
        
        output_file = f"{os.path.basename(directory_name)}.tar.bz2"
        make_tarfile(output_file, directory_name)
        
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()