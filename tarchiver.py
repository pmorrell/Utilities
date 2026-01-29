#!/usr/bin/env python3

# Peter L. Morrell - A tool for creating tarred and bzipped archives

import os
import sys
import tarfile
import subprocess
import shutil
import re

USAGE = """
Usage:
    python tarchiver.py directory_name
    
Creates a compressed tar archive (.tar.bz2) of the specified directory.
"""

def make_tarfile(output_filename, source_dir):
    """Create a compressed tar archive using pbzip2 if available."""
    # Check if pbzip2 is available for faster parallel compression
    if shutil.which('pbzip2'):
        print("Using pbzip2 for parallel compression...")
        try:
            # Use tar with pbzip2 for parallel compression
            cmd = ['tar', '-c', '-C', os.path.dirname(source_dir) or '.',
                   os.path.basename(source_dir)]
            with open(output_filename, 'wb') as f:
                tar_proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
                pbzip2_proc = subprocess.Popen(['pbzip2'],
                                               stdin=tar_proc.stdout,
                                               stdout=f)
                tar_proc.stdout.close()
                pbzip2_proc.communicate()

            if pbzip2_proc.returncode != 0:
                raise subprocess.CalledProcessError(pbzip2_proc.returncode,
                                                    'pbzip2')
            print(f"Archive created: {output_filename}")
        except Exception as e:
            print(f"pbzip2 compression failed, falling back to bzip2: {e}")
            # Fall back to standard compression
            make_tarfile_standard(output_filename, source_dir)
    else:
        print("pbzip2 not found, using standard bzip2 compression...")
        print("Tip: Install pbzip2 for faster compression:")
        print("  brew install pbzip2")
        make_tarfile_standard(output_filename, source_dir)


def make_tarfile_standard(output_filename, source_dir):
    """Create a compressed tar archive using standard bzip2."""
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
            print(f"Error: '{directory_name}' is not a directory or "
                  f"doesn't exist")
            sys.exit(1)

        # Write output to the same directory as the input
        abs_dir = os.path.abspath(directory_name)
        base = os.path.basename(abs_dir)
        # Sanitize: replace spaces and non-alphanumeric with underscores
        safe_base = re.sub(r"[^A-Za-z0-9._-]", "_", base)
        # Collapse multiple underscores
        safe_base = re.sub(r"_+", "_", safe_base).strip("_")
        output_file = os.path.join(
            os.path.dirname(abs_dir),
            f"{safe_base}.tar.bz2"
        )
        make_tarfile(output_file, directory_name)

    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
