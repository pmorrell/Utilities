#/usr/bin/env python3

# 04 October 2018, Falcon Heights, MN
# Peter L. Morrell - tired of manually creating tarred and bzipped archives
# 24 April 2019, Falcon Heights, MN
# Returned to this and got it to work!

# Needed for arguments below
import os
import sys
import tarfile

Usage = """

Usage:
python archiver.py directory_name
"""

if len(sys.argv) < 1:
    print(Usage)
    exit(1)

#####
#   Defining arguments
#####
#   A description of the program
DESCR = """A Python tool for archiving directories. This was created to archive
old projects that are taking up disk space.

Pass a directory as the only argument. The original directory is not modified,
but a new file that is tarred and compressed with bzip2 is created. The tarred
archive takes the name of the directory."""
#   Create a new argument parser

try:
    directory_name = sys.argv[1]
except (OSError) as e:
    print(f"An error occured - {e}")


def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename, "w:bz2") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))


def main:
    make_tarfile()
