#/usr/bin/env python3

# 04 October 2018, St. Paul, MN
# Peter L. Morrell - tired of manually creating tarred and bzipped archives

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
DESCR = """A Python tool for archiving directories."""
#   Create a new argument parser

try:
    directory_name = sys.argv[1]
    print(directory_name)
except:
    print('Please pass directory_name')


with tarfile.open( directory_name + ".tar.bz2", "w:bz2" ) as tar:
    name = os.listdir( directory_name )
        #print(name)
        tar.add(name)
    tar.close()
