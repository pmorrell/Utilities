#/usr/bin/env python
# 12 August 2014, St. Louis, MO
# Peter L. Morrell - my first piece of Python code
# This script is intended to split a file into smaller pieces 
# so that it can be uploaded to S3 (Amazon) storage.
# Sending the files to a specified S3 repository would be a nice addition.
# The UNIX split command has this usage:
# split -b [size]m [filename] [filename]-part-
# The line above creates files, with size in megabytes using the specified file \
# appends '-part-' and then 'aa', 'ab', etc. in sequence

# Needed for arguments below
import sys
# To get arguments
import argparse
# Needed for 'call' command to run UNIX process
from subprocess import call

Usage = """

Usage:
python splitter.py -i filename -s [size]
"""
if len(sys.argv) < 2:
    print(Usage)
    exit(1)

#####
#   Defining arguments
#####
#   A description of the program
DESCR = """A Python tool for easier file splitting."""
#   Create a new argument parser
Arguments = argparse.ArgumentParser(description=DESCR, add_help=True)
#   Add some arguments
Arguments.add_argument('-i',
                '--infile',
                metavar='INPUT_FILE',
                type=argparse.FileType('r'),
                #default=sys.stdin,
                help='Input file',
                required=True)
#   Size of split files to create
Arguments.add_argument('-s',
                '--size',
                metavar='SIZE',
                default='50',
                type=int,
                help='Size of output files in MB',
                required=False)
#   And parse them
ParsedArgs = Arguments.parse_args()

# Call split with parameters defined above
# str() coerces the variable which is an integer to a string
call(['split', '-b', str(ParsedArgs.size) + 'm', ParsedArgs.infile.name, ParsedArgs.infile.name + '-part-'])
