#/usr/bin/env python
# 12 August 2014, St. Louis, MO
# This script is intended to take a file and split it into smaller pieces
# so that it can be uploaded to S3 (Amazon) storage.
# Sending the files to a specified S3 repository would be a nice addition.

# Usage: python splitter.py -i filename

# Needed for arguments below
import sys
# Needed for python call command to run UNIX split command
from subprocess import call

# The desired size of split files
# Current set to 50 MB as small files are easier to upload over a shaky VPN
size = '50'

# Take filename from standard input
# Filename will be used as the base name for the split files
# Number of split files is determined by input file size
flag = 0
for arg in sys.argv:
	if flag:
		filename_in = arg
		break
	if arg == '-i':
		flag = 1
sys.stdin = open(filename_in, 'r')

# Call the shell utility 'split' with argument '-b' for split size
# File names are take from the command line
call(['split', '-b', size + 'm', sys.stdin.name, sys.stdin.name + '-part-'])

