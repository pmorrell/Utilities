#!/usr/bin/env python3
"""A script to rename files, especially fastq files downloaded from the \
NCBI Sequence Read Archive."""

import os
import sys

Usage = """ Provide a file name """

# Expects a filename as the only argument
if len(sys.argv) < 1:
    print (Usage)
    exit(1)

#file_data = []
file_data = {}
#   Read the file in line-by-line
try:
    with open(sys.argv[1]) as f:
        for line in f:
        #   Skip the header lines
            if line.startswith('#'):
                continue
            else:
        #   file_data.append(line.strip().split('\t'))
                names = line.split()
                file_data[names[0]] = names[1]
except IOError as ex:
	print ("File open issue: " + ex.strerror)

cwd = os.getcwd()

files = []
files = os.listdir(cwd)

re_names = {k: file_data[k] for k in files if k in file_data}


for i, j in re_names.items():
    #os.rename(str(list(re_names.keys())), str(list(re_names.values())))
    os.rename(i, re_names[i])
