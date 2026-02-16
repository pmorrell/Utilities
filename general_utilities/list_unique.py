#!/usr/bin/env python3
# 04 October 2018, St. Paul, MN
# Peter L. Morrell -

# Needed for arguments below
import sys

Usage = """

Usage:
python inclusive_file exclusive_file out_directory
"""

if len(sys.argv) < 2:
    print(Usage)
    sys.exit(1)

#####
#   Defining arguments
#####
#   A description of the program
DESCR = """Compare files and remove overlaps."""
#   Create a new argument parser

def retain(ret):
    retain_temp = []
    with open(ret, 'r') as f:
        for line in f:
            retain_temp.append(line.strip().split())
            return retain_temp
def remove(rem):
    remove_temp = []
    with open(rem, 'r') as f:
        for line in f:
            remove_temp.append(line.strip().split())
            return remove_temp

def uniq_list(retain, remove):
    for element in remove:
        retain.remove(element)

    #return list(set(retain).symmetric_difference(set(remove)))
    return retain

def main(ret,rem):
    retained = retain(ret)
    removed = remove(rem)

    #rename = ()
#        temp = os.path.basename(rem)
#        rename = os.path.splitext(temp)
#    except:
#        print('Please pass name of file that will drop lines')

#    try:
#        directory_name = sys.argv[3]
#        print(directory_name)
#    except:
#        print('Please pass directory_name')

    unique = []
    unique = uniq_list(retained, removed)
    print (unique)

main(sys.argv[1], sys.argv[2])
