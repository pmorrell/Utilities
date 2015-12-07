#!/usr/bin/env python

#   This script is designed to modify SNP names, but could be used for any \
#   type of transformation between lists of names (by changing the rows \
#	or columns it works on.) \ 
#   Alchmey genotype files http://alchemy.sourceforge.net have a standard \
#   format with four rows of comments followed by and empty comment line \ 
#   a header, and then data in 12 columns. 

# ALCHEMY $Id: alchemy.c,v 1.5 2009/12/19 16:15:12 koni Exp $
# Program copyright 2008, 2009, 2010 Mark Hamilton Wright (Koni)
# BUILD DATE: Thu Feb  9 17:45:45 CST 2012
# RUN DATE: Tue Feb 21 15:22:11 2012
#
#snp	sample	AB call	ACGT call	P(call)	P(AA)	P(AB)	P(BB)	P(NC)	A int. (adj)	B int. (adj)	A int. (raw)	B int. (raw)
#11_10090	WBDC001	AA	CC	1.0000e+00	1.0000e+00	7.0487e-12	1.5024e-10	1.1041e-09	9.20	7.12	9.37	7.15	-125.0000
#11_10090	WBDC002	AA	CC	1.0000e+00	1.0000e+00	6.0809e-12	1.6993e-10	9.5722e-10	9.31	7.09	9.48	7.07	-125.0000
#11_10090	WBDC004	AA	CC	1.0000e+00	1.0000e+00	7.4498e-12	1.6085e-10	1.0466e-09	9.29	7.18	9.39	7.17	-125.0000

#   Original version appears to be from Ana Poets, but script provenance \
#   has been lost or forgotten. \
#   Current version modified by Peter L. Morrell, St. Paul, MN 
#   29 October 2015

import sys

Usage = """
Change SNP names in an Alchemy genotype file.
Two files are required:
1)The file should include Alchemy genotype data.
2)The second file should include a list of names, with new names in the 
1st column and old names in the 2nd column.
"""

if len(sys.argv) < 3:
	print Usage
else:

	in_file = sys.argv[1]
	names_file = sys.argv[2]
	alchemy_input = open(sys.argv[1],'r')
	# NamesFile=('name_test.txt')
	names = open(names_file, 'r')

	out_file = ("renamed_"+in_file)
	outfile = open(out_file, 'w')

# Create a dictionary
	file_dic = {}
	for line in names:
		# Remove if empty key
		line_new = line.strip('\n').split('\t')
		name_id = line_new[0]
		name_long = line_new[1]
		if name_long not in file_dic: file_dic[name_long] = []
		file_dic[name_long] = name_id

	for line in alchemy_input:
		snp = line.strip('\n').split('\t')
		snp = snp[0]
		if snp in file_dic:
			new_name = file_dic[snp]
		else:
			new_name = snp
		output = line.replace(snp, new_name)

		outfile.write(output)


alchemy_input.close()
names.close()
outfile.close()
