#!/bin/env python

from Bio import SeqIO
import re
import fileinput
import os

INPUTDIR  = "."

list_of_files=[phd for phd in os.listdir(INPUTDIR) if phd.endswith(".phd.1")]

for input in list_of_files:
	print("%s %i" % (input.id, len(input)))
    SeqIO.write(input, "".join(str(input).split(".")""+"fq","fastq-sanger")

    	  #  output="".join(str(input).split(".")[:-1])+"_filtered.fasta"

#record = SeqIO.read("5Pepc_MorexR2", "abi")
#print("%s %i" % (record.id, len(record)))
#SeqIO.write(record, "5Pepc_MorexR2.fq", "fastq-sanger")
