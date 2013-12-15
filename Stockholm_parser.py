__author__ = 'Saleem'

'''
Created on: 25th November, 2013.

Comp 561 Project
Objective 1

Summary: Read a Stockholm file and extract the sequence alignment and the consensus secondary structure.
'''

import string

#setting path to file
file_path = "RF00754_seed.stockholm.txt"

#consensus secondary structure
ss_cons = ""
#sequence_alignment
sq = {}

#opening file
with open(file_path, "r") as f_sto:

    #readin all lines of the file
    lines = f_sto.readlines()

    #iterating trhough all the lines
    for line in lines:
        if line.startswith("#=GC SS_cons"):
                #storing the lines with consensus secondary structure, stripping new line and white spaces
                # and concatenating them
                ss_cons += line.split()[2]
        #seleting all the lines with sequences
        if not line.startswith("#") and not line.startswith("/") and not line.startswith("\n"):
                #concatenating all the sequences in the
                if line.split()[0] in sq.keys():
                        sq[line.split()[0]] += line.split()[1]
                else:
                        sq[line.split()[0]] = line.split()[1]


print "Consensus Secondary Structure"
print ss_cons
print "Sequence Alignment"
for key in sq.keys():
    print sq[key]+"\t\t\t"+key

