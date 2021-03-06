__author__ = 'SaleemSaddm'
'''
Created on: 15th December, 2013.

Comp 561 Project
Objective 8

Summary: Implementation of Sankoff Algorithm with penalties -1 intra group and -2 inter group and -2 for gaps.
We try to preserve the base pairing dependencies of the consensus structure by forcing the nucleotide at the closing binding index to be something that can form a pair with
corresponding opening binding index.

The tree structure is automatically generated and the taxa read from a stockholm file and the output is stored in a file compare.txt

'''
import string
from pprint import pprint
import sys
import os
import re


#####################################################################################################
#creating cost matrix

mat_list = ["A", "C", "G", "U", "."]
purines = ["A", "G"]
pyrimidines = ["C", "U"]

n = len(mat_list)
#creating the cost matrix (all entries are cost, and therefore would later be substracted.)
cost_mat = [[0,2,1,2,2],[2,0,2,1,2],[1,2,0,2,2],[2,1,2,0,2],[2,2,2,2,0]]
#the cost matrix is updated

#####################################################################################################
#giving the tree structure manually

#reading the stockholm file and getting leaf nodes

tree_taxa = {}
name_map = {}

def parse_stockholm(file_path):
    ss_cons = ""
    sq = {}
    name_map = {}
    al = list(string.ascii_letters)
    sq2 = {}
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
                        name_map[al[0]] = line.split()[0]
                        al.pop(0)
        for key in name_map.iterkeys():
            sq2[key] = sq[name_map[key]]

    num = len(sq2)

    #converting to paranthesis format
    ss_cons = re.sub("<", "(", ss_cons)
    ss_cons = re.sub(">", ")", ss_cons)
    return sq2, ss_cons, name_map, num

file_path = "RF00654_seed.stockholm.txt"


tree_taxa, ss_cons, name_map, num = parse_stockholm(file_path)

seq_length = len(tree_taxa["a"])

def parseRNAfold(seq):
    dictRNAfold = {}
    stack = []
    justlist = []
    i = 0
    for char in seq:
        if char == "(":
            stack.append(i)
            justlist.append(i)
        if char == ")":
            dictRNAfold[i] = stack.pop(-1)
            justlist.append(i)
        i+=1
    return dictRNAfold, justlist

dictRNAfold, justlistofindices = parseRNAfold(ss_cons)

#gives the indices of the binding pair

#-------------------------------------------------------------
#Calculating the number of nodes and ancestors for the tree

alphalist = list(string.ascii_letters)

treeseq = ""
for i in xrange(1,num):
    treeseq = treeseq+"("
treeseq = treeseq+alphalist[0]+","
alphalist.pop(0)
for i in xrange(1,num):
    treeseq = treeseq+alphalist[0]+"),"
    alphalist.pop(0)


treeseq = treeseq[:-1]


treeseqlist = list(treeseq)

taxa_num = 0 #leaf nodes
for char in treeseqlist:
        if char.isalpha():
                taxa_num+=1

#-------------------------------------------------------------
# converting tree sequence into hierarchy

tree_dict = {}

def parse_tree(seqlist):
    close = 0
    i = 0
    for c in seqlist:
        if c == ")":
             close = i
             break
        i+=1
    return close

def join_node(tree_dict, seqlist, close, num):
        open = close-4
        seqlist[open] = num
        tree_dict[num] = [seqlist[open+1],seqlist[open+3]]
        del seqlist[open+1: open+5]
        return seqlist, tree_dict

for i in xrange(1,taxa_num):
        close = parse_tree(treeseqlist)
        treeseqlist, tree_dict = join_node(tree_dict, treeseqlist, close, i)

#-------------------------------------------------------------


#####################################################################################################
#calculaitng ancestral sequences

ances_taxa = {}

#-------------------------------------------------------------
#initializing the leaf nodes

for taxa in tree_taxa.iterkeys(): #for all leaf nodes
        tempcal = []
        for neuclotide in tree_taxa[taxa]: #for all the neucleotides in the sequence
            allnuc = []
            for k in xrange(0,n):
                if k == mat_list.index(neuclotide):
                    allnuc.append(0)
                else:
                    allnuc.append(float("inf"))
            tempcal.append(allnuc)
        ances_taxa[taxa] = tempcal

#---------------------------------------------------------------
#the main algorithm

for i in xrange(1,taxa_num): #for each ancestor
    lc = tree_dict[i][0]
    rc = tree_dict[i][1]
    leftchild = ances_taxa[lc]
    rightchild = ances_taxa[rc]

    ances = []

    for j in xrange(0,seq_length): #for each column of the sequence
        col_cost = []
        for k in range(0, n): #for each neucloetide possible
                # min-cost change from left child
                min_left = float("inf")
                for l in range(0, n): #from each neucleotide possible
                    this_cost = cost_mat[k][l] + leftchild[j][l]
                    min_left = min(min_left, this_cost)
                # min-cost change from right child
                min_right = float("inf")
                for l in range(0, n): #from each neucleotide possible
                    this_cost = cost_mat[k][l] + rightchild[j][l]
                    min_right = min(min_right, this_cost)
                col_cost.append(min_left + min_right)
        ances.append(col_cost)
    ances_taxa[i] = ances

#---------------------------------------------------------------
#Getting the ancestral sequence

#We have to make sure of two things:
#1. There are no gaps at the indices of binding pair
#2. The closing binding pair is aligned with opening binding pair, i.e., A-U and G-C

allowed_pair = {}
allowed_pair["A"] = "U"
allowed_pair["U"] = "A"
allowed_pair["G"] = "C"
allowed_pair["C"] = "G"

ancestor = {}
for key in ances_taxa.iterkeys():
    sequence = ""
    if str(key).isalpha():
        continue
    else:
        for i in xrange(0,seq_length):
            temp = ances_taxa[key][i]
            if i in justlistofindices:#removing the possibilty of a gap at a binding index
                temp.pop(-1)
            if i in dictRNAfold.iterkeys():
                openingidx = dictRNAfold[i]
                sequence = sequence+allowed_pair[sequence[openingidx]]
            else:
                min_idx = temp.index(min(temp))
                sequence=sequence+mat_list[min_idx]
    ancestor[key] = sequence

#####################################################################################################
#Python Wrapper for RNAfold and

##Reference: Adapted from Vienna Wrappers for Python by Yann Ponty

RNAFOLDwrapper = "RNAfold"
RNADISTwrapper = "RNAdistance"

def runRNAFold(seq):
  (tmpFileOut,tmpFileIn,struct,energy) = ("tmpOut.dat","tmpIn.dat",None,None)
  inFile = open(tmpFileIn,"w")
  inFile.write(seq+"\n")
  inFile.close()
  os.system("%s > %s < %s"%(RNAFOLDwrapper,tmpFileOut,tmpFileIn))
  lineno = 0
  for l in open(tmpFileOut,"r"):
    if lineno==1:
      data = l[:-1].split()
      struct = data[0]
      energy = float(" ".join(data[1:])[1:-1])
    lineno += 1
  os.remove("rna.ps")
  os.remove(tmpFileOut)
  os.remove(tmpFileIn)
  return (struct,energy)

def runRNADistance(seq1, seq2):
  (tmpFileOut,tmpFileIn,dist) = ("tmpOut.dat","tmpIn.dat",None)
  inFile = open(tmpFileIn,"w")
  inFile.write(seq1+"\n")
  inFile.write(seq2+"\n")
  inFile.close()
  os.system("%s > %s < %s"%(RNADISTwrapper,tmpFileOut,tmpFileIn))
  with open(tmpFileOut,"r") as f_temp:
    l = f_temp.readline()
    dist = l[3:]
  os.remove(tmpFileOut)
  os.remove(tmpFileIn)
  return dist


#---------------------------------------------------------------
foldstruct = {}
diststruct = {}
#Running for every ancestor
distsum = 0
GCcont = 0
for key in ancestor.iterkeys():
    a,b = runRNAFold(ancestor[key])
    c = runRNADistance(ss_cons, a)
    foldstruct[key] = a
    diststruct[key] = c
    distsum += int(c)
    temp_seq = ancestor[key]
    for char in temp_seq:
        if char == "G" or char == "C":
            GCcont+=1


#####################################################################################################
#creating the output file
outfile = "compare.txt"
with open(outfile, "a") as f_out:
    f_out.write(file_path[:-19])
    f_out.write("\t"+str(num))
    f_out.write("\t"+str(seq_length))
    f_out.write("\t"+str(distsum))
    f_out.write("\t"+str(GCcont)+"\n")



#####################################################################################################





