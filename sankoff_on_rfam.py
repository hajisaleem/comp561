__author__ = 'SaleemSaddm'
'''
Created on: 28th November, 2013.

Comp 561 Project
Objective 4

Summary: Implementation of Sankoff Algorithm with penalties -1 intra group and -2 inter group and -2 for gaps.
The tree structure is given manually and the taxa read from a stockholm file and the output is stored in a file> family_name_out.txt

'''
from Bio import Phylo # just to display our tree
import string
from pprint import pprint

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

#Displaying the tree
tree = Phylo.read("tree_rfam.txt", "newick")
tree.rooted = True

#-------------------------------------------------------------
#Calculating the number of nodes and ancestors for the tree
with open("tree_rfam.txt", "r") as f:
    treeseq = f.readline()

treeseq = treeseq[:-2]
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
#reading the stockholm file and getting leaf nodes

tree_taxa = {}
name_map = {}

def parse_stockholm(file_path):
    ss_cons = ""
    sq = {}
    name_map = {}
    al = list(string.ascii_uppercase)
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

    return sq2, ss_cons, name_map

file_path = "RF00754_seed.stockholm.txt"
tree_taxa, ss_cons, name_map = parse_stockholm(file_path)

seq_length = len(tree_taxa["A"])

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
ancestor = {}
for key in ances_taxa.iterkeys():
    sequence = ""
    if str(key).isalpha():
        continue
    else:
        for i in xrange(0,seq_length):
            temp = ances_taxa[key][i]
            min_idx = temp.index(min(temp))
            sequence=sequence+mat_list[min_idx]
    ancestor[key] = sequence

#####################################################################################################
#creating the output file
outfile = file_path[:-19]+"out.txt"
with open(outfile, "w") as f_out:
    f_out.writelines("\tSankoff Algorithm on RNA Family "+file_path[:-19]+"\n")
    f_out.writelines("--------------------------------------------------------------------------------------------------------------------------\n")
    f_out.write("\tTree Sequence\n\n")
    f_out.writelines("\t"+treeseq+"\n")
    f_out.writelines("--------------------------------------------------------------------------------------------------------------------------\n")
    f_out.write("\tTree Structure\n")
    Phylo.draw_ascii(tree, file=f_out)
    f_out.write("\tp\t[lc, rc]\n\n")
    for key in tree_dict.iterkeys():
        f_out.write("\t"+str(key)+"\t"+str(tree_dict[key])+"\n")
    f_out.writelines("--------------------------------------------------------------------------------------------------------------------------\n")
    f_out.write("\tLeaf Nodes\n\n")
    for key in name_map.iterkeys():
        f_out.write("\t"+str(key)+"\t"+name_map[key]+"\n")
    f_out.writelines("--------------------------------------------------------------------------------------------------------------------------\n")
    f_out.write("\tLeaf Sequence\n\n")
    for key in tree_taxa.iterkeys():
        f_out.write("\t"+tree_taxa[key]+"\t"+str(key)+"\n")
    f_out.writelines("--------------------------------------------------------------------------------------------------------------------------\n")
    f_out.write("\tCost Matrix\n\n")
    pprint(mat_list, stream=f_out)
    pprint(cost_mat, stream=f_out)
    f_out.writelines("--------------------------------------------------------------------------------------------------------------------------\n")
    f_out.write("\tCalculated Ancestor Sequence\n\n")
    for key in ancestor.iterkeys():
        f_out.write("\t"+ancestor[key]+"\t"+str(key)+"\n")
    f_out.writelines("--------------------------------------------------------------------------------------------------------------------------\n")


#####################################################################################################





