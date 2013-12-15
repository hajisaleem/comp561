__author__ = 'SaleemSaddm'
'''
Created on: 27th November, 2013.

Comp 561 Project
Objective 3

Summary: Implementation of Sankoff Algorithm with gaps with penalties -1 intra group and -2 inter group and -2 for gaps.
The cost matrix is no 5x5.
The tree structure is given manually and the taxa are also given manually
'''
from Bio import Phylo # just to display our tree


mat_list = ["A", "C", "G", "U", "."]
purines = ["A", "G"]
pyrimidines = ["C", "U"]

n = len(mat_list)
#creating the cost matrix (all entries are cost, and therefore would later be substracted.)
cost_mat = [[0,2,1,2,2],[2,0,2,1,2],[1,2,0,2,2],[2,1,2,0,2],[2,2,2,2,0]]
#the cost matrix is updated

#Displaying the tree
tree = Phylo.read("tree.txt", "newick")
tree.rooted = True
Phylo.draw(tree)


#Calculating the number of nodes and ancestors for the tree
with open("tree.txt", "r") as f:
    treeseq = f.readline()

treeseq = treeseq[:-2]
print treeseq
treeseqlist = list(treeseq)

taxa_num = 0 #leaf nodes
for char in treeseqlist:
        if char.isalpha():
                taxa_num+=1

tree_dict = {}
tree_taxa = {}

#degining the sequence here
tree_taxa['A'] = "C"
tree_taxa['B'] = "A"
tree_taxa['C'] = "C"
tree_taxa['D'] = "A"
tree_taxa['E'] = "G"

seq_length = len(tree_taxa["A"])

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

'''
print tree_dict  # here we have converted the tree from a sequence to a hierarchy
print tree_taxa
print seq_length
'''
#calculaitng ancestral sequences

ances_taxa = {}

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


for i in xrange(1,taxa_num): #for each ancestor
    lc = tree_dict[i][0]
    rc = tree_dict[i][1]
    leftchild = ances_taxa[lc]
    rightchild = ances_taxa[rc]

    ances = []

    for j in xrange(0,seq_length): #for each column of the sequence
        col_cost = []
        for k in range(0, n): #for each neucloetide possible
                # Find min-cost change from left child
                min_left = float("inf")
                for l in range(0, n): #from each neucleotide possible
                    this_cost = cost_mat[k][l] + leftchild[j][l]
                    min_left = min(min_left, this_cost)
                # Find min-cost change from node on the right
                min_right = float("inf")
                for l in range(0, n):
                    this_cost = cost_mat[k][l] + rightchild[j][l]
                    min_right = min(min_right, this_cost)
                col_cost.append(min_left + min_right)
        ances.append(col_cost)

    ances_taxa[i] = ances


print ances_taxa






