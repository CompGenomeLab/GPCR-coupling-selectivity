# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 14:41:27 2019

@author: Berkay Selcuk
"""
from ete3 import PhyloTree

import argparse
ap = argparse.ArgumentParser(description="Script for ordering tree (Human at the top).")
ap.add_argument("--tree", required=True,help="The tree file in newick format.")
ap.add_argument("--out", required=True, help="The output directory and file name of the new FASTA file.")
args = vars(ap.parse_args())

tree=args["tree"]
out_file=args["out"]

sample=tree
t=open(sample,"r")
line=t.readline()
tree=PhyloTree(line)
for node in tree.traverse():
    leaves=node.get_children()
    if len(leaves)<2 or node.is_leaf():
        continue
    leaf2=leaves[1].get_leaves()
    for i in leaf2:
        leaf_node=str(i).split("|")[-1]
        #print (leaf_node)
        if leaf_node=="9606":
            node.swap_children()
            break
tree.write(outfile=out_file)
print (tree)
