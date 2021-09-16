# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 11:32:40 2020

@author: bselcuk
"""

def get_gene_clade(outgroup,tree,gene_header,delimiter,out_file):
    import ete3
    from ete3 import NCBITaxa,PhyloTree
    from statistics import stdev
    from statistics import mean
    import pickle
    from scipy.stats import ttest_ind
    from numpy import isnan

    tree=PhyloTree(tree)
    nodes_list=[]
    tree.set_outgroup(outgroup)
    target=tree.get_leaves_by_name(gene_header)
    node=target[0]
    count=0
    taxID_list=[]
    removed_list=[]
    while node!=None:
        repeating_taxID=[]
        existing_taxa=0
        count+=1
        print("Step:{}".format(count))
        node=node.up
        new_leafs=[]
        target=node.get_leaf_names()
        new_taxID=[]
        for leaf in target:
            if leaf in nodes_list or leaf in removed_list:
                continue
            else:
                nodes_list.append(leaf)
                new_leafs.append(leaf)
                leaf_list=leaf.split(delimiter)
                taxID=leaf_list[-1]
                if taxID not in taxID_list:
                    if taxID not in new_taxID:
                        new_taxID.append(taxID)
                else:
                    repeating_taxID.append(taxID)
                    if taxID=="9606":
                        print("List of repeating tax IDs:",repeating_taxID)
                        print("The search is finished!")
                        nodes_list = [x for x in nodes_list if x not in new_leafs]
                        tree.prune(nodes_list)
                        tree.write(outfile=out_file,format=2)
                        with open('outgroup.pkl', 'wb') as f:
                            pickle.dump(outgroup, f)
                        return tree
                        break
                    existing_taxa=1 
        if existing_taxa==1:
            print("List of repeating tax IDs:",repeating_taxID)
            nodes_list = [x for x in nodes_list if x not in new_leafs]
            removed_list+=new_leafs
        else:
            outgroup=new_leafs
            taxID_list+=new_taxID
            
def protein_name_obtainer(blast_result,num_human):
    import re
    blast=open(blast_result,"r")
    human_count=0
    while num_human!=human_count:
        line=blast.readline()
        line_list=line.split("\t")
        protein=line_list[1]
        fasta_info=line_list[-1]
        motif="^.+OX=([0-9]+).+$"
        motif_check=re.match(motif,fasta_info)
        if motif_check:
            taxID=motif_check.group(1)
            if taxID=="9606":
                human_count+=1
    if num_human==1:
        print("Target protein:",protein+"|9606")
    elif num_human==3:
        print ("Outgroup:",protein+"|9606")
    return "{}|9606".format(protein)

import argparse
ap = argparse.ArgumentParser(description="Script for obtaining species tree for a protein.")
ap.add_argument("--tree", required=True,help="The tree file in newick format.")
ap.add_argument("--out", required=True, help="The output directory and file name of the new tree file.")
ap.add_argument("--blastout", required=True,help="The blast result in tabular format.\
                Required for detecting the gene and the outgroup.")
args = vars(ap.parse_args())

tree=args["tree"]
blast_out=args["blastout"]
out_tree=args["out"]

print(get_gene_clade(protein_name_obtainer(blast_out,3),
                     tree,
                     protein_name_obtainer(blast_out,1),"|"
                     ,out_tree))