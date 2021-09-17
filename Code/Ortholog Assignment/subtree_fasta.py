# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import re
import sys
from ete3 import PhyloTree
def species_names(tree_file):
    header_list=[]
    tree=PhyloTree(tree_file)
    for leaf in tree:
        leaf_name=leaf.name
        header_list.append(leaf_name)
    return header_list

def fasta_obtainer(protein_list,output_fasta):
    import os
    print ("number of proteins are %i" %len(protein_list))
    protein_count=0
    new_fasta=open(output_fasta,"w")
    for protein in protein_list:
        print (protein +" is being searched..")
        protein_info=protein.split("|")
        tax_ID=protein_info[3]

	#----------------------
	#This part should be specificed based on your fasta sequences database directory.
        for file in os.listdir("/cta/users/bselcuk/GPCRA/blast_documents/fasta_files/"):
    	    if file.endswith("_"+tax_ID+".fasta"):
                file=(os.path.join("/cta/users/bselcuk/GPCRA/blast_documents/fasta_files/", file))
                fasta=open(file,"r")
                print ("Proteome file: "+file)
                break

	#----------------------
        while 1:
            line=fasta.readline()
            if line=="":
                break
            if line[0]==">":
                line_list=line.split(" ")
            if line_list[0]==">{}|{}|{}".format(protein_info[0],protein_info[1],protein_info[2]): #matching our query with the protein in the database
                new_fasta.write(">"+protein+"\n")
                line=fasta.readline()
                while line[0]!=">":
                    new_fasta.write(line)
                    line=fasta.readline()
                    if line=="":
                        break
                protein_count+=1
                print (protein +" Done!")
                print ("Protein count= "+ str(protein_count)+"\n")            
                break
        fasta.close()
    new_fasta.close() 
    return ("FASTA file is produced!\nNumber of input sequences:{} \nNumber of sequences: {} \nNumber of missing sequences: {} ".format(len(protein_list),protein_count,len(protein_list)-protein_count))         


import argparse
ap = argparse.ArgumentParser(description="Script for obtaining FASTA files from a newick tree.")
ap.add_argument("--tree", required=True,help="The tree file in newick format.")
ap.add_argument("--out", required=True, help="The output directory and file name of the new FASTA file.")
args = vars(ap.parse_args())

tree=args["tree"]
out_fasta=args["out"]

sub_tree_names=species_names(tree)
print(fasta_obtainer(sub_tree_names,out_fasta))
