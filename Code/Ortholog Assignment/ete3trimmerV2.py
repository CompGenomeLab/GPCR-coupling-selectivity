# -*- coding: utf-8 -*-
"""
Created on Sat Nov  9 23:58:01 2019

@author: Berkay Selcuk
"""

def alignment_score(x,y):
    from Bio.SubsMat import MatrixInfo as matlist
    from Bio import pairwise2
    matrix = matlist.blosum62
    alignments = pairwise2.align.globaldx(x,y,matrix,score_only=True)
    return alignments

def get_sequence(gene_of_interest,header,fasta_file):
    file=open(fasta_file,"r")
    a=0
    sequence=""
    while 1:
        line=file.readline()
        if line=="":
            return sequence
        if line.strip()==">"+header:
            a=1
            continue
        if a==1:
            if line[0]==">":
                return sequence
            sequence+=line[:-1]
    file.close()


def human_seq_find(gene_of_interest,fasta_file,t):
    for node in t.traverse():
        if node.is_leaf():   
            taxID=str(node).split("|")[-1] #Here you need to introduce TAX ID
            if taxID=="9606":
                gene=str(node).split("|")[-2]#This line and below line you need to define gene name
                if gene=="{}_HUMAN".format(gene_of_interest):
                    human_seq=get_sequence(gene_of_interest,str(node)[3:],fasta_file)
                    return human_seq
            
#To obtain all of the leaf nodes from the root
def distance_and_similarity(t,gene_of_interest):
    for node in t.traverse():
        total_leave_list=node.get_leaf_names()
        root=node
        break
    similarity_dict={}
    distance_dict={}
    common_lineage_dict={}
    for leaf in total_leave_list:
        distance=t.get_leaves_by_name(leaf)[0].get_distance(root)
        tar_seq=get_sequence(gene_of_interest,leaf,fasta_file)
        score=alignment_score(human_seq,tar_seq)
        similarity_dict[leaf]=score
        distance_dict[leaf]=distance
        taxID=leaf.split("|")[-1]
        common_taxa=common_two_lists(ncbi.get_lineage(9606),ncbi.get_lineage(int(taxID)))
        common_lineage_dict[leaf]=common_taxa
    return [distance_dict,similarity_dict,total_leave_list,common_lineage_dict]


def common_two_lists(list1,list2):
    length=min([len(list1),len(list2)])
    common=0
    for i in range(length):
        if list1[i]==list2[i]:
            common+=1
    return common

def node_analyze(node):
    children=node.get_children()
    if len(children)==1 or node.is_leaf():
        return [[],[],[],[],[]]
    leaf1=children[0].get_leaf_names()
    tax1=[]
    leaf2=children[1].get_leaf_names()
    tax2=[]
    for i in leaf1:
        if i in remove_list:
            tax1=[]
            break
        taxID=i.split("|")[-1]
        if taxID=="\n--":
            continue
        if taxID not in tax1:
            tax1.append(taxID)
            
    for a in leaf2:
        if i in remove_list:
            tax2=[]
            break
        taxID=a.split("|")[-1]
        if taxID=="\n--":
            continue
        if taxID not in tax2:
            tax2.append(taxID)
    tax_list=list(set(tax1)&set(tax2))
    return [tax_list,leaf1,leaf2,children[0],children[1]]
            
            
import argparse
ap = argparse.ArgumentParser(description="Script for obtaining FASTA files from a newick tree.")
ap.add_argument("--tree", required=True,help="The tree file in newick format.")
ap.add_argument("--out", required=True, help="The output directory and file name of the new FASTA file.")
ap.add_argument("--blastout", required=True,help="The blast result in tabular format.")
ap.add_argument("--fasta", required=True,help="The FASTA file of raw sequences.")
ap.add_argument("--protein", required=True,help="The FASTA file of raw sequences.")
args = vars(ap.parse_args())

tree=args["tree"]
blast_out=args["blastout"]
fasta_file=args["fasta"]
out_file=args["out"]
protein=args["protein"]
        
import numpy as np
from ete3 import PhyloTree
from scipy.stats import ttest_ind
from ete3 import NCBITaxa
ncbi = NCBITaxa()  
t=PhyloTree(tree,format=1)
import pickle
from numpy import isnan
with open('outgroup.pkl', 'rb') as f:
    outgroup_list = pickle.load(f)
node_list=[]
for leaf in t:
  if leaf.name not in outgroup_list:
    print(leaf.name)
    t.set_outgroup(leaf)
    print("First, control rooting")
    break

print("Second and real rooting")
if len(outgroup_list)>1:
    root =  t.get_common_ancestor(outgroup_list)
    t.set_outgroup(root)
else:
    t.set_outgroup(outgroup_list[0])
human_seq=human_seq_find(protein,fasta_file,t)
dicts=distance_and_similarity(t,protein)
similarity_dict=dicts[1]
distance_dict=dicts[0]
total_leave_list=dicts[2]
common_lineage_dict=dicts[3]
start_number=len(total_leave_list)

count=0
remove_list=[]
for node in t.traverse():
    duplicates=0
    node_results=node_analyze(node)
    tax_list=node_results[0]
    leaf1=node_results[1]
    leaf2=node_results[2]
    child1=node_results[3]
    child2=node_results[4]
    if tax_list!=[]:
        duplicates=1
    if duplicates==1:
        count+=1
        print("\n==========Comparison number: {}==========".format(count))
        print ("Number of common tax IDs between two clades:{}\nList of tax IDs:{}".format(len(tax_list), tax_list))
        child1_scores=[]
        child1_distances=[]
        child2_scores=[]
        child2_distances=[]
        child1_inter_distances=[]
        child2_inter_distances=[]
        child1_common_taxas=[]
        child2_common_taxas=[]
        for leaf in leaf1:
            taxID=leaf.split("|")[-1]
            common_taxa=common_two_lists(ncbi.get_lineage(9606),ncbi.get_lineage(int(taxID)))
            child1_common_taxas.append(common_taxa)
            score=similarity_dict[leaf]
            child1_scores.append(score)
            distance=distance_dict[leaf]
            child1_distances.append(distance)
            child1_inter_distances.append(t.get_leaves_by_name(leaf)[0].get_distance(child1))

            
        for leaf in leaf2:
            taxID=leaf.split("|")[-1]
            common_taxa=common_two_lists(ncbi.get_lineage(9606),ncbi.get_lineage(int(taxID)))
            child2_common_taxas.append(common_taxa)
            score=similarity_dict[leaf]
            child2_scores.append(score)
            distance=distance_dict[leaf]
            child2_distances.append(distance)
            child2_inter_distances.append(t.get_leaves_by_name(leaf)[0].get_distance(child2))

        if len(leaf1)<3 and len(leaf2)<3:
            print(node,"\n")
            leaf1_total_score=0
            leaf1_total_distance=0
            for leaf in leaf1:
                score=similarity_dict[leaf]
                distance=distance_dict[leaf]
                leaf1_total_score+=score
                leaf1_total_distance+=distance
            leaf1_avg_score=leaf1_total_score/len(leaf1)   
            leaf1_avg_dist=leaf1_total_distance/len(leaf1)
            print("Upper clade avg score: {}".format(leaf1_avg_score))
            print("Upper clade avg distance: {}".format(leaf1_avg_dist))
            leaf2_total_score=0
            leaf2_total_distance=0
            for leaf in leaf2:
                score=similarity_dict[leaf]
                distance=distance_dict[leaf]
                leaf2_total_score+=score
                leaf2_total_distance+=distance
            leaf2_avg_score=leaf2_total_score/len(leaf2)
            leaf2_avg_dist=leaf2_total_distance/len(leaf2)
            print("Lower clade avg score: {}".format(leaf2_avg_score))
            print("Lower clade avg distance: {}".format(leaf2_avg_dist))
            if leaf1_avg_score>leaf2_avg_score:
                for i in leaf2: 
                    remove_list.append(i)
                    print ("{} is deleted".format(i))
            elif leaf2_avg_score>leaf1_avg_score:
                for i in leaf1: 
                    remove_list.append(i)
                    print ("{} is deleted".format(i))
        else:
            rank_comparison_stat=ttest_ind(child1_common_taxas,child2_common_taxas,equal_var=False)
            score_stat=ttest_ind(child2_scores,child1_scores,equal_var=False)
            distance_stat=ttest_ind(child1_distances,child2_distances,equal_var=False)
            int_distance_stat=ttest_ind(child1_inter_distances,child2_inter_distances,equal_var=False)
            print(node)
            print()
            print("Parental branch length upper clade:",child1.dist)
            print("Parental branch length lower clade:",child2.dist)
            print("Taxa comparsion statistics: {}  ====  {}".format(rank_comparison_stat[0],rank_comparison_stat[1]))
            print("The similarity statistics: {}  ====  {}".format(score_stat[0],score_stat[1]))
            print("The distance statistics: {}  ====  {}".format(distance_stat[0],distance_stat[1]))
            print("The int-distance statistics: {}  ====  {}".format(int_distance_stat[0],int_distance_stat[1]))
            print("========")
            if score_stat[1]<=0.1:
                if score_stat[0]>0:
                    if (rank_comparison_stat[0]<0 and rank_comparison_stat[1]<0.1) or rank_comparison_stat[1]>0.01 or isnan(rank_comparison_stat[1]):
                        for i in leaf1:
                            remove_list.append(i)
                            print ("{} is deleted".format(i))
                elif score_stat[0]<0 :
                    if (rank_comparison_stat[0]>0 and rank_comparison_stat[1]<0.1) or rank_comparison_stat[1]>0.01 or isnan(rank_comparison_stat[1]):
                        for i in leaf2: 
                            remove_list.append(i)
                            print ("{} is deleted".format(i))
                
            else:
                print("Clades are not significantly different from each other.")

for leaf in remove_list:
    total_leave_list.remove(leaf)
t.prune(total_leave_list)
print("******************************")
print ("Trimming is done!\nInput number of sequences: {}\nNumber of eliminated paralogs:{}".format(start_number,start_number-len(total_leave_list)))
t.write(outfile=out_file,format=1)
print ("Ortholog tree is saved:{}".format(out_file))
