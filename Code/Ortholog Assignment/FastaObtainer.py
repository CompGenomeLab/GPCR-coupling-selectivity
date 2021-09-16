# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 16:31:40 2019

@author: Berkay Selcuk
"""
import re
import sys
import os

def protein_name_obtainer(blast_result,number_tax,tax_ID):
    conversion_dict={}
    protein_name_list=[]
    blast=open(blast_result,"r")
    human_count=0
    while number_tax!=human_count:
        line=blast.readline()
        line_list=line.split("\t")
        protein=line_list[1]
        fasta_info=line_list[-1]
        motif="^.+OX=([0-9]+).+$"
        motif_check=re.match(motif,fasta_info)
        if motif_check:
            taxID=motif_check.group(1)
            if taxID==tax_ID:
                human_count+=1
            conversion_dict[protein]=taxID
        if protein not in protein_name_list:
            protein_name_list.append(protein)
    # print (conversion_dict)
    return [protein_name_list,conversion_dict]

def fasta_obtainer(protein_list,ID_dict,output_file):
    duplicate_remove_dict={}
    print ("number of proteins are %i" %len(protein_list))
    protein_count=0
    duplicates=0
    new_fasta=open(output_file,"w")
    for protein in protein_list:
        check_string=""
        print (protein +" is being searched..")
        tax_ID=ID_dict[protein]
        
        if tax_ID not in duplicate_remove_dict:
            duplicate_remove_dict[tax_ID]=[]
            
        for file in os.listdir("/cta/users/bselcuk/GPCRA/blast_documents/fasta_files/"):
            if file.endswith("_"+tax_ID+".fasta"):
                file=(os.path.join("/cta/users/bselcuk/GPCRA/blast_documents/fasta_files/", file))
                fasta=open(file,"r")
                print ("Proteome file: "+file)
                break
        while 1:
            line=fasta.readline()
            if line=="":
                break
            if line[0]==">":
                line_list=line.split(" ")
            if line_list[0]==">"+protein: #matching our query with the protein in the database
                check_string+=">%s|%s\n"%(protein,tax_ID)
                line=fasta.readline()
                while line[0]!=">":
                    check_string+=line
                    line=fasta.readline()
                    if line=="":
                        break
                if check_string not in duplicate_remove_dict[tax_ID]:
                    duplicate_remove_dict[tax_ID].append(check_string)
                    new_fasta.write(check_string)
                    print (check_string)
                    protein_count+=1
                    print (protein +" Done!")
       	            print ("Protein count= "+ str(protein_count)+"\n")
                else:
                    duplicates+=1
                break
    fasta.close()
    new_fasta.close()
    return ("FASTA file is produced!\nNumber of input sequences:{} \nNumber of sequences: {} \nNumber of eliminated duplicates: {} ".format(len(protein_list),protein_count,duplicates))                            
                          
    
import argparse
ap = argparse.ArgumentParser(description="Script for obtaining FASTA files from a source.")
ap.add_argument("--blastout", required=True,help="The old FASTA file")
ap.add_argument("--num", required=True,help="The new FASTA file containing extra sequences")
ap.add_argument("--taxid", required=True,help="The FASTA headers to remove (newlne separated)")
ap.add_argument("--out", required=True, help="The output directory and file name of the new FASTA file.")
args = vars(ap.parse_args())

blast_out=args["blastout"]
num=args["num"]
taxid=args["taxid"]
out_fasta=args["out"]


prot_dict_list=protein_name_obtainer(blast_out,int(num),taxid)
print (fasta_obtainer(prot_dict_list[0],prot_dict_list[1],out_fasta))
