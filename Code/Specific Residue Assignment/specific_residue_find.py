# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 19:21:41 2020

@author: bselcuk
"""
from Bio.SubsMat import MatrixInfo as matlist


def blosum_score(most_freq_aa,target_aa):
    blosum80=matlist.blosum80
    key=(most_freq_aa,target_aa)
    if key in blosum80:
        return blosum80[key]
    else:
        key=(target_aa,most_freq_aa)
        return blosum80[key]

def conservation_check(place_num,align_dict,tax,similarity,min_blossum_score):
    place_num=int(place_num)
    place_count=0
    residue_count=0
    result_aa_list=[]
    result_aa_number=[]
    total=0
    dash=0
    
    for seq in align_dict:
        if seq.split("|")[-1]==tax: #Finds first seq with certain taxID
            reference_seq=align_dict[seq]
            break
    while place_count<place_num:
        amino_acid=reference_seq[place_count]
        residue_control=0
        if amino_acid!="-":
            residue_count+=1
            residue_control=1
        place_count+=1
                
    for seq in align_dict:  
        total+=1
        target_aa=align_dict[seq][place_num-1]#Fix this problemm!!!!
        if target_aa=="-":
            dash+=1
        if target_aa not in result_aa_list:
            result_aa_list.append(target_aa)
            result_aa_number.append(0)
        idx=result_aa_list.index(target_aa)
        result_aa_number[idx]+=1
    max_num=max(result_aa_number)
    idx_max=result_aa_number.index(max_num)
    aa=result_aa_list[idx_max]
    
    if aa=="-":
        return [0,"-",-1,1000,"none"]
    elif similarity:
        for aminoacid in result_aa_list:
            if aminoacid!="-" and aminoacid!=aa[0]:
                score=blosum_score(aa[0],aminoacid)
                if score>=min_blossum_score:
                    idx_similar=result_aa_list.index(aminoacid)
                    max_num+=result_aa_number[idx_similar]
                    aa+="/"+aminoacid
                    
    if residue_control==0:
        residue_count=-1
    return [max_num*100/(total-dash),aa,residue_count,100*dash/total,total-dash]

def alignment_dict(MA_file):
    result_dict={}
    #alignment_file=open(MA_file,"r")
    with open(MA_file) as alignment_file:
        for line in alignment_file:
            if line=="":
                break
            if line[0]==">":
                header=line.strip()
                result_dict[header]=""
            else:
                result_dict[header]+=line.strip()
    length=len(result_dict[header])
    alignment_file.close()
    return (result_dict,length)

def spial(alignment,groupA,groupB,specificity_threshold,conservation_threshold,
          structral_or_sequence,lower_cutoff,similarity,PDB_taxID,min_blossum_score):
    A_specific_dict={}
    B_specific_dict={}
    C_specific_dict={}
    AB_specific_dict={}
    threshold_dict={}
    cons_and_spec_dict={}
    A_conserved_dict={}
    B_conserved_dict={}
    total_group=[]
    for a in groupA:
        total_group.append(a)
    for b in groupB:
        total_group.append(b)
    
    alignment_d=alignment_dict(alignment)
    align_length=alignment_d[1]
    alignment_d=fasta_divider(total_group,alignment_d[0])
    alignmentA_dict=fasta_divider(groupA,alignment_d)
    alignmentB_dict=fasta_divider(groupB,alignment_d)
    for i in range(align_length):
        
        result_A=conservation_check(i+1,alignmentA_dict,PDB_taxID,similarity,min_blossum_score)
        conservationA=result_A[0]
        aa_A=result_A[1]
        residue_num_A=result_A[2]
        gap_percentage_A=result_A[3]
        total_num_A=result_A[4]
        gene_A=groupA[0]
        
        result_B=conservation_check(i+1,alignmentB_dict,"9606",similarity,min_blossum_score)
        conservationB=result_B[0]
        aa_B=result_B[1]
        residue_num_B=result_B[2]
        gap_percentage_B=result_B[3]
        total_num_B=result_B[4]
        gene_B=groupB[0]
        
        result_C=conservation_check(i+1,alignment_d,"9606",similarity,min_blossum_score)
        conservationC=result_C[0]
        aa_C=result_C[1]
        gap_percentage_C=result_C[3]
        total_num_C=result_C[4]

        if structral_or_sequence=="structural":
            placeA=GPCR_convert(keys[i],gene_A)
            placeB=GPCR_convert(keys[i],gene_B)
            keyA=keys[i]
            keyB=keyA
            
        elif structral_or_sequence=="sequence":
            placeA=residue_num_A
            placeB=residue_num_A
            keyA=GPCR_convert(placeA,gene_A)
            #keyB=GPCR_convert(placeB,gene_B)
            keyB=keyA
            #if you enter residue number we get GPCRdb number
        if keyA=="-":
            keyC=keyB
        else:
            keyC=keyA
        
        a=0
        b=0
        c=0
        placeC="{}".format(placeA)#,placeB)
        if aa_A=="-" and aa_B!="-":
            if conservationB>=specificity_threshold:# and gap_percentage_A<=20:
                b=1
        elif aa_A!="-" and aa_B=="-":
            if conservationA>=specificity_threshold:# and gap_percentage_A<=20:
                a=1
        elif aa_A=="-" and aa_B=="-":
            continue
        elif conservationC>=conservation_threshold and a==0 and b==0:# and gap_percentage_C<=10:
            c=1
            if conservationC not in C_specific_dict:
                C_specific_dict[conservationC]=[]
            C_specific_dict[conservationC].append([placeC,aa_C,i+1,keyC,gap_percentage_C,total_num_C])
        
        elif blosum_score(aa_A[0],aa_B[0])<min_blossum_score:
            if conservationA>=specificity_threshold:# and gap_percentage_A<=20:
                a=1
            if conservationB>=specificity_threshold:# and gap_percentage_B<=20:
                b=1
        elif blosum_score(aa_A[0],aa_B[0])>=min_blossum_score:
            if conservationA>=specificity_threshold:# and gap_percentage_A<=20:
                if conservationB<lower_cutoff:
                    a=1
                else:
                    a=2
            if conservationB>=specificity_threshold:# and gap_percentage_B<=20:
                if conservationA<lower_cutoff:
                    b=1
                else:
                    b=2
        if a==1:
            if conservationA not in A_specific_dict:
                A_specific_dict[conservationA]=[]
            A_specific_dict[conservationA].append([placeA,aa_A,i+1,keyA,gap_percentage_A,total_num_A])
        if b==1:
            if conservationB not in B_specific_dict:
                B_specific_dict[conservationB]=[]
            B_specific_dict[conservationB].append([placeB,aa_B,i+1,keyB,gap_percentage_B,total_num_B])
        if a==1 and b==1:
            if conservationB not in AB_specific_dict:
                AB_specific_dict[conservationB]=[]
            AB_specific_dict[conservationB].append([placeB,aa_B,i+1,keyB,gap_percentage_B,total_num_B])
            if conservationA not in AB_specific_dict:
                AB_specific_dict[conservationA]=[]
            AB_specific_dict[conservationA].append([placeB,aa_A,i+1,keyA,gap_percentage_A,total_num_A])
        if a==2:
            if conservationA not in A_conserved_dict:
                A_conserved_dict[conservationA]=[]
            A_conserved_dict[conservationA].append([placeA,aa_A,i+1,keyA,gap_percentage_A,total_num_A])
        if b==2:
            if conservationB not in B_conserved_dict:
                B_conserved_dict[conservationB]=[]
            B_conserved_dict[conservationB].append([placeB,aa_B,i+1,keyB,gap_percentage_B,total_num_B])
       
        #The ones that are below every threshold
        if a==0 and b==0 and c==0:
            if conservationA not in threshold_dict:
                threshold_dict[conservationA]=[]
            if conservationB not in threshold_dict:
                threshold_dict[conservationB]=[]
            threshold_dict[conservationA].append([placeA,aa_A,i+1,keyA,gap_percentage_A,total_num_A])
            threshold_dict[conservationB].append([placeB,aa_B,i+1,keyB,gap_percentage_B,total_num_B])
        if c==1 or (a==1 and b==1):
            if conservationA not in cons_and_spec_dict:
                cons_and_spec_dict[conservationA]=[]
            cons_and_spec_dict[conservationA].append([placeA,aa_A,i+1,keyA,gap_percentage_A,total_num_A])
    return [A_specific_dict,B_specific_dict,C_specific_dict,AB_specific_dict
            ,A_conserved_dict,B_conserved_dict,threshold_dict,cons_and_spec_dict]

def spial_result_parse(result_dict,structural_or_sequence):
    conservation_list=list(result_dict.keys())
    conservation_list.sort(reverse=True)
    PyMOLinput=""
    output=""
    output_list=[]
    alignment_position_list=[]
    output+="|Conservation|Residue Number|Amino Acid|MSA Position|GPCRdb Number|Gap Percentage|# of Sequences|\n"
    output+="| ------------ | ------------ | ------------ | ------------ | ------------ | ------------ | ------------ |\n"
    for value in conservation_list:
        for result in result_dict[value]:
            res_no=result[0]
            place_in_alignment=result[2]
            
            if res_no==-1:
                continue
            PyMOLinput+=str(int(res_no))+"+"
            output_list.append(int(res_no))
            alignment_position_list.append(place_in_alignment)
            if structural_or_sequence=="structural":
                output+="{:.2f}%\t{}\t{}\t{}\t{}\t{:.2f}\t{}\n".format(value,result[0],result[1],result[2],result[3],result[4],result[5])
            else:    
                output+="|{:.2f}%|{}|{}|{}|{}|{:.2f}|{}|\n".format(value,result[0],result[1],result[2],result[3],result[4],result[5])
    print ("-----------------------------")
    output+=PyMOLinput[:-1]+"\n"
    return output,output_list,alignment_position_list

def GPCR_convert(res,gene):
    import pickle
    gene=gene.lower()
    try:
        with open(r"GPCR_notations\{}_notation_dict".format(gene), 'rb') as handle:
            conv_dict = pickle.load(handle)
        
        if type(res)==str:# GPCRdb Numbering ---> Residue Numbering
            if res not in conv_dict:
                return "-"
            return conv_dict[res]
        if type(res)==int: # Residue number ---> GPCRdb Numbering
            for key in conv_dict:
                if int(conv_dict[key])==res:
                    return key
            return ("-")
    except:
        return ("-")

def fasta_divider(gene_list,align_dict):
    result_dict={}
    for gene in gene_list:
        count=0
        for seq in align_dict:
            if seq.split("|")[-2]=="{}_HUMAN".format(gene):
                result_dict[seq]=align_dict[seq]
                count=1
                continue
            if count==1:
                if seq.split("|")[-1]=="9606":
                    break
                else:
                    result_dict[seq]=align_dict[seq]
    return result_dict

def place_find(taxID,alignment_dict,residue_number):
    place_count=0
    residue_count=0
    search=0
    for seq in alignment_dict:
        if seq.split("|")[-1]==taxID: #Finds first seq with certain taxID
            reference_seq=alignment_dict[seq]
            search=1
            break
    if search==0:
        print("We couldn't find the proper place.")
    while residue_count<residue_number:
        amino_acid=reference_seq[place_count]
        if amino_acid!="-":
            residue_count+=1
        place_count+=1
    return place_count

def index_to_resno(taxID,alignment_dict,index):
    place_count=0
    residue_count=0
    search=0
    for seq in alignment_dict:
        if seq.split("|")[-1]==taxID: #Finds first seq with certain taxID
            reference_seq=alignment_dict[seq]
            search=1
            break
    if search==0:
        print("We couldn't find the proper place.")
    while place_count<index:
        amino_acid=reference_seq[place_count]
        if amino_acid!="-":
            residue_count+=1
        place_count+=1
    return residue_count

    
def pairwise_compare(list1,list2,position_object,comparison_type,MSA,specificity,consensus,lower,taxID,
                     similarity):
    residue_list=[]
    empty_list=0
    for i in list2:
        print("****",i,"*****")
        if [i]==list1:
            continue
        sequence_result=spial(MSA,
                              list1,[i]
                              ,specificity,consensus,"sequence",lower,
                              similarity,
                              PDB_taxID=taxID,
                              min_blossum_score=2)[comparison_type]
        spial_result=spial_result_parse(sequence_result,"sequence")
        if type(position_object)==dict:
            for b in spial_result[1]:
                if b not in residue_list:
                    residue_list.append(b)
            for a in spial_result[2]:
                if a not in position_object:
                    position_object[a]=1
                else:
                    position_object[a]+=1
            print(position_object)
            #print("+".join(str(x) for x in residue_list))
        if type(position_object)==list:    
            if position_object==[]:
                if empty_list==1:
                    print("No residues left.")
                    return
                empty_list=1
                position_object=spial_result[2]
                residue_list=spial_result[1]
            else:
                position_object=list(set(position_list) & set(spial_result[2]))
                residue_list=list(set(residue_list) & set(spial_result[1]))
            print(position_list)
         
    
    if type(position_object)==list:
        print("+".join(str(x) for x in residue_list))
        return position_object
    if type(position_object)==dict:
        return position_object


#You should add the following lines to your code and run it:
#create two lists to compare as couplers and non_couplers
#couplers should contain all gene names that couple to a particular g protein
#non_couplers should contain gene names that does not couple to that g protein

#Here we give an example for Gz in Inoue:

couplers=["ADA2A","ADA2B","ADA2C","DRD1","DRD2","DRD3","DRD4","DRD5","5HT1B","5HT1D","5HT1E","5HT2A","5HT2B","5HT2C","5HT4R","5HT6R","5HT7R","HRH1","HRH2","ACM3","ACM5"]
non_couplers=["ADRB1","ADRB2","ADRB3","ADA1A","ADA1B","ADA1D","5HT1A","5HT1F","HRH3","HRH4","ACM1","ACM2","ACM4"]

MSA="directory of your multiple sequence alignment"


#Specific Aproach

#First, comparing within coupler group
inside_dict={}
for i in couplers:
    inside_dict=pairwise_compare([i],couplers,inside_dict,0,MSA,90,90,70,"9606",1)

#Secondly, we compare with the non-couplers
outside_dict={}
for i in couplers:
    outside_dict=pairwise_compare([i],non_couplers,outside_dict,0,MSA,90,90,70,"9606",1)

#For a residue to be specific for this approach inside_dict value should be 0 and outside_dict value should be at leas 1

#================
#Sensitive Aproach
sensitive_dict={}
sensitive_dict=pairwise_compare([couplers],non_couplers,sensitive_dict,0,MSA,90,90,70,"9606",1)
#We get every residue from this approach.

#Last parameter controls the similarity measure from BLOSUM80 matrix. You can disable it by chaingin it to 0.

