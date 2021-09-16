# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 20:32:01 2021

@author: bselcuk
"""

def GPCR_convert(res,gene):
    import pickle
    gene=gene.lower()
    try:
        with open(r"D:\Users\suuser\Desktop\CLASS A GPCRS STUDY FILE\GPCR_notations\{}_notation_dict".format(gene), 'rb') as handle:
            conv_dict = pickle.load(handle)
        
        if type(res)==str:# GPCRdb Numbering ---> Residue Numbering
            if res not in conv_dict:
                return "-"
            return conv_dict[res]
        if type(res)==int: # Residue number ---> GPCRdb Numbering
            for key in conv_dict:
                if int(conv_dict[key])==res:
                    return key
            return (str(res))
    except:
        return ("-")

def RRCS_residue_obtain(RRCS_file,residue_list,chain,output_file,gene,convert,correction):
    file=open(RRCS_file,"r")
    output=open(output_file,"w")
    count=0
    for line in file:
        if line.count(chain+":")==2:
            line_list=line.strip().split()
            interactor1=line_list[0]
            interactor2=line_list[1]
            value=line_list[2]
            res_no1=interactor1.split(":")[1].split("_")[0]
            res_no2=interactor2.split(":")[1].split("_")[0]
            if int(res_no1)+correction not in residue_list or int(res_no2)+correction not in residue_list:
                continue
            else:
                new_line=GPCR_convert(int(res_no1)+correction,gene)+"\t"+GPCR_convert(int(res_no2)+correction,gene)+"\t"+value+"\n"
                output.write(new_line)
                count+=1
    file.close()
    output.close()
    return "We obtained",count,"of interactions in total."

def RRCS_to_dict(RRCS_file):
    out_dict={}
    file=open(RRCS_file,"r")
    for line in file:
        line=line.strip()
        line_list=line.split("\t")
        key=line_list[0]+"\t"+line_list[1]
        value=line_list[2]
        out_dict[key]=value
    return out_dict

def alignment_dict(MA_file):
    result_dict={}
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

def fasta_divider(gene_list,align_dict):
    result_dict={}
    for gene in gene_list:
        count=0
        for seq in align_dict:
            if "|" not in seq:
                continue
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

    if amino_acid!="-":    
        return residue_count
    else:
        return str(residue_count)+"-"+str(residue_count+1)

def residue_convert(gene1,gene2,alignment_dict,res_no):
    alignment_d1=fasta_divider([gene1],alignment_dict)
    alignment_d2=fasta_divider([gene2],alignment_dict)
    position=place_find("9606",alignment_d1,res_no)
    index=position
    new_residue=index_to_resno("9606", alignment_d2,index)
    return new_residue

def RRCS_convert(RRCS,gene1,gene2,align_d,outfile,threshold):
    RRCS_dict=RRCS_to_dict(RRCS)
    out_dict={}
    file=open(outfile,"w")
    for key in RRCS_dict:
        key_list=key.split()
        key1=key_list[0]
        key2=key_list[1]
        if "x" in key1:
            internal_value=residue_convert(gene1,gene2,align_d,int(GPCR_convert(key_list[0],gene1)))
            key1=GPCR_convert(int(internal_value),gene2)
        else:
            key1=GPCR_convert(residue_convert(gene1,gene2,align_d,int(key_list[0])),gene2)
        if "x" in key2:
            key2=GPCR_convert(int(residue_convert(gene1,gene2,align_d,int(GPCR_convert(key_list[1],gene1)))),gene2)
        else:
            key2=GPCR_convert(residue_convert(gene1,gene2,align_d,int(key_list[1])),gene2)
        value=RRCS_dict[key]
        if float(value)>threshold or float(value)*-1>threshold:
            out_dict[key1+"\t"+key2]=RRCS_dict[key]
    for key in out_dict:
        file.write(key+"\t"+out_dict[key]+"\n")
    file.close()
    return "DONE!!!!!"
import pandas as pd
import os

def df_pairwise_t_test(p_value,first_list,second_list,df,mode):
    from scipy.stats import ttest_ind
    from statistics import mean 
    from statistics import mode as md
    if mode==1:
        for ref in first_list:
            for tar in second_list:
                for col in df:
                    reference_list=[]
                    target_list=[]
                    idx=0
                    for value in df[col]:
                        if ref in column_list[idx]:
                            reference_list.append(value)     
                        if tar in column_list[idx]:
                            target_list.append(value)
                        idx+=1
                    # if max(reference_list)==0==min(reference_list):
                    if reference_list.count(0)>=len(reference_list)/2:
                        print("OH NO!!")
                        df.drop(col,axis='columns', inplace=True)
                        continue
                    condition=False
                    comparison_stat=ttest_ind(reference_list,target_list,equal_var=False)

                    if p_value>=comparison_stat[1] and ((mean(reference_list)>0 and mean(target_list)<=0) or (mean(reference_list)<0 and mean(target_list)>=0)):
                        condition=True
                    if condition:
                        continue
                    else:
                        df.drop(col,axis='columns', inplace=True)            
        for idx1 in range(len(first_list)-1):
            for idx2 in range(1,len(first_list)):
                for col in df:
                    reference_list=[]
                    target_list=[]
                    idx=0
                    for value in df[col]:
                        if first_list[idx1] in column_list[idx]:
                            reference_list.append(value)     
                        if first_list[idx2] in column_list[idx]:
                            target_list.append(value)
                        idx+=1
                    if max(reference_list)==0==min(reference_list):
                        df.drop(col,axis='columns', inplace=True)
                        continue
                    condition=True
                    comparison_stat=ttest_ind(reference_list,target_list,equal_var=False)
                    if p_value>=comparison_stat[1]:
                        condition=False
                    if condition:
                        continue
                    else:
                        df.drop(col,axis='columns', inplace=True)
    if mode==2:
        for col in df:
            reference_list=[]
            target_list=[]
            idx=0
            for value in df[col]:
                for ref in first_list:
                    if ref in column_list[idx]:
                        reference_list.append(value)
                        break
                for tar in second_list:    
                    if tar in column_list[idx]:
                        target_list.append(value)
                        break
                idx+=1
            if max(reference_list)==0==min(reference_list):
                df.drop(col,axis='columns', inplace=True)
                continue
            condition=False
            comparison_stat=ttest_ind(reference_list,target_list,equal_var=False)
            if p_value>=comparison_stat[1] and ((mean(reference_list)>0 and mean(target_list)<=0) or (mean(reference_list)<0 and mean(target_list)>=0)):
                condition=True
            if condition:
                print(col+"\t"+str(comparison_stat[0])+"\t"+str(comparison_stat[1]))
                continue
            else:
                df.drop(col,axis='columns', inplace=True) 
    return df

def df_to_pymol(df,gene_to_convert,Gtarget_list,separator):
    pymol_list=[]
    specs_list=[]
    count=0
    for i in df:
        count+=1
        i_list=i.replace("*","").split(separator)
        if "x" in i_list[0]:
            resno=int(GPCR_convert(i_list[0],"ADRB2"))
            first=str(residue_convert("ADRB2",gene_to_convert,alignment,resno))
        else:
            first=str(residue_convert("ADRB2",gene_to_convert,alignment,int(i_list[0])))
            
        if first not in specs_list and i_list[0] in Gtarget_list:
            specs_list.append(first)
            
        if "x" in i_list[1]:
            resno=int(GPCR_convert(i_list[1],"ADRB2"))
            second=str(residue_convert("ADRB2",gene_to_convert,alignment,resno))
        else:
            second=str(residue_convert("ADRB2",gene_to_convert,alignment,int(i_list[1])))
        if second not in specs_list and i_list[1] in Gtarget_list:
            specs_list.append(second)
        if first not in pymol_list:
            pymol_list.append(first)
        if second not in pymol_list:
            pymol_list.append(second)
    # print(df)
    print()
    print("All residues that are present in specific networks:\n")
    print("+".join(pymol_list))
    print()
    print("Specific residues that we observe in these specific networks:\n")
    print("+".join(specs_list))
    print("Total number of interactions:",count)
#%%
residue_listDRD2=[156, 79, 51, 70, 189, 87, 128, 110, 81, 434, 369, 139, 378, 
                  367, 103, 415, 441, 174, 401, 63, 75, 83, 48, 170, 55, 437, 
                  142, 204, 59, 88, 120, 89, 212, 72, 163, 54, 430, 399, 209, 
                  199, 125, 375, 418, 119, 205, 49, 390, 368, 371, 374, 147, 
                  213, 52, 65, 76, 80, 100, 107, 114, 121, 131, 132, 133, 136,
                  160, 182, 198, 201, 370, 382, 386, 388, 413, 419, 423, 426]
residue_list5HT2A=[196, 119, 91, 110, 234, 127, 169, 151, 121, 388, 319, 180, 
                   328, 317, 144, 369, 395, 215, 353, 103, 115, 123, 88, 211, 
                   95, 391, 183, 249, 99, 128, 161, 129, 257, 112, 203, 94, 
                   384, 349, 254, 244, 166, 325, 372, 160, 250, 89, 340, 318, 
                   321, 324, 188, 258, 92, 105, 116, 120, 141, 148, 155, 162, 
                   172, 173, 174, 177, 200, 227, 243, 246, 320, 332, 336, 338, 
                   367, 373, 377, 380]
residue_list5HT1B=[170, 94, 66, 85, 208, 102, 143, 125, 96, 377, 310, 154, 319, 
                   308, 118, 358, 384, 187, 344, 78, 90, 98, 63, 184, 70, 380, 
                   157, 223, 74, 103, 135, 104, 231, 87, 177, 69, 373, 340, 228, 
                   218, 140, 316, 361, 134, 224, 64, 331, 309, 312, 315, 162, 232, 
                   67, 80, 91, 95, 115, 122, 129, 136, 146, 147, 148, 151, 174, 
                   199, 217, 220, 311, 323, 327, 329, 356, 362, 366, 369]
residue_listADRB2=[154, 78, 50, 69, 199, 86, 127, 109, 80, 333, 269, 138, 278, 
                   267, 102, 315, 340, 173, 301, 62, 74, 82, 47, 168, 54, 336, 
                   141, 214, 58, 87, 119, 88, 222, 71, 161, 53, 329, 299, 219, 
                   209, 124, 275, 318, 118, 215, 48, 290, 268, 271, 274, 146, 
                   223, 51, 64, 75, 79, 99, 106, 113, 120, 130, 131, 132, 135, 
                   158, 191, 208, 211, 270, 282, 286, 288, 313, 319, 323, 326]
residue_listACM2=[144, 68, 40, 59, 186, 76, 117, 99, 70, 448, 383, 128, 392, 
                  381, 92, 429, 455, 162, 416, 52, 64, 72, 37, 158, 44, 451, 131, 201, 48, 77, 109, 78, 209, 61, 151, 43, 444, 413, 206, 196, 114, 389, 432, 108, 202, 38, 404, 382, 385, 388, 136, 210, 41, 54, 65, 69, 89, 96, 103, 110, 120, 121, 122, 125, 148, 176, 195, 198, 384, 396, 400, 402, 427, 433, 437, 440]
residue_listADA2B=[133, 57, 29, 48, 172, 65, 106, 88, 59, 434, 367, 117, 376, 
                   365, 81, 415, 441, 151, 401, 41, 53, 61, 26, 147, 33, 437, 120, 187, 37, 66, 98, 67, 195, 50, 140, 32, 430, 397, 192, 182, 103, 373, 418, 97, 188, 27, 388, 366, 369, 372, 125, 196, 30, 43, 54, 58, 78, 85, 92, 99, 109, 110, 111, 114, 137, 164, 181, 184, 368, 380, 384, 386, 413, 419, 423, 426]
residue_listDRD1=[144, 69, 40, 60, 194, 77, 117, 99, 71, 338, 268, 128, 277, 266, 93, 320, 345, 163, 307, 52, 65, 73, 37, 158, 44, 341, 131, 209, 48, 78, 109, 79, 217, 62, 151, 43, 334, 298, 214, 204, 114, 274, 323, 108, 210, 38, 289, 267, 270, 273, 136, 218, 41, 54, 66, 70, 90, 96, 103, 110, 120, 121, 122, 125, 148, 186, 203, 206, 269, 281, 285, 287, 318, 324, 328, 331]
residue_listDRD3=[154, 74, 46, 65, 188, 82, 124, 106, 76, 391, 325, 135, 334, 323, 99, 372, 398, 172, 357, 58, 70, 78, 43, 168, 50, 394, 138, 203, 54, 83, 116, 84, 211, 67, 161, 49, 387, 355, 208, 198, 121, 331, 375, 115, 204, 44, 346, 324, 327, 330, 145, 212, 47, 60, 71, 75, 96, 103, 110, 117, 127, 128, 129, 132, 158, 181, 197, 200, 326, 338, 342, 344, 370, 376, 380, 383]
residue_listHRH1=[148, 72, 44, 63, 190, 80, 121, 103, 74, 476, 411, 132, 420, 409, 96, 457, 483, 165, 443, 56, 68, 76, 41, 161, 48, 479, 135, 205, 52, 81, 113, 82, 213, 65, 155, 47, 472, 441, 210, 200, 118, 417, 460, 112, 206, 42, 432, 410, 413, 416, 140, 214, 45, 58, 69, 73, 93, 100, 107, 114, 124, 125, 126, 129, 152, 180, 199, 202, 412, 424, 428, 430, 455, 461, 465, 468]
#%%
Gs_list=["1x52","1x57","2x40","2x42","2x51","2x56","2x58","3x21","3x28",
         "3x38","3x43","34x50","34x53","4x38","4x53","173","5x39","5x48",
         "5x58","5x61","5x62","6x29","6x30","6x36","6x37","6x40","6x52","7x41","7x45","8x51","8x58"]
Gi1_list=["2x40","2x45","34x53","4x46","6x61","301","8x47","3x46","12x48","2x53","6x31","5x54","1x47"]
Gq_list=["1x46","1x47","1x52","2x40","2x42","2x45","2x49","2x51","3x28","4x46",
         "4x53","5x48","5x54","5x58","6x37","6x61","301","8x47","8x51","1x49","2x57","3x46","3x43","6x30","6x40"]
Gio_list=["2x40","2x45","34x53","5x48","5x54","3x46","7x45","5x61","5x58","6x40","8x51","3x43","6x37","4x46","8x47"]
Gi_intersection_list=["2x40","2x45","4x46","8x47","3x46","5x54"]

common_activation_mechanism=["6x48","6x44","5x51","3x39","7x45","2x50","2x46","5x55","6x41","3x43","7x49","2x46",
                             "6x40","7x50","1x49","7x52","7x53","7x54","7x55","1x53","2x43","8x50","8x51","3x46","6x37","5x58","5x57",
                             "3x49","3x50","3x51","3x53","3x54","6x33","5x61","5x62"]

layer0=["5x47","3x32","2x58","6x61","2x57","3x28","3x25","7x39","3x21","5x48","45x50","4x60","23x50","5x39","301","173","2x56","6x52","3x38","4x53"]
layer1=["4x50","1x46","2x53","6x48","7x50","1x47","2x51","2x50","7x46","1x49","7x45","2x49","1x50","3x39","6x44","7x41"]
layer2=["2x46","3x43","4x46","5x54","6x40","2x45"]
layer3=["1x53","3x46","1x57","8x54","7x53","2x40","2x42","1x52","8x51","12x50","8x58","6x37","5x58"]
layer4=["6x33","5x62","3x49","8x47","6x32","6x29","6x31","34x53","3x54","5x61","6x30","3x51","3x50","6x36","34x50","4x38","12x48"]
#%%
def get_percentages(target_list):
    target_list=list(set(target_list))
    result=[0,0,0,0,0]
    for i in target_list:
        if i in layer0:
            result[0]+=1
            continue
        if i in layer1:
            result[1]+=1
            continue
        if i in layer2:
            result[2]+=1
            continue
        if i in layer3:
            result[3]+=1
            continue
        if i in layer4:
            result[4]+=1
            continue
        else:
            print(i)
            
    for i in range(5):
        result[i]/=len(target_list)
        result[i]=result[i]*100
    return result

gs_percentages=get_percentages(Gs_list)
gi1_precentages=get_percentages(Gi1_list)
go_percentages=get_percentages(Gio_list)
gq_percentages=get_percentages(Gq_list)

data=[gi1_precentages,go_percentages,gq_percentages,gs_percentages]
#%%
import plotly.graph_objects as go

top_labels = ['Strongly<br>agree', 'Agree', 'Neutral', 'Disagree',
              'Strongly<br>disagree']

colors = ['rgba(38, 24, 74, 0.8)', 'rgba(71, 58, 131, 0.8)',
          'rgba(122, 120, 168, 0.8)', 'rgba(164, 163, 204, 0.85)',
          'rgba(190, 192, 213, 1)']

x_data = data

y_data = ['The course was effectively<br>organized',
          'The course developed my<br>abilities and skills ' +
          'for<br>the subject', 'The course developed ' +
          'my<br>ability to think critically about<br>the subject',
          'I would recommend this<br>course to a friend']

fig = go.Figure()

for i in range(0, len(x_data[0])):
    for xd, yd in zip(x_data, y_data):
        fig.add_trace(go.Bar(
            x=[xd[i]], y=[yd],
            orientation='h',
            marker=dict(
                color=colors[i],
                line=dict(color='rgb(248, 248, 249)', width=1)
            )
        ))

fig.update_layout(
    xaxis=dict(
        showgrid=False,
        showline=False,
        showticklabels=False,
        zeroline=False,
        domain=[0.15, 1]
    ),
    yaxis=dict(
        showgrid=False,
        showline=False,
        showticklabels=False,
        zeroline=False,
    ),
    barmode='stack',
    paper_bgcolor='rgb(248, 248, 255)',
    plot_bgcolor='rgb(248, 248, 255)',
    margin=dict(l=120, r=10, t=140, b=80),
    showlegend=False,
)

annotations = []

for yd, xd in zip(y_data, x_data):
    # labeling the y-axis
    annotations.append(dict(xref='paper', yref='y',
                            x=0.14, y=yd,
                            xanchor='right',
                            text=str(yd),
                            font=dict(family='Arial', size=14,
                                      color='rgb(67, 67, 67)'),
                            showarrow=False, align='right'))
    # labeling the first percentage of each bar (x_axis)
    annotations.append(dict(xref='x', yref='y',
                            x=xd[0] / 2, y=yd,
                            text=str(xd[0]) + '%',
                            font=dict(family='Arial', size=14,
                                      color='rgb(248, 248, 255)'),
                            showarrow=False))
    # labeling the first Likert scale (on the top)
    if yd == y_data[-1]:
        annotations.append(dict(xref='x', yref='paper',
                                x=xd[0] / 2, y=1.1,
                                text=top_labels[0],
                                font=dict(family='Arial', size=14,
                                          color='rgb(67, 67, 67)'),
                                showarrow=False))
    space = xd[0]
    for i in range(1, len(xd)):
            # labeling the rest of percentages for each bar (x_axis)
            annotations.append(dict(xref='x', yref='y',
                                    x=space + (xd[i]/2), y=yd,
                                    text=str(xd[i]) + '%',
                                    font=dict(family='Arial', size=14,
                                              color='rgb(248, 248, 255)'),
                                    showarrow=False))
            # labeling the Likert scale
            if yd == y_data[-1]:
                annotations.append(dict(xref='x', yref='paper',
                                        x=space + (xd[i]/2), y=1.1,
                                        text=top_labels[i],
                                        font=dict(family='Arial', size=14,
                                                  color='rgb(67, 67, 67)'),
                                        showarrow=False))
            space += xd[i]

fig.update_layout(annotations=annotations)

fig.show()
#%%
fig.write_image(r"D:\Users\suuser\Desktop\Figure2_bar.svg")

#%%
count=0
MSA=r"D:\Users\suuser\Desktop\CLASS A GPCRS STUDY FILE\aminergic_orthologs_04-07-2020_MSA_einsi.fasta"
alignment=alignment_dict(MSA)[0]
#%%
pdb_list=["7DHR","7DHI","7BZ2","6OBA","6NI3","6PS6","6PS5","6PS3","6PS1","6PS0"
          ,"6PS2","6PS4","6PRZ","6N48","6E67","6MXT","5X7D","5D6L","5JQH"
          ,"5D5A","5D5B","4QKX","4LDE","4LDO","4LDL","4GBR","3SN6","3P0G","3PDS","3NYA"
          ,"3NY8","3NY9","3KJ6","3D4S","2R4R","2R4S","2RH1"]
active_pdb_list=["3sn6","6e67","6ni3","7bz2","7dhi","7dhr"]
for active_pdb in active_pdb_list:
    if active_pdb=="6e67":
        chain="B"
    else:
        chain="R"
    print(RRCS_residue_obtain(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\ADRB2_actives_RRCS\{}.pdb.cscore".format(active_pdb.lower()),
                          residue_listADRB2,chain,
                          r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\ADRB2_active_RRCS_{}.tsv".format(active_pdb),"ADRB2","No",0))
    for pdb in  pdb_list:
        if pdb!="5JQH":
            continue
        try:
            count+=1
            print(RRCS_residue_obtain(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\ADRB2_inactives_RRCS\{}.pdb.cscore".format(pdb.lower()),
                          residue_listADRB2,"A",
                          r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\ADRB2_inactive_RRCS_{}.tsv".format(pdb),"ADRB2","No",-1000))
            inactive_dict=RRCS_to_dict((r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\ADRB2_inactive_RRCS_{}.tsv".format(pdb.lower())))
            active_dict=RRCS_to_dict((r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\ADRB2_active_RRCS_{}.tsv".format(active_pdb)))
            file=open(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\ADRB2_delta_RRCS-{}-{}.tsv".format(pdb.lower(),active_pdb),"w")
            for i in active_dict:
                if i in inactive_dict:
                    #print (i,float(active_dict[i])-float(inactive_dict[i]))
                    file.write("{}\t{}\tRepacking\t{}\t{}\n".format(i,float(active_dict[i])-float(inactive_dict[i]),float(active_dict[i]),float(inactive_dict[i])))
                else:
                    #print(i,float(active_dict[i]))
                    file.write("{}\t{}\tSwitching\t{}\t{}\n".format(i,float(active_dict[i]),float(active_dict[i]),0))
            for i in inactive_dict:
                if i not in active_dict:
                    #print (i,-1*float(inactive_dict[i]))
                    file.write("{}\t{}\tSwitching\t{}\t{}\n".format(i,-1*float(inactive_dict[i]),0,float(inactive_dict[i])))
            file.close()
                    
                    
            RRCS_file=r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\ADRB2_delta_RRCS-{}-{}.tsv".format(pdb.lower(),active_pdb)
            outfile=r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\ADRB2_delta_RRCS-{}-{}_converted_filtered.tsv".format(pdb.lower(),active_pdb)
            print(RRCS_convert(RRCS_file,"ADRB2","ADRB2",alignment,outfile,0.2))
        except FileNotFoundError:
            continue
    print(count)
#%%
pdb_list=["6cm4","6luq_corrected","7dfp"]
active_pdb_list=["6vsm","7jvr"]
for active_pdb in active_pdb_list:
    chain="R"
    print(RRCS_residue_obtain(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD2 Files\{}.pdb.cscore".format(active_pdb.lower()),
                          residue_listDRD2,chain,
                          r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD2 Files\DRD2_active_RRCS_{}.tsv".format(active_pdb),"DRD2","No",0))
    for pdb in  pdb_list:
        try:
            print(RRCS_residue_obtain(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD2 Files\{}.pdb.cscore".format(pdb.lower()),
                          residue_listDRD2,"A",r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD2 Files\DRD2_inactive_RRCS_{}.tsv".format(pdb),"DRD2","No",0))

            inactive_dict=RRCS_to_dict((r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD2 Files\DRD2_inactive_RRCS_{}.tsv".format(pdb.lower())))
            active_dict=RRCS_to_dict((r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD2 Files\DRD2_active_RRCS_{}.tsv".format(active_pdb)))
            file=open(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD2 Files\DRD2_delta_RRCS-{}-{}.tsv".format(pdb.lower(),active_pdb),"w")
            for i in active_dict:
                if i in inactive_dict:
                    #print (i,float(active_dict[i])-float(inactive_dict[i]))
                    file.write("{}\t{}\tRepacking\t{}\t{}\n".format(i,float(active_dict[i])-float(inactive_dict[i]),float(active_dict[i]),float(inactive_dict[i])))
                else:
                    #print(i,float(active_dict[i]))
                    file.write("{}\t{}\tSwitching\t{}\t{}\n".format(i,float(active_dict[i]),float(active_dict[i]),0))
            for i in inactive_dict:
                if i not in active_dict:
                    #print (i,-1*float(inactive_dict[i]))
                    file.write("{}\t{}\tSwitching\t{}\t{}\n".format(i,-1*float(inactive_dict[i]),0,float(inactive_dict[i])))
            file.close()
                    
            RRCS_file=r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD2 Files\DRD2_delta_RRCS-{}-{}.tsv".format(pdb.lower(),active_pdb)
            outfile=r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD2 Files\DRD2_delta_RRCS-{}-{}_converted_filtered.tsv".format(pdb.lower(),active_pdb)
            print(RRCS_convert(RRCS_file,"DRD2","ADRB2",alignment,outfile,0.2))
        except KeyError:
            continue
#%%

active_pdb_list=["7ckw","7ckx","7ckz","7cky","7crh","7jv5","7jvp","7jvq","7ljc","7ljd"]
pdb_list=["drd1inactivemodel"]
for active_pdb in active_pdb_list:
    chain="R"
    print(RRCS_residue_obtain(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD1 Files\{}.pdb.cscore".format(active_pdb.lower()),
                          residue_listDRD1,chain,
                          r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD1 Files\DRD1_active_RRCS_{}.tsv".format(active_pdb),"DRD1","No",0))
    for pdb in  pdb_list:
        try:
            print(RRCS_residue_obtain(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD1 Files\{}.pdb.cscore".format(pdb.lower()),
                          residue_listDRD1," ",r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD1 Files\DRD1_inactive_RRCS_{}.tsv".format(pdb),"DRD1","No",0))

            inactive_dict=RRCS_to_dict((r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD1 Files\DRD1_inactive_RRCS_{}.tsv".format(pdb.lower())))
            active_dict=RRCS_to_dict((r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD1 Files\DRD1_active_RRCS_{}.tsv".format(active_pdb)))
            file=open(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD1 Files\DRD1_delta_RRCS-{}-{}.tsv".format(pdb.lower(),active_pdb),"w")
            for i in active_dict:
                if i in inactive_dict:
                    #print (i,float(active_dict[i])-float(inactive_dict[i]))
                    file.write("{}\t{}\tRepacking\t{}\t{}\n".format(i,float(active_dict[i])-float(inactive_dict[i]),float(active_dict[i]),float(inactive_dict[i])))
                else:
                    #print(i,float(active_dict[i]))
                    file.write("{}\t{}\tSwitching\t{}\t{}\n".format(i,float(active_dict[i]),float(active_dict[i]),0))
            for i in inactive_dict:
                if i not in active_dict:
                    #print (i,-1*float(inactive_dict[i]))
                    file.write("{}\t{}\tSwitching\t{}\t{}\n".format(i,-1*float(inactive_dict[i]),0,float(inactive_dict[i])))
            file.close()
                    
            RRCS_file=r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD1 Files\DRD1_delta_RRCS-{}-{}.tsv".format(pdb.lower(),active_pdb)
            outfile=r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD1 Files\DRD1_delta_RRCS-{}-{}_converted_filtered.tsv".format(pdb.lower(),active_pdb)
            print(RRCS_convert(RRCS_file,"DRD1","ADRB2",alignment,outfile,0.2))
        except KeyError:
            continue
#%%
pdb_list=["6a93","6a94","6wh4","6wgt"]
active_pdb_list=["6wha"]
for active_pdb in active_pdb_list:
    chain="A"
    print(RRCS_residue_obtain(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\5HT2A Files\{}.pdb.cscore".format(active_pdb.lower()),
                          residue_list5HT2A,chain,
                          r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\5HT2A Files\5HT2A_active_RRCS_{}.tsv".format(active_pdb),"5HT2A","No"))
    for pdb in  pdb_list:
        if pdb=="6wh4":
            chain="C"
        elif pdb=="6wgt":
            chain="B"
        else:
            chain="A"
        print(RRCS_residue_obtain(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\5HT2A Files\{}.pdb.cscore".format(pdb.lower()),
                      residue_list5HT2A,chain,r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\5HT2A Files\5HT2A_inactive_RRCS_{}.tsv".format(pdb),"5HT2A","No"))

        inactive_dict=RRCS_to_dict((r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\5HT2A Files\5HT2A_inactive_RRCS_{}.tsv".format(pdb.lower())))
        active_dict=RRCS_to_dict((r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\5HT2A Files\5HT2A_active_RRCS_{}.tsv".format(active_pdb)))
        file=open(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\5HT2A Files\5HT2A_delta_RRCS-{}-{}.tsv".format(pdb.lower(),active_pdb),"w")
        for i in active_dict:
            if i in inactive_dict:
                #print (i,float(active_dict[i])-float(inactive_dict[i]))
                file.write("{}\t{}\tRepacking\t{}\t{}\n".format(i,float(active_dict[i])-float(inactive_dict[i]),float(active_dict[i]),float(inactive_dict[i])))
            else:
                #print(i,float(active_dict[i]))
                file.write("{}\t{}\tSwitching\t{}\t{}\n".format(i,float(active_dict[i]),float(active_dict[i]),0))
        for i in inactive_dict:
            if i not in active_dict:
                #print (i,-1*float(inactive_dict[i]))
                file.write("{}\t{}\tSwitching\t{}\t{}\n".format(i,-1*float(inactive_dict[i]),0,float(inactive_dict[i])))
        file.close()
                
        RRCS_file=r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\5HT2A Files\5HT2A_delta_RRCS-{}-{}.tsv".format(pdb.lower(),active_pdb)
        outfile=r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\5HT2A Files\5HT2A_delta_RRCS-{}-{}_converted_filtered.tsv".format(pdb.lower(),active_pdb)
        print(RRCS_convert(RRCS_file,"5HT2A","ADRB2",alignment,outfile,0.2))

        
#%%
pdb_list=["3uon","5yc8","5zk3","5zk8","5zkb","5zkc"]
active_pdb_list=["6oik"]
for active_pdb in active_pdb_list:
    chain="R"
    print(RRCS_residue_obtain(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\ACM2 Files\{}.pdb.cscore".format(active_pdb.lower()),
                          residue_listACM2,chain,
                          r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\ACM2 Files\ACM2_active_RRCS_{}.tsv".format(active_pdb),"ACM2","No"))
    for pdb in  pdb_list:

        chain="A"
        print(RRCS_residue_obtain(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\ACM2 Files\{}.pdb.cscore".format(pdb.lower()),
                      residue_listACM2,chain,r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\ACM2 Files\ACM2_inactive_RRCS_{}.tsv".format(pdb),"ACM2","No"))

        inactive_dict=RRCS_to_dict((r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\ACM2 Files\ACM2_inactive_RRCS_{}.tsv".format(pdb.lower())))
        active_dict=RRCS_to_dict((r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\ACM2 Files\ACM2_active_RRCS_{}.tsv".format(active_pdb)))
        file=open(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\ACM2 Files\ACM2_delta_RRCS-{}-{}.tsv".format(pdb.lower(),active_pdb),"w")
        for i in active_dict:
            if i in inactive_dict:
                #print (i,float(active_dict[i])-float(inactive_dict[i]))
                file.write("{}\t{}\tRepacking\t{}\t{}\n".format(i,float(active_dict[i])-float(inactive_dict[i]),float(active_dict[i]),float(inactive_dict[i])))
            else:
                #print(i,float(active_dict[i]))
                file.write("{}\t{}\tSwitching\t{}\t{}\n".format(i,float(active_dict[i]),float(active_dict[i]),0))
        for i in inactive_dict:
            if i not in active_dict:
                #print (i,-1*float(inactive_dict[i]))
                file.write("{}\t{}\tSwitching\t{}\t{}\n".format(i,-1*float(inactive_dict[i]),0,float(inactive_dict[i])))
        file.close()
                
        RRCS_file=r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\ACM2 Files\ACM2_delta_RRCS-{}-{}.tsv".format(pdb.lower(),active_pdb)
        outfile=r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\ACM2 Files\ACM2_delta_RRCS-{}-{}_converted_filtered.tsv".format(pdb.lower(),active_pdb)
        print(RRCS_convert(RRCS_file,"ACM2","ADRB2",alignment,outfile,0.2))
#%%
pdb_list=["4iaq","4iar","5v54","7c61_corrected"]
active_pdb_list=["6g79"]
for active_pdb in active_pdb_list:
    chain="S"
    print(RRCS_residue_obtain(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\5HT1B Files\{}.pdb.cscore".format(active_pdb.lower()),
                          residue_list5HT1B,chain,
                          r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\5HT1B Files\5HT1B_active_RRCS_{}.tsv".format(active_pdb),"5HT1B","No",0))
    for pdb in  pdb_list:

        chain="A"
        print(RRCS_residue_obtain(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\5HT1B Files\{}.pdb.cscore".format(pdb.lower()),
                      residue_list5HT1B,chain,r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\5HT1B Files\5HT1B_inactive_RRCS_{}.tsv".format(pdb),"5HT1B","No",0))

        inactive_dict=RRCS_to_dict((r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\5HT1B Files\5HT1B_inactive_RRCS_{}.tsv".format(pdb.lower())))
        active_dict=RRCS_to_dict((r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\5HT1B Files\5HT1B_active_RRCS_{}.tsv".format(active_pdb)))
        file=open(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\5HT1B Files\5HT1B_delta_RRCS-{}-{}.tsv".format(pdb.lower(),active_pdb),"w")
        for i in active_dict:
            if i in inactive_dict:
                #print (i,float(active_dict[i])-float(inactive_dict[i]))
                file.write("{}\t{}\tRepacking\t{}\t{}\n".format(i,float(active_dict[i])-float(inactive_dict[i]),float(active_dict[i]),float(inactive_dict[i])))
            else:
                #print(i,float(active_dict[i]))
                file.write("{}\t{}\tSwitching\t{}\t{}\n".format(i,float(active_dict[i]),float(active_dict[i]),0))
        for i in inactive_dict:
            if i not in active_dict:
                #print (i,-1*float(inactive_dict[i]))
                file.write("{}\t{}\tSwitching\t{}\t{}\n".format(i,-1*float(inactive_dict[i]),0,float(inactive_dict[i])))
        file.close()
                
        RRCS_file=r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\5HT1B Files\5HT1B_delta_RRCS-{}-{}.tsv".format(pdb.lower(),active_pdb)
        outfile=r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\5HT1B Files\5HT1B_delta_RRCS-{}-{}_converted_filtered.tsv".format(pdb.lower(),active_pdb)
        print(RRCS_convert(RRCS_file,"5HT1B","ADRB2",alignment,outfile,0.2))
#%%
pdb_list=["3pbl"]
active_pdb_list=["7cmu","7cmv"]
for active_pdb in active_pdb_list:
    chain="R"
    print(RRCS_residue_obtain(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD3 Files\{}.pdb.cscore".format(active_pdb.lower()),
                          residue_listDRD3,chain,
                          r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD3 Files\DRD3_active_RRCS_{}.tsv".format(active_pdb),"DRD3","No",0))
    for pdb in  pdb_list:

        chain="A"
        print(RRCS_residue_obtain(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD3 Files\{}.pdb.cscore".format(pdb.lower()),
                      residue_listDRD3,chain,r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD3 Files\DRD3_inactive_RRCS_{}.tsv".format(pdb),"DRD3","No",0))

        inactive_dict=RRCS_to_dict((r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD3 Files\DRD3_inactive_RRCS_{}.tsv".format(pdb.lower())))
        active_dict=RRCS_to_dict((r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD3 Files\DRD3_active_RRCS_{}.tsv".format(active_pdb)))
        file=open(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD3 Files\DRD3_delta_RRCS-{}-{}.tsv".format(pdb.lower(),active_pdb),"w")
        for i in active_dict:
            if i in inactive_dict:
                #print (i,float(active_dict[i])-float(inactive_dict[i]))
                file.write("{}\t{}\tRepacking\t{}\t{}\n".format(i,float(active_dict[i])-float(inactive_dict[i]),float(active_dict[i]),float(inactive_dict[i])))
            else:
                #print(i,float(active_dict[i]))
                file.write("{}\t{}\tSwitching\t{}\t{}\n".format(i,float(active_dict[i]),float(active_dict[i]),0))
        for i in inactive_dict:
            if i not in active_dict:
                #print (i,-1*float(inactive_dict[i]))
                file.write("{}\t{}\tSwitching\t{}\t{}\n".format(i,-1*float(inactive_dict[i]),0,float(inactive_dict[i])))
        file.close()
                
        RRCS_file=r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD3 Files\DRD3_delta_RRCS-{}-{}.tsv".format(pdb.lower(),active_pdb)
        outfile=r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD3 Files\DRD3_delta_RRCS-{}-{}_converted_filtered.tsv".format(pdb.lower(),active_pdb)
        print(RRCS_convert(RRCS_file,"DRD3","ADRB2",alignment,outfile,0.2))
        
#%%
pdb_list=["3rze"]
active_pdb_list=["7dfl"]
for active_pdb in active_pdb_list:
    chain="R"
    print(RRCS_residue_obtain(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\HRH1 Files\{}.pdb.cscore".format(active_pdb.lower()),
                          residue_listHRH1,chain,
                          r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\HRH1 Files\HRH1_active_RRCS_{}.tsv".format(active_pdb),"HRH1","No",0))
    for pdb in  pdb_list:

        chain="A"
        print(RRCS_residue_obtain(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\HRH1 Files\{}.pdb.cscore".format(pdb.lower()),
                      residue_listHRH1,chain,r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\HRH1 Files\HRH1_inactive_RRCS_{}.tsv".format(pdb),"HRH1","No",0))

        inactive_dict=RRCS_to_dict((r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\HRH1 Files\HRH1_inactive_RRCS_{}.tsv".format(pdb.lower())))
        active_dict=RRCS_to_dict((r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\HRH1 Files\HRH1_active_RRCS_{}.tsv".format(active_pdb)))
        file=open(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\HRH1 Files\HRH1_delta_RRCS-{}-{}.tsv".format(pdb.lower(),active_pdb),"w")
        for i in active_dict:
            if i in inactive_dict:
                #print (i,float(active_dict[i])-float(inactive_dict[i]))
                file.write("{}\t{}\tRepacking\t{}\t{}\n".format(i,float(active_dict[i])-float(inactive_dict[i]),float(active_dict[i]),float(inactive_dict[i])))
            else:
                #print(i,float(active_dict[i]))
                file.write("{}\t{}\tSwitching\t{}\t{}\n".format(i,float(active_dict[i]),float(active_dict[i]),0))
        for i in inactive_dict:
            if i not in active_dict:
                #print (i,-1*float(inactive_dict[i]))
                file.write("{}\t{}\tSwitching\t{}\t{}\n".format(i,-1*float(inactive_dict[i]),0,float(inactive_dict[i])))
        file.close()
                
        RRCS_file=r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\HRH1 Files\HRH1_delta_RRCS-{}-{}.tsv".format(pdb.lower(),active_pdb)
        outfile=r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\HRH1 Files\HRH1_delta_RRCS-{}-{}_converted_filtered.tsv".format(pdb.lower(),active_pdb)
        print(RRCS_convert(RRCS_file,"HRH1","ADRB2",alignment,outfile,0.2))
#%%
import os
for filename in os.listdir(os.getcwd()):
    print(filename)

#%%

directory=r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\converted_files\\"
dict_list=[]
column_list=[]
filename_dict={}
for filename in os.listdir(os.getcwd()):
    if "converted" not in filename:
        continue
    else:
        pair_list=["4iar-6g79","4iaq-6g79","5v54-6g79","6cm4-6vsm","6cm4-7jvr"
                   ,"2rh1-3sn6","2rh1-7dhi","6ps2-3sn6","6ps2-7dhi","5d5a-3sn6","5d5a-7dhi","6ps3-3sn6","6ps3-7dhi"
                   ,"3uon-6oik","5zkc-6oik","5zkb-6oik"
                   ,"5yc8-6oik","5zk8-6oik","5zk3-6oik","6wgt-6wha","6wh4-6wha"
                   ,"DRD2","5HT2A","5HT1B","ACM2","DRD1","DRD3","HRH1"]
        # pair_list=["2rh1-3sn6","2rh1-7dhi","6ps2-3sn6","6ps2-7dhi","5d5a-3sn6","5d5a-7dhi","6ps3-3sn6","6ps3-7dhi","DRD1"]
        # pair_list=["ADA2B","ACM2","5HT1B","DRD2"]
        for pair in pair_list:
            yes=0
            if pair in filename:      
                yes=1
                break
        if yes==0:
            continue
        print(filename)
        filename=str(filename)
        dict_list.append(RRCS_to_dict(directory+filename))
        filename_list=filename.split("_")
        receptor_name=filename_list[0]
        active_pdb=filename_list[-3].split("-")[2]
        passive_pdb=filename_list[-3].split("-")[1]
        label=receptor_name#+"_"+active_pdb
        if label not in filename_dict:
            filename_dict[label]=1
        else:
            filename_dict[label]+=1
        column_list.append(label)#+"_"+passive_pdb)
#%%
dataframe={}
#total_keys=list(set(list(dict1.keys())+list(dict2.keys())+list(dict3.keys())+list(dict4.keys())+list(dict5.keys())+list(dict6.keys())+list(dict7.keys())+list(dict8.keys())+list(dict9.keys())))
total_keys=[]
for d in dict_list:
    total_keys+=list(d.keys())
total_keys=list(set(total_keys))
for key in total_keys:
    dataframe[key]=[]
    for i in dict_list:
        if key in i:
            dataframe[key].append(float(i[key]))
        else:
            dataframe[key].append(0)

data=list(dataframe.values())
keys=list(dataframe.keys())
new_keys=[]
# Gtarget_list=Gs_list
# Gtarget_list=Gi_intersection_list
# Gtarget_list=Gq_list
# Gtarget_list=Gi1_list
Gtarget_list=Gio_list
# Gtarget_list+=Gi1_list
for i in keys:
    residue_list=i.split("\t")
    res1=residue_list[0]
    res2=residue_list[1]
    if res1 in Gtarget_list:
        res1=res1+"*"
    if res2 in Gtarget_list:
        res2=res2+"*"
    new_column=res1+"-"+res2
    new_keys.append(new_column)
    
df = pd.DataFrame(data,columns =column_list, index=new_keys)
data=""
keys=""
df= df.T
print(df)

# %%
for col in df:
    print(col)
    if "*" not in col:
        df.drop(col,axis='columns', inplace=True)
print(df)
# #%%
# for col in df:
#     col_list=df[col]
#     zero_count=0
#     for value in col_list:
#         if value==0:
#             zero_count+=1
#     print(col+"\t"+str(zero_count))
#     # print(count_zero)
#%%

# first_list=["DRD2","DRD3"]
# second_list=["ADRB2","DRD1","5HT2A","5HT1B","ACM2","HRH1"]
# first_list=["ADRB2","DRD1"]
# second_list=["DRD2","ACM2","5HT2A","DRD3","HRH1","5HT1B"]
# second_list=["DRD2","ACM2","5HT2A","DRD3","HRH1"]
# second_list=["DRD2","ACM2","5HT2A","DRD3"]
first_list=["ACM2","5HT1B"]
second_list=["ADRB2","DRD1","5HT2A","DRD2","DRD3","HRH1"]
# first_list=["DRD2","DRD3","ACM2","5HT1B"]
# second_list=["ADRB2","DRD1","5HT2A","HRH1"]
# first_list=["5HT2A","HRH1"]
# second_list=["DRD2","DRD1","ACM2","5HT1B","ADRB2"]
# first_list=["DRD2","5HT1B","ACM2","DRD3","5HT2A","HRH1"]
# second_list=["ADRB2","DRD1"]
df=df_pairwise_t_test(10**-2,first_list,second_list,df,2)
df_to_pymol(df,"ADRB2",Gtarget_list,"-")               
 
#%%
for i in df:
    print (i)
    if "6x61" in i:
        df.drop(i,axis='columns', inplace=True)
#%%

import seaborn as sns; sns.set(color_codes=True,font_scale = 2.5)

g = sns.clustermap(df,cmap="vlag",figsize=(60, 30),method="complete",vmax=4,vmin=-4,dendrogram_ratio=(.02, .08),cbar_pos=(0, 1, .03, .1))
for a in g.ax_row_dendrogram.collections:
    a.set_linewidth(5)

for a in g.ax_col_dendrogram.collections:
    a.set_linewidth(5)
    
g



g.savefig(r"D:\Users\suuser\Desktop\interactions_heatmap.svg")




#%%
"""
pdb_list=["6OBA","6PS6","6PS5","6PS3","6PS1","6PS0"
          ,"6PS2","6PS4","6PRZ","5X7D","5D6L","5JQH"
          ,"5D5A","5D5B","4GBR","3PDS","3NYA"
          ,"3NY8","3NY9","3KJ6","3D4S","2R4R","2R4S","2RH1"]
active_pdb_list=["3sn6","6e67","6ni3","7bz2","7dhi","7dhr"]
done=[]
for active_pdb in active_pdb_list:
    if active_pdb=="6e67":
        chain="B"
    else:
        chain="R"
    print(RRCS_residue_obtain(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\ADRB2_actives_RRCS\{}.pdb.cscore".format(active_pdb.lower()),
                          all_residues,chain,
                          r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\ADRB2_active_RRCS_{}.tsv".format(active_pdb),"ADRB2","No",0))
    for pdb in  pdb_list:
        if pdb=="5JQH":
            continue
        try:
            count+=1
            if pdb not in done:
                
                print(RRCS_residue_obtain(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\ADRB2_inactives_RRCS\{}.pdb.cscore".format(pdb.lower()),
                              all_residues,"A",
                              r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\ADRB2_inactive_RRCS_{}.tsv".format(pdb),"ADRB2","No",0))
                done.append(pdb)
            inactive_dict=RRCS_to_dict((r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\ADRB2_inactive_RRCS_{}.tsv".format(pdb.lower())))
            active_dict=RRCS_to_dict((r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\ADRB2_active_RRCS_{}.tsv".format(active_pdb)))
            file=open(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\ADRB2_delta_RRCS-{}-{}.tsv".format(pdb.lower(),active_pdb),"w")
            for i in active_dict:
                if i in inactive_dict:
                    #print (i,float(active_dict[i])-float(inactive_dict[i]))
                    file.write("{}\t{}\tRepacking\t{}\t{}\n".format(i,float(active_dict[i])-float(inactive_dict[i]),float(active_dict[i]),float(inactive_dict[i])))
                else:
                    #print(i,float(active_dict[i]))
                    file.write("{}\t{}\tSwitching\t{}\t{}\n".format(i,float(active_dict[i]),float(active_dict[i]),0))
            for i in inactive_dict:
                if i not in active_dict:
                    #print (i,-1*float(inactive_dict[i]))
                    file.write("{}\t{}\tSwitching\t{}\t{}\n".format(i,-1*float(inactive_dict[i]),0,float(inactive_dict[i])))
            file.close()
                    
                    
            RRCS_file=r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\ADRB2_delta_RRCS-{}-{}.tsv".format(pdb.lower(),active_pdb)
            outfile=r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\ADRB2_delta_RRCS-{}-{}_converted_filtered.tsv".format(pdb.lower(),active_pdb)
            print(RRCS_convert(RRCS_file,"ADRB2","ADRB2",alignment,outfile,0.2))
        except FileNotFoundError:
            continue
    print(count)

#%%
all_residues=list(range(1,444))
pdb_list=["6cm4","6luq","7dfp"]
active_pdb_list=["6vsm","7jvr"]
#pdb_list=["6luq_corrected"]
done=[]
for active_pdb in active_pdb_list:
    chain="R"
    print(RRCS_residue_obtain(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD2 Files\{}.pdb.cscore".format(active_pdb.lower()),
                          all_residues,chain,
                          r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\DRD2_active_RRCS_{}.tsv".format(active_pdb),"DRD2","No",0))
    for pdb in  pdb_list:
        try:
            if pdb not in done:
                
                print(RRCS_residue_obtain(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD2 Files\{}.pdb.cscore".format(pdb.lower()),
                          all_residues,"A",r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\\DRD2_inactive_RRCS_{}.tsv".format(pdb),"DRD2","No",0))
                done.append(pdb)
            inactive_dict=RRCS_to_dict((r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\DRD2_inactive_RRCS_{}.tsv".format(pdb.lower())))
            active_dict=RRCS_to_dict((r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\DRD2_active_RRCS_{}.tsv".format(active_pdb)))
            file=open(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\DRD2_delta_RRCS-{}-{}.tsv".format(pdb.lower(),active_pdb),"w")
            for i in active_dict:
                if i in inactive_dict:
                    #print (i,float(active_dict[i])-float(inactive_dict[i]))
                    file.write("{}\t{}\tRepacking\t{}\t{}\n".format(i,float(active_dict[i])-float(inactive_dict[i]),float(active_dict[i]),float(inactive_dict[i])))
                else:
                    #print(i,float(active_dict[i]))
                    file.write("{}\t{}\tSwitching\t{}\t{}\n".format(i,float(active_dict[i]),float(active_dict[i]),0))
            for i in inactive_dict:
                if i not in active_dict:
                    #print (i,-1*float(inactive_dict[i]))
                    file.write("{}\t{}\tSwitching\t{}\t{}\n".format(i,-1*float(inactive_dict[i]),0,float(inactive_dict[i])))
            file.close()
                    
            RRCS_file=r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\DRD2_delta_RRCS-{}-{}.tsv".format(pdb.lower(),active_pdb)
            outfile=r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\DRD2_delta_RRCS-{}-{}_converted_filtered.tsv".format(pdb.lower(),active_pdb)
            print(RRCS_convert(RRCS_file,"DRD2","ADRB2",alignment,outfile,0.2))
        except KeyError:
            continue
#%%
active_pdb_list=["7ckw","7ckx","7ckz","7cky","7crh","7jv5","7jvp","7jvq","7ljc","7ljd"]
pdb_list=["drd1inactivemodel"]
all_residues=list(range(1,447))
for active_pdb in active_pdb_list:
    chain="R"
    print(RRCS_residue_obtain(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD1 Files\{}.pdb.cscore".format(active_pdb.lower()),
                          all_residues,chain,
                          r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\DRD1_active_RRCS_{}.tsv".format(active_pdb),"DRD1","No",0))
    for pdb in  pdb_list:
        try:
            print(RRCS_residue_obtain(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\DRD1 Files\{}.pdb.cscore".format(pdb.lower()),
                          all_residues," ",r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\DRD1_inactive_RRCS_{}.tsv".format(pdb),"DRD1","No",0))

            inactive_dict=RRCS_to_dict((r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\DRD1_inactive_RRCS_{}.tsv".format(pdb.lower())))
            active_dict=RRCS_to_dict((r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\DRD1_active_RRCS_{}.tsv".format(active_pdb)))
            file=open(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\DRD1_delta_RRCS-{}-{}.tsv".format(pdb.lower(),active_pdb),"w")
            for i in active_dict:
                if i in inactive_dict:
                    #print (i,float(active_dict[i])-float(inactive_dict[i]))
                    file.write("{}\t{}\tRepacking\t{}\t{}\n".format(i,float(active_dict[i])-float(inactive_dict[i]),float(active_dict[i]),float(inactive_dict[i])))
                else:
                    #print(i,float(active_dict[i]))
                    file.write("{}\t{}\tSwitching\t{}\t{}\n".format(i,float(active_dict[i]),float(active_dict[i]),0))
            for i in inactive_dict:
                if i not in active_dict:
                    #print (i,-1*float(inactive_dict[i]))
                    file.write("{}\t{}\tSwitching\t{}\t{}\n".format(i,-1*float(inactive_dict[i]),0,float(inactive_dict[i])))
            file.close()
                    
            RRCS_file=r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\DRD1_delta_RRCS-{}-{}.tsv".format(pdb.lower(),active_pdb)
            outfile=r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\DRD1_delta_RRCS-{}-{}_converted_filtered.tsv".format(pdb.lower(),active_pdb)
            print(RRCS_convert(RRCS_file,"DRD1","ADRB2",alignment,outfile,0.2))
        except KeyError:
            continue
        
#%%
pdb_list=["6a93","6a94","6wh4","6wgt"]
active_pdb_list=["6wha"]
all_residues=list(range(1,472))
for active_pdb in active_pdb_list:
    chain="A"
    print(RRCS_residue_obtain(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\5HT2A Files\{}.pdb.cscore".format(active_pdb.lower()),
                          all_residues,chain,
                          r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\5HT2A_active_RRCS_{}.tsv".format(active_pdb),"5HT2A","No",0))
    for pdb in  pdb_list:
        if pdb=="6wh4":
            chain="C"
        elif pdb=="6wgt":
            chain="B"
        else:
            chain="A"
        print(RRCS_residue_obtain(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\5HT2A Files\{}.pdb.cscore".format(pdb.lower()),
                      all_residues,chain,r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\5HT2A_inactive_RRCS_{}.tsv".format(pdb),"5HT2A","No",0))

        inactive_dict=RRCS_to_dict((r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\5HT2A_inactive_RRCS_{}.tsv".format(pdb.lower())))
        active_dict=RRCS_to_dict((r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\5HT2A_active_RRCS_{}.tsv".format(active_pdb)))
        file=open(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\5HT2A_delta_RRCS-{}-{}.tsv".format(pdb.lower(),active_pdb),"w")
        for i in active_dict:
            if i in inactive_dict:
                #print (i,float(active_dict[i])-float(inactive_dict[i]))
                file.write("{}\t{}\tRepacking\t{}\t{}\n".format(i,float(active_dict[i])-float(inactive_dict[i]),float(active_dict[i]),float(inactive_dict[i])))
            else:
                #print(i,float(active_dict[i]))
                file.write("{}\t{}\tSwitching\t{}\t{}\n".format(i,float(active_dict[i]),float(active_dict[i]),0))
        for i in inactive_dict:
            if i not in active_dict:
                #print (i,-1*float(inactive_dict[i]))
                file.write("{}\t{}\tSwitching\t{}\t{}\n".format(i,-1*float(inactive_dict[i]),0,float(inactive_dict[i])))
        file.close()
                
        RRCS_file=r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\5HT2A_delta_RRCS-{}-{}.tsv".format(pdb.lower(),active_pdb)
        outfile=r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\5HT2A_delta_RRCS-{}-{}_converted_filtered.tsv".format(pdb.lower(),active_pdb)
        print(RRCS_convert(RRCS_file,"5HT2A","ADRB2",alignment,outfile,0.2))
#%%
pdb_list=["3uon","5yc8","5zk3","5zk8","5zkb","5zkc"]
active_pdb_list=["6oik"]
all_residues=list(range(1,467))
for active_pdb in active_pdb_list:
    chain="R"
    print(RRCS_residue_obtain(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\ACM2 Files\{}.pdb.cscore".format(active_pdb.lower()),
                          all_residues,chain,
                          r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\ACM2_active_RRCS_{}.tsv".format(active_pdb),"ACM2","No",0))
    for pdb in  pdb_list:

        chain="A"
        print(RRCS_residue_obtain(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\ACM2 Files\{}.pdb.cscore".format(pdb.lower()),
                      all_residues,chain,r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\ACM2_inactive_RRCS_{}.tsv".format(pdb),"ACM2","No",0))

        inactive_dict=RRCS_to_dict((r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\ACM2_inactive_RRCS_{}.tsv".format(pdb.lower())))
        active_dict=RRCS_to_dict((r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\ACM2_active_RRCS_{}.tsv".format(active_pdb)))
        file=open(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\ACM2_delta_RRCS-{}-{}.tsv".format(pdb.lower(),active_pdb),"w")
        for i in active_dict:
            if i in inactive_dict:
                #print (i,float(active_dict[i])-float(inactive_dict[i]))
                file.write("{}\t{}\tRepacking\t{}\t{}\n".format(i,float(active_dict[i])-float(inactive_dict[i]),float(active_dict[i]),float(inactive_dict[i])))
            else:
                #print(i,float(active_dict[i]))
                file.write("{}\t{}\tSwitching\t{}\t{}\n".format(i,float(active_dict[i]),float(active_dict[i]),0))
        for i in inactive_dict:
            if i not in active_dict:
                #print (i,-1*float(inactive_dict[i]))
                file.write("{}\t{}\tSwitching\t{}\t{}\n".format(i,-1*float(inactive_dict[i]),0,float(inactive_dict[i])))
        file.close()
                
        RRCS_file=r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\ACM2_delta_RRCS-{}-{}.tsv".format(pdb.lower(),active_pdb)
        outfile=r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\ACM2_delta_RRCS-{}-{}_converted_filtered.tsv".format(pdb.lower(),active_pdb)
        print(RRCS_convert(RRCS_file,"ACM2","ADRB2",alignment,outfile,0.2))
        
        
#%%
pdb_list=["4iaq","4iar","5v54","7c61_corrected"]
active_pdb_list=["6g79"]
all_residues=list(range(1,391))
for active_pdb in active_pdb_list:
    chain="S"
    print(RRCS_residue_obtain(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\5HT1B Files\{}.pdb.cscore".format(active_pdb.lower()),
                          all_residues,chain,
                          r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\5HT1B_active_RRCS_{}.tsv".format(active_pdb),"5HT1B","No",0))
    for pdb in  pdb_list:

        chain="A"
        print(RRCS_residue_obtain(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\5HT1B Files\{}.pdb.cscore".format(pdb.lower()),
                      all_residues,chain,r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\5HT1B_inactive_RRCS_{}.tsv".format(pdb),"5HT1B","No",0))

        inactive_dict=RRCS_to_dict((r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\5HT1B_inactive_RRCS_{}.tsv".format(pdb.lower())))
        active_dict=RRCS_to_dict((r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\5HT1B_active_RRCS_{}.tsv".format(active_pdb)))
        file=open(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\5HT1B_delta_RRCS-{}-{}.tsv".format(pdb.lower(),active_pdb),"w")
        for i in active_dict:
            if i in inactive_dict:
                #print (i,float(active_dict[i])-float(inactive_dict[i]))
                file.write("{}\t{}\tRepacking\t{}\t{}\n".format(i,float(active_dict[i])-float(inactive_dict[i]),float(active_dict[i]),float(inactive_dict[i])))
            else:
                #print(i,float(active_dict[i]))
                file.write("{}\t{}\tSwitching\t{}\t{}\n".format(i,float(active_dict[i]),float(active_dict[i]),0))
        for i in inactive_dict:
            if i not in active_dict:
                #print (i,-1*float(inactive_dict[i]))
                file.write("{}\t{}\tSwitching\t{}\t{}\n".format(i,-1*float(inactive_dict[i]),0,float(inactive_dict[i])))
        file.close()
                
        RRCS_file=r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\5HT1B_delta_RRCS-{}-{}.tsv".format(pdb.lower(),active_pdb)
        outfile=r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\5HT1B_delta_RRCS-{}-{}_converted_filtered.tsv".format(pdb.lower(),active_pdb)
        print(RRCS_convert(RRCS_file,"5HT1B","ADRB2",alignment,outfile,0.2))
"""
#%%
actives_directory=r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\Only Actives\\"
active_pdb_list=["6k41","6k42"]
chain="R"
all_residues=list(range(1,451))
for active_pdb in active_pdb_list:
    print(RRCS_residue_obtain(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\original scores\{}.pdb.cscore".format(active_pdb.lower()),
                          all_residues,chain,
                          r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\RRCS_tables\All Residues\ADA2B_active_RRCS_{}.tsv".format(active_pdb),"ADA2B","No",0))

#%%
import os
for filename in os.listdir(os.getcwd()):
    if "_active" not in filename:
        continue
    gene_name=filename.split("_")[0]
    print(gene_name)
    RRCS_file=actives_directory+filename
    outfile=actives_directory+filename.split(".")[0]+"-ADRB2_converted.tsv"
    print(RRCS_convert(RRCS_file,gene_name,"ADRB2",alignment,outfile,0.2))