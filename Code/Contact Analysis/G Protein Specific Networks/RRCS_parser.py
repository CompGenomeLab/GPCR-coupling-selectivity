# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 20:32:01 2021

@author: bselcuk
"""


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


def df_pairwise_t_test(p_value,first_list,second_list,df,mode):
    from scipy.stats import ttest_ind
    from statistics import mean 
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



#%%
#directory should be adjusted 
directory=r"D:\Users\suuser\Desktop\gitdocs\GPCR-coupling-selectivity\Code\Contact Analysis\G Protein Specific Networks\RRCS_tables\all_converted_files\\"
directory=r"D:\Users\suuser\Desktop\G protein selectivity manuscript\false_negative\\"
dict_list=[]
column_list=[]
filename_dict={}
for filename in os.listdir(directory):
    if "converted" not in filename:
        continue
    else:
        pair_list=["ADRB2","DRD2","5HT2A","5HT1B","ACM2","DRD1","DRD3","HRH1"]
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
        label=receptor_name
        if label not in filename_dict:
            filename_dict[label]=1
        else:
            filename_dict[label]+=1
        column_list.append(label)
        
#%%
Gs_list=["1x52","1x57","2x40","2x42","2x51","2x56","2x58","3x21","3x28",
         "3x38","3x43","34x50","34x53","4x38","4x53","173","5x39","5x48",
         "5x58","5x61","5x62","6x29","6x30","6x36","6x37","6x40","6x52","7x41","7x45","8x51","8x58"]
Gi1_list=["2x40","2x45","34x53","4x46","6x61","301","8x47","3x46","12x48","2x53","6x31","5x54","1x47"]
Gq_list=["1x46","1x47","1x52","2x40","2x42","2x45","2x49","2x51","3x28","4x46",
         "4x53","5x48","5x54","5x58","6x37","6x61","301","8x47","8x51","1x49","2x57","3x46","3x43","6x30","6x40"]
Gio_list=["2x40","2x45","34x53","5x48","5x54","3x46","7x45","5x61","5x58","6x40","8x51","3x43","6x37","4x46","8x47"]
Gi_Gio_intersection_list=["2x40","2x45","4x46","8x47","3x46","5x54"]
Gi_Gio_Gq_list=Gi_Gio_intersection_list.copy()

#%%
dataframe={}
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
#%%
#Choose which one to use based on your comparison

Gtarget_list=Gs_list
# Gtarget_list=Gi_Gio_intersection_list
# Gtarget_list=Gq_list
# Gtarget_list=Gi1_list
# Gtarget_list=Gio_list
# Gtarget_list=Gq_list
# Gtarget_list=list(set(Gi_list+Gq_list+Gio_list))
# 
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
#%%
new_keys=keys
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
#%%
#Gi1 vs others
# first_list=["DRD2","DRD3"]
# second_list=["ADRB2","DRD1","5HT2A","5HT1B","ACM2","HRH1"]


#Gs vs others
# first_list=["ADRB2","DRD1"]
# second_list=["DRD2","ACM2","5HT2A","DRD3","HRH1","5HT1B"]

#Gio vs others
# first_list=["ACM2","5HT1B"]
# second_list=["ADRB2","DRD1","5HT2A","DRD2","DRD3","HRH1"]

#Gi+Gio vs Gq+Gs
first_list=["DRD2","DRD3","ACM2","5HT1B"]
second_list=["ADRB2","DRD1","5HT2A","HRH1"]

#Gq vs others
# first_list=["5HT2A","HRH1"]
# second_list=["DRD2","DRD1","ACM2","5HT1B","ADRB2"]

#Gi1+Gio+Gq vs Gs
# first_list=["DRD2","5HT1B","ACM2","DRD3","5HT2A","HRH1"]
# second_list=["ADRB2","DRD1"]


df=df_pairwise_t_test(10**-2,first_list,second_list,df,2)
             
 
