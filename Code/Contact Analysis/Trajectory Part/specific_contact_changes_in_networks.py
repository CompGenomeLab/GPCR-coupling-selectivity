# -*- coding: utf-8 -*-
"""
Created on Tue May 25 13:23:59 2021

@author: bselcuk
"""

def RRCS_residue_obtain(RRCS_file,residue_list,chain):
    file=open(RRCS_file,"r")
    for line in file:
        if line.count(chain+":")==2:
            line_list=line.strip().split()
            interactor1=line_list[0]
            interactor2=line_list[1]
            value=line_list[2]
            res_no1=interactor1.split(":")[1].split("_")[0]
            res_no2=interactor2.split(":")[1].split("_")[0]
            if int(res_no1) not in residue_list or int(res_no2) not in residue_list:
                continue            
            else:
                return value
    file.close()
    return 0

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

def obtain_simulation_contacts(result_dict,simulation_name,residue_list,step,contact_directory):
    result_dict[simulation_name]=[]
    for frame_no in range(0,1001,step):
        directory=r"{}\{}_{}.pdb.cscore".format(contact_directory,simulation_name,frame_no)
        contact_score=RRCS_residue_obtain(directory,residue_list,"A")
        result_dict[simulation_name].append(float(contact_score))
    return result_dict


def df_pairwise_t_test(residue_pair,p_value,first_list,second_list):
    from scipy.stats import ttest_ind
    comparison_stat=ttest_ind(first_list,second_list,equal_var=False)
    if p_value>=comparison_stat[1]:
        print (str(residue_pair[0])+"-"+str(residue_pair[1])+"\t"+str(comparison_stat[0])+"\t"+str(comparison_stat[1])+"\tSignificant")
    else:
        print (str(residue_pair[0])+"-"+str(residue_pair[1])+"\t"+str(comparison_stat[0])+"\t"+str(comparison_stat[1])+"\tNot Significant")
    return [comparison_stat[0],comparison_stat[1]]

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
    
def mutation_test(residue_list,mutation,replicate_num,step,contact_directory,active_or_inactive):
    output_file_name="G315"+mutation+"_"+str(int(1000/step)+1)+".tsv"
    impact_file=open(r"{}\{}".format(contact_directory,output_file_name),"w")
    WT_replicates=[]
    MT_replicates=[]
    for i in range(1,replicate_num+1):
        if active_or_inactive=="inactive":
            WT_replicates.append("WT_I_"+str(i))
            MT_replicates.append("G315"+mutation+"_I_"+str(i))
        else:
            WT_replicates.append("WT_"+str(i))
            MT_replicates.append("G315"+mutation+"_"+str(i))
    simulation_list=[WT_replicates,MT_replicates]
    for pair in residue_list:
        print("------{}-{}-----".format(pair[0],pair[1]))
        if type(pair[0])==str:
            pair[0]=int(GPCR_convert(pair[0],"ADRB2"))
        if type(pair[1])==str:
            pair[1]=int(GPCR_convert(pair[1],"ADRB2"))
        WT_contacts=[]
        MT_contacts=[]
        result_dict_wt={}
        result_dict_mt={}

        for simulation in simulation_list[0]:
            result_dict_wt=obtain_simulation_contacts(result_dict_wt,simulation,pair,step,contact_directory)
        for key in result_dict_wt:
            WT_contacts+=result_dict_wt[key]

        for simulation in simulation_list[1]:
            result_dict_wt=obtain_simulation_contacts(result_dict_mt,simulation,pair,step,contact_directory)
        for key in result_dict_mt:
            MT_contacts+=result_dict_mt[key]
        print("Comparing contact scores...")
        statistical_test=df_pairwise_t_test(pair,0.05,MT_contacts,WT_contacts)
        line="{}\t{}\t{}\t{}\t".format(GPCR_convert(pair[0],"ADRB2"),GPCR_convert(pair[1],"ADRB2"),statistical_test[0],statistical_test[1])
        if statistical_test[1]<0.05:    
            impact_file.write(line+"Significant\n")
        else:
            impact_file.write(line+"Non-significant\n")
    impact_file.close()
    
    
def get_all_residue_pairs(simulation_list,step,chain,contact_directory):
    residue_pair_list=[]
    for simulation_name in simulation_list:
        for frame_no in range(0,1001,step):
            directory=r"{}\{}_{}.pdb.cscore".format(contact_directory,simulation_name,frame_no)
            file=open(directory,"r")
            for line in file:
                if line.count(chain+":")==2:
                    line_list=line.strip().split()
                    interactor1=line_list[0]
                    interactor2=line_list[1]
                    value=line_list[2]
                    res_no1=interactor1.split(":")[1].split("_")[0]
                    res_no1=int(res_no1)
                    res_no2=interactor2.split(":")[1].split("_")[0]
                    res_no2=int(res_no2)
                    if [res_no1,res_no2] not in residue_pair_list:
                        residue_pair_list.append([res_no1,res_no2])
            file.close()
    return residue_pair_list
    

#%% Active Part
simulation_list=["G315C_1","G315C_2","G315C_3","G315C_4","G315C_5","G315C_6","G315C_7",
                 "G315L_1","G315L_2","G315L_3","G315L_4","G315L_5","G315L_6","G315L_7",
                 "G315Q_1","G315Q_2","G315Q_3", "G315Q_4","G315Q_5","G315Q_6","G315Q_7",
                 "WT_1","WT_2","WT_3", "WT_4","WT_5","WT_6","WT_7"]

directory=r"D:\Users\suuser\Desktop\gitdocs\GPCR-coupling-selectivity\Code\Contact Analysis\Trajectory Part\Active Trajectories"

total_res_list=get_all_residue_pairs(simulation_list,100,"A",directory)
#Here we obtain every contact that was observed in active-state simılations.
#%%
mutation_test(total_res_list,"C",7,100,directory,"active")
# 7 indicates the number of replicates
# mutation_test(total_res_list,"L",7,100,"active")
# mutation_test(total_res_list,"Q",7,100,"active")

#%% Inactive Part

simulation_list=['G315Q_I_1', 'G315C_I_1', 'G315L_I_1', 'WT_I_1', 'G315Q_I_2', 'G315C_I_2', 'G315L_I_2', 'WT_I_2']
directory=r"D:\Users\suuser\Desktop\gitdocs\GPCR-coupling-selectivity\Code\Contact Analysis\Trajectory Part\Inactive Trajectories"
total_res_list=get_all_residue_pairs(simulation_list,100,"A",directory)
#Here we obtain every contact that was observed in active-state simılations.
#%%
mutation_test(total_res_list,"C",2,100,directory,"inactive")  
# 2 indicates the number of replicates
# mutation_test(total_res_list,"L",2,100,directory,"inactive")
# mutation_test(total_res_list,"Q",2,100,directory,"inactive")