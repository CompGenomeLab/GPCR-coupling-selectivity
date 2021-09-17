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
    #print ("We obtained",count,"of interactions in total.")
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

def obtain_simulation_contacts(result_dict,simulation_name,residue_list,step):
    result_dict[simulation_name]=[]
    for frame_no in range(0,1001,step):
        directory=r"D:\Users\suuser\Documents\All_Monomer_Systems\{}_PDB\{}_{}.pdb.cscore".format(simulation_name,simulation_name,frame_no)
        contact_score=RRCS_residue_obtain(directory,residue_list,"A")
        result_dict[simulation_name].append(float(contact_score))
    return result_dict

def obtain_cluster_contacts(simulation_name,residue_list,repli_no):
    result_list=[]
    for repli in range(1,repli_no+1):
        directory=r"D:\Users\suuser\Desktop\Cluster_Structures\Contact Scores\{}_cluster{}.pdb.cscore".format(simulation_name,repli)
        contact_score=RRCS_residue_obtain(directory,residue_list,"A")
        result_list.append(float(contact_score))
    return result_list

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
    
def mutation_test(residue_list,mutation,replicate_num,step):
    output_file_name="G315"+mutation+"_"+str(int(1000/step)+1)+".tsv"
    impact_file=open(r"D:\Users\suuser\Documents\All_Monomer_Systems\{}".format(output_file_name),"w")
    WT_replicates=[]
    MT_replicates=[]
    for i in range(1,replicate_num+1):
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
            result_dict_wt=obtain_simulation_contacts(result_dict_wt,simulation,pair,step)
        for key in result_dict_wt:
            WT_contacts+=result_dict_wt[key]

        for simulation in simulation_list[1]:
            result_dict_wt=obtain_simulation_contacts(result_dict_mt,simulation,pair,step)
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
    
def mutation_cluster_test(residue_list,mutation,replicate_num,output_file):
    impact_file=open(output_file,"w")
    for pair in residue_list:
        print("------{}-{}-----".format(pair[0],pair[1]))
        pair[0]=int(GPCR_convert(pair[0],"ADRB2"))
        pair[1]=int(GPCR_convert(pair[1],"ADRB2"))
        WT_contacts=obtain_cluster_contacts("WT",pair,7)
        MT_contacts=obtain_cluster_contacts("G315"+mutation,pair,7)
        print("Comparing contact scores...")
        print(MT_contacts)
        print(WT_contacts)
        statistical_test=df_pairwise_t_test(pair,0.05,MT_contacts,WT_contacts)
        line="{}\t{}\t{}\t{}\n".format(GPCR_convert(pair[0],"ADRB2"),GPCR_convert(pair[1],"ADRB2"),statistical_test[0],statistical_test[1])
        impact_file.write(line)
        print(line)
    impact_file.close()
    
def get_all_residue_pairs(simulation_list,step,chain):
    residue_pair_list=[]
    for simulation_name in simulation_list:
        for frame_no in range(0,1001,step):
            directory=r"D:\Users\suuser\Documents\All_Monomer_Systems\{}_PDB\{}_{}.pdb.cscore".format(simulation_name,simulation_name,frame_no)
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
    
#%%
result_dict={}
result2_dict={}
simulation_list=["G315C_1","G315C_2","G315C_3","G315C_4","G315C_5","G315C_6","G315C_7",
                 "G315L_1","G315L_2","G315L_3","G315L_4","G315L_5","G315L_6","G315L_7",
                 "G315Q_1","G315Q_2","G315Q_3", "G315Q_4","G315Q_5","G315Q_6","G315Q_7",
                 "WT_1","WT_2","WT_3", "WT_4","WT_5","WT_6","WT_7"]
residue_list=[332,326]
for sim in simulation_list:
    result_dict=obtain_simulation_contacts(result_dict,sim,residue_list,100)
    
for sim in result_dict:
    state=sim.split("_")[0]
    if state not in result2_dict:
        result2_dict[state]=result_dict[sim]
    else:
        result2_dict[state]+=result_dict[sim]
#%%
import pandas as pd
df=pd.DataFrame.from_dict(result2_dict)
import plotly.express as px
fig = px.box(df,y=["G315C","G315L","G315Q","WT"],width=1600, height=1000,points="all")
fig.show()
fig.write_html(r"D:\Users\suuser\Desktop\Simulation_contact_score_analysis\7x53-8x50_box_77vs77.html")
#%%
simulation_list=[["WT_1","WT_2","WT_3"],["G315Q_1","G315Q_2","G315Q_3"]]
impact_file=open("Monomer_Leucine_GsNetwork.tsv","w")
for pair in residue_pair_list:
    print("------{}-{}-----".format(GPCR_convert(pair[0],"ADRB2"),GPCR_convert(pair[1],"ADRB2")))
    WT_contacts=[]
    MT_contacts=[]
    result_dict_wt={}
    result_dict_mt={}
    print("Obtaining WT contact scores...")
    for simulation in simulation_list[0]:
        result_dict_wt=obtain_simulation_contacts(result_dict_wt,simulation,pair)
    for key in result_dict_wt:
        WT_contacts+=result_dict_wt[key]
    print("WT contacts scores are ready!")
    print("Obtaining MT contact scores...")
    for simulation in simulation_list[1]:
        result_dict_wt=obtain_simulation_contacts(result_dict_mt,simulation,pair)
    for key in result_dict_mt:
        MT_contacts+=result_dict_mt[key]
    print("Comparing contact scores...")
    statistical_test=df_pairwise_t_test(pair,0.05,MT_contacts,WT_contacts)
    line="{}\t{}\t{}\t{}\n".format(GPCR_convert(pair[0],"ADRB2"),GPCR_convert(pair[1],"ADRB2"),statistical_test[0],statistical_test[1])
    impact_file.write(line)
impact_file.close()
    

#%%
common_activation_layer1=[["6x48","3x40"],["6x48","7x45"],["3x39","2x50"],["5x51","6x44"],["6x44","6x48"]]
common_activation_layer1and2=[["7x45","7x49"],["7x45","6x44"],["2x50","2x46"],["2x50","7x49"]]
common_activation_layer2=[["5x55","6x41"],["6x41","3x43"],["3x43","6x40"],
                          ["3x43","7x49"],["7x50","1x49"]]
common_activation_layer2and3=[["3x43","7x53"],["7x50","7x55"]]
common_activation_layer3=[["6x37","3x46"],["3x46","7x53"],["7x52","7x53"],["7x53","1x53"],["7x53","2x43"]
                          ,["7x53","8x50"],["1x53","7x54"],["8x50","7x54"],["7x54","8x51"]]
common_activation_layer3and4=[["5x57","3x51"],["6x37","5x62"],["6x37","3x50"],["3x46","3x50"],["7x53","3x50"]]
common_activation_layer4=[["3x49","3x50"],["3x50","3x53"]]
gpcrdb=[["2x41","6x38"],["3x44","7x53"]]
common_activation_list=common_activation_layer1+common_activation_layer1and2+common_activation_layer2+common_activation_layer2and3+common_activation_layer3+common_activation_layer3and4+common_activation_layer4

# common_activation_list=[["3x43","6x40"]]
#mutation_cluster_test(common_activation_list,"C",7,"Cysteine_cluster_common_activation.tsv")
mutation_test(common_activation_list,"L",7,100)

#%%
residue_pair_list=[[286,315],[286,290],[318,120],[278,326],[274,271],[271,275],
                   [271,223],[268,269],[268,141],[132,222],[219,275],[219,222],
                   [124,127],[118,161],[119,161],[161,158],[173,168],[80,51],
                   [336,340],[336,53],[54,58],[58,333],[333,64],[64,69],[269,141]]

mutation_test(residue_pair_list,"Q",7,200)

#%%
simulation_list=["G315C_1","G315C_2","G315C_3","G315C_4","G315C_5","G315C_6","G315C_7",
                 "G315L_1","G315L_2","G315L_3","G315L_4","G315L_5","G315L_6","G315L_7",
                 "G315Q_1","G315Q_2","G315Q_3", "G315Q_4","G315Q_5","G315Q_6","G315Q_7",
                 "WT_1","WT_2","WT_3", "WT_4","WT_5","WT_6","WT_7"]
total_res_list=get_all_residue_pairs(simulation_list,100,"A")
#%%
mutation_test(total_res_list,"C",7,100)

#%% Inactive Part
def get_all_residue_pairs(simulation_list,step,chain):
    residue_pair_list=[]
    for simulation_name in simulation_list:
        for frame_no in range(0,1001,step):
            directory=r"D:\Users\suuser\Desktop\Inactive Trajectories\\{}_{}.pdb.cscore".format(simulation_name,frame_no)
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
simulation_list=['G315Q_I_1', 'G315C_I_1', 'G315L_I_1', 'WT_I_1', 'G315Q_I_2', 'G315C_I_2', 'G315L_I_2', 'WT_I_2']
total_res_list=get_all_residue_pairs(simulation_list,100,"A")
#%%
def obtain_simulation_contacts(result_dict,simulation_name,residue_list,step):
    result_dict[simulation_name]=[]
    for frame_no in range(0,1001,step):
        directory=r"D:\Users\suuser\Desktop\Inactive Trajectories\{}_{}.pdb.cscore".format(simulation_name,frame_no)
        contact_score=RRCS_residue_obtain(directory,residue_list,"A")
        result_dict[simulation_name].append(float(contact_score))
    return result_dict

def mutation_test(residue_list,mutation,replicate_num,step):
    output_file_name="G315"+mutation+"_"+str(int(1000/step)+1)+".tsv"
    impact_file=open(r"D:\Users\suuser\Desktop\Inactive Trajectories\{}.tsv".format(output_file_name),"w")
    WT_replicates=[]
    MT_replicates=[]
    for i in range(1,replicate_num+1):
        WT_replicates.append("WT_I_"+str(i))
        MT_replicates.append("G315"+mutation+"_I_"+str(i))
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
        simulation_list=[WT_replicates,MT_replicates]
        for simulation in simulation_list[0]:
            result_dict_wt=obtain_simulation_contacts(result_dict_wt,simulation,pair,step)
        for key in result_dict_wt:
            WT_contacts+=result_dict_wt[key]

        for simulation in simulation_list[1]:
            result_dict_wt=obtain_simulation_contacts(result_dict_mt,simulation,pair,step)
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
    
mutation_test(total_res_list,"L",2,100)