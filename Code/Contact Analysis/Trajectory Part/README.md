# Trajectory Analysis

In this part you can find the scripts that we used to analyze molecular dynamics simulation trajectories. We calculated GPCRdb 
distances for each simulation and concataneted all of the replicates. GPCRdb distances are given in 'GPCRdb_distances.txt'. 

**Contact analysis of active simulations**
```python
#Active Part
simulation_list=["G315C_1","G315C_2","G315C_3","G315C_4","G315C_5","G315C_6","G315C_7",
                 "G315L_1","G315L_2","G315L_3","G315L_4","G315L_5","G315L_6","G315L_7",
                 "G315Q_1","G315Q_2","G315Q_3", "G315Q_4","G315Q_5","G315Q_6","G315Q_7",
                 "WT_1","WT_2","WT_3", "WT_4","WT_5","WT_6","WT_7"]

# You should modify the directory within the script to run it.
directory=r"D:\Users\suuser\Desktop\gitdocs\GPCR-coupling-selectivity\Code\Contact Analysis\Trajectory Part\Active Trajectories"

total_res_list=get_all_residue_pairs(simulation_list,100,"A",directory) #Here we obtain every contact that was observed in active-state simılations.


mutation_test(total_res_list,"C",7,100,directory,"active") # 7 indicates the number of replicates

# mutation_test(total_res_list,"L",7,100,"active")
# mutation_test(total_res_list,"Q",7,100,"active")
```

**Contact analysis of inactive simulations**
```python
#%% Inactive Part

simulation_list=['G315Q_I_1', 'G315C_I_1', 'G315L_I_1', 'WT_I_1', 'G315Q_I_2', 'G315C_I_2', 'G315L_I_2', 'WT_I_2']
# You should modify the directory.
directory=r"D:\Users\suuser\Desktop\gitdocs\GPCR-coupling-selectivity\Code\Contact Analysis\Trajectory Part\Inactive Trajectories"
total_res_list=get_all_residue_pairs(simulation_list,100,"A",directory)
#Here we obtain every contact that was observed in active-state simılations.

mutation_test(total_res_list,"C",2,100,directory,"inactive")  #2 indicates the number of replicates

# mutation_test(total_res_list,"L",2,100,directory,"inactive")
# mutation_test(total_res_list,"Q",2,100,directory,"inactive")
```