# G protein specific networks

**RRCS_parser.py usage: **\
Here we are choosing to which G protein to analyze:
```python
#You should comment out the G protein list you want to use

Gtarget_list=Gs_list
# Gtarget_list=Gi_Gio_intersection_list
# Gtarget_list=Gq_list
# Gtarget_list=Gi1_list
# Gtarget_list=Gio_list
# Gtarget_list=Gq_list
# Gtarget_list=list(set(Gi_list+Gq_list+Gio_list))
```

In this part we are choosing structures that contain our G protein of interest:
```python
#You should comment out gene name lists based on your previous choice

#Gi1 vs others
# first_list=["DRD2","DRD3"]
# second_list=["ADRB2","DRD1","5HT2A","5HT1B","ACM2","HRH1"]


#Gs vs others
first_list=["ADRB2","DRD1"]
second_list=["DRD2","ACM2","5HT2A","DRD3","HRH1","5HT1B"]

#Gio vs others
# first_list=["ACM2","5HT1B"]
# second_list=["ADRB2","DRD1","5HT2A","DRD2","DRD3","HRH1"]

#Gi+Gio vs Gq+Gs
# first_list=["DRD2","DRD3","ACM2","5HT1B"]
# second_list=["ADRB2","DRD1","5HT2A","HRH1"]

#Gq vs others
# first_list=["5HT2A","HRH1"]
# second_list=["DRD2","DRD1","ACM2","5HT1B","ADRB2"]

#Gi1+Gio+Gq vs Gs
# first_list=["DRD2","5HT1B","ACM2","DRD3","5HT2A","HRH1"]
# second_list=["ADRB2","DRD1"]


df=df_pairwise_t_test(10**-2,first_list,second_list,df,2)
#10**-2 is the threshold we use for the p value. It is adjusted to 0.1 when analyzing Gi1 contacts.
```