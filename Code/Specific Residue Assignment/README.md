# Specific Residue Assignment

Requirements:
- A multiple sequence alignment containing orthologs of the receptors you want to compare
- Header format should be `">sp or tr|UniprotID|GeneName_SpeciesName|SpeciesTaxID"` e.g. `"\>sp|P08908|5HT1A_HUMAN|9606"`
- Sequences should be ordered as following: 
\>sp|P08908|5HT1A_HUMAN|9606 
\>5HT1A orthologs 
\>sp|P07550|ADRB2_HUMAN|9606 
\>ADRB2 orthologs 
\>sp|P21728|DRD1_HUMAN|9606 
\>DRD1 orthologs 

Usage of specific_residue_find.py:

```python
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
```

