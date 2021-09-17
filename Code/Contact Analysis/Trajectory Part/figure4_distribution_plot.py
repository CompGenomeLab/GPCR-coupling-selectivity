# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 15:54:00 2021

@author: bselcuk
"""

x1=[]
x2=[]
x3=[]
x4=[]

file=open(r"D:\Users\suuser\Desktop\GPCRdb_activity.txt","r")

count=0
for line in file:
    if count==0:
        count+=1
        continue
    if count%10!=0:
        count+=1
        continue
    count+=1
    line=line.strip()
    line_list=line.split("\t")
    x1.append(float(line_list[0]))
    x2.append(float(line_list[2]))
    x3.append(float(line_list[1]))
    x4.append(float(line_list[3]))
file.close()

#%%

import plotly.figure_factory as ff
group_labels = ['Group 1', 'Group 2', 'Group 3', 'Group 4']
hist_data=[x1,x2,x3,x4]
# Create distplot with custom bin_size
fig = ff.create_distplot(hist_data, group_labels, show_hist=False)

# Add title
fig.update_layout(title_text='Curve and Rug Plot')
fig.show()
    