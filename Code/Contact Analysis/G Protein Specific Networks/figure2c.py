# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 14:10:28 2021

@author: bselcuk
"""

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
