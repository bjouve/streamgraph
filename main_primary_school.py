import numpy as np
import scipy as sc 
import matplotlib.pyplot as plt
import random as random
import math as math
from heapq import heappush, heappop
from itertools import count
import networkx as nx
import random
import operator 
import pandas as pd
import fct

dir='/Users/mehdidjellabi/Desktop/code_article/'
mon_fichier=open(dir+'primary_school_data/primaryschool.txt', "r")
contenu=mon_fichier.readlines()
mon_fichier.close()


########Data importation ########################
temp_edge=[]
classe=dict()
for l in range(len(contenu)):
    ligne = contenu[l]
    chaine=(ligne.split("\t"))
    u=int(chaine[1])
    v=int(chaine[2])
    edge=(int(chaine[0]),u,v,)
    temp_edge.append(edge)
    classe[u]=chaine[3]
    classe[v]=chaine[4][:-1]

################################################
day_min=54940
day=dict()
i=0
for d in range(1,3):
    j=[]
    while i < len(temp_edge)-1 and temp_edge[i+1][0]-temp_edge[i][0]<day_min:
        j.append(temp_edge[i])
        i+=1
    day[d]=j
    i+=1

t0={i : day[i][0][0] for i in range(1,3)}

day_0={d:[(day[d][i][0]-t0[d],day[d][i][1],day[d][i][2]) for i in range(len(day[d]))]
    for d in range(1,3)}

T_d={d : set([t for t,u,v in day_0[d]]) for d in range(1,3)}



####################################################################
####################################################################
####################################################################
####################################################################
####################################################################
####################################################################
####################################################################
####################################################################
####################################################################
####################################################################










##chose the day over which to carry the calculation
d=1
######chose the window size######
largeur=300
#chose weather there is overlap or not if slide = largeure, it means there is no overlap between windows
slide=largeur


N_null = 100    ##The number of null models to average the higher the more precise and slower
smoothing=1     ##No smoothing when value is 1
seed=random.randint(0,1000)









####################################################################
####################################################################
####################################################################
####################################################################
####################################################################
####################################################################
####################################################################
####################################################################
####################################################################
####################################################################



#########            Code execution          #######################
day_0_graph=dict()
for d in range(d,d+1):
    G_t=dict()
    t_max=max(T_d[d])
    t1=0
    t2=largeur
    while t2<=t_max:
        G=nx.Graph()
        edge_list=[(edge[1],edge[2]) for edge in day_0[d] if t2>edge[0]>=t1]
        G.add_edges_from(edge_list)
        G_t[t1]=G
        t1+=slide
        t2+=slide
    G=nx.Graph()
    edge_list=[(edge[1],edge[2]) for edge in day_0[d] if t_max>=edge[0]>=t1]
    G.add_edges_from(edge_list)
    G_t[t1]=G
    day_0_graph[d]=G_t
####################################################################
####################################################################
T_uv={d:{(u,v):set() for u in classe for v in classe} for d in range(d,d+1)}
for d in range(d,d+1):
    for t,u,v in day_0[d]:
        T_uv[d][(u,v)].add(t)
        T_uv[d][(v,u)].add(t)

T_uv_inter=dict()

for d in range(d,d+1):
    T_uv_inter[d]=dict()
    for e in T_uv[d]:
        ls=sorted(list(T_uv[d][e]))
        if len(ls)>0:
            dure_liens=fct.duree(ls,20)
            i0=0
            intervals=[]
            for l in dure_liens:
                intervals.append((ls[i0],ls[i0+l-1]+20))
                i0+=l
            
            T_uv_inter[d][e]=intervals
####################################################################
####################################################################

segmente_d_t_temp=dict()
W_t_d=dict()

for d in range(d,d+1):
    segmente_d_t_temp[d]=dict()
    W_t_d[d]=dict()
    for t in day_0_graph[d]:
        t1=t
        t2=t+largeur
        W=fct.edge_constr_dynamic(T_uv[d],t1,t2)
        W_t_d[d][t]=W
        f = open(dir+'results/Input_W_primary_school/W_day_%d_t0_%d.txt'%(d,t), "w")
        f.write('day %d \t t_0= %d \n' %(d,t))
        for e in W.edges():
            u=e[0]
            v=e[1]
            f.write(str(u)+'\t'+str(v)+'\t'+str(W.get_edge_data(u,v)['weight'])+'\n')
        f.close()

        segmente_temp=fct.segmentation_temporal(W,N_null,smoothing,seed)
        segmente_d_t_temp[d][t]=segmente_temp

        f = open(dir+'results/RC_primary_school/segmentation_day_%d_t0_%d.txt'%(d,t), "w")
        f.write('t %d \n' %(t))
        if segmente_d_t_temp[d][t] != None:
            segmente=segmente_d_t_temp[d][t]
            for i in range(len(segmente)):
                seg=segmente[i]
                f.write('RC ')
                for u in list(seg[0].nodes()):
                    f.write('%d '%(u))
                f.write('\n')
                f.write('Q %f %f \n'%(seg[5],seg[6]))

        else :
            f.write('RC \n')
            f.write('Q \n')
        f.close()


