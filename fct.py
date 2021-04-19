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

def weight_shuffle(G,weights):
    random.shuffle(weights)
    null=nx.Graph()
    edges=list(G.edges())
    for i in range(len(edges)):
        edge=edges[i]
        poids=weights[i]
        null.add_edge(*edge,weight=poids)
    return null

def strength(G,node):
    force=0.
    neigh=G.neighbors(node)
    for vois in neigh:
        force+= G.get_edge_data(node,vois)['weight']
    return force

def delta(G,node):
    return strength(G,node)*1./(len(G)-1)

def harmonic_mean(a,b):
    if a ==0 and b==0:
        return 0.
    else:
        return 2*a*b*1./(a+b)

def overlap(A,B):
    return len(set(A).intersection(set(B)))*1./min([len(set(A)),len(set(B))])

def kronecker(a,b):
    if a==b :
        return 1
    else:
        return 0 

def lissage(Lx,Ly,p):
    Lxout=[]
    Lyout=[]
    for i in range(p,len(Lx)-p):   
        Lxout.append(Lx[i])
    for i in range(p,len(Ly)-p):
        val=0
        for k in range(2*p):
            val+=Ly[i-p+k]
        Lyout.append(val/2/p)
            
    return (Lxout,Lyout)

def rand_top_wight(randG,weights):
    randW=nx.Graph()
    randW.add_nodes_from(list(randG.nodes()))
    w0=weights.copy()
    for edge in randG.edges():
        w=random.choice(w0)
        randW.add_edge(*edge,weight=w)
        w0.remove(w)
    return randW

def edge_constr(G):
    deg=nx.degree(G)
    N=len(G)
    W=nx.Graph()
    topologiecal_poids={}
    if len(G)>2:
        for edge in G.edges():
            i=edge[0]
            j=edge[1]
            topologiecal_poids[edge]=len(set(G.neighbors(i)).intersection(set(G.neighbors(j))))*harmonic_mean(G.degree(i),G.degree(j))
    for edge in topologiecal_poids.keys():
        W.add_edge(min(edge[0], edge[1]),max(edge[0], edge[1]), weight=topologiecal_poids[edge])
    return W

def indicatrice(A,e):
    if e in A : 
        return 1.
    else:
        return 0.
        
def omega(i,j,involved_nodes,T_uv,T_reduit):
    t1,t2=min(T_reduit),max(T_reduit)
    w=0.
    i_k=set([k for k in involved_nodes  if k!=i and len(T_uv[(i,k)])>0])
    j_k=set([k for k in involved_nodes  if k!=j and len(T_uv[(j,k)])>0])
    i_l=set([l for l in involved_nodes  if i!= l and len(T_uv[(i,l)])>0])
    j_m=set([m for m in involved_nodes  if j!= m and len(T_uv[(j,m)])>0])
    triples=set([(k,l,m) for k in i_k.intersection(j_k) for l in i_l  for m in j_m if len(set.intersection(T_uv[(i,k)],T_uv[(j,k)],T_uv[(i,l)],T_uv[(j,m)]))>0])
    for t in list(T_reduit):
        num=0.
        denom=0.
        for k,l,m in triples: 
            num+=indicatrice(set.intersection(T_uv[(i,k)], T_uv[(j,k)], T_uv[(i,l)],T_uv[(j,m)]),t)    
        if num>0.:
            k_list=set([k for k in involved_nodes if len(T_uv[(i,k)].union(T_uv[(j,k)]))>0 ])
            for k in k_list:
                denom+=indicatrice(T_uv[(i,k)].union(T_uv[(j,k)]),t)+indicatrice(T_uv[(i,k)].intersection(T_uv[(j,k)]),t)
            if denom>0.:
                w+=num*2./denom 
    return w*1./len(T_reduit)

def edge_constr_dynamic(T_uv,t1,t2):
    involved_edges=list(set([e for e in T_uv for t in T_uv[e] if t2>t>=t1 ]))
    involved_nodes=list(set([elem for e in involved_edges for elem in e]))
    T=list(set().union(*list(T_uv.values())))
    T_reduit=sorted([t for t in T if t2>t>=t1])
    W_t=nx.Graph()
    for i in range(len(involved_nodes)):
        ni=involved_nodes[i]
        for j in range(i+1,len(involved_nodes)):
            nj=involved_nodes[j]
            if (ni,nj) in involved_edges:
                w=omega(ni,nj,involved_nodes,T_uv,T_reduit)
                W_t.add_edge(ni, nj, weight=w)
    return W_t

def segmentation_temporal(W,Nrep,wid,seed):
    sortie=dict()
    scales={}
    u=0
    W0=W.copy()
    Wtst=W.copy()
    compare_qual=[]
    minstrength=math.inf
    S=W.size(weight='weight')
    if S==0:
        return None 
    else:   
        while W.size(weight='weight')>0.:
            weights=[W[u][v]['weight'] for u,v in W.edges()]
            S=W.size(weight='weight')
            strengthlist={node : strength(W,node) for node in W.nodes()}
            trille=sorted(strengthlist.items(), key=operator.itemgetter(1),reverse=True) 
            res=[trille[i][0] for i in range(len(trille))]   
            tot=[]
            A=[]
            if S==0:
                break 
            for node in (res):
                A.append(node)
                tot.append((W.subgraph(A).size(weight='weight'))*1./S)
            mean=[]
            for i in range(Nrep):
                degree_sequence=[]
                for node in W.nodes():
                    degree_sequence.append(W.degree(node))
                randW=nx.configuration_model(degree_sequence, create_using=W.copy())
                randW=rand_top_wight(randW,weights)
                RS=randW.size(weight='weight')
                randstrengthlist={node : strength(randW,node) for node in randW.nodes()}      
                randtrille=sorted(randstrengthlist.items(), key=operator.itemgetter(1),reverse=True)
                randres=[randtrille[i][0] for i in range(len(randtrille))]
                randA=[]
                randtot_i=[]
                for node in (randres):
                    randA.append(node)
                    if RS==0: 
                        randtot_i.append(0)
                    else:
                        randtot_i.append((randW.subgraph(randA).size(weight='weight'))*1./RS)
                mean.append(randtot_i)
            Q_list=[sum([(tot[j]-mean[i][j])*1./len(tot) for j in range(len(tot))]) for i in range(Nrep)]
            Q_mean=sum(Q_list)*1./len(Q_list)
            Q_std_dev=sum([(q-Q_mean)**2 for q in Q_list])*1./(len(Q_list)-1)
            randtot_mean=np.zeros(len(tot))
            for i in range(Nrep):
                act=mean[i]    
                for j in range(len(act)):
                    randtot_mean[j]+=act[j]*1./Nrep
            randtot=list(randtot_mean)
    
            std_dev=[]
            for i in range(len(W)):
                std_i=0
                for iteration in range(Nrep):
                    std_i+=(mean[iteration][i])**2
                std_dev.append(math.sqrt(round(abs((1./Nrep)*(std_i)-(randtot[i])**2),10)))        
            diff=[tot[i]-randtot[i] for i in range(len(tot))]
            liss=diff
            if max(diff)>0.0:

                minode=np.argmax(liss)+wid-1
                ND=[]                      
                for i in range(minode+1,len(res)):
                    ND.append(res[i])

                DW=W.subgraph(ND)
                Dense=[x for x in W if x not in ND]
                mini=min([strength(W,node) for node in Dense])
                for node in ND:
                    if strength(W,node)==mini:
                        Dense.append(node)
                scales[u]=Dense
                Quality=Q_mean
                sortie[u]=[W0.subgraph(Dense),tot,randtot,std_dev,liss,Quality,math.sqrt(Q_std_dev)]
                if len(Dense)>0:
                    W=DW
                    u+=1
                else:
                    break
            else:
                break
        return sortie

def segmentation(G0,Nrep,wid,seed):
    G=G0.copy()
    sortie=dict()
    scales={}
    u=0
    W=edge_constr(G)    
    W0=W.copy()
    Wtst=W.copy()
    compare_qual=[]
    minstrength=math.inf
    S=W.size(weight='weight')
    if S==0:
        return None 
    else:   
        while Wtst.size(weight='weight')>0.:
            weights=[W[u][v]['weight'] for u,v in W.edges()]
            S=W.size(weight='weight')
            strengthlist={node : strength(W,node) for node in W.nodes()}
            trille=sorted(strengthlist.items(), key=operator.itemgetter(1),reverse=True)
            res=[trille[i][0] for i in range(len(trille))]   
            tot=[]
            A=[]
            if S==0:
                break 
            for node in (res):
                A.append(node)
                tot.append((W.subgraph(A).size(weight='weight'))*1./S)
            mean=[]
            for i in range(Nrep):
                degree_sequence=[]
                for node in G.nodes():
                    degree_sequence.append(G.degree(node))
                randG=nx.configuration_model(degree_sequence, create_using=G.copy())
                randW=rand_top_wight(randG,weights)
                RS=randW.size(weight='weight')
                randstrengthlist={node : strength(randW,node) for node in randW.nodes()}      
                randtrille=sorted(randstrengthlist.items(), key=operator.itemgetter(1),reverse=True)
                randres=[randtrille[i][0] for i in range(len(randtrille))]
                randA=[]
                randtot_i=[]
                for node in (randres):
                    randA.append(node)
                    if RS==0: 
                        randtot_i.append(0)
                    else:
                        randtot_i.append((randW.subgraph(randA).size(weight='weight'))*1./RS)
                mean.append(randtot_i)
            Q_list=[sum([(tot[j]-mean[i][j])*1./len(tot) for j in range(len(tot))]) for i in range(Nrep)]
            Q_mean=sum(Q_list)*1./len(Q_list)
            Q_std_dev=sum([(q-Q_mean)**2 for q in Q_list])*1./(len(Q_list)-1)
            randtot_mean=np.zeros(len(tot))
            for i in range(Nrep):
                act=mean[i]    
                for j in range(len(act)):
                    randtot_mean[j]+=act[j]*1./Nrep
            randtot=list(randtot_mean)

            std_dev=[]
            for i in range(len(W)):
                std_i=0
                for iteration in range(Nrep):
                    std_i+=(mean[iteration][i])**2
                std_dev.append(math.sqrt(round(abs((1./Nrep)*(std_i)-(randtot[i])**2),10)))
            diff=[tot[i]-randtot[i] for i in range(len(tot))]
            liss=diff
            if max(diff)>0.0:
                minode=np.argmax(liss)+wid-1
                ND=[]                      
                for i in range(minode+1,len(res)):
                    ND.append(res[i])
                DW=W.subgraph(ND)
                tout=G.nodes()
                t=[x for x in tout if x not in ND]
                scales[u]=t
                Quality=Q_mean
                sortie[u]=[G0.subgraph(t),tot,randtot,std_dev,liss,Quality,math.sqrt(Q_std_dev)]
                print(len(t))
                if len(t)>0:
                    W=DW
                    G=G.subgraph(W)
                    u+=1
                else:
                    break
            else:
                break
        print(" ")
        return sortie
def duree(ls,dt):
    len_seq=[]
    if len(ls)>0:
        i=0
        while i<len(ls)-1:
            j=0
            while i+j+1<=len(ls)-1 and ls[i+j+1]-ls[i+j]==dt :
                j+=1
            if j==0:
                i+=1
                len_seq.append(1)
            else :
                i+=j+1
                len_seq.append((j+1))
    if sum(len_seq)<len(ls):
        len_seq.append(1)
    return len_seq

def duree(ls,dt):
    len_seq=[]
    if len(ls)>0:
        i=0
        while i<len(ls)-1:
            j=0
            while i+j+1<=len(ls)-1 and ls[i+j+1]-ls[i+j]==dt :
                j+=1
            if j==0:
                i+=1
                len_seq.append(1)
            else : 
                i+=j+1
                len_seq.append((j+1))
        if sum(len_seq)<len(ls):
            len_seq.append(1)
    return len_seq



