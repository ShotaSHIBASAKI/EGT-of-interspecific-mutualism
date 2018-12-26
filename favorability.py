#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 16 07:22:15 2018
@author: shibasakishota

Calculate the probability that each species becomes selfish given M
"""

import numpy as np
from scipy.stats import norm
from scipy.stats import uniform 
import matplotlib.pyplot as plt
import math

def main (M):
    k_array=np.array([0.5,1.5])
    species=np.linspace(0,M-1, M)
    for i in range (np.size(k_array)):
        k=k_array[i]
        m=math.ceil((M-1)/(3-k)) # maximum number of generous species
        print(str('maximum value of generous is %d' %(m)))
        prob_f1=np.zeros([m+1])#p(f1)
        prob_f1[0]=1#initial condition
        q_fav=np.zeros([M])# favorability or probability to be selfish
        q_cond=-np.zeros([m+1])#conditional probability to be selfish given f1
        for j in range(M-2):
            for f1 in range(m+1):
                #calculate conditional probability to be selfish
                Tf=((M-1)/(3-k)-f1)/(M-j-1)
                sd=(1/(12*(M-j-1)))**0.5 
                q_cond[f1]=1-norm.cdf(Tf,0.5,sd)# conditional prob to be selfish
                q_fav[j]+=prob_f1[f1]*q_cond[f1]
            prob_new=np.zeros([m+1])
            for f1 in range (m+1):
                #calculate p(f1) for the next species
               if f1==0:
                   prob_new[f1]=prob_f1[f1]*q_cond[f1]
               else:
                   prob_new[f1]=prob_f1[f1-1]*(1-q_cond[f1-1])+prob_f1[f1]*q_cond[f1] 
            prob_f1=prob_new#update
            
        #second slowest (M-1)species
        for f1 in range(m+1):
            Tf=((M-1)/(3-k)-f1)  
            q_cond[f1]=1-uniform.cdf(Tf,0,1)                
            q_fav[M-2]+=prob_f1[f1]*q_cond[f1]
            #print(q_fav[M-2])
        prob_new=np.zeros([m+1])    
        for f1 in range (m+1):
                #calculate p(f1) for the next species
               if f1==0:
                   prob_new[f1]=prob_f1[f1]*q_cond[f1]
               else:
                   prob_new[f1]=prob_f1[f1-1]*(1-q_cond[f1-1])+prob_f1[f1]*q_cond[f1] 
        prob_f1=prob_new#update
        #slowest M species
        q_fav[M-1]=prob_f1[m]
        #plotthe results
        #print(q_fav)
        if i==0:
            plt.plot(species, q_fav[::-1], color="c",marker="o",markersize=10)
        else:
            plt.plot(species, q_fav[::-1], color="b",marker="D",markersize=10)
    plt.xticks([0,M-1],["slowest","fastest"],fontsize=12)
    plt.xlabel("evolutionary rate",fontsize=12)
    plt.ylim(0.0,1)
    plt.ylabel("favorablity",fontsize=12)
    fname=str('favorability_M%d'%(M))
    plt.savefig(fname+'.pdf')
    plt.show()
         
            
            
        