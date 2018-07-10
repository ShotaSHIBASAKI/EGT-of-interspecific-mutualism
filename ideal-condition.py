#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 26 06:26:34 2018
compute foavorability under the ideal conditions, see inequality (9)
@author: shibasakishota
"""

import numpy as np
from scipy.stats import norm
from scipy.stats import uniform 
import matplotlib.pyplot as plt
import math

def Image():
    M=3
    f0=0
    f1=0
    R=M-(f0+f1)
    k=0.5
    Tf=1/(R-1)*((M-1)/(3-k)-f1)
    fname=str('ideal_condition.pdf')
    print(Tf)
    mu=0.5
    sd=(1/(12*(R-1)))**0.5
    x=np.arange(0,1,0.01)
    y=norm.pdf(x,mu,sd) #pdf (mean sd)
    y_cdf=norm.cdf(x,mu,sd)#cdf
    plt.plot(x,y,'r')
    plt.show()
    plt.plot(x,y_cdf,'b')
    plt.show()
    generous_range=np.arange(0,Tf+0.01,0.01)
    y1=norm.pdf(generous_range,mu,sd)
    selfish_range=np.arange(Tf,1,0.01)
    y2=norm.pdf(selfish_range,mu,sd)
    generous_prob=norm.cdf(Tf,mu,sd)
    
    plt.fill_between(generous_range,y1,np.zeros(np.size(generous_range)),facecolor='r',alpha=0.5)
    plt.fill_between(selfish_range,y2,np.zeros(np.size(selfish_range)),facecolor='b',alpha=0.5)
    plt.xlim(0,1)
    plt.xlabel('$x$',fontsize=16)
    plt.text(0.2,0.5,'generous' ,fontsize=16)
    plt.text(0.8,0.5,'selfish' ,fontsize=16)
    #plt.savefig(fname)
    plt.show()
    
def Analyze (M,k):
    fav=np.zeros(M)#favorability
    #M number of species in the community
    # k degree of conflict
    
    #inital condition
    comb=1
    f0=np.zeros(comb)
    f1=np.zeros(comb)
    prob=np.ones(comb)
    gen_prob=np.zeros(comb)
    #l.h.s.'s distribution is approximated by normal distiribution whose mean and sd are
    
    
    for i in range(comb):
        
        R=M-(f0[i]+f1[i])
        Tf=((M-1)/(3-k)-f1[i])/(R-1)
        mu=0.5
        sd=(1/(12*(R-1)))**0.5   
        #g_range=np.arange(0,Tf+0.01,0.01)
        gen_prob[i]=norm.cdf(Tf,mu,sd)
        fav[0]+=(1-gen_prob[i])*prob[i]
    #print(fav[0])
    #calculate sequentially
    if M>3:
        for m in range(1, M-2):
            #print(m)
            # We do not have to calculate the slowest species
            #set the condition
            comb=2*comb
            prob=np.concatenate((prob,prob),axis=0)
            gen_prob=np.concatenate((gen_prob,gen_prob),axis=0)
            f0=np.concatenate((f0,f0),axis=0)
            f1=np.concatenate((f1,f1),axis=0)
            for j in range (comb):
                if j<comb/2:#become selfish
                    f0[j]+=1
                    prob[j]=prob[j]*(1-gen_prob[j])
                else:#become generous
                    f1[j]+=1
                    prob[j]=prob[j]*gen_prob[j]
            
            gen_prob=np.zeros(comb)
            for i in range(comb): 
                R=M-(f0[i]+f1[i])
                Tf=((M-1)/(3-k)-f1[i])/(R-1)
                if Tf<0:
                    gen_prob[i]=0
                elif Tf>1:
                    gen_prob[i]=1
                else:
                    mu=0.5
                    sd=(1/(12*(R-1)))**0.5   
                    gen_prob[i]=norm.cdf(Tf,mu,sd)
                fav[m]+=(1-gen_prob[i])*prob[i]
            print(fav[m])
    #calculate the 2nd slowest
    #notice that we can use uniform distribution in this case
    comb=2*comb
    prob=np.concatenate((prob,prob),axis=0)
    gen_prob=np.concatenate((gen_prob,gen_prob),axis=0)
    f0=np.concatenate((f0,f0),axis=0)
    f1=np.concatenate((f1,f1),axis=0)
    for j in range (comb):
        if j<comb/2:#become selfish
            f0[j]+=1
            prob[j]=prob[j]*(1-gen_prob[j])
        else:#become generous
            f1[j]+=1
            prob[j]=prob[j]*gen_prob[j]
        
    gen_prob=np.zeros(comb)
    for i in range(comb): 
        R=M-(f0[i]+f1[i])
        Tf=((M-1)/(3-k)-f1[i])/(R-1)
        gen_prob[i]=uniform.cdf(Tf,0,1)#using uniform
        fav[M-2]+=(1-gen_prob[i])*prob[i]
    
    #the slowest species
    comb=2*comb
    prob=np.concatenate((prob,prob),axis=0)
    gen_prob=np.concatenate((gen_prob,gen_prob),axis=0)
    f0=np.concatenate((f0,f0),axis=0)
    f1=np.concatenate((f1,f1),axis=0)
    for j in range (comb):
        if j<comb/2:#become selfish
            f0[j]+=1
            prob[j]=prob[j]*(1-gen_prob[j])
        else:#become generous
            f1[j]+=1
            prob[j]=prob[j]*gen_prob[j]  
    total_gen=math.floor((M-k+2)/(3-k))
    for i in range (comb):
        if f1[i]>=total_gen:
            fav[M-1]+=prob[i]        
    return(fav)

def Analyze2():
    # test code for calculating the favorability in the three and four species community
    case=4
    K=np.array([0.5,1.5,0.5,1.6])
    M=np.array([3,  3,  4,  4])
    result3=np.zeros([2,3])
    result4=np.zeros([2,4])
    for i in range(case):
        if i<2:
            result3[i,:]=Analyze(M[i],K[i])
        else:
            result4[i-2,:]=Analyze(M[i],K[i])
    #plot the figure
    plt.figure(figsize=(5,7))
    three_species=np.linspace(0,2,3)
    plt.subplot(2,1,1)
    plt.plot(three_species,result3[0,::-1], color="c", marker = "o", markersize=10, label="small $k, (M=3)$")
    plt.plot(three_species,result3[1,::-1], color="b",marker="D",markersize=10, label="large $k, (M=3)$")
    plt.xticks([],[])
    plt.ylim(0.0,1)
    plt.yticks([0,0.5,1],[0,0.5,1], fontsize=12)
    plt.ylabel("favorabililty",fontsize=16)
    plt.legend(loc='lower center',ncol=2,fontsize=12)
    
    four_species=np.linspace(0,3,4)
    plt.subplot(2,1,2)
    plt.plot(four_species,result4[0,::-1],color="c",marker="o",markersize=10, label="small $k, (M=4)$")
    plt.plot(four_species,result4[1,::-1],color="b",marker="D",markersize=10,label="large $k, (M=4)$")
    plt.xticks([0,3],["slowest","fastest"],fontsize=12)
    plt.xlabel("evolutionary rate",fontsize=12)
    plt.ylim(0.0,1)
    plt.yticks([0,0.5,1],[0,0.5,1],fontsize=12)
    #plt.ylabel("favorable",fontsize=14)
    plt.legend(loc='upper center',ncol=2,fontsize=12)
    plt.tight_layout()
    plt.savefig("ideal_favorability.pdf")
    plt.show()
    print(result3)
    print(result4)

        
        
    
