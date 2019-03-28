#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 13:50:56 2019

@author: shibasakishota
"""

import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt

def func(t, x, r, k, M):
    SumX=sum(x)
    f=[r*x*(1-x)*(1+(k-3)*(SumX-x)/(M-1))]
    return f

def MC_ODE(M):
    if M==4:
        k_array=np.array([0.5, 1.6])
    else:
        k_array=np.array([0.5, 1.5])
    
    #M nunmber of species
    
    #r_array array for evolutionary rates
    num_sample=1000 #number of sample in eacnreplicate
    num_replicate=10 #number of replicate
    CSV=np.zeros([num_sample,3*M])#save 
    """
    When M=4,
    CSV= init             results     1st fixed 2nd fixed, 3rd fixed, 4th fixed
    0.1, 0.2, 0.3, 0.4 | 1, 0, 0, 1 |  1           2           3          4
    """
    t0=0#init time
    T=5000#end time
    t=np.linspace(t0, T, T*10)
    threshold =10**(-3)#thresholdwhere selfish is fixed
    if M==3:
        r=np.array([1/8,1,8])#evolutionary rate for 3 species model
    elif M==4:
        r=np.array([1/8,1/2,2,8])#evolutionary rate for 4 species model
    else:
        r=np.zeros([M])
        for i in range (M):
            r[i]=2**(M-1)
    for a in range(np.size(k_array)):
        k=k_array[a] #k is a parameter for payoff
        print(str('k=%.1f'%(k)))
        for rep in range(num_replicate):
            print(str('replicate %d' %(rep)))
            for i in range (num_sample):
                #initial condition is given by  Uniform dist
                init=np.random.uniform(0,1,M)
                #flags for fixation
                flag=0
                list_already_fixed=[]
                
                #solve ODE
                solver=ode(func)
                solver.set_integrator('dopri5')#set integrator
                solver.set_f_params(r, k, M) #set parameter
                solver.set_initial_value(init, t0) # set initial conditions freq and time
                sol=np.empty((T*10,M))
                sol[0]=init
                l=1
                fixation_list=np.zeros([M])
                #calculate the dynamics
                while solver.successful() and solver.t < T:
                    solver.integrate(t[l])
                    sol[l] = solver.y
                    if M-1>flag:
                        #if flag >=M-1, the all speceis has already fixed the strategy
                        for n in range(M):
                            if str(n) not in list_already_fixed:
                                #species n fixes the strategy either generous or selfish in flagth order
                                if sol[l,n]<threshold or sol[l,n]>1-threshold:
                                    list_already_fixed.append(str(n))
                                    fixation_list[flag]=n
                                    flag+=1#next fixation
                    l += 1
                #x=1 generous, x=0 selfish
                CSV[i,:]=np.concatenate([init, sol[-1,:], fixation_list])#save the init and results
            #write csv file 
            fname=str('MC-results-k%.1f-M%d-replicate%d.csv' %(k,M, rep))
            np.savetxt(fname,CSV,delimiter=',',fmt="%.5f",comments="")    
    #plot        
    Fav_Plot(num_replicate,num_sample,M)
    
def Fav_Plot(num_replicate, num_sample, M):
    if M==4:
        k_array=np.array([0.5, 1.6])
    else:
        k_array=np.array([0.5, 1.5])
    for m in range (np.size(k_array)):
        k=k_array[m]
        Fav=np.zeros([num_replicate, M])
        for n in range(num_replicate):
            fname=str('MC-results-k%.1f-M%d-replicate%d.csv' %(k,M, n))
            D=np.loadtxt(fname, delimiter=',')
            for i in range(num_sample):
                for j in range (M):
                    if D[i, j+M]<10**(-3):
                        Fav[n, j]+=1
            Fav[n,:]=Fav[n,:]/num_sample#normalize
        Mean_Fav=np.mean(Fav,0)
        SE_Fav=np.std(Fav,0)/(num_replicate)**0.5
        X=np.linspace(1,M,M)
        if k<1:    
              plt.plot(X, Mean_Fav,color='c', marker='o',markersize=10,label="small $k$")
              plt.errorbar(X, Mean_Fav, SE_Fav,  color='c')
                    
        else:
            plt.plot(X, Mean_Fav,color='b', marker='D',markersize=10,label="large $k$ ")
            plt.errorbar(X,Mean_Fav, SE_Fav,  color='b')
    plt.ylim(0,0.8)
    #plt.xticks([1, M], ['', ''], fontsize=14)
    plt.xticks([1, M], ['slowest', 'fastest'], fontsize=14)
    plt.yticks([0.0, 0.4, 0.8], ['0.0', '0.4', '0.8'], fontsize=14)
    gname=str('MCfav-%dspecies.pdf' %(M))
    plt.legend(loc='lower center', ncol=2, fontsize=14)
    plt.savefig(gname)
    plt.show()
    
                
def Figure1(num_replicate,num_sample):
    #plot figure 1
    #image 
    plt.figure(figsize=(6,9))
    plt.subplot(311)
    x=np.linspace(0,0.8,801)
    plt.plot(x, -x+0.8, color='r', linestyle='--', label='Red King')
    plt.plot(x, x, color='r', linestyle='-', label='Red Queen')
    plt.legend(loc='lower center', ncol=2, fontsize=14)
    plt.xticks([1,], [ ''], fontsize=14)
    plt.yticks([0.0, 0.4, 0.8], ['0.0', '0.4', '0.8'], fontsize=14)
    M_array=np.array([3,4])
    for m in range (np.size(M_array)):
        M=M_array[m]
        if M==4:
            k_array=np.array([0.5, 1.6])
            plt.subplot(313)
        else:
            k_array=np.array([0.5, 1.5])
            plt.subplot(312)
        Fav=np.zeros([num_replicate, M])
        for l in range (np.size(k_array)):
            k=k_array[l]
            for n in range(num_replicate):
                fname=str('MC-results-k%.1f-M%d-replicate%d.csv' %(k,M, n))
                D=np.loadtxt(fname, delimiter=',')
                for i in range(num_sample):
                    for j in range (M):
                        if D[i, j+M]<10**(-3):
                            Fav[n, j]+=1
                Fav[n,:]=Fav[n,:]/num_sample#normalize
            Mean_Fav=np.mean(Fav,0)
            SE_Fav=np.std(Fav,0)/(num_replicate)**0.5
            X=np.linspace(1,M,M)
            if k<1:    
              plt.plot(X, Mean_Fav,color='c', marker='o',markersize=10,label="small $k$")
              plt.errorbar(X, Mean_Fav, SE_Fav,  color='c')
                    
            else:
                plt.plot(X, Mean_Fav,color='b', marker='D',markersize=10,label="large $k$ ")
                plt.errorbar(X,Mean_Fav, SE_Fav,  color='b')
        plt.ylim(0,0.8)
        if M==3:
            plt.xticks([1, M], ['', ''], fontsize=14)
            plt.ylabel('favorability', fontsize=16)
        else:
            plt.xticks([1, M], ['slowest', 'fastest'], fontsize=14)
            plt.xlabel('evolutionary rate', fontsize=16)
        plt.yticks([0.0, 0.4, 0.8], ['0.0', '0.4', '0.8'], fontsize=14)
        plt.legend(loc='lower center', ncol=2, fontsize=14)
    plt.savefig('figure1.pdf')
    plt.show()
        
def RateFix(M):
    #analyzing the relationship between the evolutionary rate and the order of the fixation
    if M==4:
        k_array=np.array([0.5, 1.6])
    else:
        k_array=np.array([0.5, 1.5])
    
    #M nunmber of species
    num_replicate=10 #number of replicate
    for i in range (np.size(k_array)):
        k=k_array[i]
        #Data=[]#including all 
        for rep in range(num_replicate):
            fname=str('MC-results-k%.1f-M%d-replicate%d.csv' %(k,M, rep))
            #read csv
            data=np.loadtxt(fname, delimiter=',')
            D=data[:,2*M:]
            if rep==0:
                Data=D
            else:
                Data=np.concatenate([Data,D])
        #analyze the data
        Matrix=np.zeros([M,M])
        for m in range(M):#order
            for i in range(np.size(Data,0)):#sample
                for a in range (M):#species
                    if Data[i,m]> a-0.5 and Data[i,m]<a+0.5:
                        Matrix[m,a]+=1 #frequency that species a gets fixed at mth order
        for a in range(M):
            #normalizing according to the order
            Matrix[:,a]= Matrix[:,a]/sum( Matrix[:,a])
        #plot example when M=4
        fig, ax =plt.subplots()
        heatmap = ax.pcolor(Matrix, cmap=plt.cm.Blues)
        ax.invert_yaxis()
        plt.yticks(np.arange(0.5, M, 1), ['1st', '2nd', '3rd', '4th'],\
                   fontsize=14)
        plt.xticks(np.linspace(0.5, M-0.5, 2), ['slowest', 'fastest'],\
                   fontsize=14)
        fig.colorbar(heatmap,ticks=[0,0.25,0.5,0.75,1.0])
        if k<1:
            
            plt.title('small $k$',fontsize=18)
            plt.savefig('fixation-order-smallK.pdf')
        else:
            plt.title('large $k$',fontsize=18)
            plt.savefig('fixation-order-largeK.pdf')
        plt.show()
        print(Matrix)
            
                        
            
            
        
    