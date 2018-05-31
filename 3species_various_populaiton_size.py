#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 11:51:45 2018

@author: shibasakishota
In this code, there exist 3 spiceis with populaiton sizes
"""
from scipy.integrate import ode
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#fixed paraemters
M=3 # number of species
T=25000 # end time
t0=0 # init time
threshold =10**(-3)#thresholdwhere selfish is fixed

def func(t, x, r, k, M, W):
    """
    replicator dynamics in mutualistic communities
    t time
    x frequency of generous in each species, vector
    r evolutionary rates of each species , vector
    k payoff when the two palyers are generous
    M number of species
    W wieght matrix, diagonal elements are zeros
    """
    weight=np.dot(W,x)#weighted effect from other species vector
    f=[r*x*(1-x)*(1+(k-3)*weight)]
    return f

def main (model):
    #deceide the model anf values of paraemters
    if 1==model:
        #model1 baseline
        k=0.5
        N=np.array([3,3,4])#relative population size
        r=np.array([1/8,1,8])#evolutionary rate
        fname="RKBQ7_model1"
    elif 2==model:
        #change the evolutionary rates
        k=0.5
        N=np.array([3,3,4])#relative population size
        r=np.array([8,1,1/8])#evolutionary rate
        fname="RKBQ7_model2"
    elif 3==model:
        #change the value of k
        k=1.5
        N=np.array([3,3,4])#relative population size
        r=np.array([1/8,1,8])#evolutionary rate
        fname="RKBQ7_model3"
        
    elif 4 == model:
        #change the balance of reralize population size
        k=0.5
        N=np.array([1,1,5])#relative population size
        r=np.array([1/8,1,8])#evolutionary rate
        fname="RKBQ7_model4"
    else:
        print("Error in model veriation.")
        return 1
    #prodcing wieght matrix W given N
    W=np.zeros([M,M])
    for i in range(M):
        N_i=sum(N)-N[i]
        for j in range(M):
            if i==j:
                #diagonal
                W[i,j]=0
            else:
                W[i,j]=N[j]/N_i
                
    #give the initial condition
    init_low=0.01#lowest freq of generous in the initial condition
    init_high=0.99#highest freq
    init_step=0.05
    size=(init_high-init_low)//init_step+1
    x0=np.ones(M)*init_low
    #set init values 
    counter=0
    counter_end=int (size**M)
    result=np.zeros([counter_end,M*3+2])#the result listss
    """
    in result
    init cond x3, final condition x3, cond after first fix x3, 1st fix,equilibrium
    therefore counter_end x 10 matrix
    """
    while counter<counter_end:
        solver=ode(func)#ordinary differential equations
        solver.set_integrator('dopri5')#set integrator
        solver.set_f_params(r, k, M, W) #set parameter
        solver.set_initial_value(x0, t0) # set initial conditions freq and time
        t=np.linspace(t0, T, T*10)
        sol=np.empty((T*10,M))
        sol[0]=x0
        l=1
        #calculate the dynamics
        flag=0#flag reset
        while solver.successful() and solver.t < T:
            solver.integrate(t[l])
            sol[l] = solver.y
            for i in range(M):
                #species i fixes the strategy either generous or selfish
                if 0==flag:
                #if flag !=0, the fixation condition is already saved
                    if sol[l,i]<threshold or sol[l,i]>1-threshold:
                        fixation=np.concatenate((sol[l,:],np.array([i+1])),axis=0)
                        #print(fixation)
                        flag=1
                        break
            l += 1
        #check the equilibrium type
        equ_type=-1#reset
        if sol[-1,0]<threshold and sol[-1,1]<threshold and sol[-1,2]<threshold:
            #type 0 (0,0,0)
            equ_type=0
        elif sol[-1,0]<threshold and sol[-1,1]<threshold and sol[-1,2]>1-threshold:
            #type 1 (0,0,1)
            equ_type=1
        elif sol[-1,0]<threshold and sol[-1,1]>1-threshold and sol[-1,2]<threshold:
            #type 2 (0,1,0)
            equ_type=2
        elif sol[-1,0]>1-threshold and sol[-1,1]<threshold and sol[-1,2]<threshold:
            #type 3 (1,0,0)
            equ_type=3
        elif sol[-1,0]<threshold and sol[-1,1]>1-threshold and sol[-1,2]>1-threshold:
            #type 4 (0,1,1)
            equ_type=4
        elif sol[-1,0]>1-threshold and sol[-1,1]<threshold and sol[-1,2]>1-threshold:
            #type 5 (1,0,1)
            equ_type=5
        elif sol[-1,0]>1-threshold and sol[-1,1]>1-threshold and sol[-1,2]<threshold:
            #type 6 (1,1,0)
            equ_type=6
        elif sol[-1,0]>1-threshold and sol[-1,1]>1-threshold and sol[-1,2]>1-threshold:
            #type 7 (1,1,1)
            equ_type=7
        #save the result file
        result[counter,:]=np.hstack((x0,sol[-1,:],fixation,equ_type))
        #change the initial condition
        for i in range (M-1,-1,-1):
            if x0[i]+init_step>init_high:
                x0[i]=init_low
            else:
                x0[i]+=init_step
                break # end changin the init condition
        counter+=1
        print(counter)
        
    #end while loop of all patterns of initial conditions
    #save the data and/or plot the figure
    #output the csv file
    head='init1, init2, init3, final1,final2, final3, ff1,ff2,ff3,ff, equilibrium'
    np.savetxt(fname+".csv",result,delimiter=',',fmt="%.4f",comments="",header=head)
    
    #plot the initial conditions and the fate of dynamics
    fig = plt.figure(figsize=(10,7.5))
    plt.rcParams["font.size"]=10
    ax = Axes3D(fig)
    #ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel("Species 1")
    ax.set_ylabel("Species 2")
    ax.set_zlabel("Species 3")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_zlim(0, 1)
    ax.invert_xaxis()
    number_eq=np.zeros([1,8])
    equib0=result[result[:,-1]==0]
    number_eq[0,0]=np.size(equib0,0)/counter_end
    equib1=result[result[:,-1]==1]
    number_eq[0,1]=np.size(equib1,0)/counter_end
    equib2=result[result[:,-1]==2]
    number_eq[0,2]=np.size(equib2,0)/counter_end
    equib3=result[result[:,-1]==3]
    number_eq[0,3]=np.size(equib3,0)/counter_end
    equib4=result[result[:,-1]==4]
    number_eq[0,4]=np.size(equib4,0)/counter_end
    equib5=result[result[:,-1]==5]
    number_eq[0,5]=np.size(equib5,0)/counter_end
    equib6=result[result[:,-1]==6]
    number_eq[0,6]=np.size(equib6,0)/counter_end
    equib7=result[result[:,-1]==7]
    number_eq[0,7]=np.size(equib7,0)/counter_end
    ax.plot(equib0[:,0], equib0[:,1], equib0[:,2],"o", color="k",\
            markeredgecolor='k',label='$(0,0,0)$')
    ax.plot(equib1[:,0], equib1[:,1], equib1[:,2],"o", color="b",\
            markeredgecolor='k',label='$(0,0,1)$')
    ax.plot(equib2[:,0], equib2[:,1], equib2[:,2],"o", color="g",\
            markeredgecolor='k',label='$(0,1,0)$')
    ax.plot(equib3[:,0], equib3[:,1], equib3[:,2],"o", color="r",\
            markeredgecolor='k',label='$(1,0,0)$')
    ax.plot(equib4[:,0], equib4[:,1], equib4[:,2],"o", color="c",\
            markeredgecolor='k',label='$(0,1,1)$')
    ax.plot(equib5[:,0], equib5[:,1], equib5[:,2],"o", color="m",\
            markeredgecolor='k',label='$(1,0,1)$')
    ax.plot(equib6[:,0], equib6[:,1], equib6[:,2],"o", color="y",\
            markeredgecolor='k',label='$(1,1,0)$')
    ax.plot(equib7[:,0], equib7[:,1], equib7[:,2],"o", color="w",\
            markeredgecolor='k',label='$(1,1,1)$')
    ax.legend(loc="upper left")   
    plt.savefig(fname+"_Basin_Atractor_RKBQ7.pdf",transparent=True)
    fig.show()  
    
    #show the volume of each steady state
    head='000,001,010,100,011,101,110,111'
    np.savetxt(fname+"_number_equilibrium.csv",number_eq,delimiter=',',fmt="%.4f",comments="",header=head)
    
    #show the table of selfish of each species
    selfish=np.zeros([1,3])
    selfish[0,0]=number_eq[0,0]+number_eq[0,1]+number_eq[0,2]+number_eq[0,4]
    selfish[0,1]=number_eq[0,0]+number_eq[0,1]+number_eq[0,3]+number_eq[0,5]
    selfish[0,2]=number_eq[0,0]+number_eq[0,2]+number_eq[0,3]+number_eq[0,6]
    head='selfish Sp1, selfish Sp2, selfish Sp3'
    np.savetxt(fname+"_number_selfish.csv",selfish,delimiter=',',fmt="%.4f",comments="",header=head)
    
    
    #plot the conditions of first fix and the fate
    
    ff1=result[result[:,3*M]==1] #first fix is specie1
    ff2=result[result[:,3*M]==2]
    ff3=result[result[:,3*M]==3]
    ff_size=np.zeros([1,M])
    ff_size[0,0]=np.size(ff1,0)/counter_end
    ff_size[0,1]=np.size(ff2,0)/counter_end
    ff_size[0,2]=np.size(ff3,0)/counter_end
    head = "species 1, species 2, speceis 3"
    np.savetxt(fname+"_FirstFix.csv",ff_size,delimiter=',',fmt="%.4f",comments="",header=head)
    
    #plot when speceis 1 is FF
    if(len(ff1)!=0):
        fig = plt.figure()
        equib0=ff1[ff1[:,-1]==0]
        equib1=ff1[ff1[:,-1]==1]
        equib2=ff1[ff1[:,-1]==2]
        equib3=ff1[ff1[:,-1]==3]
        equib4=ff1[ff1[:,-1]==4]
        equib5=ff1[ff1[:,-1]==5]
        equib6=ff1[ff1[:,-1]==6]
        equib7=ff1[ff1[:,-1]==7]      
        ax = fig.add_subplot(1,1,1)
        ax.scatter(equib0[i:,2*M+1], equib0[:,2*M+2],marker="o", color="k",\
                edgecolors='k',label='$(0,0,0)$')
        ax.scatter(equib1[:,2*M+1], equib1[:,2*M+2],marker="o", color="b",\
            edgecolors='k',label='$(0,0,1)$')
        ax.scatter(equib2[:,2*M+1], equib2[:,2*M+2],marker="o", color="g",\
            edgecolors='k',label='$(0,1,0)$')
        ax.scatter(equib3[:,2*M+1], equib3[:,2*M+2],marker="o", color="r",\
            edgecolors='k',label='$(1,0,0)$')
        ax.scatter(equib4[:,2*M+1], equib4[:,2*M+2],marker="o", color="c",\
            edgecolors='k',label='$(0,1,1)$')
        ax.scatter(equib5[:,2*M+1], equib5[:,2*M+2],marker="o", color="m",\
            edgecolors='k',label='$(1,0,1)$')
        ax.scatter(equib6[:,2*M+1], equib6[:,2*M+2],marker="o", color="y",\
            edgecolors='k',label='$(1,1,0)$')
        ax.scatter(equib7[:,2*M+1], equib7[:,2*M+2],marker="o", color="w",\
            edgecolors='k',label='$(1,1,1)$')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)   
        if model !=2:
            ax.set_xlabel('generous in normal', fontsize=14)
            ax.set_ylabel('generous in fast',fontsize=14)
        else:
            ax.set_xlabel('generous in normal', fontsize=14)
            ax.set_ylabel('generous in slow',fontsize=14)
        plt.xlim(0,1)
        plt.ylim(0,1)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.savefig(fname+"_FirstFix_sp1.pdf",transparent=True)
        fig.show()
    #sp2 is FF
    if(len(ff2)!=0):
        fig = plt.figure()
        equib0=ff2[ff2[:,-1]==0]
        equib1=ff2[ff2[:,-1]==1]
        equib2=ff2[ff2[:,-1]==2]
        equib3=ff2[ff2[:,-1]==3]
        equib4=ff2[ff2[:,-1]==4]
        equib5=ff2[ff2[:,-1]==5]
        equib6=ff2[ff2[:,-1]==6]
        equib7=ff2[ff2[:,-1]==7]      
        ax = fig.add_subplot(1,1,1)
        ax.scatter(equib0[:,2*M+0], equib0[:,2*M+2],marker="o", color="k",\
                   edgecolors='k',label='$(0,0,0)$')
        ax.scatter(equib1[:,2*M+0], equib1[:,2*M+2],marker="o", color="b",\
            edgecolors='k',label='$(0,0,1)$')
        ax.scatter(equib2[:,2*M+0], equib2[:,2*M+2],marker="o", color="g",\
            edgecolors='k',label='$(0,1,0)$')
        ax.scatter(equib3[:,2*M+0], equib3[:,2*M+2],marker="o", color="r",\
            edgecolors='k',label='$(1,0,0)$')
        ax.scatter(equib4[:,2*M+0], equib4[:,2*M+2],marker="o", color="c",\
            edgecolors='k',label='$(0,1,1)$')
        ax.scatter(equib5[:,2*M+0], equib5[:,2*M+2],marker="o", color="m",\
                   edgecolors='k',label='$(1,0,1)$')
        ax.scatter(equib6[:,2*M+0], equib6[:,2*M+2],marker="o", color="y",\
            edgecolors='k',label='$(1,1,0)$')
        ax.scatter(equib7[:,2*M+0], equib7[:,2*M+2],marker="o", color="w",\
            edgecolors='k',label='$(1,1,1)$')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)   
        if model !=2:
            ax.set_xlabel('generous in slow', fontsize=14)
            ax.set_ylabel('generous in fast',fontsize=14)
        else:
            ax.set_xlabel('generous in fast', fontsize=14)
            ax.set_ylabel('generous in slow',fontsize=14)
        plt.xlim(0,1)
        plt.ylim(0,1)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.savefig(fname+"_FirstFix_sp2.pdf",transparent=True)
        fig.show()
        
        #sp3 is FF
    if(len(ff3)!=0):
        fig = plt.figure()
        equib0=ff3[ff3[:,-1]==0]
        equib1=ff3[ff3[:,-1]==1]
        equib2=ff3[ff3[:,-1]==2]
        equib3=ff3[ff3[:,-1]==3]
        equib4=ff3[ff3[:,-1]==4]
        equib5=ff3[ff3[:,-1]==5]
        equib6=ff3[ff3[:,-1]==6]
        equib7=ff3[ff3[:,-1]==7]      
        ax = fig.add_subplot(1,1,1)
        ax.scatter(equib0[:,2*M+0], equib0[:,2*M+1],marker="o", color="k",\
            edgecolors='k',label='$(0,0,0)$')
        ax.scatter(equib1[:,2*M+0], equib1[:,2*M+1],marker="o", color="b",\
            edgecolors='k',label='$(0,0,1)$')
        ax.scatter(equib2[:,2*M+0], equib2[:,2*M+1],marker="o", color="g",\
            edgecolors='k',label='$(0,1,0)$')
        ax.scatter(equib3[:,2*M+0], equib3[:,2*M+1],marker="o", color="r",\
            edgecolors='k',label='$(1,0,0)$')
        ax.scatter(equib4[:,2*M+0], equib4[:,2*M+1],marker="o", color="c",\
            edgecolors='k',label='$(0,1,1)$')
        ax.scatter(equib5[:,2*M+0], equib5[:,2*M+1],marker="o", color="m",\
            edgecolors='k',label='$(1,0,1)$')
        ax.scatter(equib6[:,2*M+0], equib6[:,2*M+1],marker="o", color="y",\
            edgecolors='k',label='$(1,1,0)$')
        ax.scatter(equib7[:,2*M+0], equib7[:,2*M+1],marker="o", color="w",\
            edgecolors='k',label='$(1,1,1)$')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)   
        if model !=2:
            ax.set_xlabel('generous in slow', fontsize=14)
            ax.set_ylabel('generous in normal',fontsize=14)
        else:
            ax.set_xlabel('generous in fast', fontsize=14)
            ax.set_ylabel('generous in normal',fontsize=14)
        plt.xlim(0,1)
        plt.ylim(0,1)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.savefig(fname+"_FirstFix_sp3.pdf",transparent=True)
        fig.show()
    