#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 12:13:14 2018
1.read csv file
2. check which species is more likey to be fixed first
3. classfy the SS and the other two species states when one sp is fixed 
@author: shibasakishota
"""

from numpy import genfromtxt
import numpy as np
import matplotlib.pyplot as plt
import csv

#parameters
M=3#number of spe

def Save_Csv(file_name, array):
    with open(file_name, 'w') as f:
             writer=csv.writer(f)
             writer.writerows(array)

def main(k):
    int(M)
    if k<1:
        #read csv file and save as data
        fname="k_Small_3species_infinite_Converge_FirstFixation.csv"
        data=genfromtxt(fname,delimiter=',',skip_header=1)
        tail="_smallK"
    else:
        fname="k_Large__3species_infinite_Converge_FirstFixation.csv"
        data=genfromtxt(fname,delimiter=',',skip_header=1)
        tail="_largeK"
    #check which species fixes strategies first
    #1: slow, 2: normal, 3: fast...
    first_fix_data=data[:,2*M]#data of which sp fixes first (ff)
    first_fix_count=np.zeros([2, M])#counter of ff speceis and check selfish or generous
    result_fix1=np.empty((0,M+1))#the other speceis strategy and result
    result_fix2=np.empty((0,M+1))
    result_fix3=np.empty((0,M+1))
    for i in range(len(first_fix_data)):
        sp=first_fix_data[i]#speceis which fixes the st first in simulation i
        sp=int(sp)
        state=np.hstack((data[i,M:2*M],data[i,-1]))
        state=np.array([state])
         #state at one sp fixes
         #data[i,M:2*M-1] is the state of all species
         #adta[i,-1] is the steady state the dynamics converege
        if 1==sp:
            result_fix1=np.append(result_fix1, state,axis=0)
        elif 2==sp:
            result_fix2=np.append(result_fix2, state,axis=0)
        elif 3==sp:
            result_fix3=np.append(result_fix3, state,axis=0)
        else:
            print("EROOR in classification of steady states")
            return(1)
        
        sp=int(sp-1)
        #check wheter ff species is selfish or generous
        if data[i, sp+3]<0.01:
            #selfish
            first_fix_count[0,sp]+=1
        else:
            #generous
            first_fix_count[1,sp]+=1
    np.savetxt("first_fixed_species"+tail+".csv",first_fix_count,delimiter=',',\
               header = "slow,normal,fast,1st row: selfish, 2nd row: generous",fmt="%d")
        
    #plot the state
    if k<1:
        #in this case, two species are selfish
        if(len(result_fix1)!=0):
            fig = plt.figure()
            equib1=result_fix1[result_fix1[:,-1]==1]#slow and normal are selfish
            equib2=result_fix1[result_fix1[:,-1]==2]#slow and fast are selfish
            equib3=result_fix1[result_fix1[:,-1]==3]#normal and fast are selfish
            ax = fig.add_subplot(1,1,1)
            ax.scatter(equib1[:,1],equib1[:,2], color='c',label='s & n')
            ax.scatter(equib2[:,1],equib2[:,2], color='y',label='s & f')
            ax.scatter(equib3[:,1],equib3[:,2], color='m',label='n & f')
            ax.set_xlabel('generous in normal', fontsize=14)
            ax.set_ylabel('generous in fast',fontsize=14)
            #ax.set_title("The slow fixes first",fontsize=14)
            ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=3, mode="expand", borderaxespad=0,fontsize=12)
            plt.xlim(0,1)
            plt.ylim(0,1)
            plt.gca().set_aspect('equal', adjustable='box')
            plt.savefig("Result_slow"+tail+".pdf",transparent=True)
            fig.show()
    
        if (len(result_fix2)!=0):
            fig = plt.figure()
            equib1=result_fix2[result_fix2[:,-1]==1]#slow and normal are selfish
            equib2=result_fix2[result_fix2[:,-1]==2]#slow and fast are selfish
            equib3=result_fix2[result_fix2[:,-1]==3]#normal and fast are selfish
            ax = fig.add_subplot(1,1,1)
            ax.scatter(equib1[:,0],equib1[:,2], color='c',label='s & n')
            ax.scatter(equib2[:,0],equib2[:,2], color='y',label='s & f')
            ax.scatter(equib3[:,0],equib3[:,2], color='m',label='n & f')
            ax.set_xlabel('generous in slow', fontsize=14)
            ax.set_ylabel('generous in fast',fontsize=14)
            #ax.set_title("The slow fixes first",fontsize=14)
            ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=3, mode="expand", borderaxespad=0,fontsize=12)
            plt.xlim(0,1)
            plt.ylim(0,1)
            plt.gca().set_aspect('equal', adjustable='box')
            plt.savefig("Result_normal"+tail+".pdf",transparent=True)
            fig.show()
            
        if (len(result_fix3)!=0):
            fig = plt.figure()
            equib1=result_fix3[result_fix3[:,-1]==1]#slow and normal are selfish
            equib2=result_fix3[result_fix3[:,-1]==2]#slow and fast are selfish
            equib3=result_fix3[result_fix3[:,-1]==3]#normal and fast are selfish
            ax = fig.add_subplot(1,1,1)
            ax.scatter(equib1[:,0],equib1[:,1], color='c',label='s & n')
            ax.scatter(equib2[:,0],equib2[:,1], color='y',label='s & f')
            ax.scatter(equib3[:,0],equib3[:,1], color='m',label='n & f')
            ax.set_xlabel('generous in slow', fontsize=14)
            ax.set_ylabel('generous in normal',fontsize=14)
            #ax.set_title("The slow fixes first",fontsize=14)
            ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=3, mode="expand", borderaxespad=0,fontsize=12)
            plt.xlim(0,1)
            plt.ylim(0,1)
            plt.gca().set_aspect('equal', adjustable='box')
            plt.savefig("Result_faster"+tail+".pdf",transparent=True)
            fig.show()
    else:
        #in this case,  only one speceis is selfish
        if(len(result_fix1)!=0):
            fig = plt.figure()
            equib1=result_fix1[result_fix1[:,-1]==1]#slow is selfish
            equib2=result_fix1[result_fix1[:,-1]==2]#normal is selfish
            equib3=result_fix1[result_fix1[:,-1]==3]#fast is selfish
            ax = fig.add_subplot(1,1,1)
            ax.scatter(equib1[:,1],equib1[:,2], color='c',label='s')
            ax.scatter(equib2[:,1],equib2[:,2], color='y',label='n')
            ax.scatter(equib3[:,1],equib3[:,2], color='m',label='f')
            ax.set_xlabel('generous in normal', fontsize=14)
            ax.set_ylabel('generous in fast',fontsize=14)
            #ax.set_title("The slow fixes first",fontsize=14)
            ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=3, mode="expand", borderaxespad=0,fontsize=12)
            plt.xlim(0,1)
            plt.ylim(0,1)
            plt.gca().set_aspect('equal', adjustable='box')
            plt.savefig("Result_slow"+tail+".pdf",transparent=True)
            fig.show()
    
        if (len(result_fix2)!=0):
            fig = plt.figure()
            equib1=result_fix2[result_fix2[:,-1]==1]#slow isselfish
            equib2=result_fix2[result_fix2[:,-1]==2]#normal selfish
            equib3=result_fix2[result_fix2[:,-1]==3]#fast is selfish
            ax = fig.add_subplot(1,1,1)
            ax.scatter(equib1[:,0],equib1[:,2], color='c',label='s')
            ax.scatter(equib2[:,0],equib2[:,2], color='y',label='n')
            ax.scatter(equib3[:,0],equib3[:,2], color='m',label='f')
            ax.set_xlabel('generous in slow', fontsize=14)
            ax.set_ylabel('generous in fast',fontsize=14)
            #ax.set_title("The slow fixes first",fontsize=14)
            ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=3, mode="expand", borderaxespad=0,fontsize=12)
            plt.xlim(0,1)
            plt.ylim(0,1)
            plt.gca().set_aspect('equal', adjustable='box')
            plt.savefig("Result_normal"+tail+".pdf",transparent=True)
            fig.show()
            
        if (len(result_fix3)!=0):
            fig = plt.figure()
            equib1=result_fix3[result_fix3[:,-1]==1]#slowis is selfish
            equib2=result_fix3[result_fix3[:,-1]==2]#normal is selfish
            equib3=result_fix3[result_fix3[:,-1]==3]#fast is selfish
            ax = fig.add_subplot(1,1,1)
            ax.scatter(equib1[:,0],equib1[:,1], color='c',label='s')
            ax.scatter(equib2[:,0],equib2[:,1], color='y',label='n')
            ax.scatter(equib3[:,0],equib3[:,1], color='m',label='f')
            ax.set_xlabel('generous in slow', fontsize=14)
            ax.set_ylabel('generous in normal',fontsize=14)
            #ax.set_title("The slow fixes first",fontsize=14)
            ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=3, mode="expand", borderaxespad=0,fontsize=12)
            plt.xlim(0,1)
            plt.ylim(0,1)
            plt.gca().set_aspect('equal', adjustable='box')
            plt.savefig("Result_faster"+tail+".pdf",transparent=True)
            fig.show()