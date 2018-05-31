#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 08:20:58 2018

@author: shibasakishota
4 speceis interspecific mutualism
"""
import numpy as np
import scipy.misc as scm
from scipy.integrate import ode
import csv

#parameter
M=4 #numper of species
T=5000 # end time
t0=0 # init time
threshold =10**(-3)#thresholdwhere selfish is fixed
def func(t, x, r, k, M):
    """
    replicator dynamics of the bi-matrix snow drift game
    given number of species M
    """
    SumX=sum(x)
    f=[r*x*(1-x)*(1+(k-3)*(SumX-x)/(M-1))]
    return f

def Save_Csv(file_name, array):
    with open(file_name, 'w') as f:
             writer=csv.writer(f)
             writer.writerows(array)

def main(state):
    if 0==state:
        fname="4species_2generous_2selfish"
        k=0.5
    elif 1==state:
        fname="4species_3generous_1selfish" 
        k=1.6
    else:
        return (1)#ignore k=1
    interior=np.ones(M)/(3-k)#interior equilibrium
    
    #note total number of cases is size**M
    r=np.array([1/8,1/2,2,8]) # evolutionary rate
    
    #x0 =np.ones(M)*0.999
    # set initial condition
    init_low=0.01
    init_high=0.99
    init_step=0.05
    size=(init_high-init_low)//init_step+1
    x0=np.ones(M)*init_low
    #set init values 
    counter=0
    counter_end=int (size**M)    
    init_list=np.zeros([counter_end,int(scm.comb(M,1))])
    fixation_list=np.zeros([counter_end,int(scm.comb(M,1))])
    direction_list=np.zeros([counter_end,int(scm.comb(M,1))])
    while counter<counter_end:
        #x0=np.empty((1,M)) 
        """
        for i in range(M):
            num=init_counter[i]
            print(num)
            #num=int(num)
            x0[0,i]=x[num]  
        """
        #print(x0)
        #solution
        solver=ode(func)
        solver.set_integrator('dopri5')#set integrator
        solver.set_f_params(r, k, M) #set parameter
        solver.set_initial_value(x0, t0) # set initial conditions freq and time
        t=np.linspace(t0, T, T*10)
        sol=np.empty((T*10,M))
        sol[0]=x0
        init_list[counter,:]=x0#save the initial condition
        l=1
        #calculate the dynamics
        flag1=0#flag one species has already fixed the strategy 
        flag2=0#flag two speceis have already fixed the strategy 
        flag3=0#flag three speceis have already fixed the strategy 
        flag4=0#flag four species have already fixed the strategy
        list_already_fixed=[]
        while solver.successful() and solver.t < T:
            solver.integrate(t[l])
            sol[l] = solver.y
            if 0==flag4:
                #if flag4 !=0, the all speceis has already fixed the strategy
                for i in range(M):
                    if str(i) not in list_already_fixed:
                        #species i fixes the strategy either generous or selfish
                        if sol[l,i]<threshold or sol[l,i]>1-threshold:
                            if sol[l,i]<threshold:
                                    direction_list[counter,i]=0#selfish
                            else:
                                direction_list[counter,i]=1#generous
                            list_already_fixed.append(str(i))
                            if 0==flag1:
                                #first fixed speceis
                                fixation_list[counter,0]=i+1
                                flag1=1
                            elif 0==flag2:
                                #second fixed species
                                fixation_list[counter,1]=i+1
                                flag2=1
                            elif 0==flag3:
                                #third fixed speceis
                                fixation_list[counter,2]=i+1
                                flag3=1
                            else:
                                fixation_list[counter,3]=i+1
                                flag4=1

            l += 1
        
        #change the initial condition
        for i in range (M-1,-1,-1):
            if x0[i]+init_step>init_high:
                x0[i]=init_low
            else:
                x0[i]+=init_step
                break # end changin the init condition
        counter+=1
        print(counter)
        
    #end while loop
    #output the csv file
    head='init_sp1,init_sp2,init_sp3,init_sp4'\
    +'1st_fix,2nd_fix,3rd_fix,4th_fix'\
    +'sp1_dir,sp2_dir,sp3_dir,sp4_dir'
    result=np.hstack((init_list,fixation_list,direction_list))
    np.savetxt(fname+".csv",result,delimiter=',',fmt="%.2f",comments="",header=head)
    
    #analyze which species are more favorable
    Fav=np.zeros([1,M])
    for i in range(M):
        Fav[0,i]=(counter-sum(direction_list[:,i]))/counter
    head="sp1,sp2,sp3,sp4"
    np.savetxt(fname+"_favorable.csv",Fav,delimiter=',',fmt="%.5f",comments="",header=head)
    
    #analyze the order of fixation and
    #first fixed species and the evolutionary direction
    first_fix=np.zeros([M,1])
    first_fix_generous=np.zeros([M,1])
    first_fix_selfish=np.zeros([M,1])
    for i in range (M):
        first_fix[i,0]=np.size((fixation_list[fixation_list[:,0]==i+1]),0)
        for j in range(counter):
            if result[j,M]==i+1 and result[j,2*M+i]==1:
                first_fix_generous[i,0]+=1
        first_fix_selfish[i,0]=first_fix[i,0]-first_fix_generous[i,0]
    first_result=np.hstack((first_fix,first_fix_generous))
    first_result=np.hstack((first_result,first_fix_selfish))
    head="total_first_fix,generous,selfish"
    np.savetxt(fname+"_first_fix.csv",first_result,delimiter=',',fmt="%d",comments="",header=head)
    #second fixed species and the evolutionary direction
    second_fix=np.zeros([M,1])
    second_fix_generous=np.zeros([M,1])
    second_fix_selfish=np.zeros([M,1])
    for i in range (M):
        second_fix[i,0]=np.size((fixation_list[fixation_list[:,1]==i+1]),0)
        for j in range(counter):
            if result[j,M+1]==i+1 and result[j,2*M+i]==1:
                second_fix_generous[i,0]+=1
        second_fix_selfish[i,0]=second_fix[i,0]-second_fix_generous[i,0]
    second_result=np.hstack((second_fix,second_fix_generous))
    second_result=np.hstack((second_result,second_fix_selfish))
    head="total_second_fix,generous,selfish"
    np.savetxt(fname+"_second_fix.csv",second_result,delimiter=',',fmt="%d",comments="",header=head)

def third_fix(m):
    m= int(m)
    M=4#number of species
    # m is the number of generous species
    if m==2 or m==3:
        print("Good value of m")
        print(m)
    else:
        print("Error in vale of m")
        return 1
    
    data_name=str("%dspecies_%dgenerous_%dselfish.csv" %(M,m,M-m))
    fname=str("%dspecies_%dgenerous_%dselfish" %(M,m,M-m))
    #read csv file
    data=np.loadtxt(data_name,delimiter=',',skiprows=1)
    #third fixed species and the evolutionary direction
    third_fix=np.zeros([M,1])
    third_fix_generous=np.zeros([M,1])
    third_fix_selfish=np.zeros([M,1])
    for i in range (M):
        third_fix[i,0]=np.size((data[data[:,M+2]==i+1]),0)
        for j in range(np.size(data,0)):
            if data[j,M+2]==i+1 and data[j,2*M+i]==1:
                third_fix_generous[i,0]+=1
        third_fix_selfish[i,0]=third_fix[i,0]-third_fix_generous[i,0]
    third_result=np.hstack((third_fix,third_fix_generous))
    third_result=np.hstack((third_result,third_fix_selfish))
    head="total_third_fix,generous,selfish"
    np.savetxt(fname+"_third_fix.csv",third_result,delimiter=',',fmt="%d",comments="",header=head)
        
        