#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 15:00:16 2018

@author: shibasakishota
"""

"""
Created on Wed Jan 17 15:00:42 2018
@author: shibasakishota
Investigate which species fix the strategy first
and test the condition of other species at that time

"""

import numpy as np
import scipy.misc as scm
from scipy.integrate import ode
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv

#parameter
M=3 #numper of species
k=1.5#benefit of both generous, large k
#k=0.5#small k
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

def main():
    if k<1:
        fname="k_Small_3species_infinite"
        Fav=np.zeros([1,int(scm.comb(M,2)+1)])#two speceis selfish+interipor equi
    elif k>1:
        fname="k_Large_Redfav_3species_infinite"
        Fav=np.zeros([1,int(scm.comb(M,1)+1)])#only one speceis selfish+interior equi
    else:
        return (1)#ignore k=1
    interior=np.ones(3)/(3-k)#interior equilibrium
    
    #note total number of cases is size**M
    r=np.array([1/8,1,8]) # evolutionary rate
    
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
    if k<1:
        init_list=np.zeros([counter_end,int(scm.comb(M,2)+5)])#the last col is result
    else:
        init_list=np.zeros([counter_end,int(scm.comb(M,1)+5)])
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
        
        #calculate to which equilibrium the dynamics reach
        if k<1:#small k
            if sol[-1,0]<threshold and sol[-1,1]<threshold:
                #slow and normal are selfish
                Fav[0,0]+=1
                init_list[counter,:]=np.hstack((x0,fixation,1))
            elif sol[-1,0]<threshold and sol[-1,2]<threshold:
                 #slow and fast are selfish
                Fav[0,1]+=1
                init_list[counter,:]=np.hstack((x0,fixation,2))
            elif sol[-1,1]<threshold and sol[-1,2]<threshold:
                 #normal and fast are selfish
                Fav[0,2]+=1
                init_list[counter,:]=np.hstack((x0,fixation,3))
            """
            elif sum(abs(sol[-1,0:M]-interior))<M*threshold:
                # interior equilbrium
                Fav[0,3]+=1 
                init_list[counter,:]=np.hstack((x0,4))
            """
        else: #large k
            if sol[-1,0]<threshold:
                #slower is selfish
                Fav[0,0]+=1
                init_list[counter,:]=np.hstack((x0,fixation,1))
            elif sol[-1,1]<threshold:
                #normal is selfish
                Fav[0,1]+=1
                init_list[counter,:]=np.hstack((x0,fixation,2))
            elif sol[-1,2]<threshold:
                #faster is selfish
                Fav[0,2]+=1
                init_list[counter,:]=np.hstack((x0,fixation,3))
            """
            elif sum(abs(sol[-1,0:M]-interior))<M*threshold:
                # interior equilbrium
                Fav[0,3]+=1 
                init_list[counter,:]=np.hstack((x0,4))
            """
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
    
    #plot the 3D and shows the attractor
    fig = plt.figure(figsize=(10,7.5))
    plt.rcParams["font.size"]=10
    ax = Axes3D(fig)
    #ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel("Slow")
    ax.set_ylabel("Intermidiate")
    ax.set_zlabel("fast")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_zlim(0, 1)
    ax.invert_xaxis()
    #print(init_list)
    equib1=init_list[init_list[:,-1]==1]
    equib2=init_list[init_list[:,-1]==2]
    equib3=init_list[init_list[:,-1]==3]
    equib4=init_list[init_list[:,-1]==4]
    ax.plot(equib1[:,0], equib1[i:,1], equib1[i:,2],"o", color="c")
    ax.plot(equib2[:,0], equib2[i:,1], equib2[i:,2],"o", color="y")
    ax.plot(equib3[:,0], equib3[i:,1], equib3[i:,2],"o", color="m")
    ax.plot(equib4[:,0], equib4[i:,1], equib4[i:,2],"o", color="k")
    plt.savefig(fname+"_equilibrium.pdf",transparent=True)
    fig.show()  
    #save as csv file
    #initial condition and evolutionary fate
    head='init_slow,init_normal,init_fast,@fixed_slow,@fixed_normal,@fixed_fast,first_fixation,equilibrium type'
    np.savetxt(fname+"_Converge_FirstFixation.csv",init_list,delimiter=',',fmt="%.5f",comments="",header=head)
    #volume of basin of attraction
    for i in range(M+1):
        Fav[0,i]=Fav[0,i]/(counter_end)  
    if k<1:
        head='slow-normal,slow-fast,normal-fast,interior'
    else:
        head='slow,normal,fast,interior'
    #np.savetxt(fname+"_Bequilibrium",Fav,delimiter=',',fmt="%.5f",comments="",header=head)
    
if __name__ == '__main__':
    main()    