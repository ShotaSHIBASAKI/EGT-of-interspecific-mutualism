#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 10:52:40 2018
read csv and classfy the Steady stetes
@author: shibasakishota
"""

from numpy import genfromtxt
import numpy as np

def main(k):
    if k<1:
        #read csv file and save as data
        fname="k_Small_3species_infinite_Converge_FirstFixation.csv"
        data=genfromtxt(fname,delimiter=',',skip_header=1)
        tail="_smallK"
    else:
        fname="k_Large_Redfav_3species_infinite_Converge_FirstFixation.csv"
        data=genfromtxt(fname,delimiter=',',skip_header=1)
        tail="_largeK"
    
    ss=data[:,-1]#read the steady state
    #all possible ss and their number
    snumber=np.zeros(3)
    for i in range(np.size(ss)):
        if ss[i]<1.1:
            snumber[0]=snumber[0]+1
        elif ss[i]<2.1:
            snumber[1]=snumber[1]+1
        else:
            snumber[2]=snumber[2]+1
    print(snumber)
    snumber=snumber/sum(snumber)#normalize
    print(snumber)
    np.savetxt("Steady_state"+tail+".csv",snumber,delimiter=',',fmt="%.5f")
    
    
