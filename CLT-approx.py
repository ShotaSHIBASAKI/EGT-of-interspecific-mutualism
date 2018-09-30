#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 30 09:51:14 2018
#CLT approximation
@author: shibasakishota
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm


def main ():
    M_array=np.array([2,3,4,5,10,15])#array of species number
    x=np.linspace(0,1,500) # array of average generous fraction
    num_sample=10**6 #number of sampling for calculation of pdf
    fig,axes = plt.subplots(nrows=2,ncols=3,figsize=(10,8))
    for i in range (np.size(M_array)):
        sigma=(1/(12*M_array[i]))**(1/2)
        CLT=norm.pdf(x, 0.5,sigma)
        pdf=np.zeros([num_sample])
        tex=str('M=%d'%(M_array[i]))
        for j in range(num_sample):
            pdf[j]=sum(np.random.rand(M_array[i]))/M_array[i]
        if i<3:
            r=0
            c=i
        else:
            r=1
            c=i-3
        
        axes[r,c].plot(x,CLT,color='red', label='approximation')
        axes[r,c].hist(pdf,bins=100,density=True)
        axes[r,c].text(0.7,1.8,tex, fontsize=14)
    plt.savefig('approximation.pdf')
    plt.show()
    
