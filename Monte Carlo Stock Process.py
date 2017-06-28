# -*- coding: utf-8 -*-
"""
Created on Fri Jan 16 02:41:55 2015

@author: stevengoldschmidt
"""

import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import csv
from scipy.stats import poisson

def monte_carlo_stock_process():
    print()
    print("Monte Carlo Stock Process and Option Valuation")
    print()
    
    n = 500 #time steps
    Trials = 20000 #simulations per run
    
    print("Time steps = %s" %n)
    print("Simulations = %s" %Trials)
    print()
    
    #Black-Scholes Pricing Inputs
    S_0=118.79
    K=120.00
    sigma=0.4615
    r=0.0066
    T=240.0/365.0
    h = T/n
    divYield = 0.0000
    
    print("S_0 = %s" %S_0)
    print("K = %s" %K)
    print("sigma = %s" %sigma)
    print("r = %s" %r)
    print("T = %s" %T)
    print("divYield = %s" %divYield)
    print()
        
    #nname='Monte_Carlo'
    #name=nname+'.csv'

    #file = open(name,'w',newline='')
    #writer = csv.writer(file, quoting=csv.QUOTE_ALL)
    #writer.writerow(['Monte Carlo Simulation']) #file header

    Expiration_Assets = np.zeros(Trials,dtype=np.float32)
    Call_Exercises = np.zeros(Trials,dtype=np.float32)
    Put_Exercises = np.zeros(Trials,dtype=np.float32)
    
    trial=1
   
    while trial < Trials + 1:
    
        S=np.zeros(n+1,dtype=np.float32)
        N=np.zeros(n+1,dtype=np.float32)
    
        i=0 #time step counter
    
        for i in range(0,n+1):
            N[i] = i*(T/n)
        
        
        S[0] = S_0
    
        k=1 #time step counter
    
        for k in range(1,n+1):
            
            Drift = math.exp((r - divYield - 0.5*sigma**2)*h)
            Noise = math.exp(sigma*math.sqrt(h)*np.random.normal(0,1))
            
            S[k] = S[k-1]*Drift*Noise
            
            if k == n:
                Expiration_Assets[trial-1] = S[k]

        #file = open(name,'a',newline='') #open the file to append
        #writer = csv.writer(file, quoting=csv.QUOTE_ALL)
        #writer.writerow(['Trial', trial]) #adjusts because the counter starts at 0
    
        #for values in N,S: #write the trial's values
                    #writer.writerow(values)
    
        print(trial)
        plt.plot(N,S,'b--')
        trial+=1
        
    print()
    
    Average_Expiration_Asset_Price = np.mean(Expiration_Assets)
    
    print("Average Asset Price = %s" %Average_Expiration_Asset_Price)
    print()
    
    i=0
    
    for i in range(0,Trials):
            Call_Exercises[i] = max(Expiration_Assets[i] - K,0)
            Put_Exercises[i] = max(K - Expiration_Assets[i],0)
    
    
    Call_Price = math.exp(-r*T)*np.mean(Call_Exercises)
    Put_Price = math.exp(-r*T)*np.mean(Put_Exercises)
    
    Call_Price_STDev = np.std(Call_Exercises)*(1/math.sqrt(Trials))
    Put_Price_STDev = np.std(Put_Exercises)*(1/math.sqrt(Trials))
    
    print("Call Price = %s" %Call_Price)
    print("Call Price STDev = %s" %Call_Price_STDev)
    print()
    
    print("Put Price = %s" %Put_Price)
    print("Put Price STDev = %s" %Put_Price_STDev)
    print()
    

# main program starts here
monte_carlo_stock_process()