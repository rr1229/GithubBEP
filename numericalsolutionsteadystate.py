# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 17:58:21 2021

@author: ronar
"""
import steadystatefunctions as ssf
import parameters as p
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm 

'number of steps:'
M=1000
steps=np.arange(0,M,1)

T_numerical=np.zeros([M,p.N])
T_wall_numerical=np.zeros([M,p.N])
V_numerical=np.zeros(M)

T_numerical[0]=p.T_steadystate0
T_wall_numerical[0]=p.Tb_steadystate0
V_numerical[0]=p.v_steadystate0

for j in tqdm(np.arange(1,M,1)):
    T_numerical[j]=ssf.temperature(V_numerical[j-1], T_numerical[j-1] , T_wall_numerical[j-1])
    T_wall_numerical[j]=ssf.Twall(V_numerical[j-1],T_numerical[j-1])
    V_numerical[j]=ssf.velocity(V_numerical[j-1], T_numerical[j-1])
    




if True:
    plt.plot(steps,V_numerical) 
    plt.title('velocity of fluid')
    plt.ylabel('v (m/s)')
    plt.xlabel('steps')    
    
    fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)
    
    ax1.plot(steps,T_numerical[:,int(3*p.N/8)])
    ax1.title.set_text('plot of temperature in cilinder 1')
    ax1.set_ylabel('T (k)')
    ax1.set_xlabel('steps')
    
    ax2.plot(steps,T_numerical[:,int(3*p.N/8)])
    ax2.title.set_text('plot of temperature in cilinder 2')
    ax2.set_ylabel('T (k)')
    ax2.set_xlabel('steps')
   
    ax3.plot(steps,T_numerical[:,int(5*p.N/8)])
    ax3.title.set_text('plot of temperature in cilinder 3')
    ax3.set_ylabel('T (k)')
    ax3.set_xlabel('steps')
    
    ax4.plot(steps,T_numerical[:,int(7*p.N/8)])
    ax4.title.set_text('plot of temperature in cilinder 4')
    ax4.set_ylabel('T (k)')
    ax4.set_xlabel('steps')
    
    plt.suptitle('Temperature of fluid')
    
    
    
    fig,((bx1,bx2),(bx3,bx4)) = plt.subplots(2,2)
    
    bx1.plot(steps,T_wall_numerical[:,int(p.N/8)])
    bx1.title.set_text('plot of wall temperature in cilinder 1')
    bx1.set_ylabel('T (k)')
    bx1.set_xlabel('steps')
    
    bx2.plot(steps,T_wall_numerical[:,int(3*p.N/8)])
    bx2.title.set_text('plot of wall temperature in cilinder 2')
    bx2.set_ylabel('T (k)')
    bx2.set_xlabel('steps')
   
    bx3.plot(steps,T_wall_numerical[:,int(5*p.N/8)])
    bx3.title.set_text('plot of wall temperature in cilinder 3')
    bx3.set_ylabel('T (k)')
    bx3.set_xlabel('steps')
    
    bx4.plot(steps,T_wall_numerical[:,int(7*p.N/8)])
    bx4.title.set_text('plot of wall temperature in cilinder 4')
    bx4.set_ylabel('T (k)')
    bx4.set_xlabel('steps')
    
    plt.suptitle('Temperature of wall')
    
     
    