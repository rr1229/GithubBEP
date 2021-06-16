# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 15:03:45 2021

@author: ronar
"""
import numpy as np
import parameters as p
import functions as f
import matplotlib.pyplot as plt
from tqdm import tqdm
import rungakutamethod as RK4

'comparing the different methods'
# fig,((ex1,ex2)) = plt.subplots(1,2)
# plt.suptitle('End result of Runga kutta method and Newton Raphson method, \n with RK: v= %.3e and NR: v=%.3e' %(vend[steps-1],v) )
# ex1.plot(lengthvector,Tend[:,steps-1])
# ex1.plot(lengthvector,T)
# ex1.set_title('Temperature of fluid')
# ex1.set_ylabel('$T$')
# ex1.set_xlabel('$l$')

# ex2.plot(lengthvector,Tbend[:,steps-1])
# ex2.plot(lengthvector,Tb)
# ex2.set_title('Temperature of wall')
# ex2.set_ylabel('$T_B$')
# ex2.set_xlabel('$l$')


# ex2.legend(['Runga kutta method','Newton Rahpson method'])



'Varieren van diktes van de buis'







'Varieren van Initial condition '
if True:
    steps=RK4.steps
    fig,((cx1,cx2),(cx3,cx4)) = plt.subplots(2,2)
    plt.suptitle('plots of mean Temperatures of fluid in different cilinders, with different starting temperatures') 
    fig,((dx1,dx2),(dx3,dx4)) = plt.subplots(2,2)
    plt.suptitle(' plots of mean Temperatures of wall in different cilinders, with different starting temperatures') 
    fig,((nx1)) = plt.subplots(1,1)
    vvari=0.00000001 #np.arange(0.02,0.201,0.02)
    T0=np.arange(10+273.15,70+273.15,10) #p.T0
    Tb0=np.arange(30+273.15,90+273.15,10) #p.Tb0
    k=0
    for T in T0:
        T_0=np.ones(p.N)*T
        Tb_0=np.ones(p.N)*Tb0[k]
        RK=RK4.RK4solver(vvari, T_0, Tb_0)
        vend=RK[0]
        Tend=RK[1]
        Tbend=RK[2]
        tend=RK[3]
        k=k+1
        nx1.plot(tend,vend) 
        
        #'plot of mean values inside cilinders, temperature of fluid' 
        TendC1=np.zeros(int(steps+1))
        TendC2=np.zeros(int(steps+1))
        TendC3=np.zeros(int(steps+1))
        TendC4=np.zeros(int(steps+1))
        for i in np.arange(0,steps+1,1):
            TendC1[i]=np.mean((Tend[:,i])[0:int(p.N/4)])
            TendC2[i]=np.mean((Tend[:,i])[int(p.N/4):int(p.N/2)])
            TendC3[i]=np.mean((Tend[:,i])[int(p.N/2):int(3*p.N/4)])
            TendC4[i]=np.mean((Tend[:,i])[int(3*p.N/4):int(p.N-1)])
    
            
        cx1.plot(TendC1)
        cx2.plot(TendC2)
        cx3.plot(TendC3)
        cx4.plot(TendC4)
        #'plot of mean values inside cilinders, temperature of wall'    
        
        
        TbendC1=np.zeros(int(steps+1))
        TbendC2=np.zeros(int(steps+1))
        TbendC3=np.zeros(int(steps+1))
        TbendC4=np.zeros(int(steps+1))
        for i in np.arange(0,steps+1,1):
            TbendC1[i]=np.mean((Tbend[:,i])[0:int(p.N/4)])
            TbendC2[i]=np.mean((Tbend[:,i])[int(p.N/4):int(p.N/2)])
            TbendC3[i]=np.mean((Tbend[:,i])[int(p.N/2):int(3*p.N/4)])
            TbendC4[i]=np.mean((Tbend[:,i])[int(3*p.N/4):int(p.N-1)])
        
        dx1.plot(tend,TbendC1)
        dx2.plot(tend,TbendC2)
        dx3.plot(tend,TbendC3)
        dx4.plot(tend,TbendC4)
        
        
        
        
        
    nx1.title.set_text('velocity of fluid, with different starting temperatures')
    nx1.set_ylabel('v (m/s)')
    nx1.set_xlabel('t (s)')  
    legend=[]
    z=0
    for Tint in T0:
        legend.append('T_0 = %.3f and Tb_0= %.3f ' %(Tint,Tb0[z]))
        z=z+1
    nx1.legend(legend)
        
    cx1.title.set_text('cilinder 1')
    cx1.set_ylabel('T (k)')
    cx1.set_xlabel('t (s)')
    
    cx2.plot(TendC2)
    cx2.title.set_text('cilinder 2')
    cx2.set_ylabel('T (k)')
    cx2.set_xlabel('t (s)')
    
    cx3.title.set_text('cilinder 3')
    cx3.set_ylabel('T (k)')
    cx3.set_xlabel('t (s)')
    
    cx4.title.set_text('cilinder 4')
    cx4.set_ylabel('T (k)')
    cx4.set_xlabel('t (s)')
    
    
    #plt.suptitle('plots of mean Temperatures of fluid in different cilinders')  
        
    
    dx1.title.set_text('cilinder 1')
    dx1.set_ylabel('$T_B (k)$')
    dx1.set_xlabel('t (s)')
    
    dx2.title.set_text('cilinder 2')
    dx2.set_ylabel('$T_B (k)$')
    dx2.set_xlabel('t (s)')
    
    dx3.title.set_text('cilinder 3')
    dx3.set_ylabel('$T_B (k)$')
    dx3.set_xlabel('t (s)')
    
    dx4.title.set_text('cilinder 4')
    dx4.set_ylabel('$T_B (k)$')
    dx4.set_xlabel('t (s)')
    
    #plt.suptitle(' plots of mean Temperatures of wall in different cilinders') 
    