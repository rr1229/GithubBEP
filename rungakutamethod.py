# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 20:09:25 2021

@author: ronar
"""
import numpy as np
import parameters as p
import functions as f
import matplotlib.pyplot as plt
from tqdm import tqdm 
plt.close('all')
def euler(R,h):
    v=R[0]
    T=0
    k1v=h*f.f1v(v, T)
    v=v+k1v
    return [v,0,0,k1v,0,0]

def rungakutta(R,h):
    k1Tn=np.zeros(p.N)
    k2Tn=np.zeros(p.N)
    k3Tn=np.zeros(p.N)
    k4Tn=np.zeros(p.N)
    k1Tb=np.zeros(p.N)
    k2Tb=np.zeros(p.N)
    k3Tb=np.zeros(p.N)
    k4Tb=np.zeros(p.N)
    dTndt=np.zeros(p.N)
    dTbdt=np.zeros(p.N)
    v=R[0]
    T=R[1]
    Tb=R[2]
    
    k1v=h*f.f1v(v,T)
    dvdt=f.f1v(v, T)
    for i in np.arange(0,p.N,1):
        k1Tn[i]=h*f.f2Tn(v,T,i,Tb)
        k1Tb[i]=h*f.f3Tb(T,i,Tb,v)
        dTndt[i]=f.f2Tn(v, T, i, Tb)
        dTbdt[i]=f.f3Tb(T,i,Tb, v)
        
    k2v=h*f.f1v(v+0.5*k1v,T+0.5*k1Tn)
    for i in np.arange(0,p.N,1):
        k2Tn[i]=h*f.f2Tn(v+0.5*k1v,T+0.5*k1Tn,i,Tb+0.5*k1Tb)
        k2Tb[i]=h*f.f3Tb(T+0.5*k1Tn,i,Tb+0.5*k1Tb, v+0.5*k1v)
    
    k3v=h*f.f1v(v+0.5*k2v,T+0.5*k2Tn)
    for i in np.arange(0,p.N,1):
        k3Tn[i]=h*f.f2Tn(v+0.5*k2v,T+0.5*k2Tn,i,Tb+0.5*k2Tb)
        k3Tb[i]=h*f.f3Tb(T+0.5*k2Tn,i,Tb+0.5*k2Tb, v+0.5*k2v)
        
    k4v=h*f.f1v(v+k3v,T+k3Tn)
    for i in np.arange(0,p.N,1):
        k4Tn[i]=h*f.f2Tn(v+k3v,T+k3Tn,i,Tb+k3Tb)
        k4Tb[i]=h*f.f3Tb(T+k3Tn,i,Tb+k3Tb, v+k3v)
    
    v=v+(1/6)*(k1v+2*k2v+2*k3v+k4v)
    T=T+(1/6)*(k1Tn+2*k2Tn+2*k3Tn+k4Tn)
    Tb=Tb+(1/6)*(k1Tb+2*k2Tb+2*k3Tb+k4Tb)
    RK=[v,T,Tb,dvdt,dTndt,dTbdt]
    return RK

#steps=5000
#def RK4solver(v0,T0,Tb0):
initial=[p.v0,p.T0,p.Tb0]
#    initial=[v0,T0,Tb0]
t=0
h=0.5
hmax=2*h
steps=5000#400000
RK=initial

vend=np.zeros(steps+1)
vend[0]=RK[0]
tend=np.zeros(steps+1)
tend[0]=t
Tend=np.zeros([p.N,steps+1])
Tbend=np.zeros([p.N,steps+1])
Tend[:,0]=RK[1]
Tbend[:,0]=RK[2]
dvdtend=np.zeros(steps+1)
dTndtend=np.zeros([p.N,steps+1])
dTbdtend=np.zeros([p.N,steps+1])


for tn in tqdm(np.arange(0,steps,1)):
    #RK=euler(RK,h)
    RK=rungakutta(RK,h)
    #print(tn)
    Tend[:,tn+1]=RK[1]
    Tbend[:,tn+1]=RK[2]
    vend[tn+1]=RK[0]
    t=h+t
    tend[tn+1]=t
    dvdtend[tn+1]=RK[3]
    dTndtend[:,tn+1]=RK[4]
    dTbdtend[:,tn+1]=RK[5]

    #return vend, Tend, Tbend, tend
    



if True:
#if False:
    plt.plot(tend,vend)    
    #tend3=tend
    #vend3=vend
    #'plot of velocity'
    plt.plot(tend,vend) 
    plt.title('velocity of fluid')
    plt.ylabel('v (m/s)')
    plt.xlabel('t (s)')    
    
    #'plot of points inside cilinder, temperature of fluid'
    fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)

    ax1.plot(tend,Tend[int(p.N/8)])
    ax1.title.set_text('cilinder 1')
    ax1.set_ylabel('T (k)')
    ax1.set_xlabel('t (s)')
    
    ax2.plot(tend,Tend[int(3*p.N/8)])
    ax2.title.set_text('cilinder 2')
    ax2.set_ylabel('T (k)')
    ax2.set_xlabel('t (s)')
   
    ax3.plot(tend,Tend[int(5*p.N/8)])
    ax3.title.set_text('cilinder 3')
    ax3.set_ylabel('T (k)')
    ax3.set_xlabel('t (s)')
    
    ax4.plot(tend,Tend[int(7*p.N/8)])
    ax4.title.set_text('cilinder 4')
    ax4.set_ylabel('T (k)')
    ax4.set_xlabel('t (s)')
    
    plt.suptitle('Temperature of fluid in middle point of each cilinder') 
    
    #'plot of points inside cilinder, temperature of wall'
    fig,((bx1,bx2),(bx3,bx4)) = plt.subplots(2,2)
    
    bx1.plot(tend,Tbend[int(p.N/8)])
    bx1.title.set_text('cilinder 1')
    bx1.set_ylabel('$T_B (k)$')
    bx1.set_xlabel('t (s)')
    
    bx2.plot(tend,Tbend[int(3*p.N/8)])
    bx2.title.set_text('cilinder 2')
    bx2.set_ylabel('$T_B (k)$')
    bx2.set_xlabel('t (s)')
   
    bx3.plot(tend,Tbend[int(5*p.N/8)])
    bx3.title.set_text('cilinder 3')
    bx3.set_ylabel('$T_B (k)$')
    bx3.set_xlabel('t (s)')
    
    bx4.plot(tend,Tbend[int(7*p.N/8)])
    bx4.title.set_text('cilinder 4')
    bx4.set_ylabel('$T_B (k)$')
    bx4.set_xlabel('t (s)')
    
    plt.suptitle('Temperature of wall in middle point of each cilinder')
    
    
    #'plot of mean values inside cilinders, temperature of fluid'
    fig,((cx1,cx2),(cx3,cx4)) = plt.subplots(2,2)
    
    TendC1=np.zeros(int(steps+1))
    TendC2=np.zeros(int(steps+1))
    TendC3=np.zeros(int(steps+1))
    TendC4=np.zeros(int(steps+1))
    #Tendall3=np.zeros((int(steps+1)))
    for i in np.arange(0,steps+1,1):
        TendC1[i]=np.mean((Tend[:,i])[0:int(p.N/4)])
        TendC2[i]=np.mean((Tend[:,i])[int(p.N/4):int(p.N/2)])
        TendC3[i]=np.mean((Tend[:,i])[int(p.N/2):int(3*p.N/4)])
        TendC4[i]=np.mean((Tend[:,i])[int(3*p.N/4):int(p.N-1)])
     #   Tendall3[i]=np.mean((Tend[:,i])[0:int(p.N-1)])

        
    cx1.plot(tend,TendC1)
    cx1.title.set_text('cilinder 1')
    cx1.set_ylabel('T (k)')
    cx1.set_xlabel('t (s)')
    
    cx2.plot(tend,TendC2)
    cx2.title.set_text('cilinder 2')
    cx2.set_ylabel('T (k)')
    cx2.set_xlabel('t (s)')

    cx3.plot(tend,TendC3)
    cx3.title.set_text('cilinder 3')
    cx3.set_ylabel('T (k)')
    cx3.set_xlabel('t (s)')
    
    cx4.plot(tend,TendC4)
    cx4.title.set_text('cilinder 4')
    cx4.set_ylabel('T (k)')
    cx4.set_xlabel('t (s)')
    
    plt.suptitle('plots of mean Temperatures of fluid in different cilinders')   
    
    
    #'plot of mean values inside cilinders, temperature of wall'    
    fig,((dx1,dx2),(dx3,dx4)) = plt.subplots(2,2)
    
    TbendC1=np.zeros(int(steps+1))
    TbendC2=np.zeros(int(steps+1))
    TbendC3=np.zeros(int(steps+1))
    TbendC4=np.zeros(int(steps+1))
    #Tbendall3=np.zeros(int(steps+1))
    for i in np.arange(0,steps+1,1):
        TbendC1[i]=np.mean((Tbend[:,i])[0:int(p.N/4)])
        TbendC2[i]=np.mean((Tbend[:,i])[int(p.N/4):int(p.N/2)])
        TbendC3[i]=np.mean((Tbend[:,i])[int(p.N/2):int(3*p.N/4)])
        TbendC4[i]=np.mean((Tbend[:,i])[int(3*p.N/4):int(p.N-1)])
     #   Tbendall3[i]=np.mean((Tbend[:,i])[0:int(p.N-1)])
    
    dx1.plot(tend,TbendC1)
    dx1.title.set_text('cilinder 1')
    dx1.set_ylabel('$T_B (k)$')
    dx1.set_xlabel('t (s)')
    
    dx2.plot(tend,TbendC2)
    dx2.title.set_text('cilinder 2')
    dx2.set_ylabel('$T_B (k)$')
    dx2.set_xlabel('t (s)')
   
    dx3.plot(tend,TbendC3)
    dx3.title.set_text('cilinder 3')
    dx3.set_ylabel('$T_B (k)$')
    dx3.set_xlabel('t (s)')
    
    dx4.plot(tend,TbendC4)
    dx4.title.set_text('cilinder 4')
    dx4.set_ylabel('$T_B (k)$')
    dx4.set_xlabel('t (s)')
    
    plt.suptitle(' plots of mean Temperatures of wall in different cilinders') 
    
    
    
    #plot of end temperature of wall and fluid
    lengthvector=np.linspace(0,p.length,p.N)
    fig,((ex1,ex2)) = plt.subplots(1,2)
    plt.suptitle('End result of Runga kutta method, v= %.3e' %(vend[steps-1]))
    ex1.plot(lengthvector,Tend[:,steps-1])
    ex1.set_title('Temperature of fluid')
    ex1.set_ylabel('$T$')
    ex1.set_xlabel('$l$')

    ex2.plot(lengthvector,Tbend[:,steps-1])
    ex2.set_title('Temperature of wall')
    ex2.set_ylabel('$T_B$')
    ex2.set_xlabel('$l$')
    
    

    # fig,((te1,te2,te3)) = plt.subplots(1,3)
    # plt.suptitle('Runga kutta method with different time steps')
    # te1.plot(tend2_9,vend2_9)
    # tb1.plot(tend3,vend3)
    # te1.set_title('velocity, with $\Delta t = 2.9$s')
    # tb1.set_title('velocity, with $\Delta t = 3$s')
    # te1.set_ylabel('$v$')
    # te1.set_xlabel('$t$')
    # tb1.set_ylabel('$v$')
    # tb1.set_xlabel('$t$')
    
    
    # te2.plot(tend2_9,Tendall2_9)
    # tb2.plot(tend3,Tendall3)
    # te2.set_title('mean temperature of the fluid, with $\Delta t = 2.9$s')
    # tb2.set_title('mean temperature of the fluid, with $\Delta t = 3$s')
    # te2.set_ylabel('$T$')
    # te2.set_xlabel('$t$')
    # tb2.set_ylabel('$T$')
    # tb2.set_xlabel('$t$')
    
    
    # te3.plot(tend2_9,Tbendall2_9)
    # tb3.plot(tend3,Tbendall3)
    # te3.set_title('mean temperature of the wall, with $\Delta t = 2.9$s')
    # tb3.set_title('mean temperature of the wall, with $\Delta t = 3$s')
    # te3.set_ylabel('$T$')
    # te3.set_xlabel('$t$')
    # tb3.set_ylabel('$T$')
    # tb3.set_xlabel('$t$')
    
    
    