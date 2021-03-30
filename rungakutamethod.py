# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 20:09:25 2021

@author: ronar
"""
import numpy as np
import parameters as p
import functions as f
import matplotlib.pyplot as plt



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
        k1Tb[i]=h*f.f3Tb(T,i,Tb)
        dTndt[i]=f.f2Tn(v, T, i, Tb)
        dTbdt[i]=f.f3Tb(T,i,Tb)
        
    k2v=h*f.f1v(v+0.5*k1v,T+0.5*k1Tn)
    for i in np.arange(0,p.N,1):
        k2Tn[i]=h*f.f2Tn(v+0.5*k1v,T+0.5*k1Tn,i,Tb+0.5*k1Tb)
        k2Tb[i]=h*f.f3Tb(T+0.5*k1Tn,i,Tb+0.5*k1Tb)
    
    k3v=h*f.f1v(v+0.5*k2v,T+0.5*k2Tn)
    for i in np.arange(0,p.N,1):
        k3Tn[i]=h*f.f2Tn(v+0.5*k2v,T+0.5*k2Tn,i,Tb+0.5*k2Tb)
        k3Tb[i]=h*f.f3Tb(T+0.5*k2Tn,i,Tb+0.5*k2Tb)
        
    k4v=h*f.f1v(v+k3v,T+k3Tn)
    for i in np.arange(0,p.N,1):
        k4Tn[i]=h*f.f2Tn(v+k3v,T+k3Tn,i,Tb+k3Tb)
        k4Tb[i]=h*f.f3Tb(T+k3Tn,i,Tb+k3Tb)
    
    v=v+(1/6)*(k1v+2*k2v+2*k3v+k4v)
    T=T+(1/6)*(k1Tn+2*k2Tn+2*k3Tn+k4Tn)
    Tb=Tb+(1/6)*(k1Tb+2*k2Tb+2*k3Tb+k4Tb)
    RK=[v,T,Tb,dvdt,dTndt,dTbdt]
    return RK

initial=[p.v0,p.T0,p.Tb0]

t=0
h=5
steps=1000
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

for tn in np.arange(0,steps,1):
    RK=rungakutta(RK,h)
    Tend[:,tn+1]=RK[1]
    Tbend[:,tn+1]=RK[2]
    vend[tn+1]=RK[0]
    t=tn*h+t
    tend[tn+1]=t
    dvdtend[tn+1]=RK[3]
    dTndtend[:,tn+1]=RK[4]
    dTbdtend[:,tn+1]=RK[5]
        
        

# def k1f1(h,v,T):
#     k=h*f.f1v(v,T)
#     return k

# def k1f2(h,v,T,n,Tb):
#     k=h*f.f2Tn(v,T,n,Tb)
#     return k

# def k1f3(h,Tb,n,T):
#     k=h*f.f3Tb(T,n,Tb)
#     return k

# T=p.T0
# Tb=p.Tb0
# v=p.v0
# T1=np.zeros(p.N)
# Tb1=np.zeros(p.N)
# t=0

# vend=np.zeros(10001)
# vend[0]=v
# tend=np.zeros(10001)
# tend[0]=t
# Tend=np.zeros([p.N,10001])
# Tbend=np.zeros([p.N,10001])
# Tend[:,0]=p.T0
# Tbend[:,0]=p.Tb0

# for tn in np.arange(0,10000,1):
#     v1=v+0.5*k1f1(h,v,T)
#     for i in np.arange(0,p.N,1):
#         T1[i]=T[i]+0.5*k1f2(h, v, T, i, Tb)
#         Tb1[i]=Tb[i]+0.5*k1f3(h, Tb, i, T)
#     Tend[:,tn+1]=T1
#     Tbend[:,tn+1]=Tb1
#     T=T1
#     Tb=Tb1
#     vend[tn+1]=v1
#     v=v1
#     t=tn*h+t
#     tend[tn+1]=t
    
plt.plot(tend,vend)    
    
    
    
    
    