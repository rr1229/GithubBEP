# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 15:14:12 2021

@author: ronar
"""
import parameters as p
import numpy as np
import matplotlib.pyplot as plt



def density(T):
    rho=p.rho_0*(1-p.beta*(T-p.T_0))
    return rho 

def gammaheating(dr):
    Vwall=(p.l/p.N)*np.pi*((p.r+dr)**2-p.r**2)
    P=p.u*p.rho_wall*Vwall
    return P

def gravity(T):
    grav=0
    for j in np.arange(0,p.N,1) :
        grav=grav+(density(T[j])*np.sin(p.angle[j]))
    return grav

"ODE`s"
#ODE from momentum balance: (dv/dt)
def f1v(v,T):
    f=(1/(p.tVsys*p.rho_0))*(-(v**2)*np.pi*p.r*p.rho_0*(p.f*p.l+p.r*(p.kw1+p.kw2))+np.pi*p.g*p.l*(p.r**2)*gravity(T)) 
    return f

#ODE from heat balance over fluid: (dTn/dt)
def f2Tn(v,T,n,Tb):
    deltal=p.l/p.N
    if (n>0 and n<0.25*p.N) or n==0 or (n>0.5*p.N and n<p.N) or n==0.5*p.N:
        phiq=p.h_AB*p.Opp_wallAB*(T[n]-Tb[n])
    elif (n>0.25*p.N and n<0.5*p.N) or n==0.25*p.N:
        phiq=0
        
    if n>0 and n<p.N:
        deltaT=(T[n-1]-T[n])
    elif n==0:
        deltaT=(T[p.N-1]-T[n])
    elif n==p.N or n>p.N:
        print('error in deltaT from function f2Tn n: n is exceeding bounds')
    
    f=(1/(p.rho_0*p.C_pfluid*p.Vsys_segm))*((p.rho_0*p.C_pfluid*np.pi*p.r**2)*v*deltaT+phiq) # nog even goed controleren!!!!
    return f

#ODE from heat balance over wall: (dTb/dt)
def f3Tb(T,n,Tb):
    if (n>0 and n<0.25*p.N) or n==0 or (n>0.5*p.N and n<p.N) or n==0.5*p.N:
        phiq=p.h_AB*p.Opp_wallAB*(T[n]-Tb[n])
    elif (n>0.25*p.N and n<0.5*p.N) or n==0.25*p.N:
        phiq=0    
        
    if n<0.25*p.N:
        P=gammaheating(p.dr1)
    elif (n>0.25*p.N or n==0.25*p.N) and (n<0.5*p.N) :
        P=0
    elif (n>0.5*p.N or n==0.5*p.N) and (n<0.75*p.N) :
        P=gammaheating(p.dr2)
    elif (n>0.75*p.N or n==0.75*p.N) and (n<p.N) :
        P=gammaheating(p.dr2)
    elif n>p.N:
        print('error in gammaheating calculation in f3Tb')
    
    f=(1/(p.rho_wall*p.C_pwall*p.Vwall[n]))*(p.h_BC*p.Opp_wallBC[n]*(p.T_c-Tb[n])+phiq+P)   # nog goed checken of dit klopt!!!!
    return f
    


