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
    Vwall=(p.length/p.N)*np.pi*((p.r+dr)**2-p.r**2)
    P=p.u*p.rho_wall*Vwall
    return P

#def gravity(T):
#    #grav=0
#    gravitysegment=(density(T)*np.sin(p.angle))
#    grav=-np.sum(gravitysegment)
#   # for j in np.arange(0,p.N,1) :
#   #     grav=grav+(density(T[j])*np.sin(p.angle[j]))
#    return grav

def gravity(v,T):
    #grav=0
    gravitysegment=(p.rho_0-p.rho_0*p.beta*(T-p.T_0))*np.sin(p.angle) 
    grav=np.sum(gravitysegment)
    if v>0:
        grav1=grav
    else:
        grav1=-grav
   # for j in np.arange(0,p.N,1) :
   #     grav=grav+(density(T[j])*np.sin(p.angle[j]))
    return grav1

def darcyfriction(v,T):
    Re=p.Reynolds(v, T)
    a=1/(1+((Re/2712)**8.4))
    b=1/(1+((Re/(150*2*p.r/p.eff))**1.8))
    df=((64/Re)**a)*((0.75*np.log(Re/5.37))**(2*(a-1)*b))*((0.88*np.log(6.82*2*p.r/p.eff))**(2*(a-1)*(1-b)))
    f=np.sum(df)
    return f

"ODE`s"
#ODE from momentum balance: (dv/dt)
def f1v(v,T):
    friction=darcyfriction(abs(v), T)/p.N
    if v>0:
        f=((-2*friction*v**2)/p.r)-(((p.kw1+p.kw2)*v**2)/p.length)+((p.g/(p.N*p.rho_0))*gravity(v,T))
    else:
        f=((+2*friction*v**2)/p.r)+(((p.kw1+p.kw2)*v**2)/p.length)+((p.g/(p.N*p.rho_0))*(-1)*gravity(v,T))
    #f=(1/(p.tVsys*p.rho_0))*(-(v**2)*np.pi*p.r*p.rho_0*(friction*(p.length/p.N)+p.r*(p.kw1+p.kw2))+np.pi*p.g*(p.length/p.N)*(p.r**2)*gravity(T)) 
    return f

#ODE from heat balance over fluid: (dTn/dt)
def f2Tn(v,T,n,Tb):
    if (n>0 and n<0.25*p.N) or n==0 or (n>0.5*p.N and n<p.N) or n==0.5*p.N:
        #phiq_coef=2*p.h_AB(n,abs(v),T[n])*(1/(p.r*p.rho_0*p.C_pfluid))*(T[n]-Tb[n])
        phiq=p.h_AB(n,v,T[n])*p.Opp_wallAB*(T[n]-Tb[n])
    elif (n>0.25*p.N and n<0.5*p.N) or n==0.25*p.N:
        #phiq_coef=0 #p.h_AB(n,v,T[n])*p.Opp_wallAB*(T[n]-Tb[n])
        phiq=p.h_AB(n,v,T[n])*p.Opp_wallAB*(T[n]-Tb[n])
        
    if v>0:
        if n>0 and n<p.N:
            deltaT=(T[n-1]-T[n])
        elif n==0:
            deltaT=(T[p.N-1]-T[n])
        elif n==p.N or n>p.N:
            print('error in deltaT from function f2Tn n: n is exceeding bounds')
    else:
        if n>=0 and n<(p.N-1):
            deltaT=(T[n+1]-T[n])
        elif n==p.N-1:
            deltaT=(T[0]-T[p.N-1])
        elif n==p.N or n>p.N:
            print('error in deltaT from function f2Tn n: n is exceeding bounds')
            
    #f=-phiq_coef+(p.N/p.length)*v*deltaT
    f=(1/(p.rho_0*p.C_pfluid*p.Vsys_segm))*((p.rho_0*p.C_pfluid*np.pi*p.r**2)*v*deltaT-phiq) # nog even goed controleren!!!!
    return f

#ODE from heat balance over wall: (dTb/dt)
def f3Tb(T,n,Tb,v):
    
    if n<0.25*p.N:
        #dr=p.dr1
        P=gammaheating(p.dr1)
    elif (n>0.25*p.N or n==0.25*p.N) and (n<0.5*p.N) :
        P=gammaheating(p.dr1)
        #dr=p.dr1
    elif (n>0.5*p.N or n==0.5*p.N) and (n<0.75*p.N) :
        P=gammaheating(p.dr2)
        #dr=p.dr2
    elif (n>0.75*p.N or n==0.75*p.N) and (n<p.N) :
        P=gammaheating(p.dr2)
        #dr=p.dr1
    elif n>p.N:
        
        print('error in gammaheating calculation in f3Tb')
    if (n>0 and n<0.25*p.N) or n==0 or (n>0.5*p.N and n<p.N) or n==0.5*p.N:
        #phiq_coef=2*p.h_AB(n, abs(v), T[n])*p.r*(1/(2*p.r*dr+dr**2))*(T[n]-Tb[n])
        phiqAB=p.h_AB(n,v,T[n])*p.Opp_wallAB*(T[n]-Tb[n])
    elif (n>0.25*p.N and n<0.5*p.N) or n==0.25*p.N:
        #phiq_coef=0 #2*p.h_AB(n, v, T[n])*p.r*(1/(2*p.r*dr+dr**2))*(T[n]-Tb[n])
        phiqAB=p.h_AB(n,v,T[n])*p.Opp_wallAB*(T[n]-Tb[n])    
        
    phiqBC=p.h_BC(n)*p.Opp_wallBC[n]*(Tb[n]-p.T_c)
    #f=(p.u/p.C_pwall)+ (phiq_coef/(p.C_pwall*p.rho_wall))-(2*p.h_BC(n)*(p.r+dr)*(Tb[n]-p.T_c)/(p.C_pwall*p.rho_wall*(2*p.r*dr+dr**2)))
    f=(1/(p.rho_wall*p.C_pwall*p.Vwall[n]))*(-phiqBC+phiqAB+P)   # nog goed checken of dit klopt!!!!
    return f
    


