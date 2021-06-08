# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 16:06:37 2021

@author: ronar
"""

import parameters as p
import numpy as np
import matplotlib.pyplot as plt


def phim(v):
    phi=v*p.rho_0*np.pi*p.r**2    #Mass flow of fluid through the pipe
    return phi

def density(T):
    rho=p.rho_0*(1-p.beta*(T-p.T_0)) # density dependend on temperature
    return rho 

def gravity(T):                     # the total gravity force in the system on the fluid in the direction of the flow
    #grav=0
    gravitysegment=(density(T)*np.sin(p.angle))
    grav=np.sum(gravitysegment)
   # for j in np.arange(0,p.N,1) :
   #     grav=grav+(density(T[j])*np.sin(p.angle[j]))
    return abs(grav)

def darcyfriction(v,T):
    Re=p.Reynolds(v, T)
    a=1/(1+((Re/2712)**8.4))
    b=1/(1+((Re/(150*2*p.r/p.eff))**1.8))
    df=((64/Re)**a)*((0.75*np.log(Re/5.37))**(2*(a-1)*b))*((0.88*np.log(6.82*2*p.r/p.eff))**(2*(a-1)*(1-b)))
    f=np.sum(df)
    return f

def phiAB(v,T,Tb):           #flow of heat trough wall exiting the fluid
    phiAB=np.zeros(p.N)
    for n in np.arange(0,p.N,1):
        if n<0.25*p.N:
            phiAB[n]=p.h_BC(n)*p.Opp_wallBC[n]*(Tb[n]-p.T_c)-p.Pgamma[n]
        elif (n>0.25*p.N or n==0.25*p.N) and (n<0.5*p.N) :
            phiAB[n]=0
        elif (n>0.5*p.N or n==0.5*p.N) and (n<0.75*p.N) :
            phiAB[n]=p.h_BC(n)*p.Opp_wallBC[n]*(Tb[n]-p.T_c)-p.Pgamma[n]
        elif (n>0.75*p.N or n==0.75*p.N) and (n<p.N) :
            phiAB[n]=p.h_BC(n)*p.Opp_wallBC[n]*(Tb[n]-p.T_c)-p.Pgamma[n]
        elif n>p.N:
                    print('error in temperature making function')
    return phiAB 

def Twall(v,T):               #temperature of the wall
    Tb=np.zeros(p.N)
    for n in np.arange(0,p.N,1):
        if n<p.N-1:
            Ta=T[n+1] # could later be changed: now we consider outflowing temperature to be temperature of segement
        elif n==p.N-1:
            Ta=T[0]
            
            
        if n<0.25*p.N:
            dr=p.dr1
            Tb[n]=((((p.h_AB(n, v, T[n])*p.r/(p.r+dr)))+p.h_BC(n))**(-1))*( (p.rho_wall*p.u*(2*p.r+dr)*dr/(2*p.r+2*dr)) + (p.h_AB(n, v, T[n])*p.r*T[n]/(p.r+dr)) + p.h_BC(n)*p.T_c)
            #Tb[n]=(p.Pgamma[n]+p.h_BC(n)*p.Opp_wallBC[n]*p.T_c+p.h_AB(n,v,T[n])*p.Opp_wallAB*Ta)/(p.h_AB(n,v,T[n])*p.Opp_wallAB+p.h_BC(n)*p.Opp_wallBC[n])
        elif (n>0.25*p.N or n==0.25*p.N) and (n<0.5*p.N) :
            dr=p.dr1
            Tb[n]=( (p.rho_wall*p.u*(2*p.r+dr)*dr/(p.h_BC(n)*(2*p.r+2*dr))) +  p.T_c)
            
            #Tb[n]=((((p.h_AB(n, v, T)*p.r/(p.r+dr)))+p.h_BC(n))**(-1))*( (p.rho_wall*p.u*(2*p.r+dr)*dr/(2*p.r+2*dr)) + (p.h_AB(n, v, T)*p.r*T[n]/(p.r+dr)) + p.h_BC(n)*p.T_c)
            
            #Tb[n]=(p.Pgamma[n]+p.h_BC(n)*p.Opp_wallBC[n]*p.T_c)/(p.h_BC(n)*p.Opp_wallBC[n])
        elif (n>0.5*p.N or n==0.5*p.N) and (n<0.75*p.N) :
            dr=p.dr2
            Tb[n]=((((p.h_AB(n, v, T[n])*p.r/(p.r+dr)))+p.h_BC(n))**(-1))*( (p.rho_wall*p.u*(2*p.r+dr)*dr/(2*p.r+2*dr)) + (p.h_AB(n, v, T[n])*p.r*T[n]/(p.r+dr)) + p.h_BC(n)*p.T_c)
                        
            #Tb[n]=(p.Pgamma[n]+p.h_BC(n)*p.Opp_wallBC[n]*p.T_c+p.h_AB(n,v,T[n])*p.Opp_wallAB*Ta)/(p.h_AB(n,v,T[n])*p.Opp_wallAB+p.h_BC(n)*p.Opp_wallBC[n])
        elif (n>0.75*p.N or n==0.75*p.N) and (n<p.N) :
            dr=p.dr2
            Tb[n]=((((p.h_AB(n, v, T[n])*p.r/(p.r+dr)))+p.h_BC(n))**(-1))*( (p.rho_wall*p.u*(2*p.r+dr)*dr/(2*p.r+2*dr)) + (p.h_AB(n, v, T[n])*p.r*T[n]/(p.r+dr)) + p.h_BC(n)*p.T_c)
            
            #Tb[n]=(p.Pgamma[n]+p.h_BC(n)*p.Opp_wallBC[n]*p.T_c+p.h_AB(n,v,T[n])*p.Opp_wallAB*Ta)/(p.h_AB(n,v,T[n])*p.Opp_wallAB+p.h_BC(n)*p.Opp_wallBC[n])
        elif n>p.N:
                    print('error in temperature wall making function')
    return Tb 
     
def velocity(v,T):                    #velocity of the fluid in the pipe
    friction=darcyfriction(v, T)
    #v=(gravity(T)/((4*friction*p.rho_0/(2*p.r))+(p.kw1+p.kw2)*p.rho_0))**(1/2)
    v=((p.g/(p.N*p.rho_0))*gravity(T)/((4*friction/(2*p.r))+(p.kw1+p.kw2)*(1/p.length)))**(1/2)
    return v

def temperature(v,T,Tb):                #temperature of the fluid within the pipe
    Tnew=np.zeros(p.N)
    #PhiAB=phiAB(v, T, Tb)
    #phimass=phim(v)
    #Tnew[p.N-1]=PhiAB[p.N-1]/(phimass*p.C_pfluid)+T[0]
    Tnew[0]=(((2*p.h_AB(0, v, T[0])/(p.r*p.rho_0*p.C_pfluid))+(p.N/p.length))**(-1))*((2*p.h_AB(0, v, T[0])*Tb[0]/(p.r*p.rho_0*p.C_pfluid))+(p.N*T[p.N-1]*v/p.length))
    for n in np.arange(1,p.N,1):
        Tnew[n]=(((2*p.h_AB(n, v, T[n])/(p.r*p.rho_0*p.C_pfluid))+(p.N/p.length))**(-1))*((2*p.h_AB(n, v, T[n])*Tb[n]/(p.r*p.rho_0*p.C_pfluid))+(p.N*T[n-1]*v/p.length))
        #Tnew[n]=PhiAB[n]/(phimass*p.C_pfluid)+T[n+1]
    return Tnew
    
