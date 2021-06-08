# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 16:57:00 2021

@author: ronar
"""

import parameters as p
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

def darcyfriction(v,T,Tb=0):
    Re=p.Reynolds(v, T)
    a=1/(1+((Re/2712)**8.4))
    b=1/(1+((Re/(150*2*p.r/p.eff))**1.8))
    ''' KLOPT NIET '''
    df=64/Re
    #df=((64/Re)**a)*((0.75*np.log(Re/5.37))**(2*(a-1)*b))*((0.88*np.log(6.82*2*p.r/p.eff))**(2*(a-1)*(1-b)))
    f=np.sum(df)
    return f

def gravity(T):
    # grav=np.zeros(p.N)
    # for n in np.arange(0,p.N,1):
    #     grav[n]=(p.rho_0-p.rho_0*p.beta*(T[n]-p.T_0))*np.sin(p.angle[n])        
    # grav1=sum(grav)
    grav1=np.sum((p.rho_0-p.rho_0*p.beta*(T-p.T_0))*np.sin(p.angle))
    return grav1


def phiAB(v,T,Tb):
    ''' energy transport form the fluid to the wall  '''
    phi=np.zeros(p.N)

    for n in range(phi.size()):
        phi[n]=p.h_AB(n, v, T[n])*p.Opp_wallAB*(T[n]-Tb[n])
    return phi

def phiBC(v,T,Tb):
    ''' energy transport from the wall to the surrounding water '''
    phi=np.zeros(p.N)

    for n in np.arange(0,p.N,1):
        phi[n]=p.h_BC(n)*p.Opp_wallBC[n]*(Tb[n]-p.T_c)
    return phi

def Twallzero(v,T,Tb):
    twallzero=np.zeros(p.N)
    PhiAB=phiAB(v, T, Tb)
    PhiBC=phiBC(v, T, Tb)
    for n in np.arange(0,p.N,1):
        twallzero[n]=Tb[n]-p.Vwall[n]*p.u*p.rho_wall - PhiAB[n] + PhiBC[n]
    return twallzero
        
def Tfluidzero(v,T,Tb):
    tfluidzero=np.zeros(p.N)
    PhiAB=phiAB(v, T, Tb)
    for n in np.arange(0,p.N,1):
        tfluidzero[n]=PhiAB[n]+v*p.rho_0*p.C_pfluid*np.pi*p.r**2
    return tfluidzero

def velocityzero(v,T,Tb):
    f_D=darcyfriction(v, T)
    Velocityzero=-2*f_D*(1/p.r)*v**2 - ((p.kw1+p.kw2)*v**2)/p.length + (p.g/(p.rho_0*p.N))*gravity(T)
    return Velocityzero

def system(U):
    v=U[0]
    T,Tb=np.array_split(U[1:],2)
    print('v=', v ,'\n T=',T[0], '\n Tb=',Tb)
    vzero=velocityzero(v, T[0], Tb[0])
    Tnzero=Tfluidzero(v, T[0], Tb[0])
    Tbzero=Twallzero(v, T[0], Tb[0])
    return np.concatenate((vzero),(Tnzero),(Tbzero))


initialguess=(p.v_steadystate0,p.T_steadystate0,p.Tb_steadystate0)

answer=optimize.fsolve(system,initialguess)























