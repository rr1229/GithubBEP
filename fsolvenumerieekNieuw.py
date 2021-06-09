# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 16:57:00 2021

@author: ronar
"""

import parameters as p
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
plt.close('all')


def darcyfriction(v,T,Tb=0):
    Re=p.Reynolds(v, T)
    a=1/(1+((Re/2712)**8.4))
    b=1/(1+((Re/(150*2*p.r/p.eff))**1.8))
    #df=64/Re
    df=((64/Re)**a)*((0.75*np.log(Re/5.37))**(2*(a-1)*b))*((0.88*np.log(6.82*2*p.r/p.eff))**(2*(a-1)*(1-b)))
    f=np.sum(df)
    if v>0:
        f1=f
    else:
        f1=-f
    return f1

def gravity(v,T):
    grav=np.zeros(p.N)
    for n in np.arange(0,p.N,1):
        grav[n]=(p.rho_0-p.rho_0*p.beta*(T[n]-p.T_0))*np.sin(p.angle[n])        
    #grav1=np.sum((p.rho_0-p.rho_0*p.beta*(T-p.T_0))*np.sin(p.angle))
    if v>0:
        grav1=sum(grav)
    else:
        grav1=-sum(grav)
    return grav1


def phiAB(v,T,Tb):
    ''' energy transport form the fluid to the wall  '''
    phi=np.zeros(p.N)
    for n in range(len(phi)):
        phi[n]=p.h_AB(n, v, T[n])*p.Opp_wallAB*(T[n]-Tb[n])
    return phi

def phiBC(v,T,Tb):
    ''' energy transport from the wall to the surrounding water '''
    phi=np.zeros(p.N)

    for n in np.arange(0,p.N,1):
        phi[n]=p.h_BC(n,v,T,Tb)*p.Opp_wallBC[n]*(Tb[n]-p.T_c)
    return phi

def Twallzero(v,T,Tb):
    twallzero=np.zeros(p.N)
    PhiAB=phiAB(v, T, Tb)
    PhiBC=phiBC(v, T, Tb)
    for n in np.arange(0,p.N,1):
        twallzero[n]=-PhiAB[n] + PhiBC[n]-p.Vwall[n]*p.u*p.rho_wall 
    return twallzero
        
def Tfluidzero(v,T,Tb):
    
    tfluidzero=np.zeros(p.N)
    PhiAB=phiAB(v, T, Tb)
    
    Constatns=p.rho_0*p.C_pfluid*np.pi*p.r**2
    
    if v>0:
        tfluidzero[0]=PhiAB[0]+v*Constatns*(T[0]-T[p.N-1])
        for n in np.arange(0,p.N,1):
            if n > 0:
                tfluidzero[n]=PhiAB[n]+v*Constatns*(T[n]-T[n-1])
                
    if v<=0:
        tfluidzero[p.N-1]=PhiAB[p.N-1]+v*Constatns*(T[p.N-1]-T[0])
        for n in np.arange(0,p.N,1):
            if n < p.N-1:
                tfluidzero[n]=PhiAB[n]+v*Constatns*(T[n]-T[n+1])
            #tfluidzero[n]=(-2*p.h_AB(n, v, T[n])*(T[n]-Tb[n])/(p.r*p.rho_0*p.C_pfluid))-(p.N*(T[n]-T[n-1])*v/p.length)
    #tfluidzero[0]=(-2*p.h_AB(n, v, T[0])*(T[0]-Tb[0])/(p.r*p.rho_0*p.C_pfluid))-(p.N*(T[0]-T[p.N-1])*v/p.length)
        #tfluidzero[n]=PhiAB[n]+v*p.rho_0*p.C_pfluid*np.pi*p.r**2*(T[n]-T[n-1])
    return tfluidzero

def velocityzero(v,T,Tb):
    f_D=darcyfriction(v, T)/p.N
    if v>0:
        Velocityzero=-2*f_D*(1/p.r)*v**2 - ((p.kw1+p.kw2)*v**2)/p.length + (p.g/(p.rho_0*p.N))*gravity(v,T)
    else:
        Velocityzero=-2*f_D*(1/p.r)*v**2 + ((p.kw1+p.kw2)*v**2)/p.length + (p.g/(p.rho_0*p.N))*gravity(v,T)
    #Velocityzero=-2*f_D*v**2/p.r + (p.g/(p.rho_0*p.N))*gravity(v,T) - ((p.kw1+p.kw2)*v**2)/p.length 
    #f_D=darcyfriction(v, T)
    #Velocityzero=-2*f_D*(1/p.r)*v**2 + (p.g/(p.rho_0))*gravity(v,T)
    return Velocityzero

def system(U):
    v=U[0]
    T,Tb=np.array_split(U[1:],2)
    
    #print('v=', v ,'\n T=',T, '\n Tb=',Tb)
    
    vzero=velocityzero(v, T, Tb)
    Tnzero=Tfluidzero(v, T, Tb)
    Tbzero=Twallzero(v, T, Tb)
    return np.concatenate([np.array([vzero]),Tnzero,Tbzero])

def con(v,T,Tb):
    return np.concatenate([np.array([v]),T,Tb])

def split(U):
    v=U[0]
    T,Tb=np.array_split(U[1:],2)
    return v, T, Tb


initialguess=np.concatenate([np.array([p.v_steadystate0]),p.T_steadystate0,p.Tb_steadystate0])

if True:
    answer=optimize.fsolve(system,initialguess)
else:
    answer=initialguess

v=answer[0]
T,Tb=np.array_split(answer[1:],2)


print('- - - - - - - \n \n \t velocity \n \t v = %.3e m/s \n \n mean temperature \n \t first quadrant \t T = %3.f K | Tb = %3.f \n \t second quadrant \t T = %3.f K | Tb = %3.f \n \t third quadrant \t T = %3.f K | Tb = %3.f \n \t fourth quadrant \t T = %3.f K | Tb = %3.f \n - - - - - - -  ' %(v,np.mean(T[0:int(p.N/4)]),np.mean(Tb[0:int(p.N/4)]),np.mean(T[int(p.N/4)+1:int(p.N/2)]),np.mean(Tb[int(p.N/4)+1:int(p.N/2)]),np.mean(T[int(p.N/2)+1:int(p.N*3/4)]),np.mean(Tb[int(p.N/2)+1:int(p.N*3/4)]),np.mean(T[int(p.N*3/4)+1:int(p.N-1)]),np.mean(Tb[int(p.N*3/4)+1:int(p.N-1)])) )



fig, (ax1,ax2) = plt.subplots(1,2)
plt.suptitle('Temperature')
ax1.plot(T,'k')
ax1.set_xlabel('$l$')
ax1.set_ylabel('$T$')
ax1.set_title('Temperature of fluid')


ax2.plot(Tb,'k')
ax2.set_xlabel('$l$')
ax2.set_ylabel('$T_B$')
ax2.set_title('Temperature of wall')


### TEST FUNCTIONS ###
if True:
    fig, (ax1,ax2) = plt.subplots(1,2)
    plt.suptitle('Transport')
    ax1.plot(phiAB(v,T,Tb),'k')
    ax1.set_ylabel('$\phi$')
    ax1.set_xlabel('l')
    ax1.set_title('$\phi_{AB}$')    
    
    ax2.plot(phiBC(v,T,Tb),'k')
    ax2.set_ylabel('$\phi$')
    ax2.set_xlabel('l')
    ax2.set_title('$\phi_{BC}$')

    fig, (ax1) = plt.subplots(1,1)
    
    Y=np.zeros(p.N)
    for i in range(len(Y)):
        Y[i]=p.h_fluid(v, T[i])
    plt.suptitle('Transport')
    ax1.plot(Y,'k')
    ax1.set_ylabel('$h_{fluid}$')
    ax1.set_xlabel('l')
    #ax1.set_title('$\phi_{AB}$')    
    
    fig, (ax1,ax2,ax3) = plt.subplots(1,3)    
    ax1.plot(p.Reynolds(v, T),'k')
    ax1.set_ylabel('$Reynodls$')
    ax1.set_xlabel('l')
    ax1.set_title('$\phi_{BC}$')

    ax2.plot(p.Greatz(v),'k')
    ax2.set_ylabel('$Greatz$')
    ax2.set_xlabel('l')
    ax2.set_title('$\phi_{BC}$')
    
    ax3.plot(p.Prandtl(T),'k')
    ax3.set_ylabel('$Prandtl$')
    ax3.set_xlabel('l')
    ax3.set_title('$\phi_{BC}$')
    
    fig, (gx1)=plt.subplots(1,1)
    h_c=np.zeros(p.N)
    for i in np.arange(0,p.N,1):
        h_c[i]=p.h_outside(i,v,T,Tb,i)
        
    gx1.plot(h_c,'k')
    gx1.set_ylabel('$h_C$')
    gx1.set_xlabel('$l$')
    gx1.set_title('Heat transfer coefficient of surrounding water')

    fig, (ax1,ax2) = plt.subplots(1,2)
    ax1.plot(velocityzero(v,T,Tb))
    ax2.plot(Tfluidzero(v,T,Tb),'r')
    ax2.plot(Twallzero(v,T,Tb),'b')









