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
    df=((64/Re)**a)*((0.75*np.log(Re/5.37))**(2*(a-1)*b))*((0.88*np.log(3.41*2*p.r/p.eff))**(2*(a-1)*(1-b)))
    #df=((64/Re)**a)*((0.75*np.log(Re/5.37))**(2*(a-1)*b))*((0.88*np.log(6.82*2*p.r/p.eff))**(2*(a-1)*(1-b)))
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
        Velocityzero=-f_D*(1/(4*p.r))*v**2 - ((p.kw1+p.kw2)*v**2)/p.length + (p.g/(p.rho_0*p.N))*gravity(v,T)
    else:
        Velocityzero=-f_D*(1/(4*p.r))*v**2 + ((p.kw1+p.kw2)*v**2)/p.length + (p.g/(p.rho_0*p.N))*gravity(v,T)
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
#initialguess=np.concatenate([np.array([vend[len(vend)-1]]),p.T_steadystate0,p.Tb_steadystate0])

if True:
    #answer=optimize.newton_krylov(system,initialguess) # doet het nu niet 
    #answer=optimize.newton(system,initialguess) # doet het alleeen als we de uitkomst van rungakutta pakken als initial condition
    answer=optimize.fsolve(system,initialguess)

else:
    answer=initialguess

v=answer[0]
T,Tb=np.array_split(answer[1:],2)


print('- - - - - - - \n \n \t velocity \n \t v = %.3e m/s \n \n mean temperature \n \t first quadrant \t T = %3.f K | Tb = %3.f \n \t second quadrant \t T = %3.f K | Tb = %3.f \n \t third quadrant \t T = %3.f K | Tb = %3.f \n \t fourth quadrant \t T = %3.f K | Tb = %3.f \n - - - - - - -  ' %(v,np.mean(T[0:int(p.N/4)]),np.mean(Tb[0:int(p.N/4)]),np.mean(T[int(p.N/4)+1:int(p.N/2)]),np.mean(Tb[int(p.N/4)+1:int(p.N/2)]),np.mean(T[int(p.N/2)+1:int(p.N*3/4)]),np.mean(Tb[int(p.N/2)+1:int(p.N*3/4)]),np.mean(T[int(p.N*3/4)+1:int(p.N-1)]),np.mean(Tb[int(p.N*3/4)+1:int(p.N-1)])) )


fig, (ax1,ax2) = plt.subplots(1,2)
plt.suptitle('Temperature')
ax1.plot(T,'k')
#ax1.plot(p.T_steadystate0,'r')
#ax1.plot(Trk,'b')
ax1.set_xlabel('$l$')
ax1.set_ylabel('$T$')
ax1.set_title('Temperature of fluid')


ax2.plot(Tb,'k')
ax2.set_xlabel('$l$')
ax2.set_ylabel('$T_B$')
ax2.set_title('Temperature of wall')


### TEST FUNCTIONS ###
if False:
    fig, (arx1,arx2) = plt.subplots(1,2)
    plt.suptitle('Transport')
    arx1.plot(phiAB(v,T,Tb),'k')
    arx1.set_ylabel('$\phi_{AB}$')
    arx1.set_xlabel('l')
    arx1.set_title('Heat transport from the \n fluid to the wall')    
    
    arx2.plot(phiBC(v,T,Tb),'k')
    arx2.set_ylabel('$\phi_{BC}$')
    arx2.set_xlabel('l')
    arx2.set_title('Heat transport from the \n wall to the surrounding water')


    fig, px1 = plt.subplots(1,1)
    Y=np.zeros(p.N)
    Y2=np.zeros(p.N)
    for i in range(len(Y)):
        Y[i]=p.h_fluid(v, T[i])
        Y2[i]=p.h_outside(i, v, T, Tb, i)
    plt.suptitle('Transport')
    px1.plot(Y,'k')
    px1.set_ylabel('$h_{fluid}$')
    px1.set_xlabel('l')
    px1.set_title('heat transfer coefficient of the fluid')    
    
    
    fig, (qx1,qx3) = plt.subplots(1,2)   
    plt.suptitle('Dimensionless numbers,\n  We have a Greatz number of: Gr= %.3e' %(p.Greatz(v)))
    qx1.plot(p.Reynolds(v, T),'k')
    qx1.set_ylabel('$Re$')
    qx1.set_xlabel('l')
    qx1.set_title('Reynolds Number')
    
    qx3.plot(p.Prandtl(T),'k')
    qx3.set_ylabel('$Pr$')
    qx3.set_xlabel('l')
    qx3.set_title('Prandtl')
    
    
    fig, (gx1)=plt.subplots(1,1)
    h_c=np.zeros(p.N)
    for i in np.arange(0,p.N,1):
        h_c[i]=p.h_outside(i,v,T,Tb,i)
        
    gx1.plot(h_c,'k')
    gx1.set_ylabel('$h_C$')
    gx1.set_xlabel('$l$')
    gx1.set_title('Heat transfer coefficient of surrounding water')


    fig, (abx2) = plt.subplots(1,1)
    plt.suptitle('Order of Error of result of numerical method, \n  $\epsilon (v) =$ %.3e' %(velocityzero(v, T, Tb)))
    abx2.plot(Tfluidzero(v,T,Tb),'r')
    abx2.plot(Twallzero(v,T,Tb),'b')
    abx2.set_ylabel('$\epsilon(T)$ & $\epsilon(T_B)$')
    abx2.set_xlabel('l')
    









