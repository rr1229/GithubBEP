# -*- coding: utf-8 -*-
"""
Created on Sat May 15 16:23:30 2021

@author: ronar
"""
import parameters as p
import numpy as np
import matplotlib.pyplot as plt

blocks=100 #number of parts we use for numerical integration, would work best in multiples of 4

#initial condition of Temperature:


def drwall(i):
    if i<0.25*p.N:
        dr_wall=p.dr1
    elif (i>0.25*p.N or i==0.25*p.N) and (i<0.5*p.N) :
        dr_wall=p.dr1
    elif (i>0.5*p.N or i==0.5*p.N) and (i<0.75*p.N) :
        dr_wall=p.dr2
    elif (i>0.75*p.N or i==0.75*p.N) and (i<p.N) :
        dr_wall=p.dr2
    elif i>p.N:
        print('error in dr_wall')
        
    return dr_wall

def darcyfriction(v,T):
    Re=p.Reynolds(v, T)
    a=1/(1+((Re/2712)**8.4))
    b=1/(1+((Re/(150*2*p.r/p.eff))**1.8))
    df=((64/Re)**a)*((0.75*np.log(Re/5.37))**(2*(a-1)*b))*((0.88*np.log(6.82*2*p.r/p.eff))**(2*(a-1)*(1-b)))
    f=np.sum(df)
    return f
        
def constantA(n,T,v):
    ans1=p.h_BC(n)*p.h_AB(n, v, T)*2*(p.r+drwall(n))
    ans2=v * p.r * p.rho_0 * p.C_pfluid * (p.h_AB(n, v, T)*p.r + p.h_BC(n)*(p.r+drwall(n)))
    An=ans1/ans2
    #An=(p.h_AB(n,v,T)*p.h_BC(n)*4*np.pi*((p.r**2)+p.r*drwall(n)))/(p.h_BC(n)*2*np.pi*(p.r+drwall(n))+p.h_AB(n,v,T)*2*np.pi*p.r)
    return An

def constantB(n,T,v):
    ans1=p.h_AB(n, v, T)*drwall(n)*(2*p.r+drwall(n))*p.u*p.rho_wall
    ans2=2*p.h_AB(n, v, T)*p.h_BC(n)*(p.r+drwall(n))*p.T_c
    ans3=v * p.r * p.rho_0 * p.C_pfluid * (p.h_AB(n, v, T)*p.r + p.h_BC(n)*(p.r+drwall(n))) 
    Bn=(ans1+ans2)/ans3
    #ans1=(np.pi*p.C_pfluid*p.rho_0*p.r**2)**-1
    #ans2=((-2*p.rho_wall*p.mu_20*(np.pi**2)*p.h_AB(n,v,T)*((2*drwall(n)*p.r**2)+p.r*drwall(n)**2))
    #      -p.h_AB(n,v,T)*p.h_BC(n)*4*p.r*(p.r+drwall(n))*p.T_c*np.pi**2)
    #ans3=(p.h_BC(n)*2*np.pi*(p.r+drwall(n))+p.h_AB(n,v,T)*2*np.pi*p.r)**(-1)
    #Bn=ans1*ans2*ans3
    return Bn

def constantC(n,v,T):
    friction=darcyfriction(v, T)/blocks
    Cn=(-p.g*p.beta*np.sin(p.angle[n]))/(p.kw1+p.kw2+friction*p.length*(p.r**-1))
    return Cn

#def constant2C(n,v,T):                          Darcy friction is actually temperature and so space dependend!!!!
#    Cn=(-p.g*p.beta*np.sin(p.angle[n]))/(p.kw1+p.kw2+darcyfrinction(v,T[n])*p.l*(p.r**-1))
#    return Cn


def T_l(l,v,Told):
    #constantes defenieren van de vergelijking voor de temperatuur
    A1=constantA(round(p.N*0.10),Told[round(p.N*0.10)],v)
    A2=constantA(round(p.N*0.40),Told[round(p.N*0.40)],v)
    A3=constantA(round(p.N*0.60),Told[round(p.N*0.60)],v)
    A4=constantA(round(p.N*0.90),Told[round(p.N*0.90)],v)
    B1=constantB(round(p.N*0.10),Told[round(p.N*0.10)],v)
    B2=constantB(round(p.N*0.40),Told[round(p.N*0.40)],v)
    B3=constantB(round(p.N*0.60),Told[round(p.N*0.60)],v)
    B4=constantB(round(p.N*0.90),Told[round(p.N*0.90)],v)
    F1=A1/B1
    F2=A2/B2
    F3=A3/B3
    F4=A4/B4
    
    #lengtes definieren op verschillende punten  in de buis (randvoorwaarden op hoeken)
    l1=0.25*p.length
    l2=0.5*p.length
    l3=0.75*p.length
    l4=p.length
    
    #Constantes bepalen in vergelijking wat volgt uit de ODE en randvoorwaarden
    ans1=(1-np.e**((v**-1)*(A4*l4+A1*l1-A2*l1+A2*l2-A3*l2+A3*l3-A4*l3)))**-1
    ans2=(F4-F3)*(np.e**((v**-1)*(-A4*l3)))
    +(F3-F2)*np.e**((v**-1)*(-A3*l2+A3*l3-A4*l3))
    +(F2-F1)*np.e**((v**-1)*(-A2*l1+A2*l2-A3*l2+A3*l3-A4*l3))
    +(F1-F4)*np.e**((v**-1)*(A1*l1-A2*l1+A2*l2-A3*l2+A3*l3-A4*l3))
    
    constant4=(ans1*ans2)   
    constant1=F1-F4+constant4*np.e**(A4*l4/v)
    constant2=(F2-F1+constant1*np.e**(A1*l1/v))*np.e**(-A2*l1/v)   
    constant3=(F3-F2+constant2*np.e**(A2*l2/v))*np.e**(-A3*l2/v)   
    
    #Bepalen in welk deel van de buis we zitten, dus welke constantes we moeten hebben
    if l<0.25*p.length:
        n=round(p.N*0.10)
        constantn=constant1
    elif (l>0.25*p.length or l==0.25*p.length) and (l<0.5*p.length) :
        n=round(p.N*0.40)
        constantn=constant2
    elif (l>0.5*p.length or l==0.5*p.length) and (l<0.75*p.length) :
        n=round(p.N*0.60)
        constantn=constant3
    elif (l>0.75*p.length or l==0.75*p.length) and (l<p.length) :
        n=round(p.N*0.90)
        constantn=constant4
    elif l>p.length:
        print('error in dr_wall')
        
    #het in elkaar zetten van de vergelijking voor de Temperatuur.
    B=constantB(n,Told[n],v)
    A=constantA(n,Told[n],v)
    T_l=-B/A+constantn*np.e**(A*l)
    
    return T_l

def velocity(v_i,Told):
    #constantes bereken in de vergelijking voor de snelheid
    C1=constantC(round(p.N*0.10),v_i,Told)
    C2=constantC(round(p.N*0.40),v_i,Told)
    C3=constantC(round(p.N*0.60),v_i,Told)
    C4=constantC(round(p.N*0.90),v_i,Told)
    
    #lengtes definieren in buis in elke deel van de buis
    values1=np.linspace(0,0.25*p.length-(0.25*p.length/25),round(blocks*0.25))
    values2=np.linspace(0.25*p.length,0.5*p.length-(0.25*p.length/25),round(blocks*0.25))
    values3=np.linspace(0.5*p.length,0.75*p.length-(0.25*p.length/25),round(blocks*0.25))
    values4=np.linspace(0.75*p.length,p.length-(0.25*p.length/25),round(blocks*0.25))
    
    #Lege vectoren maken voor het berekenen voor elke Temperatuur in de verschillende delen 
    T_l_values1=np.zeros(round(blocks*0.25))
    T_l_values2=np.zeros(round(blocks*0.25))
    T_l_values3=np.zeros(round(blocks*0.25))
    T_l_values4=np.zeros(round(blocks*0.25))
    
    #Temperatuur berekenen deel 1 buis
    i=0
    for l in values1:
        T_l_values1[i]=T_l(l,v_i,Told)
        i=i+1
    #Temperatuur berekenen deel 2 buis   
    i=0
    for l in values2:
        T_l_values2[i]=T_l(l,v_i,Told)
        i=i+1
    #Temperatuur berekenen deel 3 buis  
    i=0
    for l in values3:
        T_l_values3[i]=T_l(l,v_i,Told)
        i=i+1
    #Temperatuur berekenen deel 4 buis    
    i=0
    for l in values4:
        T_l_values4[i]=T_l(l,v_i,Told)
        i=i+1
    T1new=np.append(T_l_values1,T_l_values2)
    T2new=np.append(T_l_values3,T_l_values4)
    Tnew=np.append(T1new,T2new)
    #uiteindelijk in elkaar zetten van de vergelijking voor de snelheid
    ans1=C1*sum(T_l_values1)*(p.length/blocks)  +C1*(p.length/blocks)*(p.T_0+p.beta**-1)
    ans2=C2*sum(T_l_values2)*(p.length/blocks)  +C2*(p.length/blocks)*(p.T_0+p.beta**-1)
    ans3=C3*sum(T_l_values3)*(p.length/blocks)  +C3*(p.length/blocks)*(p.T_0+p.beta**-1)
    ans4=C4*sum(T_l_values4)*(p.length/blocks)  +C4*(p.length/blocks)*(p.T_0+p.beta**-1)
    v=np.sqrt(abs(ans1+ans2+ans3+ans4))
    return v, Tnew


#initial condition of Temperature:
T_0=np.zeros(blocks)
for n in np.arange(0,blocks,1):
    if n<0.25*blocks:
        T_0[n]=(70+273.15)
    elif (n>0.25*blocks or n==0.25*blocks) and (n<0.5*blocks) :
        T_0[n]=(50+273.15)
    elif (n>0.5*blocks or n==0.5*blocks) and (n<0.75*blocks) :
        T_0[n]=(55+273.15)
    elif (n>0.75*blocks or n==0.75*blocks) and (n<blocks) :
        T_0[n]=(55+273.15)
    elif n>blocks:
        print('error in T_0 creation')

#berekenen van initial snelheid aan de hand van the initial condition van de temperature. 
C1=constantC(round(p.N*0.10),0.03,T_0) #!!!!!!!!!! hier nog even naar de velocity kijken 
C2=constantC(round(p.N*0.40),0.03,T_0)
C3=constantC(round(p.N*0.60),0.03,T_0)
C4=constantC(round(p.N*0.90),0.03,T_0)
ans1=C1*sum(T_0[round(0):round(0.25*blocks)])*(p.length/blocks)
ans2=C2*sum(T_0[round(0.25*blocks):round(0.5*blocks)])*(p.length/blocks)
ans3=C3*sum(T_0[round(0.5*blocks):round(0.75*blocks)])*(p.length/blocks)
ans4=C4*sum(T_0[round(0.75*blocks):round(blocks)])*(p.length/blocks)
v_0=np.sqrt(abs(ans1+ans2+ans3+ans4))

# array maken voor velocity 
Velocityend=np.array([v_0])
i=0
vi=Velocityend[i]
vi_1= 0.00205150   #velocity(vi,T_0)[0]    #
Velocityend=np.append(Velocityend,vi_1)
Tnew=T_0
# Iteratieve methode om de snelheid te bepalen.
while abs(vi-vi_1)>0.0000003 and i<15:
    i=i+1
    vi=Velocityend[i]
    vi_1=velocity(vi,Tnew)[0]
    Tnew=velocity(vi,Tnew)[1]
    Velocityend=np.append(Velocityend,vi_1)
    


