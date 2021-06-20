# -*- coding: utf-8 -*-
"""
Created on Sun Jun 20 15:36:18 2021

@author: ronar
"""

import parameters as p
import numpy as np
import matplotlib.pyplot as plt
import fsolvenumerieekNieuw as fs
from tqdm import tqdm 

blocks=p.N 
vo=fs.answer[0]
To,Tbo=np.array_split(fs.answer[1:],2)

def h_BC(n,v,T,Tb,m):
    h_BC= ((1/p.h_wall(n))+(1/p.h_outside(n, v, T, Tb,m)))**-1 
    return h_BC

def Ch_BC(n,l):
    hbc=np.zeros(p.N)
    for i in np.arange(0,p.N,1):
        hbc[i]=h_BC(n, vo, To, Tbo, i)
        
    if l >= 0 and l<p.length/4:
        hbcC1=np.mean(hbc[0:int(p.N/4)])
    elif l >= p.length/4 and l<p.length/2:
        hbcC1=np.mean(hbc[int(p.N/4):int(p.N/2)])
    elif l >= p.length/2 and l<3*p.length/4:
        hbcC1=np.mean(hbc[int(p.N/2):int(3*p.N/4)])
    elif l >= 3*p.length/4 and l<=p.length:
        hbcC1=np.mean(hbc[int(3*p.N/4):int(p.N-1)])  
    if l<p.length:
        hbcC2=h_BC(n,vo,To,Tbo,int(l*p.N/p.length))
    else:
        hbcC2=h_BC(n,vo,To,Tbo,int(p.N-1))
    
    return hbcC2

def drwall(i):
    if i<0.25*p.N:
        dr_wall=p.dr[i]
    elif (i>0.25*p.N or i==0.25*p.N) and (i<p.deel*p.N) :
        dr_wall=p.dr[i]
    elif (i>p.deel*p.N or i==p.deel*p.N) and (i<0.75*p.N) :
        dr_wall=p.dr[i]
    elif (i>0.75*p.N or i==0.75*p.N) and (i<p.N) :
        dr_wall=p.dr[i]
    elif i>p.N:
        print('error in dr_wall')
    return dr_wall
        
def angle(l):
    if l<0.25*p.length:
        angle=-p.anglepipe
    elif (l>0.25*p.length or l==0.25*p.length) and (l<p.length*0.5) :
        angle=-0.5*np.pi
    elif (l>p.length*0.5 or l==p.length*0.5) and (l<0.75*p.length) :
        angle=p.anglepipe
    elif (l>0.75*p.length or l==0.75*p.length) and (l<p.length) :
        angle=0.5*np.pi
    elif l>p.length:
        print('error in angle')    
    return angle

def darcyfriction(v,T):
    Re=p.Reynolds(v, T)
    a=1/(1+((Re/2712)**8.4))
    b=1/(1+((Re/(150*2*p.r/p.eff))**1.8))
    df=((64/Re)**a)*((0.75*np.log(Re/5.37))**(2*(a-1)*b))*((0.88*np.log(3.41*2*p.r/p.eff))**(2*(a-1)*(1-b)))
    #df=((64/Re)**a)*((0.75*np.log(Re/5.37))**(2*(a-1)*b))*((0.88*np.log(6.82*2*p.r/p.eff))**(2*(a-1)*(1-b)))
    f=np.sum(df)
    return f
        
def constantA(n,T,v,l):
    ans1=-Ch_BC(n,l)*p.h_AB(n, v, T)*2*(p.r+drwall(n))
    ans2=v * p.r * p.rho_0 * p.C_pfluid * (p.h_AB(n, v, T)*p.r + Ch_BC(n,l)*(p.r+drwall(n)))
    An=ans1/ans2
    #An=(p.h_AB(n,v,T)*p.h_BC(n)*4*np.pi*((p.r**2)+p.r*drwall(n)))/(p.h_BC(n)*2*np.pi*(p.r+drwall(n))+p.h_AB(n,v,T)*2*np.pi*p.r)
    return An

def constantB(n,T,v,l):
    ans1=p.h_AB(n, v, T)*drwall(n)*(2*p.r+drwall(n))*p.u*p.rho_wall
    ans2=2*p.h_AB(n, v, T)*Ch_BC(n,l)*(p.r+drwall(n))*p.T_c
    ans3=v * p.r * p.rho_0 * p.C_pfluid * (p.h_AB(n, v, T)*p.r + Ch_BC(n,l)*(p.r+drwall(n))) 
    Bn=(ans1+ans2)/ans3
    #ans1=(np.pi*p.C_pfluid*p.rho_0*p.r**2)**-1
    #ans2=((-2*p.rho_wall*p.mu_20*(np.pi**2)*p.h_AB(n,v,T)*((2*drwall(n)*p.r**2)+p.r*drwall(n)**2))
    #      -p.h_AB(n,v,T)*p.h_BC(n)*4*p.r*(p.r+drwall(n))*p.T_c*np.pi**2)
    #ans3=(p.h_BC(n)*2*np.pi*(p.r+drwall(n))+p.h_AB(n,v,T)*2*np.pi*p.r)**(-1)
    #Bn=ans1*ans2*ans3
    return Bn

def constantC(n,v,T,l):
    friction=darcyfriction(v, T)/blocks
    Cn=(-p.g*p.beta*np.sin(angle(l)))/(p.kw1+p.kw2+friction*p.length*((4*p.r)**-1))
    return Cn



def T_l(l,v,Told):
    #constantes defenieren van de vergelijking voor de temperatuur
    Avector=np.zeros(blocks)
    Bvector=np.zeros(blocks)
    Fvector=np.zeros(blocks)
    for n in np.arange(0,blocks,1):
        if n<0.25*p.N or n==0.25*p.N:
            length_for_n=n*p.length/p.N
        if n>0.25*p.N and n<0.5*p.N or n==0.5*p.N:
            length_for_n=p.length*0.25 + (p.deel-0.25)*p.length*(n-0.25*p.N)/(0.25*p.N)
        if n>0.5*p.N and n<0.75*p.N or n==0.75*p.N:
            length_for_n=p.length*0.25 + (p.deel-0.25)*p.length + (n-0.5*p.N)*(0.75-p.deel)*(p.length)/(p.N*0.25)
        if n>0.75*p.N:
            length_for_n=p.length*0.25 + (p.deel-0.25)*p.length + (0.75-p.deel)*(p.length) + (n-0.75*p.N)*p.length/p.N
        
        Avector[n]=constantA(n, Told[n], v, length_for_n)
        Bvector[n]=constantB(n, Told[n], v, length_for_n)
        Fvector[n]=Bvector[n]/Avector[n]
        
    F1=np.mean(Fvector[0:round(0.25*p.N)])
    F2=np.mean(Fvector[round(0.25*p.N):round(0.5*p.N)])
    F3=np.mean(Fvector[round(0.5*p.N):round(0.75*p.N)])
    F4=np.mean(Fvector[round(0.75*p.N):round(p.N)])
    
    A1=np.mean(Avector[0:round(0.25*p.N)])
    A2=np.mean(Avector[round(0.25*p.N):round(0.5*p.N)])
    A3=np.mean(Avector[round(0.5*p.N):round(0.75*p.N)])
    A4=np.mean(Avector[round(0.75*p.N):round(p.N)])
    
    B1=np.mean(Bvector[0:round(0.25*p.N)])
    B2=np.mean(Bvector[round(0.25*p.N):round(0.5*p.N)])
    B3=np.mean(Bvector[round(0.5*p.N):round(0.75*p.N)])
    B4=np.mean(Bvector[round(0.75*p.N):round(p.N)])
    
    
    
    deel2=((p.deel-0.25)/2)+0.25
    # A1=constantA(round(p.N*0.125),Told[round(p.N*0.125)],v,0.125*p.length)
    # A2=constantA(round(p.N*deel2),Told[round(p.N*deel2)],v,deel2*p.length)
    # A3=constantA(round(p.N*0.625),Told[round(p.N*0.625)],v,0.625*p.length)
    # A4=constantA(round(p.N*0.875),Told[round(p.N*0.875)],v,0.875*p.length)
    # B1=constantB(round(p.N*0.125),Told[round(p.N*0.125)],v,0.125*p.length)
    # B2=constantB(round(p.N*deel2),Told[round(p.N*deel2)],v,deel2*p.length)
    # B3=constantB(round(p.N*0.625),Told[round(p.N*0.625)],v,0.625*p.length)
    # B4=constantB(round(p.N*0.875),Told[round(p.N*0.875)],v,0.875*p.length)
    # A1=constantA(round(p.N*0),Told[round(p.N*0)],v,0*p.length)
    # A2=constantA(round(p.N*0.25),Told[round(p.N*0.25)],v,0.25*p.length)
    # A3=constantA(round(p.N*p.deel),Told[round(p.N*p.deel)],v,p.deel*p.length)
    # A4=constantA(round(p.N*0.75),Told[round(p.N*0.75)],v,0.75*p.length)
    # B1=constantB(round(p.N*0),Told[round(p.N*0)],v,0*p.length)
    # B2=constantB(round(p.N*0.25),Told[round(p.N*0.25)],v,0.25*p.length)
    # B3=constantB(round(p.N*p.deel),Told[round(p.N*p.deel)],v,p.deel*p.length)
    # B4=constantB(round(p.N*0.75),Told[round(p.N*0.75)],v,0.75*p.length)
    

    # F1=B1/A1 #A1/B1
    # F2=B2/A2#A2/B2
    # F3=B3/A3#A3/B3
    # F4=B4/A4#A4/B4
    # C=np.array((F1,F2,F3,F4))
    
    #lengtes definieren op verschillende punten  in de buis (randvoorwaarden op hoeken)
    l1=0.25*p.length
    l2=p.deel*p.length
    l3=0.75*p.length
    l4=p.length
    
    #Constantes bepalen in vergelijking wat volgt uit de ODE en randvoorwaarden
    ans1=(1-np.e**((-A4*l4-A1*l1+A2*l1-A2*l2+A3*l2-A3*l3+A4*l3)))**-1
    ans2=(F1-F2)*(np.e**((-A1*l1)))
    +(F2-F3)*np.e**((-A2*l2+A2*l1-A1*l1))
    +(F3-F4)*np.e**((-A3*l3+A3*l2-A2*l2+A2*l1-A1*l1))
    +(F4-F1)*np.e**((-A4*l4+A4*l3-A3*l3+A3*l2-A2*l2+A2*l1-A1*l1))

    constant1=(ans1*ans2)   
    constant4=(-F1+F4+constant1)*np.e**(-A4*l4)
    constant3=(-F4+F3+constant4*np.e**(A4*l3))*np.e**(-A3*l3)   
    constant2=(-F3+F2+constant3*np.e**(A3*l2))*np.e**(-A2*l2)   
    
    #Bepalen in welk deel van de buis we zitten, dus welke constantes we moeten hebben
    if l<0.25*p.length:
        n=round(p.N*0.10)
        constantn=constant1
        A=A1
        B=B1
        sign=1
    elif (l>0.25*p.length or l==0.25*p.length) and (l<p.deel*p.length) :
        n=round(p.N*0.40)
        constantn=constant2
        A=A2
        B=B2
        sign=1
    elif (l>p.deel*p.length or l==p.deel*p.length) and (l<0.75*p.length) :
        n=round(p.N*deel2)
        constantn=constant3
        sign=-1
        A=A3
        B=B3
    elif (l>0.75*p.length or l==0.75*p.length) and (l<p.length) :
        n=round(p.N*0.90)
        constantn=constant4
        A=A4
        B=B4
        sign=-1
    elif l>p.length:
        print('error in dr_wall')
    n=round(l/p.length) 
    #het in elkaar zetten van de vergelijking voor de Temperatuur.
    #B=constantB(n,Told[n],v,l)
    #A=constantA(n,Told[n],v,l)
    
    T_l=(-B/A+constantn*np.e**(A*l))
    
    return T_l, constantn, A, B

def velocity(v_i,Told):
    #constantes bereken in de vergelijking voor de snelheid
    deel2=((p.deel-0.25)/2)+0.25
    lengthpart2=p.deel-0.25
    lengthpart3=0.75-p.deel
    #C1=constantC(round(p.N*0.10),v_i,Told)
    #C2=constantC(round(p.N*deel2),v_i,Told)
    #C3=constantC(round(p.N*0.60),v_i,Told)
    #C4=constantC(round(p.N*0.90),v_i,Told)
    
    #lengtes definieren in buis in elke deel van de buis
    step1=(0.25*p.length/round(blocks*0.25))
    step2=(lengthpart2*p.length/round(blocks*0.25))
    step3=(lengthpart3*p.length/round(blocks*0.25))
    step4=(0.25*p.length/round(blocks*0.25))
    
    # values1=np.linspace(0,0.25*p.length-step1,round(blocks*0.25))
    # values2=np.linspace(0.25*p.length,p.deel*p.length-step2,round(blocks*0.25))
    # values3=np.linspace(p.deel*p.length,0.75*p.length-step3,round(blocks*0.25))
    # values4=np.linspace(0.75*p.length,p.length-step4,round(blocks*0.25))
    
    values1=np.linspace(0+0.5*step1,0.25*p.length-0.5*step1,round(blocks*0.25))
    values2=np.linspace(0.25*p.length+0.5*step2,p.deel*p.length-0.1*step2,round(blocks*0.25))
    values3=np.linspace(p.deel*p.length+0.1*step3,0.75*p.length-0.5*step3,round(blocks*0.25))
    values4=np.linspace(0.75*p.length+0.5*step4,p.length-0.5*step4,round(blocks*0.25))
    #Lege vectoren maken voor het berekenen voor elke Temperatuur in de verschillende delen 
    T_l_values1=np.zeros(round(blocks*0.25))
    T_l_values2=np.zeros(round(blocks*0.25))
    T_l_values3=np.zeros(round(blocks*0.25))
    T_l_values4=np.zeros(round(blocks*0.25))
    
    #Lege vectoren maken voor berkenen van constante C in elk deel
    C1=np.zeros(round(blocks*0.25))
    C2=np.zeros(round(blocks*0.25))
    C3=np.zeros(round(blocks*0.25))
    C4=np.zeros(round(blocks*0.25))
    
    #Temperatuur berekenen deel 1 buis
    i=0
    j=0
    Co=np.zeros(blocks)
    C_A=np.zeros(blocks)
    C_B=np.zeros(blocks)
    for l in values1:
        valuesofT=T_l(l,v_i,Told)
        T_l_values1[i]=valuesofT[0]
        Co[j]=valuesofT[1]
        C_A[j]=valuesofT[2]
        C_B[j]=valuesofT[3]
        C1[i]=constantC(j,v_i,Told,l)
        j=j+1
        i=i+1
    #Temperatuur berekenen deel 2 buis   
    i=0
    for l in values2:
        valuesofT=T_l(l,v_i,Told)
        T_l_values2[i]=valuesofT[0]
        Co[j]=valuesofT[1]
        C_A[j]=valuesofT[2]
        C_B[j]=valuesofT[3]
        #T_l_values2[i]=T_l(l,v_i,Told)[0]
        #Co[j]=T_l(l,v_i,Told)[1]
        C2[i]=constantC(j,v_i,Told,l)
        j=j+1
        i=i+1
    #Temperatuur berekenen deel 3 buis  
    i=0
    for l in values3:
        valuesofT=T_l(l,v_i,Told)
        T_l_values3[i]=valuesofT[0]
        Co[j]=valuesofT[1]
        C_A[j]=valuesofT[2]
        C_B[j]=valuesofT[3]
        # T_l_values3[i]=T_l(l,v_i,Told)[0]
        # Co[j]=T_l(l,v_i,Told)[1]
        C3[i]=constantC(j,v_i,Told,l)
        j=j+1
        i=i+1
    #Temperatuur berekenen deel 4 buis    
    i=0
    for l in values4:
        valuesofT=T_l(l,v_i,Told)
        T_l_values4[i]=valuesofT[0]
        Co[j]=valuesofT[1]
        C_A[j]=valuesofT[2]
        C_B[j]=valuesofT[3]
        #T_l_values4[i]=T_l(l,v_i,Told)[0]
        #Co[j]=T_l(l,v_i,Told)[1]
        C4[i]=constantC(j,v_i,Told,l)
        j=j+1
        i=i+1
    T1new=np.append(T_l_values1,T_l_values2)
    T2new=np.append(T_l_values3,T_l_values4)
    Tnew=np.append(T1new,T2new)
    length12=np.append(values1,values2)
    length34=np.append(values3,values4)
    length=np.append(length12,length34)
    #uiteindelijk in elkaar zetten van de vergelijking voor de snelheid
    ans1=sum(C1*T_l_values1)*step1  #+C1*step1*(p.T_0+p.beta**-1)
    ans2=sum(C2*T_l_values2)*step2  #+C2*step2*(p.T_0+p.beta**-1)
    ans3=sum(C3*T_l_values3)*step3  #+C3*step3*(p.T_0+p.beta**-1)
    ans4=sum(C4*T_l_values4)*step4  #+C4*step4*(p.T_0+p.beta**-1)
    v=np.sqrt(abs(ans1+ans2+ans3+ans4))
    return v, Tnew, Co, length, C_A, C_B


#initial condition of Temperature:
#T_0=np.zeros(blocks)
#for n in np.arange(0,blocks,1):
#    if n<0.25*blocks:
#        T_0[n]=(70+273.15)
#    elif (n>0.25*blocks or n==0.25*blocks) and (n<0.5*blocks) :
#        T_0[n]=(50+273.15)
#    elif (n>0.5*blocks or n==0.5*blocks) and (n<0.75*blocks) :
#        T_0[n]=(55+273.15)
#    elif (n>0.75*blocks or n==0.75*blocks) and (n<blocks) :
#        T_0[n]=(55+273.15)
#    elif n>blocks:
#        print('error in T_0 creation')
T_0=To





'''berekenen van initial snelheid aan de hand van the initial condition van de temperature. '''
#deel2=((p.deel-0.25)/2)+0.25
#C1=constantC(round(p.N*0.10),vo,T_0,p.length*0.1) #!!!!!!!!!! hier nog even naar de velocity kijken 
#C2=constantC(round(p.N*deel2),vo,T_0,p.length*0.3)
#C3=constantC(round(p.N*0.60),vo,T_0,p.length*0.6)
#C4=constantC(round(p.N*0.90),vo,T_0,p.length*0.9)

#lengtes definieren in buis in elke deel van de buis
values1=np.linspace(0,0.25*p.length-p.length/p.N,round(blocks*0.25))
values2=np.linspace(0.25*p.length,p.deel*p.length-p.length/p.N,round(blocks*0.25))
values3=np.linspace(p.deel*p.length,0.75*p.length-p.length/p.N,round(blocks*0.25))
values4=np.linspace(0.75*p.length,p.length-p.length/p.N,round(blocks*0.25))
C1=np.zeros(round(blocks*0.25))
C2=np.zeros(round(blocks*0.25))
C3=np.zeros(round(blocks*0.25))
C4=np.zeros(round(blocks*0.25))

#constantes berekenen buis
i=0
for l in values1:
    C1[i]=constantC(0+i,vo,T_0,values1[i])
    C2[i]=constantC(round(0.25*p.N+i),vo,T_0,values2[i])
    C3[i]=constantC(round(0.5*p.N+i),vo,T_0,values3[i])
    C4[i]=constantC(round(0.75*p.N+i),vo,T_0,values4[i])
    i=i+1



ans1=sum(C1*T_0[round(0):round(0.25*blocks)])*(p.length/blocks)
ans2=sum(C2*T_0[round(0.25*blocks):round(p.deel*blocks)])*(p.length/blocks)
ans3=sum(C3*T_0[round(p.deel*blocks):round(0.75*blocks)])*(p.length/blocks)
ans4=sum(C4*T_0[round(0.75*blocks):round(blocks)])*(p.length/blocks)
v_0=np.sqrt(abs(ans1+ans2+ans3+ans4))
'''_________________________________________________________________'''






'''array maken voor velocity '''
Velocityend1=np.array([v_0])
i=0
vi=Velocityend1[i]
vi_1= velocity(vi,T_0)[0]    #
Velocityend1=np.append(Velocityend1,vi_1)
Tnew1=T_0

''' Iteratieve methode om de snelheid te bepalen.'''
while abs(vi-vi_1)>0.00000003 and i<40:
#for iterations in tqdm(np.arange(0,8,1)):
    i=i+1
    vi=Velocityend1[i]
    veli=velocity(vi,Tnew1)
    vi_1=veli[0]
    Tnew1=veli[1]
    C0=veli[2]
    length=veli[3]
    Velocityend1=np.append(Velocityend1,vi_1)
    
    
    
    
    
"de plots"    
fig,((kx1,kx2)) = plt.subplots(1,2)

kx1.plot(Velocityend1)
kx1.set_title('Velocity over iterations')
kx1.set_ylabel('$v$')
kx1.set_xlabel('$i$')

lengthvector=np.linspace(0,p.length,p.N)
kx2.plot(length,Tnew1)
kx2.set_title('Temperature of fluid')
kx2.set_ylabel('$T$')
kx2.set_xlabel('$l$')
kx2.set_xlim(0, p.length)



fig,((ex1)) = plt.subplots(1,1)
plt.suptitle('End result of Newton Raphson method and the analytical method, \n with NR: v= %.3e and Analytisch: v=%.3e' %(vo,Velocityend1[-1]) )
ex1.plot(lengthvector,To)
ex1.plot(length,Tnew1)
ex1.set_title('Temperature of fluid')
ex1.set_ylabel('$T$')
ex1.set_xlabel('$l$')
ex1.legend(['Newton Rahpson method','Analytical method'])









