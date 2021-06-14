# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 22:38:30 2021

@author: ronar
"""
import parameters as p
import numpy as np
import matplotlib.pyplot as plt
import fsolvenumerieekNieuw as fs
from tqdm import tqdm 

blocks=p.N #number of parts we use for numerical integration, would work best in multiples of 4

#initial condition of Temperature:
vo=fs.answer[0]
To,Tbo=np.array_split(fs.answer[1:],2)
# vo=0.0010461944706357197
# To=np.array([325.10297393, 325.14567897, 325.18621626, 325.22469582,
#        325.26122209, 325.29589421, 325.32880629, 325.36004766,
#        325.38970312, 325.41785317, 326.388868  , 327.31859889,
#        328.20878404, 329.06109018, 329.87711519, 330.65839084,
#        331.40638538, 332.12250607, 332.80810163, 333.46446455,
#        332.61074395, 331.79916799, 331.02759573, 330.29399865,
#        329.5964546 , 328.93314206, 328.30233466, 327.70239603,
#        327.131775  , 326.58900095, 326.41285568, 326.24228654,
#        326.07711264, 325.91715904, 325.76225671, 325.61224231,
#        325.466958  , 325.3262512 , 325.18997447, 325.05798528])

# Tbo=np.array([325.85214267, 325.85681948, 325.86125857, 325.86547205,
#        325.8694714 , 325.8732675 , 325.8768707 , 325.8802908 ,
#        325.88353712, 325.8866185 , 342.558573  , 342.80082607,
#        343.032481  , 343.25401174, 343.46586999, 343.66848637,
#        343.86227151, 344.04761707, 344.22489665, 344.39446675,
#        319.7886221 , 319.6100216 , 319.439269  , 319.27602511,
#        319.1199654 , 318.97077937, 318.82816987, 318.69185258,
#        318.56155538, 318.43701784, 323.76730879, 323.68049039,
#        323.59634803, 323.51479868, 323.43576193, 323.35915991,
#        323.28491716, 323.21296063, 323.14321953, 323.07562532])

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

def darcyfriction(v,T):
    Re=p.Reynolds(v, T)
    a=1/(1+((Re/2712)**8.4))
    b=1/(1+((Re/(150*2*p.r/p.eff))**1.8))
    #a=1/(1+((Re/2.712)**8.4))
    #df=((64/Re)**a)*((0.75*np.log(Re/5.37))**(2*(a-1)*b))*((0.88*np.log(3.41*2*p.r/p.eff))**(2*(a-1)*(1-b)))
    df=((64/Re)**a)*((0.75*np.log(Re/5.37))**(2*(a-1)*b))*((0.88*np.log(6.82*2*p.r/p.eff))**(2*(a-1)*(1-b)))
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

def constantC(n,v,T):
    friction=darcyfriction(v, T)/blocks
    Cn=(-p.g*p.beta*np.sin(p.angle[n]))/(p.kw1+p.kw2+friction*p.length*((4*p.r)**-1))
    return Cn

#def constant2C(n,v,T):                          Darcy friction is actually temperature and so space dependend!!!!
#    Cn=(-p.g*p.beta*np.sin(p.angle[n]))/(p.kw1+p.kw2+darcyfrinction(v,T[n])*p.l*(p.r**-1))
#    return Cn


def T_l(l,v,Told):
    #constantes defenieren van de vergelijking voor de temperatuur
    deel2=((p.deel-0.25)/2)+0.25
    A1=constantA(round(p.N*0.10),Told[round(p.N*0.10)],v,0.10*p.length)
    #A2=constantA(round(p.N*deel2),Told[round(p.N*deel2)],v,deel2*p.length)
    A2=constantA(round(p.N*0.60),Told[round(p.N*0.60)],v,0.60*p.length)
    #A4=constantA(round(p.N*0.90),Told[round(p.N*0.90)],v,0.90*p.length)
    B1=constantB(round(p.N*0.10),Told[round(p.N*0.10)],v,0.10*p.length)
    #B2=constantB(round(p.N*deel2),Told[round(p.N*deel2)],v,deel2*p.length)
    B2=constantB(round(p.N*0.60),Told[round(p.N*0.60)],v,0.60*p.length)
    #B4=constantB(round(p.N*0.90),Told[round(p.N*0.90)],v,0.90*p.length)
    

    F1=B1/A1 #A1/B1
    #F2=B2/A2#A2/B2
    F2=B2/A2#A3/B3
    #F4=B4/A4#A4/B4
    #C=np.array((F1,F2,F3,F4))
    
    #lengtes definieren op verschillende punten  in de buis (randvoorwaarden op hoeken)
    l4=p.length
    
    #Constantes bepalen in vergelijking wat volgt uit de ODE en randvoorwaarden
    ans1=(np.e**(p.deel*A2*l4)-np.e**((A2*l4+p.deel*A1*l4)))**-1
    ans2=F2*(1-np.e**(p.deel*A1*l4)) + F1*(np.e**(p.deel*A1*l4)-1)

    constant2=(ans1*ans2)   
    constant1=-F2+F1+constant2*np.e**(A2*l4)
    
    #Bepalen in welk deel van de buis we zitten, dus welke constantes we moeten hebben
    if l<p.deel*p.length:
        n=round(p.N*0.10)
        constantn=constant1
        A=A1
        B=B1
        sign=1
    elif (l>p.deel*p.length or l==p.deel*p.length) and (l<p.length) :
        n=round(p.N*0.90)
        constantn=constant2
        A=A2
        B=B2
        sign=1
    elif l>p.length:
        print('error in dr_wall')
    n=round(l/p.length) 
    #het in elkaar zetten van de vergelijking voor de Temperatuur.
    #B=constantB(n,Told[n],v,l)
    #A=constantA(n,Told[n],v,l)
    
    T_l=-B/A+constantn*np.e**(A*l)
    
    return T_l, constantn

def velocity(v_i,Told):
    # #constantes bereken in de vergelijking voor de snelheid
    # deel2=((p.deel-0.25)/2)+0.25
    # C1=constantC(round(p.N*0.10),v_i,Told)
    # C2=constantC(round(p.N*deel2),v_i,Told)
    # C3=constantC(round(p.N*0.60),v_i,Told)
    # C4=constantC(round(p.N*0.90),v_i,Told)
    
    # #lengtes definieren in buis in elke deel van de buis
    # values1=np.linspace(0,0.25*p.length-(0.25*p.length/25),round(blocks*0.25))
    # values2=np.linspace(0.25*p.length,p.deel*p.length-(0.25*p.length/25),round(blocks*0.25))
    # values3=np.linspace(p.deel*p.length,0.75*p.length-(0.25*p.length/25),round(blocks*0.25))
    # values4=np.linspace(0.75*p.length,p.length-(0.25*p.length/25),round(blocks*0.25))
    
    deel2=((p.deel-0.25)/2)+0.25
    lengthpart2=p.deel-0.25
    lengthpart3=0.75-p.deel
    C1=constantC(round(p.N*0.10),v_i,Told)
    C2=constantC(round(p.N*deel2),v_i,Told)
    C3=constantC(round(p.N*0.60),v_i,Told)
    C4=constantC(round(p.N*0.90),v_i,Told)
    
    #lengtes definieren in buis in elke deel van de buis
    step1=(0.25*p.length/round(blocks*0.25))
    step2=(lengthpart2*p.length/round(blocks*0.25))
    step3=(lengthpart3*p.length/round(blocks*0.25))
    step4=(0.25*p.length/round(blocks*0.25))
    
    values1=np.linspace(0,0.25*p.length-step1,round(blocks*0.25))
    values2=np.linspace(0.25*p.length,p.deel*p.length-step2,round(blocks*0.25))
    values3=np.linspace(p.deel*p.length,0.75*p.length-step3,round(blocks*0.25))
    values4=np.linspace(0.75*p.length,p.length-step4,round(blocks*0.25))
    
    #Lege vectoren maken voor het berekenen voor elke Temperatuur in de verschillende delen 
    T_l_values1=np.zeros(round(blocks*0.25))
    T_l_values2=np.zeros(round(blocks*0.25))
    T_l_values3=np.zeros(round(blocks*0.25))
    T_l_values4=np.zeros(round(blocks*0.25))
    
    #Temperatuur berekenen deel 1 buis
    i=0
    j=0
    Co=np.zeros(blocks)
    for l in values1:
        T_l_values1[i]=T_l(l,v_i,Told)[0]
        Co[j]=T_l(l,v_i,Told)[1]
        j=j+1
        i=i+1
    #Temperatuur berekenen deel 2 buis   
    i=0
    for l in values2:
        T_l_values2[i]=T_l(l,v_i,Told)[0]
        Co[j]=T_l(l,v_i,Told)[1]
        j=j+1
        i=i+1
    #Temperatuur berekenen deel 3 buis  
    i=0
    for l in values3:
        T_l_values3[i]=T_l(l,v_i,Told)[0]
        Co[j]=T_l(l,v_i,Told)[1]
        j=j+1
        i=i+1
    #Temperatuur berekenen deel 4 buis    
    i=0
    for l in values4:
        T_l_values4[i]=T_l(l,v_i,Told)[0]
        Co[j]=T_l(l,v_i,Told)[1]
        j=j+1
        i=i+1
    T1new=np.append(T_l_values1,T_l_values2)
    T2new=np.append(T_l_values3,T_l_values4)
    Tnew=np.append(T1new,T2new)
    length12=np.append(values1,values2)
    length34=np.append(values3,values4)
    length=np.append(length12,length34)
    #uiteindelijk in elkaar zetten van de vergelijking voor de snelheid
    ans1=C1*sum(T_l_values1)*step1 #(p.length/blocks)  #+C1*(p.length/blocks)*(p.T_0+p.beta**-1)
    ans2=C2*sum(T_l_values2)*step2 #(p.length/blocks)  #+C2*(p.length/blocks)*(p.T_0+p.beta**-1)
    ans3=C3*sum(T_l_values3)*step3 #(p.length/blocks)  #+C3*(p.length/blocks)*(p.T_0+p.beta**-1)
    ans4=C4*sum(T_l_values4)*step4 #(p.length/blocks)  #+C4*(p.length/blocks)*(p.T_0+p.beta**-1)
    v=np.sqrt(abs(ans1+ans2+ans3+ans4))
    return v, Tnew, Co, length


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
#berekenen van initial snelheid aan de hand van the initial condition van de temperature. 
deel2=((p.deel-0.25)/2)+0.25
C1=constantC(round(p.N*0.10),vo,T_0) #!!!!!!!!!! hier nog even naar de velocity kijken 
C2=constantC(round(p.N*deel2),vo,T_0)
C3=constantC(round(p.N*0.60),vo,T_0)
C4=constantC(round(p.N*0.90),vo,T_0)
ans1=C1*sum(T_0[round(0):round(0.25*blocks)])*(p.length/blocks)
ans2=C2*sum(T_0[round(0.25*blocks):round(p.deel*blocks)])*(p.length/blocks)
ans3=C3*sum(T_0[round(p.deel*blocks):round(0.75*blocks)])*(p.length/blocks)
ans4=C4*sum(T_0[round(0.75*blocks):round(blocks)])*(p.length/blocks)
v_0=np.sqrt(abs(ans1+ans2+ans3+ans4))

# array maken voor velocity 
Velocityend=np.array([v_0])
i=0
vi=Velocityend[i]
vi_1= velocity(vi,T_0)[0]    #
Velocityend=np.append(Velocityend,vi_1)
Tnew=T_0
# Iteratieve methode om de snelheid te bepalen.
while abs(vi-vi_1)>0.000000003 and i<15:
#for iterations in tqdm(np.arange(0,15,1)):
    i=i+1
    vi=Velocityend[i]
    veli=velocity(vi,Tnew)
    vi_1=veli[0]
    Tnew=veli[1]
    C0=veli[2]
    length=veli[3]
    Velocityend=np.append(Velocityend,vi_1)
    
    
fig,((kx1,kx2)) = plt.subplots(1,2)

kx1.plot(Velocityend)
kx1.set_title('Velocity over iterations')
kx1.set_ylabel('$v$')
kx1.set_xlabel('$i$')

lengthvector=np.linspace(0,p.length,p.N)
kx2.plot(length,Tnew)
kx2.set_title('Temperature of fluid')
kx2.set_ylabel('$T$')
kx2.set_xlabel('$l$')
kx2.set_xlim(0, p.length)










