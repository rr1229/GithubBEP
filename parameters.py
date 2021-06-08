

# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 17:03:47 2021

@author: ronar
"""
import numpy as np
import matplotlib.pyplot as plt

#Like Dresen & Haffmans we will use zircaloy (zirconium alloy) for the wall.

#we will first test the code with water properties as fluid inside the pipe


"-----------------------Variables & Parameters----------------------------"

g=9.81  #gravitational accelaration in the netherlands
C_pfluid= 41875 # (J/(kg K)) specific heat capacity of fluid (now taken water at a pressure 10^5 Pa at 60 degrees)
C_pwall= 285 # (J/(kg K)) specific heat capacity of wall (taken from Haffmans)
T_c= 40+273.15   # temperature outside water in Kelvin

'Geometrie' 
N= 40 #number of segments, has to be multiple of 4
length= 4*110.4*10**-3 #length of entire system (each segment is 1/4 length)
anglepipe=np.pi*2*5/360  # (radialen)  angle of horizontal pipes in system with horizontal, now 5 degrees 
r= 3*10**-3 #3*10**-3 (meters) radius of inner tube
dr1= 6*10**(-3)#7*10**-3  # (meters) thickness of wall of tube for part 1
dr2= 2*10**-3 #1*10**-3 (meters) thickness of wall of tube for part 2

eff=1.5*(10**(-6)) #effective roughness of tube, 

tVsys= length*np.pi*r**2 #volume of whole system within the pipe
Vsys_segm=(length/N)*np.pi*r**2
Vwallsegm1=(length/N)*np.pi*((r+dr1)**2-r**2) # volume of wall of segment of system in part 1
Vwallsegm2=(length/N)*np.pi*((r+dr2)**2-r**2) # volume of wall of segment of system in part 2
Opp_wallAB=(length/N)*r*2*np.pi # area of inner wall of segment of system
Opp_wallBC1=(length/N)*2*np.pi*((r+dr1)) #area of outer wall of segment of system
Opp_wallBC2=(length/N)*2*np.pi*((r+dr2)) #area of outer wall of segment of system

kw1=1.30   # pressure loss term for bend 1, for now a schatting from a sharp bend
kw2=1.30   # pressure loss term for bend 2, for now a schatting from a sharp bend
f= 0.2#0.05#0.2  #friction term, ordergroote

angle=np.zeros(N)
for i in np.arange(0,N,1):
    if i<0.25*N:
        angle[i]=-anglepipe
    elif (i>0.25*N or i==0.25*N) and (i<0.5*N) :
        angle[i]=-0.5*np.pi
    elif (i>0.5*N or i==0.5*N) and (i<0.75*N) :
        angle[i]=anglepipe
    elif (i>0.75*N or i==0.75*N) and (i<N) :
        angle[i]=0.5*np.pi
    elif i>N:
        print('error in angle creation')
        
Opp_wallBC=np.zeros(N)
for i in np.arange(0,N,1):
    if i<0.25*N:
        Opp_wallBC[i]=Opp_wallBC1
    elif (i>0.25*N or i==0.25*N) and (i<0.5*N):
        Opp_wallBC[i]=Opp_wallBC1
    elif (i>0.5*N or i==0.5*N):
        Opp_wallBC[i]=Opp_wallBC2        
    elif i>N:
        print('error in Opp_wallBC creation')
        
Vwall=np.zeros(N)
for i in np.arange(0,N,1):
    if i<0.25*N:
        Vwall[i]=Vwallsegm1
    elif (i>0.25*N or i==0.25*N) and (i<0.5*N):
        Vwall[i]=Vwallsegm1
    elif (i>0.5*N or i==0.5*N):
        Vwall[i]=Vwallsegm2
    elif i>N:
        print('error in Vwall creation')


'material properties'
rho_0= 983.23  # (kg/m^3) reference density of fluid in pipe at temperature T_0, taken at atmospheric pressure (10^5 Pa)
T_0= 60+273.15 # reference temperature for fluid, now taken at 60 degrees
beta= 4.57*10**-4  # (K**-1) ORDERGROOTESCHATTING Thermal expansion coefficient of fluid (constant for determining density at different temperature)
lambda_fluid=0.665 #thermal conductivity of reference fluid, water at 60 degrees
mu_20=1.002*10**-3 #dynamic viscosity of reference fluid, water, at 20 degrees

lambda_wall=21.5

def mu_fluid(T):             #dynamic viscosity of reference fluid, water at temperature T
    mu=mu_20*np.e**((1.1709*(20+273.15-T)-0.001827*(T-20)**2)/(T+89.93))
    return mu 

a_fluid=lambda_fluid/(rho_0*C_pfluid)     #thermal diffusivity of the fluid, now estimated for water now

def nu_fluid(T):
    nu=mu_fluid(T)/rho_0
    return nu

h_A= 382.7  #ORDERGROOTESHCATTING # heat transfer coefficient of fluid
h_B= 13115  #ORDERGROOTESHCATTING    # heat transfer coefficient of wall
h_C= 382.7  #ORDERGROOTESHCATTING   # heat transfer coefficient of water

h_AB= ((1/h_A)+(1/h_B))**-1  # heat transfer coefficient of fluid to wall
h_BC= ((1/h_B)+(1/h_C))**-1  # heat transfer coefficient of wall to water

rho_wall= 6.55*10**3 # (kg/m^3) density of wall, taken from Haffmans
u= 300    # (W/kg) production term of gamma- heating by wall, Taken from Haffmans

'heat producion in pipe'

Pgamma=np.zeros(N)
for n in np.arange(0,N,1):
    Pgamma[n]=u*rho_wall*Vwall[n]        # heat produced in each wall segment of the pipe


'initial conditions'
v0=0.000000000000000000000001#0.00000000000000000001
T0=(20+273.15)*np.ones(N)
Tb0=(20+273.15)*np.ones(N)

'dimensionless numbers'

def Reynolds(v,T):
    mu=mu_20*np.e**((1.1709*(20+273.15-T)-0.001827*(T-20)**2)/(T+89.93))
    re=rho_0*v*2*r*(1/mu)
    return abs(re) #re

def Greatz(v):
    Gz=a_fluid*length/(abs(v)*((2*r)**2))
    return Gz

def Prandtl(T):
    Pr=nu_fluid(T)/a_fluid
    return Pr

def h_wall(n):
    if n<0.25*N:
        d=dr1
    elif (n>0.25*N or n==0.25*N) and (n<0.5*N) :
        d=dr1
    elif (n>0.5*N or n==0.5*N) and (n<0.75*N) :
        d=dr2
    elif (n>0.75*N or n==0.75*N) and (n<N) :
        d=dr2
    h=lambda_wall/d
    return h

def h_fluid(v,T):
    Re=Reynolds(v, T)
    Gz=Greatz(v)
    Pr=Prandtl(T)
    h= 382.7   # hier nog naar kijken wat te doen als geen van deze voldaan worden
    #h= np.zeros(Re.shape)

    
    if True:
    #for i in range(np.size(h)):
        if Re>10**4 and Pr>=0.7:
            h=0.027*(Re**0.8)*(Pr**0.33)*lambda_fluid/(2*r)
        elif Gz<0.05:
            h=1.62*(1/np.cbrt(Gz))*lambda_fluid/(2*r) #np.cbrt = 3e machtswortel
        elif Gz>0.1:
            h=3.66*lambda_fluid/(2*r)
    #else:
        #print('error in creation of h_fluid')
        
    return h

def h_AB(n,v,T):
    h_AB= ((1/h_fluid(v, T))+(1/h_wall(n)))**-1 
    return h_AB

def h_BC(n):
    h_BC= ((1/h_wall(n))+(1/h_C))**-1 
    return h_BC

v_steadystate0=0.00115#1.0*10**-8#7.07*10**-13#0.03150#5*10**-5
T_steadystate0=np.zeros(N)
Tb_steadystate0=np.zeros(N)
# for n in np.arange(0,N,1):
#     if n<0.25*N:
#         T_steadystate0[n]=(70+273.15)
#         Tb_steadystate0[n]=(72+273.15)
#     elif (n>0.25*N or n==0.25*N) and (n<0.5*N) :
#         T_steadystate0[n]=(60+273.15)
#         Tb_steadystate0[n]=(62+273.15)
#     elif (n>0.5*N or n==0.5*N) and (n<0.75*N) :
#         T_steadystate0[n]=(59+273.15)
#         Tb_steadystate0[n]=(55+273.15)
#     elif (n>0.75*N or n==0.75*N) and (n<N) :
#         T_steadystate0[n]=(57+273.15)
#         Tb_steadystate0[n]= (55+273.15)
#     elif n>N:
#         print('error in angle creation')


for n in np.arange(0,N,1):
    if n<0.25*N:
        T_steadystate0[n]=(55+273.15)
        Tb_steadystate0[n]=(60+273.15)
    elif (n>0.25*N or n==0.25*N) and (n<0.5*N) :
        T_steadystate0[n]=(58+273.15)
        Tb_steadystate0[n]=(60+273.15)
    elif (n>0.5*N or n==0.5*N) and (n<0.75*N) :
        T_steadystate0[n]=(57+273.15)
        Tb_steadystate0[n]=(52+273.15)
    elif (n>0.75*N or n==0.75*N) and (n<N) :
        T_steadystate0[n]=(54+273.15)
        Tb_steadystate0[n]= (50+273.15)
    elif n>N:
        print('error in temp_0 creation')



















