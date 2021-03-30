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
l= 4*110.4*10**-3 #length of entire system (each segment is 1/4 length)
anglepipe=np.pi*2*5/360  # (radialen)  angle of horizontal pipes in system with horizontal, now 5 degrees 
r= 3*10**-3 # (meters) radius of inner tube
dr1= 7*10**-3  # (meters) thickness of wall of tube for part 1
dr2= 2*10**-3 # (meters) thickness of wall of tube for part 2

tVsys= l*np.pi*r**2 #volume of whole system within the pipe
Vsys_segm=(l/N)*np.pi*r**2
Vwallsegm1=(l/N)*np.pi*((r+dr1)**2-r**2) # volume of wall of segment of system in part 1
Vwallsegm2=(l/N)*np.pi*((r+dr2)**2-r**2) # volume of wall of segment of system in part 2
Opp_wallAB=(l/N)*r*2*np.pi # area of inner wall of segment of system
Opp_wallBC1=(l/N)*2*np.pi*((r+dr1)) #area of outer wall of segment of system
Opp_wallBC2=(l/N)*2*np.pi*((r+dr2)) #area of outer wall of segment of system

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
    elif (i>0.25*N or i==0.25*N):
        Opp_wallBC[i]=Opp_wallBC2
    elif i>N:
        print('error in Opp_wallBC creation')
        
Vwall=np.zeros(N)
for i in np.arange(0,N,1):
    if i<0.25*N:
        Vwall[i]=Vwallsegm1
    elif (i>0.25*N or i==0.25*N):
        Vwall[i]=Vwallsegm2
    elif i>N:
        print('error in Vwall creation')


'material properties'
rho_0= 983.23  # (kg/m^3) reference density of fluid in pipe at temperature T_0, taken at atmospheric pressure (10^5 Pa)
T_0= 60+273.15 # reference temperature for fluid, now taken at 60 degrees
beta= 4.57*10**-4  # (K**-1) ORDERGROOTESCHATTING Thermal expansion coefficient of fluid (constant for determining density at different temperature)

h_A= 382.7  #ORDERGROOTESHCATTING # heat transfer coefficient of fluid
h_B= 13115  #ORDERGROOTESHCATTING    # heat transfer coefficient of wall
h_C= 382.7  #ORDERGROOTESHCATTING   # heat transfer coefficient of water

h_AB= ((1/h_A)+(1/h_B))**-1  # heat transfer coefficient of fluid to wall
h_BC= ((1/h_B)+(1/h_C))**-1  # heat transfer coefficient of wall to water

rho_wall= 6.55*10**3 # (kg/m^3) density of wall, taken from Haffmans
u= 300    # (W/kg) production term of gamma- heating by wall, Taken from Haffmans



'initial conditions'
v0=0
T0=(20+273.15)*np.ones(N)
Tb0=(20+273.15)*np.ones(N)























