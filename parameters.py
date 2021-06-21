

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
'=====================changeble parameters========================'
N= 40 #number of segments, has to be multiple of 4
anglepipe=5 * np.pi*2/360  # (radialen)  angle of horizontal pipes in system with horizontal, now 5 degrees 
r= 3*10**-3 #15*10**-3 # (meters) radius of inner tube
dr1= 7*10**-3 #20*10**-3 # (meters) thickness of wall of tube for part 1
dr2= 7*10**-3#20*10**-3 #  (meters) thickness of wall of tube for part 4
dr3= 1*10**-3 #1*10**-3 # (meters) thickness of wall of tube for part 3
dr4= 1*10**-3 #1*10**-3 # (meters) thickness of wall of tube for part 4
deel=0.5 #0.37 #0.5#where in the tube the smaller diameter will appear, works between 0.25 and 0.75



'general parameters'
g=9.81  #gravitational accelaration in the netherlands

'.......This Reactor........'
flux=3.5*10**16 # m^-2 s^-1 neutron flux in reactor, form Laurens haffmans

'......parameters of fluid.......'
C_pfluid=3327# (J/(kg C)) specific heat capacity of fluid now taken for NaCl solution at 52 degrees
rho_0=1147 #(kg/m^3) reference density of NaCl solution at 52 C
T_0=52+273.15 # K reference temp. for NaCl solutuion taken 52 degrees
lambda_fluid=0.530 #(W/(m C)) thermal conductivity of NaCl solution at 52 degrees

# C_pfluid= 41875 # (J/(kg K)) specific heat capacity of fluid (now taken water at a pressure 10^5 Pa at 60 degrees)
# rho_0= 983.23  # (kg/m^3) reference density of WATER in pipe at temperature T_0, taken at atmospheric pressure (10^5 Pa)
# T_0= 60+273.15 # reference temperature for fluid, now taken at 60 degrees
# lambda_fluid=0.6506 #0.665 W/(m K)   #thermal conductivity of reference fluid, WATER at 60 degrees

beta= 4.57*10**-4  # (K**-1) Thermal expansion coefficient of WATER 
mu_20=1.002*10**-3 #Pa s  dynamic viscosity of reference fluid, water, at 20 degrees
a_fluid=lambda_fluid/(rho_0*C_pfluid)     #thermal diffusivity of the fluid, now estimated for water now

def mu_fluid(T):             #dynamic viscosity of reference fluid, NaCl dissolved at temperature T
    mu_water=mu_20*np.e**((1.1709*(20+273.15-T)-0.001827*(T-20)**2)/(T+89.93))
    Molair=58.44 #g/mol molair mass of NaCl
    m=solution/Molair
    A=0.008-(0.005/60)*(T-293.15)
    B=0.06+(0.06/60)*(T-293.15)
    mu_rel=1+(A*m**0.5)+B*m  # see bron Viscosity of Aqueous Solutions of Sodium Chloride
                            #A. A. Aleksandrov, E. V. Dzhuraeva, and V. F. Utenkov
    mu=mu_water*mu_rel
    return mu 

def nu_fluid(T):
    nu=mu_fluid(T)/rho_0
    return nu
'production Mo99'
#solution=36/1000 #mol/kg maximum oplosbaarheid of salt in water
solution=0.84 #kg/L = 84 g/100ml from pubchem solubility for Sodium molybdate
Molair_mass_mosalt=207.864601 *10**3 # kg/mol
MolMo_98=solution/Molair_mass_mosalt #mol/L Mo98 in water
Na=6.022045*10**23 # avogadro's number (6.022×1023 atoms= 1 mol)
Sol_Mo98=Na*MolMo_98*10**3 # atoms/m^3 = Na*MolMo_98 atoms/L, max solution Mo98 in water
Sol_salt=Na*solution/rho_0 #atoms/m^3 maximum oplosbaarheid of salt in water
Cross_section_b=130*10**-3 #b (barn=1^-28) neutron cross section of molybdenum98   
Cross_section=Cross_section_b *10**-28 #m^-2
Molair_Mo99=98.907707 #g/mol


'..........parameters of wall.........'
C_pwall= 285 # (J/(kg K)) specific heat capacity of wall 
rho_wall= 6.55*10**3 # (kg/m^3) density of wall
u= 300    # (W/kg) production term of gamma- heating by wall
lambda_wall=21.5 #thermal conductivity of wall

'..........parameters of surrounding water..........'
T_c= 40+273.15   # temperature outside water in Kelvin
mu_40=0.653*10**-3 #Pa s  dynamic viscosity of surrounding water at temperature 40 degrees
rho_water=992.25 #kg/m^3  density of surrounding water a temperature 40 degrees
lambda_water= 0.627 # W/mK aprosimated thermal conductivity of water at 40 degrees 
Cp_water=4.1816*10**3 # J/(kg K) specific heat of water at temperature 40 degrees
nu_water=mu_40/rho_water
a_water=lambda_water/(rho_water*Cp_water)

'..........Geometrie.............' 
Ltube=140*10**-3 #width of outer tube in reactor
length=(4*(1+np.sin(anglepipe))**-1)*(Ltube-dr1-dr3-2*r) #hier gaan we ervan uit dat in de 1e buis de dikte altijd dr1 is en in de 3e buis altijd dr2
dr=np.zeros(N)
for i in np.arange(0,N,1):
    if i<0.25*N:
        dr[i]=dr1
    elif (i>0.25*N or i==0.25*N) and (i<deel*N) :
        dr[i]=dr2
    elif (i>deel*N or i==deel*N) and (i<0.75*N) :
        dr[i]=dr3
    elif (i>0.75*N or i==0.75*N) and (i<N) :
        dr[i]=dr4


eff=1.5*(10**(-6)) #effective roughness of tube, 

tVsys= length*np.pi*r**2 #volume of whole system within the pipe
Vsys_segm=(length/N)*np.pi*r**2 # volume of segment within the pipe
Opp_wallAB=(length/N)*r*2*np.pi # area of inner wall of segment of system

kw1=1.30   # pressure loss term for bend 1, for now a schatting from a sharp bend
kw2=1.30   # pressure loss term for bend 2, for now a schatting from a sharp bend


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
        
Opp_wallBC=np.zeros(N)  #area of outer wall of segment of system
for i in np.arange(0,N,1):
    Opp_wallBC[i]=(length/N)*2*np.pi*(r+dr[i])    

        
Vwall=np.zeros(N)
for i in np.arange(0,N,1):# volume of wall of segment of system 
    Vwall[i]=(length/N)*np.pi*((r+dr[i])**2-r**2)


'____________heat producion in pipe___________'

Pgamma=np.zeros(N)
for n in np.arange(0,N,1):
    Pgamma[n]=u*rho_wall*Vwall[n]        # heat produced in each wall segment of the pipe
    
'___________________Mo98 production_______________________'
def Reaction():
    N= Sol_Mo98#atoms Mo98 per volume
    R=Cross_section*flux*N  #atoms mo99 per volume
    R_mass=R*Molair_Mo99/Na #g per volume
    return R_mass*tVsys #g mo99 per second


'_______________dimensionless numbers________________'
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

def Grashof(n,Tb):
    if N/4<=n<N/2 or 3*N/4<=n<N: #vertical pipe
        ans1=(g*beta*abs(Tb[n]-T_c)*(length/4)**3)/(nu_water**2)
    if 0<=n<N/4 or N/2<=n<3*N/4: #horizontal pipe
        ans1=(g*beta*abs(Tb[n]-T_c)*(2*(r+dr[n]))**3)/(nu_water**2)
    return ans1

'_________________Heat transfer coefficients_________________________'
def h_wall(n):
    if n<0.25*N:
        d=dr[n]
    elif (n>0.25*N or n==0.25*N) and (n<deel*N) :
        d=dr[n]
    elif (n>deel*N or n==deel*N) and (n<0.75*N) :
        d=dr[n]
    elif (n>0.75*N or n==0.75*N) and (n<N) :
        d=dr[n]
    #h=lambda_wall/d
    h=4.36*lambda_wall/d
    return h

def h_fluid(v,T):
    Re=Reynolds(v, T)
    Gz=Greatz(v)
    Pr=Prandtl(T)  
    
    if True:
    #for i in range(np.size(h)):
        if Re>10**4 and Pr>=0.7:
            h=0.027*(Re**0.8)*(Pr**0.33)*lambda_fluid/(2*r)
        elif Gz<0.05:
            h=1.62*(1/np.cbrt(Gz))*lambda_fluid/(2*r) #np.cbrt = 3e machtswortel
        elif Gz>0.1:
            h=3.66*lambda_fluid/(2*r)
        else:
            h1=1.62*(1/np.cbrt(0.05))*lambda_fluid/(2*r)
            h2=3.66*lambda_fluid/(2*r)
            a=(h2-h1)/0.05
            h=a*(Gz-0.05)+h1
            
    return h

def h_outside(n,v,T,Tb,m):
    Gr=Grashof(m,Tb)
    Pr=nu_water/a_water
    Ra=Gr*Pr
    if N/4<=n<N/2 or 3*N/4<=n<N: #vertical pipe
        Nu_ans1=(4/3)* (Ra**0.25) * (7*Pr/(100+105*Pr))**0.25 
        Nu_ans2=(4/35)* ((272+315*Pr)/(64+63*Pr))*(length/(4*2*(r+dr[n])))
        Nu=Nu_ans1+Nu_ans2
        
        if Gr<10**8 and Ra>10**4 and Ra<10**9: 
            Nu=0.48*Ra**0.25

        h_c=Nu*lambda_water/(length/4)
    
    elif 0<=n<N/4 or N/2<=n<3*N/4: #horizontal pipe
        Nu_ans1=0.387*Ra**(1/6)
        Nu_ans2=( 1+ ((0.559**(9/16))/Pr) )**(8/27)
        Nu=(0.6+(Nu_ans1/Nu_ans2))**2
        
        h_c=Nu*lambda_water/(2*(r+dr[n]))
        #h_c=
        if Ra>10**12 or Ra<10**-5:
            print('value Ra not in right regime within h_outside horizontal')
    else:
        print(n)
    return h_c
    

def h_AB(n,v,T):
    h_AB= ((1/h_fluid(v, T))+(1/h_wall(n)))**-1 
    return h_AB

def h_BC(n,v,T,Tb):
    h_BC= ((1/h_wall(n))+(1/h_outside(n, v, T, Tb, n)))**-1 
    return h_BC



'____________initial conditions____________'
v0=0.000000000000000000000001
T0=(20+273.15)*np.ones(N)
Tb0=(20+273.15)*np.ones(N)


'_______________Initial guess_______________________'
v_steadystate0=0.00115
T_steadystate0=np.zeros(N)
Tb_steadystate0=np.zeros(N)

for n in np.arange(0,N,1):
    if n<0.25*N:
        T_steadystate0[n]=(50+273.15)
        Tb_steadystate0[n]=(60+273.15)
    elif (n>0.25*N or n==0.25*N) and (n<deel*N) :
        T_steadystate0[n]=(53+273.15)
        Tb_steadystate0[n]=(60+273.15)
    elif (n>deel*N or n==deel*N) and (n<0.75*N) :
        T_steadystate0[n]=(52+273.15)
        Tb_steadystate0[n]=(52+273.15)
    elif (n>0.75*N or n==0.75*N) and (n<N) :
        T_steadystate0[n]=(49+273.15)
        Tb_steadystate0[n]= (50+273.15)
    elif n>N:
        print('error in temp_0 creation')

# v_steadystate0=0.003221123130911058 
# T_steadystate0=np.array([323.67259372, 323.8151801 , 323.955052  , 324.09226461,
#         324.22686818, 324.35891259, 324.48844662, 324.61551796,
#         324.74017337, 324.86245877, 325.02576036, 325.18601826,
#         325.34328905, 325.49762833, 325.64909067, 325.79772964,
#         325.9435978 , 326.08674675, 326.2272271 , 326.36508856,
#         326.17199886, 325.98251249, 325.79656068, 325.61407551,
#         325.43499034, 325.25923991, 325.08676024, 324.91748865,
#         324.75136369, 324.58832514, 324.47001233, 324.35507973,
#         324.24263452, 324.13323156, 324.02720411, 323.92300349,
#         323.82029242, 323.71989032, 323.62251058, 323.52738638])
# Tb_steadystate0=np.array([330.4392481 , 330.45319331, 330.46687147, 330.48028774,
#         330.49344716, 330.50635471, 330.51901521, 330.53143344,
#         330.54361405, 330.55556162, 332.49484976, 332.51589233,
#         332.53653916, 332.55679776, 332.57667548, 332.59617954,
#         332.61531701, 332.63409483, 332.65251982, 332.67059864,
#         318.6038022 , 318.55557156, 318.50815763, 318.46154677,
#         318.4157256 , 318.37068094, 318.32639987, 318.28286963,
#         318.24007771, 318.1980118 , 319.83781147, 319.84949599,
#         319.86942196, 319.89561226, 319.90501333, 319.79924521,
#         319.76140237, 319.81372857, 319.76409232, 319.77831997])
















