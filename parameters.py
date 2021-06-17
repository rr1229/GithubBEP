

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
r= 3*10**-3 #3*10**-3 (meters) radius of inner tube
dr1= 0.6*10**(-3)#7*10**-3  # (meters) thickness of wall of tube for part 1
dr2= 0.6*10**-3 #1*10**-3 (meters) thickness of wall of tube for part 4
dr3= 0.9*10**-3 #1*10**-3 (meters) thickness of wall of tube for part 3
dr4= 0.9*10**-3 #1*10**-3 (meters) thickness of wall of tube for part 4
deel=0.5 #where in the tube the smaller diameter will appear, works between 0.25 and 0.75



'general parameters'
g=9.81  #gravitational accelaration in the netherlands

'.......This Reactor........'
flux=3.5*10**16 # m^-2 s^-1 neutron flux in reactor, form Laurens haffmans

'......parameters of fluid.......'
C_pfluid= 41875 # (J/(kg K)) specific heat capacity of fluid (now taken water at a pressure 10^5 Pa at 60 degrees)
rho_0= 983.23  # (kg/m^3) reference density of WATER in pipe at temperature T_0, taken at atmospheric pressure (10^5 Pa)
T_0= 60+273.15 # reference temperature for fluid, now taken at 60 degrees
beta= 4.57*10**-4  # (K**-1) ORDERGROOTESCHATTING Thermal expansion coefficient of WATER (data compagnion) (constant for determining density at different temperature)
lambda_fluid=0.6506 #0.665 W/(m K)   #thermal conductivity of reference fluid, WATER at 60 degrees
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
solution=36/1000 #mol/kg maximum oplosbaarheid of salt in water
Na=6.022045*10**23 # avogadro's number (6.022×1023 atoms= 1 mol)
Sol_salt=Na*solution/rho_0 #atoms/m^3 maximum oplosbaarheid of salt in water
Cross_section_b=130*10**-3 #b (barn=1^-28) neutron cross section of molybdenum98   van bron: Can Enriched Molybdenum-98 Replace Enriched Uranium? Mushtaq Ahmad 
Cross_section=Cross_section_b *10**-28 #m^-2
Molair_Mo99=98.907707 #g/mol


'..........parameters of wall.........'
C_pwall= 285 # (J/(kg K)) specific heat capacity of wall (taken from Haffmans)
rho_wall= 6.55*10**3 # (kg/m^3) density of wall, taken from Haffmans
u= 300    # (W/kg) production term of gamma- heating by wall, Taken from Haffmans
lambda_wall=21.5 #thermal conductivity of wall

'..........parameters of surrounding water..........'
T_c= 40+273.15   # temperature outside water in Kelvin
mu_40=0.653*10**-3 #Pa s  dynamic viscosity of surrounding water at temperature 40 degrees
rho_water=992.25 #kg/m^3  density of surrounding water a temperature 40 degrees
lambda_water= 0.627 # W/mK aprosimated thermal conductivity of water at 40 degrees (from data companion)
Cp_water=4.1816*10**3 # J/(kg K) specific heat of water at temperature 40 degrees
nu_water=mu_40/rho_water
a_water=lambda_water/(rho_water*Cp_water)

'..........Geometrie.............' 
Ltube=140*10**-3 #width of outer tube in reactor, from Laurens haffmans
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
        
Opp_wallBC=np.zeros(N)  #area of outer wall of segment of system
for i in np.arange(0,N,1):
    Opp_wallBC[i]=(length/N)*2*np.pi*(r+dr[i])    

        
Vwall=np.zeros(N)
for i in np.arange(0,N,1):# volume of wall of segment of system 
    Vwall[i]=(length/N)*np.pi*((r+dr[i])**2-r**2)

#h_A= 382.7  #ORDERGROOTESHCATTING # heat transfer coefficient of fluid
#h_B= 13115  #ORDERGROOTESHCATTING    # heat transfer coefficient of wall
#h_C= 382.7  #ORDERGROOTESHCATTING   # heat transfer coefficient of water

#h_AB= ((1/h_A)+(1/h_B))**-1  # heat transfer coefficient of fluid to wall
#h_BC= ((1/h_B)+(1/h_C))**-1  # heat transfer coefficient of wall to water



'____________heat producion in pipe___________'

Pgamma=np.zeros(N)
for n in np.arange(0,N,1):
    Pgamma[n]=u*rho_wall*Vwall[n]        # heat produced in each wall segment of the pipe
    
'___________________Mo98 production_______________________'
def Reaction():
    N= Sol_salt#atoms Mo98 per volume
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
            #h=270   #approximation for this regime
            #print('error in creation hfluid: Re= %.3e , Gz= %.3e and Pr= %.3e' %(Re, Gz, Pr))
        
    return h

def h_outside(n,v,T,Tb,m):
    #Re=Reynolds(v, T[n])
    Gr=Grashof(m,Tb)
    Pr=nu_water/a_water
    Ra=Gr*Pr
    if N/4<=n<N/2 or 3*N/4<=n<N: #vertical pipe
        Nu_ans1=(4/3)* (Ra**0.25) * (7*Pr/(100+105*Pr))**0.25 
        Nu_ans2=(4/35)* ((272+315*Pr)/(64+63*Pr))*(length/(4*2*(r+dr[n])))
        Nu=Nu_ans1+Nu_ans2
        
        if Gr<10**8 and Ra>10**4 and Ra<10**9: #:#Nureth Lin Xiana, b, Guangming Jianga, b, Hongxing Yua, b
            Nu=0.48*Ra**0.25 #!!!!!!!!!!!!!hier nog even goed naar kijken!!!
            
        #elif Gr>4*10**9 and Ra<10**11 and Ra>10**10:# :     #komt niet in deze region
        #    Nu=0.148*Ra**0.333
            
            #print('value Gr is smaller than 10^8, so falls in wrong regime Gr= %.3e' %(Gr))
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
v0=0.000000000000000000000001#0.00000000000000000001
T0=(20+273.15)*np.ones(N)
Tb0=(20+273.15)*np.ones(N)


'_______________Initial guess_______________________'
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

v_steadystate0=0.0024382357215700025
T_steadystate0=np.array([321.04991726, 321.135062  , 321.21787652, 321.29842567,
       321.37677203, 321.45297481, 321.52709019, 321.59917234,
       321.66927426, 321.73744847, 321.89634732, 322.05104885,
       322.20166411, 322.34830161, 322.49106717, 322.63006372,
       322.76539118, 322.8971466 , 323.02542426, 323.1503158 ,
       322.97663356, 322.80709401, 322.64159672, 322.48004332,
       322.3223371 , 322.1683835 , 322.01809034, 321.87136784,
       321.72812865, 321.5882878 , 321.51622683, 321.4464261 ,
       321.38004648, 321.3156411 , 321.25309461, 321.19245877,
       321.13330484, 321.07527462, 321.01829334, 320.96237686])
Tb_steadystate0=np.array([323.82182495, 323.83092064, 323.83976568, 323.84836698,
       323.85673115, 323.86486446, 323.87277304, 323.88046289,
       323.88793994, 323.89521008, 326.92687129, 326.94882399,
       326.97019043, 326.99098651, 327.01122771, 327.03092901,
       327.05010499, 327.06876976, 327.08693707, 327.10462026,
       317.62167487, 317.5798693 , 317.53898203, 317.49899338,
       317.4598839 , 317.42163451, 317.38422655, 317.34764178,
       317.31186237, 317.27687093, 319.34279385, 319.26365335,
       319.32451786, 319.3333531 , 319.35077446, 319.331713  ,
       319.30782417, 319.28434534, 319.26128467, 319.23864859])

#v_steadystate0=0.0009422762947661532
#T_steadystate0=np.array([320.57192048, 320.88700595, 321.18361936, 321.4628414 ,
#         321.7256897 , 321.97312268, 322.20604312, 322.42530148,
#         322.63169895, 322.82599018, 323.21060477, 323.57339125,
#         323.91558769, 324.23836189, 324.54281527, 324.82998664,
#         325.10085583, 325.35634708, 325.59733232, 325.82463412,
#         325.27060295, 324.74768182, 324.2540814 , 323.78811856,
#         323.34820993, 322.93286583, 322.54068453, 322.17034689,
#         321.82061129, 321.49030897, 321.33514334, 321.18941691,
#         321.05028279, 320.91743668, 320.79059261, 320.66947574,
#         320.55382395, 320.44338716, 320.33792679, 320.23721512])
#Tb_steadystate0=np.array([325.35397596, 325.38874528, 325.42145844, 325.45223781,
#         325.48119839, 325.50844828, 325.53408912, 325.55821651,
#         325.58092034, 325.60228515, 328.70600763, 328.75693842,
#         328.80495324, 328.85022033, 328.89289804, 328.93313543,
#         328.97107286, 329.00684249, 329.04056882, 329.07236911,
#         318.1319643 , 318.00992085, 317.89414657, 317.78432665,
#         317.68016179, 317.58136751, 317.48767335, 317.39882224,
#         317.31456979, 317.23468367, 319.3480959 , 319.31210228,
#         319.25790354, 319.20608794, 319.15655161, 319.10919501,
#         319.06392284, 319.02064387, 318.97927069, 318.93971961])


















