# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 16:47:31 2021

@author: ronar
"""

import parameters as p
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import tqdm as tqdm
plt.close('all')



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


'Varieren van diktes van de buis'
#begindr1=1
#enddr1=90.01
#begindr3=0.5
#enddr3=40.1
#step1=2
#step3=2

#listnamev=['v_(dr1_0=%.2f,dr1_end=%.2f,step1=%.2f,dr3_0=%.2f,dr3_end=%.2f,step3=%.2f).npy' %(begindr1,enddr1,step1,begindr3,enddr3,step3)]
#listnameTn=['Tn_(dr1_0=%.2f,dr1_end=%.2f,step1=%.2f,dr3_0=%.2f,dr3_end=%.2f,step3=%.2f).npy' %(begindr1,enddr1,step1,begindr3,enddr3,step3)]
#listnameTbn=['Tbn_(dr1_0=%.2f,dr1_end=%.2f,step1=%.2f,dr3_0=%.2f,dr3_end=%.2f,step3=%.2f).npy' %(begindr1,enddr1,step1,begindr3,enddr3,step3)]
#listnamemaxT=['maxtemperature_(dr1_0=%.2f,dr1_end=%.2f,step1=%.2f,dr3_0=%.2f,dr3_end=%.2f,step3=%.2f).npy' %(begindr1,enddr1,step1,begindr3,enddr3,step3)]

#r=3*10**-3
#deel=0.5
#vnsend=np.zeros((len(np.arange(begindr1,enddr1,step1)),len(np.arange(begindr3,enddr3,step3))))
#Tnsend=np.zeros((p.N,len(np.arange(begindr1,enddr1,step1)),len(np.arange(begindr3,enddr3,step3))))
#Tnsbend=np.zeros((p.N,len(np.arange(begindr1,enddr1,step1)),len(np.arange(begindr3,enddr3,step3))))
#ite=0
#drbelow=np.arange(begindr1,enddr1,step1)
#drup=np.arange(begindr3,enddr3,step3)



'variation of R and dr1'
begindr1=10
enddr1=30.01
step1=0.5
beginr=0.25
endr=0.5
stepr=0.005



#dr3=1*10**-3
#dr4=1*10**-3
#deel=0.5
listnamev=['v_(dr1_0=%.2f,dr1_end=%.2f,step1=%.2f,deel_0=%.2f,deel_end=%.2f,stepdeel=%.2f).npy' %(begindr1,enddr1,step1,beginr,endr,stepr)]
listnameTn=['Tn_(dr1_0=%.2f,dr1_end=%.2f,step1=%.2f,deel_0=%.2f,deel_end=%.2f,stepdeel=%.2f).npy' %(begindr1,enddr1,step1,beginr,endr,stepr)]
listnameTbn=['Tbn_(dr1_0=%.2f,dr1_end=%.2f,step1=%.2f,deel_0=%.2f,deel_end=%.2f,stepdeel=%.2f).npy' %(begindr1,enddr1,step1,beginr,endr,stepr)]
listnamemaxT=['maxtemperature_(dr1_0=%.2f,dr1_end=%.2f,step1=%.2f,deel_0=%.2f,deel_end=%.2f,stepdeel=%.2f).npy' %(begindr1,enddr1,step1,beginr,endr,stepr)]

vnsend=np.zeros((len(np.arange(begindr1,enddr1,step1)),len(np.arange(beginr,endr,stepr))))
Tnsend=np.zeros((p.N,len(np.arange(begindr1,enddr1,step1)),len(np.arange(beginr,endr,stepr))))
Tnsbend=np.zeros((p.N,len(np.arange(begindr1,enddr1,step1)),len(np.arange(beginr,endr,stepr))))
ite=0
drbelow=np.arange(begindr1,enddr1,step1)
rtube=np.arange(beginr,endr,stepr)


for thick1 in drbelow:
    dr1=thick1*10**-3
    dr2=thick1*10**-3
    j=0
    for thickr in rtube:
        deel=thickr
        
        
        "-----------------------Variables & Parameters----------------------------"
        '=====================changeble parameters========================'
        N= 40 #number of segments, has to be multiple of 4
        anglepipe=5 * np.pi*2/360  # (radialen)  angle of horizontal pipes in system with horizontal, now 5 degrees 
        r= 15*10**-3#3*10**-3 #3*10**-3 (meters) radius of inner tube
        dr1= dr1#20*10**-3#7*10**(-3)#7*10**-3  # (meters) thickness of wall of tube for part 1
        dr2= dr2 #20*10**-3#7*10**-3 #1*10**-3 (meters) thickness of wall of tube for part 4
        dr3= 1*10**-3 #1*10**-3 (meters) thickness of wall of tube for part 3
        dr4= 1*10**-3 #1*10**-3 (meters) thickness of wall of tube for part 4
        #deel=0.5 #where in the tube the smaller diameter will appear, works between 0.25 and 0.75
        
        
        
        'general parameters'
        g=9.81  #gravitational accelaration in the netherlands
        
        '.......This Reactor........'
        flux=3.5*10**16 # m^-2 s^-1 neutron flux in reactor, form Laurens haffmans
        
        '......parameters of fluid.......'
        #C_pfluid=3327# (J/(kg C)) specific heat capacity of fluid now taken for NaCl solution at 52 degrees
        #rho_0=1147 #(kg/m^3) reference density of NaCl solution at 52 C
        #T_0=52+273.15 # K reference temp. for NaCl solutuion taken 52 degrees
        #lambda_fluid=0.530 #(W/(m C)) thermal conductivity of NaCl solution at 52 degrees
        C_pfluid= 41875 # (J/(kg K)) specific heat capacity of fluid (now taken water at a pressure 10^5 Pa at 60 degrees)
        rho_0= 983.23  # (kg/m^3) reference density of WATER in pipe at temperature T_0, taken at atmospheric pressure (10^5 Pa)
        T_0= 60+273.15 # reference temperature for fluid, now taken at 60 degrees
        lambda_fluid=0.6506 #0.665 W/(m K)   #thermal conductivity of reference fluid, WATER at 60 degrees
        beta= 4.57*10**-4  # (K**-1) ORDERGROOTESCHATTING Thermal expansion coefficient of WATER (data compagnion) (constant for determining density at different temperature)
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
        #Sol_salt=Na*solution/rho_0 #atoms/m^3 maximum oplosbaarheid of salt in water
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
        
        
    
        
        
        
        
        
        
        
        
        
        
        def gravity(v,T):
            grav=np.zeros(p.N)
            for n in np.arange(0,p.N,1):
                grav[n]=(rho_0-rho_0*beta*(T[n]-p.T_0))*np.sin(angle[n])        
            #grav1=np.sum((p.rho_0-p.rho_0*p.beta*(T-p.T_0))*np.sin(p.angle))
            if v>0:
                grav1=sum(grav)
            else:
                grav1=-sum(grav)
            return grav1
                
        
        def darcyfriction(v,T,Tb=0):
            Re=p.Reynolds(v, T)
            a=1/(1+((Re/2712)**8.4))
            b=1/(1+((Re/(150*2*r/p.eff))**1.8))
            #df=64/Re
            df=((64/Re)**a)*((0.75*np.log(Re/5.37))**(2*(a-1)*b))*((0.88*np.log(3.41*2*r/p.eff))**(2*(a-1)*(1-b)))
            #df=((64/Re)**a)*((0.75*np.log(Re/5.37))**(2*(a-1)*b))*((0.88*np.log(6.82*2*p.r/p.eff))**(2*(a-1)*(1-b)))
            f=np.sum(df)
            if v>0:
                f1=f
            else:
                f1=-f
            return f1        
                
        
        def phiAB(v,T,Tb):
            ''' energy transport form the fluid to the wall  '''
            phi=np.zeros(p.N)
            for n in range(len(phi)):
                phi[n]=h_AB(n, v, T[n])*Opp_wallAB*(T[n]-Tb[n])
            return phi
        
        def phiBC(v,T,Tb):
            ''' energy transport from the wall to the surrounding water '''
            phi=np.zeros(p.N)
        
            for n in np.arange(0,p.N,1):
                phi[n]=h_BC(n,v,T,Tb)*Opp_wallBC[n]*(Tb[n]-p.T_c)
            return phi
        
        def Twallzero(v,T,Tb):
            twallzero=np.zeros(p.N)
            PhiAB=phiAB(v, T, Tb)
            PhiBC=phiBC(v, T, Tb)
            for n in np.arange(0,p.N,1):
                twallzero[n]=-PhiAB[n] + PhiBC[n]-Vwall[n]*p.u*p.rho_wall 
            return twallzero
                
        def Tfluidzero(v,T,Tb):
            
            tfluidzero=np.zeros(p.N)
            PhiAB=phiAB(v, T, Tb)
            
            Constatns=p.rho_0*p.C_pfluid*np.pi*r**2
            
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
                Velocityzero=-f_D*(1/(4*r))*v**2 - ((kw1+kw2)*v**2)/length + (g/(rho_0*p.N))*gravity(v,T)
            else:
                Velocityzero=-f_D*(1/(4*r))*v**2 + ((kw1+kw2)*v**2)/length + (g/(rho_0*p.N))*gravity(v,T)
            #Velocityzero=-2*f_D*v**2/p.r + (p.g/(p.rho_0*p.N))*gravity(v,T) - ((p.kw1+p.kw2)*v**2)/p.length 
            #f_D=darcyfriction(v, T)
            #Velocityzero=-2*f_D*(1/p.r)*v**2 + (p.g/(p.rho_0))*gravity(v,T)
            return Velocityzero        
        
        
        
        answer=optimize.fsolve(system,initialguess)
        
        v2=answer[0]
        T2,Tb2=np.array_split(answer[1:],2)

        vnsend[ite,j]=v2
        Tnsend[:,ite,j]=T2
        Tnsbend[:,ite,j]=Tb2
        
        #initialguess=np.concatenate([np.array([v2]),T2,Tb2])
        #initialguess=np.concatenate([np.array([p.v_steadystate0]),p.T_steadystate0,p.Tb_steadystate0])
        j=j+1
    ite=ite+1

mx=0
Maxtemperature=np.zeros((len(drbelow),len(rtube)))
interval=0
for nTL in np.arange(0,p.N,1):
    for DRBTL in np.arange(0,len(drbelow),1):
        for DRUTL in np.arange(0,len(rtube),1):
            mx=max(Tnsend[nTL,DRBTL,DRUTL],mx)
            Maxtemperature[DRBTL,DRUTL]=max(Maxtemperature[DRBTL,DRUTL],Tnsend[nTL,DRBTL,DRUTL])
            if Tnsend[nTL,DRBTL,DRUTL]>=273.15+90 and interval==0:
                print('Temperature of fluid goes above 90 degrees, at segment %.f2 in the loop with dr1=%.3f and verticaldeel=%.3f, with temperature %.3f K' %(nTL, drbelow[DRBTL],rtube[DRUTL],Tnsend[nTL,DRBTL,DRUTL]))
                interval=1
mxb=0
intervalb=0
for nTbL in np.arange(0,p.N,1):
    for DRBTbL in np.arange(0,len(drbelow),1):
        for DRUTbL in np.arange(0,len(rtube),1):
            mxb=max(Tnsbend[nTL,DRBTbL,DRUTbL],mxb)
            if Tnsend[nTbL,DRBTbL,DRUTbL]>=273.15+150 and intervalb==0:
                print('Temperature of wall goes above 150 degrees, at segment %.f2 in the loop with dr1=%.3f and verticaldeel=%.3f' %(nTL, drbelow[DRBTL],rtube[DRUTL]))
                intervalb=1

plt.figure()
np.save(listnameTn[0],Tnsend)
np.save(listnameTbn[0],Tnsbend)
np.save(listnamev[0],vnsend)
plt.imshow(vnsend,extent=[rtube[0],rtube[-1],drbelow[0],drbelow[-1]],origin='lower', aspect='auto')
plt.colorbar()

plt.xlabel('location in loop when thickness changes')
plt.ylabel('thickness of heating wall')
plt.title('velocity at different thickness profile of the wall')

plt.figure()
np.save(listnamemaxT[0],Maxtemperature)

plt.imshow(Maxtemperature,extent=[rtube[0],rtube[-1],drbelow[0],drbelow[-1]],origin='lower', aspect='auto')
plt.colorbar()

plt.xlabel('location in loop when thickness changes')
plt.ylabel('thickness of heating wall')
plt.title('maximum temperature inside the loop \n at different thickness profile of the wall')

'determining maxtemperature and safety region '
if True:
    mx=0
    Maxtemperature=np.zeros((len(drbelow),len(rtube)))
    geldigregimeofniet=np.zeros((len(drbelow),len(rtube)))
    interval=0
    for nTL in np.arange(0,p.N,1):
        for DRBTL in np.arange(0,len(drbelow),1):
            for DRUTL in np.arange(0,len(rtube),1):
                mx=max(Tnsend[nTL,DRBTL,DRUTL],mx)
                Maxtemperature[DRBTL,DRUTL]=max(Maxtemperature[DRBTL,DRUTL],Tnsend[nTL,DRBTL,DRUTL])
                if Maxtemperature[DRBTL,DRUTL]>=273.15+90:
                    geldigregimeofniet[DRBTL,DRUTL]=1
    
    plt.figure()
    plt.imshow(Maxtemperature,extent=[rtube[0],rtube[-1],drbelow[0],drbelow[-1]],origin='lower', aspect='auto')
    plt.colorbar()
    plt.xlabel('location in loop when thickness changes')
    plt.ylabel('thickness of heating wall')
    plt.title('Maximum temperature in the loop at \n different thickness profile of the wall')
    
    plt.figure()
    plt.imshow(geldigregimeofniet,extent=[rtube[0],rtube[-1],drbelow[0],drbelow[-1]],origin='lower', aspect='auto')
    plt.colorbar()
    plt.xlabel('location in loop when thickness changes')
    plt.ylabel('thickness of heating wall')
    plt.title('safety region inside the loop with temperature below 90$^o$ C \n of the loop at different thickness profile of the wall')



'determining maximum velocity'
if True:
    mxv=0
    for DRBTL in np.arange(0,len(drbelow),1):
        for DRUTL in np.arange(0,len(rtube),1):
            mxv=max(vnsend[DRBTL,DRUTL],mxv)
            #Maxtemperature[DRBTL,DRUTL]=max(Maxtemperature[DRBTL,DRUTL],Tnsend[nTL,DRBTL,DRUTL])
            if vnsend[DRBTL,DRUTL]>=mxv:
                Vplaceofdr=drbelow[DRBTL]
                Vplaceofr=rtube[DRUTL]
                VNdr=DRBTL
                VNr=DRUTL


'determining length of sytem and feasible region'
if False:
    mx=0
    lengthofsystem=np.zeros((len(drbelow),len(rtube)))
    geldigregimeofnietvoorlength=np.zeros((len(drbelow),len(rtube)))
    interval=0
    for DRBTL in np.arange(0,len(drbelow),1):
        dr1=drbelow[DRBTL]*10**-3
        for DRUTL in np.arange(0,len(rtube),1):
            r=rtube[DRUTL]*10**-3
            #mx=max(Tnsend[nTL,DRBTL,DRUTL],mx)
            lengthofsystem[DRBTL,DRUTL]=(4*(1+np.sin(anglepipe))**-1)*(Ltube-dr1-dr3-2*r)
            if lengthofsystem[DRBTL,DRUTL]<=0:
                    geldigregimeofnietvoorlength[DRBTL,DRUTL]=1           
    
    plt.figure()
    plt.imshow(lengthofsystem,extent=[rtube[0],rtube[-1],drbelow[0],drbelow[-1]],origin='lower', aspect='auto')
    plt.colorbar()
    plt.xlabel('location in loop when thickness changes')
    plt.ylabel('thickness of heating wall')
    plt.title('length of the loop at different thickness profile of the wall')
    
    plt.figure()
    plt.imshow(geldigregimeofnietvoorlength,extent=[rtube[0],rtube[-1],drbelow[0],drbelow[-1]],origin='lower', aspect='auto')
    plt.colorbar()
    plt.xlabel('location in loop when thickness changes')
    plt.ylabel('thickness of heating wall')
    plt.title('correct region of length of the loop at /n different thickness profile of the wall')
    

'determining production and volume of system'
if False:
    mx=0
    Natoms= p.Sol_Mo98#atoms Mo98 per volume
    Ratoms=p.Cross_section*p.flux*Natoms #atoms mo99 per volume
    R_mass=Ratoms*p.Molair_Mo99/p.Na #g per volume
    volumeofsystem=np.zeros((len(drbelow),len(rtube)))
    Reactionpersecond=np.zeros((len(drbelow),len(rtube))) #in gram 
    interval=0
    for DRBTL in np.arange(0,len(drbelow),1):
        dr1=drbelow[DRBTL]*10**-3
        for DRUTL in np.arange(0,len(rtube),1):
            r=rtube[DRUTL]*10**-3
            #mx=max(Tnsend[nTL,DRBTL,DRUTL],mx)
            lengthofsystem=(4*(1+np.sin(anglepipe))**-1)*(Ltube-dr1-dr3-2*r)
            volumeofsystem[DRBTL,DRUTL]=lengthofsystem*np.pi*r**2
            Reactionpersecond[DRBTL,DRUTL]=R_mass*volumeofsystem[DRBTL,DRUTL]
    
    plt.figure()
    plt.imshow(volumeofsystem,extent=[rtube[0],rtube[-1],drbelow[0],drbelow[-1]],origin='lower', aspect='auto')
    plt.colorbar()
    plt.xlabel('location in loop when thickness changes')
    plt.ylabel('thickness of heating wall')
    plt.title('volume of the loop at different thickness profile of the wall')
    
    plt.figure()
    plt.imshow(Reactionpersecond,extent=[rtube[0],rtube[-1],drbelow[0],drbelow[-1]],origin='lower',aspect='auto')
    plt.colorbar()
    plt.xlabel('location in loop when thickness changes')
    plt.ylabel('thickness of heating wall')
    plt.title('Reaction Rate the loop at different thickness profile of the wall')
    
    'determining max production'
    if True:
        mxP=0
        for DRBTL in np.arange(0,len(drbelow),1):
            for DRUTL in np.arange(0,len(rtube),1):
                mxP=max(Reactionpersecond[DRBTL,DRUTL],mxv)
                #Maxtemperature[DRBTL,DRUTL]=max(Maxtemperature[DRBTL,DRUTL],Tnsend[nTL,DRBTL,DRUTL])
                if Reactionpersecond[DRBTL,DRUTL]>=mxP:
                    Rplaceofdr=drbelow[DRBTL]
                    Rplaceofr=rtube[DRUTL]
                    RNdr=DRBTL
                    RNr=DRUTL
    
    
    
    
    
    
    