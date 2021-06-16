# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 20:06:52 2021

@author: ronar
"""

import parameters as p
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import tqdm as tqdm
plt.close('all')




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
begindr1=0.5
enddr1=6.1
begindr3=begindr1
enddr3=enddr1
step1=0.3
step3=step1


r=3*10**-3
deel=0.5
vnsend=np.zeros((len(np.arange(begindr1,enddr1,step3)),len(np.arange(begindr3,enddr3,step3))))
Tnsend=np.zeros((p.N,len(np.arange(begindr1,enddr1,step3)),len(np.arange(begindr3,enddr3,step3))))
Tnsbend=np.zeros((p.N,len(np.arange(begindr1,enddr1,step3)),len(np.arange(begindr3,enddr3,step3))))
ite=0
drbelow=np.arange(begindr1,enddr1,step3)
drup=np.arange(begindr3,enddr3,step3)
for thick1 in drbelow:
    dr1=thick1*10**-3
    dr2=thick1*10**-3
    j=0
    for thick2 in drup:
        dr3=thick2*10**-3
        dr4=thick2*10**-3
        
        
        
        
        
        anglepipe=p.anglepipe
        Ltube=p.Ltube
        N=p.N
        length=(4*(1+np.sin(anglepipe))**-1)*(Ltube-dr1-dr3-2*r) #hier gaan we ervan uit dat in de 1e buis de dikte altijd dr1 is en in de 3e buis altijd dr2
        dr=np.zeros(N)
        u=p.u
        rho_wall=p.rho_wall
        Sol_salt=p.Sol_salt
        Cross_section=p.Cross_section
        flux=p.flux
        Molair_Mo99=p.Molair_Mo99
        Na=p.Na
        a_fluid=p.a_fluid
        g=p.g
        T_c=p.T_c
        beta=p.beta
        nu_water=p.nu_water
        lambda_wall=p.lambda_wall
        lambda_fluid=p.lambda_fluid
        a_water=p.a_water
        lambda_water=p.lambda_water
        for i in np.arange(0,N,1):
            if i<0.25*N:
                dr[i]=dr1
            elif (i>0.25*N or i==0.25*N) and (i<deel*N) :
                dr[i]=dr2
            elif (i>deel*N or i==deel*N) and (i<0.75*N) :
                dr[i]=dr3
            elif (i>0.75*N or i==0.75*N) and (i<N) :
                dr[i]=dr4
        
        tVsys= length*np.pi*r**2 #volume of whole system within the pipe
        Vsys_segm=(length/N)*np.pi*r**2 # volume of segment within the pipe
        Opp_wallAB=(length/N)*r*2*np.pi
        Opp_wallBC=np.zeros(N)  #area of outer wall of segment of system
        for i in np.arange(0,N,1):
            Opp_wallBC[i]=(length/N)*2*np.pi*(r+dr[i])    
        Vwall=np.zeros(N)
        for i in np.arange(0,N,1):# volume of wall of segment of system 
            Vwall[i]=(length/N)*np.pi*((r+dr[i])**2-r**2)
            
        Pgamma=np.zeros(N)
        for n in np.arange(0,N,1):
            Pgamma[n]=u*rho_wall*Vwall[n]        # heat produced in each wall segment of the pipe
            
        '___________________Mo98 production_______________________'
        def Reaction():
            N= Sol_salt#atoms Mo98 per volume
            R=Cross_section*flux*N  #atoms mo99 per volume
            R_mass=R*Molair_Mo99/Na #g per volume
            return R_mass*tVsys
        def Greatz(v):
            Gz=a_fluid*length/(abs(v)*((2*r)**2))
            return Gz
        
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
            Re=p.Reynolds(v, T)
            Gz=Greatz(v)
            Pr=p.Prandtl(T)  
            
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
                Velocityzero=-f_D*(1/(4*r))*v**2 - ((p.kw1+p.kw2)*v**2)/length + (p.g/(p.rho_0*p.N))*gravity(v,T)
            else:
                Velocityzero=-f_D*(1/(4*r))*v**2 + ((p.kw1+p.kw2)*v**2)/length + (p.g/(p.rho_0*p.N))*gravity(v,T)
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
        
        j=j+1
    ite=ite+1


plt.imshow(vnsend,extent=[drup[0],drup[-1],drbelow[0],drup[-1]])
plt.colorbar()

plt.xlabel('thickness of cooling wall')
plt.ylabel('thickness of heating wall')
plt.title('velocity at different thicknesses')
