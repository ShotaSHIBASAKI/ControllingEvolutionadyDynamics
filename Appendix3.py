#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 08:56:46 2019
Code for Appendix 3
Jacobian matrix analysis
and Figs A.6 and A.7
@author: shota
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import sympy as sym
from sympy import init_printing

#fixed parameters
dmax=1 #maximum mortality rate caused by toxin
fmax=0.5 #maximum degradation
Kd=0.2 # median-effet density for the toxic degradation
rmax=1

def AnalyzeA3():
    init_printing()
    """
    Analyzing the 3x3 Jacobian matrix for linear stability
    """
    T0=sym.Symbol('T_{in}')#Toxin concentration flow into chemostat
    T=sym.Symbol('T') #toxin concentration
    ri=sym.Symbol('r_{i}')#intrisic growth rate of i
    rj=sym.Symbol('r_{j}')#intrisic growth rate of j
    xi=sym.Symbol('x_{i}')#density of i
    xj=sym.Symbol('x_{j}')#density of j
    fx=sym.Symbol('f(x)')#total detoxification efficiency
    df=sym.Symbol('df')#df/dx_i where i=sCo, rCo 
    dfi=sym.Symbol('df/dx_{i}')
    dfj=sym.Symbol('df/dx_{j}')
    a=sym.Symbol('a')#dilution rate, instead of alpha
    di=sym.Symbol('d_{i}(T)')#death rate of i
    dj=sym.Symbol('d_{j}(T)')#death rate of j
    ddi=sym.Symbol("d'_{i}(T)")#d(T)_i/dT
    ddj=sym.Symbol("d'_{j}(T)")#d(T)_j/dT
    #Define 3x3 Jacobian matrix in general
    
    J=sym.Matrix([[-(a+fx),-T*dfi, -T*dfj],
                  [-xi*ddi, ri*(1-2*xi-xj)-di-a, -ri*xi],
                  [-xj*ddj,-rj*xj, rj*(1-xi-2*xj)-dj-a]])
    print("The Jacobian matirx for the invasion analysis is  as shown in Eq (A.32),")
    display(J)
    """
    From Routh-Hurwiz criteria, the equilibrium is stable iff
    trJ<0
    detJ<0
    M11+M22+M33>0 (sum of minodr determinant)
    """
    trJ=J.trace()
    detJ=J.det()
    M11=J[1,1]*J[2,2]-J[1,2]*J[2,1] #minor determinant M11
    M22=J[0,0]*J[2,2]-J[0,2]*J[2,0] #minor determinanr M22
    M33=J[0,0]*J[1,1]-J[0,1]*J[1,0] #minor determinant M33
    
    #Mono-culture of sCo rCo
    print('\n Case 1: Mono-culture of i=sCo (or rCo) while j= rCh (sCh) goes extinct')
    J1=J
    J1[1,1]=-ri*xi #ri(1-xi)-di-a=0 while xj=0
    J1=sym.simplify(J1.subs([(dfi,df), (dfj,0), (xj,0)]))
    print("The Jacobian matrix is then")
    display(J1)
    print("Eqs (A.35a) - (A.35e) are checked as below")
    M11=J1[1,1]*J1[2,2]-J1[1,2]*J1[2,1] #minor determinant M11
    M22=J1[0,0]*J1[2,2]-J1[0,2]*J1[2,0] #minor determinanr M22
    M33=J1[0,0]*J1[1,1]-J1[0,1]*J1[1,0] #minor determinant M33
    display(J1.trace(), J1.det(), M11, M22, M33)
    
    #Coexistence of sCo (sCh) with rCo (rCh)
    print("\n Case 2: Coexistence of sCo (sCh) with rCo (rCh)")
    print("Proof of Eq( A.40)")
    print("When sCo coexists with rCo, the Jacobiam matrix is")
    J_sCorCo=J
    """
    Notice that ri(1-xi-xj)-di-a=rj(1-xi-xj)-dj-a=0
    <=> ri(1-2xi-xj)-di-a=-rixi, 
        rj(1-xi-2xj)-dj-a=-rjxj
    
    """
    J_sCorCo[1,1]=-ri*xi
    J_sCorCo[2,2]=-rj*xj
    J_sCorCo=sym.simplify(J_sCorCo.subs([(dfi,df),(dfj, df)]))
    #J_sCorCo[1,1]=J_sCorCo[1,1].subs(xi, 1-xj-(di+a)/ri)
    display(J_sCorCo)
    print("and its determinat is")
    #display(sym.simplify(J_sCorCo[1,1].subs(xi, 1-xj-(di+a)/ri)))
    display(J_sCorCo.det())
    print('meaning that this equilibrium should not be linearly stable.')
    print("On the other hand, when sCh coexist with rCh")
    J_sChrCh=J_sCorCo.subs(df,0)
    print("the Jacobian matrix is")
    display(J_sChrCh)
    print('and the determinant is again')
    display(J_sChrCh.det())
    print('meaning that this coexistence should not be linearly stable.')
    
    #Coexistence of i=sCo (rCo) with j=rCh (sCh)
    """
    In this case the trace is negative and the analysis of determinant is not difficult
    Notice that ri(1-xi-xj)-di-a=rj(1-xi-xj)-dj-a=0
    <=> ri(1-2xi-xj)-di-a=-rixi, 
        rj(1-xi-2xj)-dj-a=-rjxj
    """
    J2=J
    J2[1,1]=-ri*xi
    J2[2,2]=-rj*xj
    J2=sym.simplify(J2.subs([(dfi,df),(dfj, 0)]))
    print("\n Case 3: Coexistence of sCo (rCo) with rCh (sCh)")
    print('In this case, the Jacobian matrix is')
    display(J2)
    print('where the trace of the Jacobian is obviously negative')
    print("The determinant if the Jacobian is then,")
    display(sym.factor(J2.det())<0)
    print("and the following analysis is not difficult as shown in Eqs (A.41) - (A.46).")
    print("Now let consider the sum of the minor determinant.")
    print("M_11 is ovbiously zero. So, we should consider M_22 and M_33:")
    M22=J2[0,0]*J2[2,2]-J2[0,2]*J2[2,0] #minor determinanr M22
    M33=J2[0,0]*J2[1,1]-J2[0,1]*J2[1,0]
    M=sym.factor(M22+M33)
    print("M_22 + M_33 equals to")
    display(M)
    M_RHS=T*ddi*df*xi
    M_LHS=sym.factor(M+M_RHS)
    display(M_LHS>M_RHS)
    #simplifiing inequality
    M_RHS=M_RHS/ddi/xi/(a+fx)*a*T0/(a+fx)/T
    M_LHS=M_LHS/ddi/xi/(a+fx)
    print("Simplifying the inequality, we obtain inequality (A.47)")
    display(M_LHS>M_RHS)
"""
Drawing Fig A.6 eco-evo dynamics in phase plains
Run FigA6

""" 
def DynamicsCoop(x, t, T0, alpha, cost, K, n, r):
    """
    mono-culture of cooperators to obtain the T^* and x_1^* 
    where the cheaters are excluded
    """
    
    #defune the degradation effect
    deg=fmax*x[1]/(x[1]+Kd)
    D=dmax*x[0]**n/(x[0]**n+K**n)
    #define ODEs
    y=np.zeros([np.size(x)])
    y[0]=alpha*T0-deg*x[0]-alpha*x[0] #dt/dt
    y[1]=x[1]*(r*(1-cost)*(1-x[1])-D-alpha)#d Co/dt
    
    return y

def DynamicsCheat(x, t, T0, alpha, cost, K, n, r):
    """
    mono-culture of cheaterss to obtain the T^* and x_1^* 
    """
    D=dmax*x[0]**n/(x[0]**n+K**n)
    #defune the degradation effect
    deg=0
    #define ODEs
    y=np.zeros([np.size(x)])
    y[0]=alpha*T0-deg*x[0]-alpha*x[0] #dt/dt
    y[1]=x[1]*(r*(1-cost)*(1-x[1])-D-alpha)#d Co/dt
    
    return y

def DynamicsCo(x, t, T0, alpha, cost_co, cost_ch, K_co, K_ch, n, r):
    """
    Co-culture of the coopertors and the cheaters
    with the environmental feedback thorough changing the toxic concentration
    """
    y=np.zeros([np.size(x)])
    D=np.zeros([2])   
    #define fitnss
    D[0]=dmax*x[0]**n/(K_co**n+x[0]**n) #cooperator
    D[1]=dmax*x[0]**n/(K_ch**n+x[0]**n) #cheater   
    #degradation
    deg=fmax*x[1]/(x[1]+Kd)   
    #ODE of eco-evo dynamics
    y[0]=alpha*T0-deg*x[0]-alpha*x[0] #dt/dt
    y[1]=x[1]*(r*(1-cost_co)*(1-x[1]-x[2])-D[0]-alpha)#d Co/dt
    y[2]=x[2]*(r*(1-cost_ch)*(1-x[1]-x[2])-D[1]-alpha) #d Ch/dt
       
    return y
 
def FigA6(T0, alpha, cd, cr, Ks, Kr, n, case):
    """
    T0: toxin flowing into the chemostat
    alpha: dilution rate
    cd: cost of cooperation
    cr: cost for resistance
    Ks: median-effect toxin concentration for sensitive
    Kr: median-effect toxin concentration for resistant
    n: Hill coefficient
    
    """
    r=rmax#maximum intrinsic growth rate
    # check whether the cooperator can grow in mono-culture
    if case==0:
        W_T0=r*(1-cd-cr)-dmax*T0**n/(Kr**n+T0**n)
    elif case==1:
        W_T0=r*(1-cd)-dmax*T0**n/(Ks**n+T0**n)
    else:
        print("Error in case; case should be 0 or 1.")
    if W_T0>alpha:
        print('cooperators can grow in mono-culture')
        init=np.array([T0, 10**(-6)])
        time=np.linspace(0,100000,100001)
        sol=odeint(DynamicsCheat, init, time, args=(T0, alpha, cr, Kr, n, r))
        if case==0:
            sol=odeint(DynamicsCoop, init, time, args=(T0, alpha, cr+cd, Kr, n, r))
        else:
            sol=odeint(DynamicsCoop, init, time, args=(T0, alpha, cd, Ks, n, r))
        
        plt.plot(time, sol[:,0])
        plt.ylim(0, T0+0.1)
        plt.xlabel('time', fontsize=14)
        plt.ylabel('toxic concentration', fontsize=14)
        plt.show()
        print(sol[-1,0])
        
    
    else:
        print('Error: cooperator cannot grow in mono-culture')
        print(r*W_T0)
        return 1
    #plot the initial condition mono-culture of cheaters
    init=np.array([T0, 10**(-6)])
    time=np.linspace(0,100000,500001)
    if case==0:
        #cheater is sCh
        sol=odeint(DynamicsCheat, init, time, args=(T0, alpha, 0, Ks, n, r))
    else:
        #cheater is rCh
        sol=odeint(DynamicsCheat, init, time, args=(T0, alpha, cr, Kr, n, r))
    x2s=sol[-1,1]    
    init=np.array([sol[-1,0], 10**(-6), x2s-10**(-6),])#mutation
    if sol[-1,1]<10**(-6):
        print('Cheaters died')
        return 1
        
    if case==0:
        #cooperator is rCo
        sol=odeint(DynamicsCo, init, time, args=(T0, alpha, cd+cr ,0, Kr, Ks, n, r))
    else:
        #cooperator is sCo
        sol=odeint(DynamicsCo, init, time, args=(T0, alpha, cd ,cr, Ks, Kr, n, r))
    plt.plot(time, sol[:,1], label='cooperator')
    plt.plot(time, sol[:,2], label='cheater')
    plt.xlabel('time', fontsize=14)
    plt.ylabel('density', fontsize=14)
    plt.ylim(0,1)
    plt.legend(loc='best')
    plt.show()
    
    plt.plot(sol[:,1], sol[:,2])
    plt.xlabel('cooperator', fontsize=14)
    plt.ylabel('cheater', fontsize=14)
    #plt.xlim(0,1)
    #plt.ylim(0,1)
    plt.savefig('phaseplain-Co-Ch.pdf')
    plt.show()
    
    plt.plot(sol[:,0], sol[:,1])
    plt.ylabel('cooperator', fontsize=14)
    plt.xlabel('toxin', fontsize=14)
    #plt.xlim(0,1)
    #plt.ylim(0,1)
    plt.savefig('phaseplain-Tox-Co.pdf')
    plt.show()
    
    plt.plot(sol[:,0], sol[:,2])
    plt.ylabel('cheater', fontsize=14)
    plt.xlabel('toxin', fontsize=14)
    #plt.xlim(0,1)
    #plt.ylim(0,1)
    plt.savefig('phaseplain-Tox-Ch.pdf')
    plt.show()
    fname='parameter-list'
    parm=np.array([[case,T0, alpha, cd, cr, Ks, Kr, n, r]]).T
    np.savetxt(fname+'.csv',parm.T, delimiter=',', header='case, T0, alpha, cd, cr, Ks, Kr, Hill, rmax', fmt='%.6f') 
    
"""
Fig A.7 Monte Calro simulation to estimate the space
where two equilibria are stable, respectively.
Run Fig A7
"""    
def DyCoop(x, t, T0, r1, K, alpha, n):
    """
    mono-culture of cooperators to obtain the T^* and x_1^* 
    where the cheaters are excluded
    """
    
    #defune the degradation effect
    deg=fmax*x[1]/(x[1]+Kd)
    D=dmax*x[0]**n/(x[0]**n+K**n)
    #define ODEs
    y=np.zeros([np.size(x)])
    y[0]=alpha*T0-deg*x[0]-alpha*x[0] #dt/dt
    y[1]=x[1]*(r1*(1-x[1])-D-alpha)#d Co/dt
    
    return y

def DyCheat(x, t, T0, r2, K, alpha, n):
    """
    mono-culture of cheaterss to obtain the T^* and x_1^* 
    """
    D=dmax*x[0]**n/(x[0]**n+K**n)
    #defune the degradation effect
    deg=0
    #define ODEs
    y=np.zeros([np.size(x)])
    y[0]=alpha*T0-deg*x[0]-alpha*x[0] #dt/dt
    y[1]=x[1]*(r2*(1-x[1])-D-alpha)#d Co/dt
    
    return y

def DyEvo(x, t, T0, r1, r2, K_co, K_ch, alpha, n):
    """
    Co-culture of the coopertors and the cheaters
    with the environmental feedback thorough changing the toxic concentration
    """
    y=np.zeros([np.size(x)])
    D=np.zeros([2])   
    #define fitnss
    D[0]=dmax*x[0]**n/(K_co**n+x[0]**n) #cooperator
    D[1]=dmax*x[0]**n/(K_ch**n+x[0]**n) #cheater   
    #degradation
    deg=fmax*x[1]/(x[1]+Kd)   
    #ODE of eco-evo dynamics
    y[0]=alpha*T0-deg*x[0]-alpha*x[0] #dt/dt
    y[1]=x[1]*(r1*(1-x[1]-x[2])-D[0]-alpha)#d Co/dt
    y[2]=x[2]*(r2*(1-x[1]-x[2])-D[1]-alpha) #d Ch/dt
       
    return y
def Quad(case,alpha,c1,c2,K1, K2,n):
    A=(dmax+alpha)*(c2-c1)
    B=dmax*(1-c1)*K1**n-dmax*(1-c2)*K2**n+alpha*(c2-c1)*(K1**n+K2**n)
    C=alpha*(c2-c1)*K1**n*K2**n

    
    if case==0:
        solp=(-B+(B**2-4*A*C)**(1/2))/(2*A)
        solm=(-B-(B**2-4*A*C)**(1/2))/(2*A)
        # then there exist real roots
        return min(solp,solm)
    else:
        if B**2-4*A*C<0:
            # no real roots and therefore cheaters are excluded
            solp=1
            solm=1# with any value of T0, cheaters go extinct
        else:
            solp=max((-B+(B**2-4*A*C)**(1/2))/(2*A),(-B-(B**2-4*A*C)**(1/2))/(2*A))
            solm=min((-B+(B**2-4*A*C)**(1/2))/(2*A),(-B-(B**2-4*A*C)**(1/2))/(2*A))
        return [solp,solm]
    
def Stab_cond(alpha, T0, T,x1,x2, r1,r2,n, K):
    #rhs=-alpha*T0*fmax*Kd/(alpha*(x1+Kd)+fmax*x1)**2
    rhs=-alpha*T0*(fmax*Kd/(x1+Kd)**2)/((alpha+fmax*x1/(x1+Kd))**2)
    delta_p=dmax*n*K**n*T**(n-1)/(T**n+K**n)**2
    lhs=-(r1*x1+r2*x2)/(x1*delta_p)
    
    if lhs<rhs:
        return 0
    else:
         return 1
     
def FigA7(case):
    """
    Classify the evolutionary fate given n and case
    case: 0 rCo vs Sch, 1 sCO vs rCo
    n: Hill coeff.
    
    Parameter sweep with below
    Kr K of resistant
    Ks K of sensitive (Ks<Kr)
    cd cost of detoxification
    cr cost of resistance (cr>cd)
    alpha dilution rate
    T0 initial toxin conc. 
    """
    
    #set the parameter, arrays
   
    n_array=np.array([1,2,3])

    #set the result arrays
    if case==0:
        class_number=5
    elif case==1:
        class_number=6
    fate=np.zeros([class_number])#number of evolutionary fate
    fate_matrix=np.zeros([np.size(n_array),np.size(fate)])
    
    time=np.linspace(0,100000, 1000000)
    loop=10**6
    """
    0 Co and/or Ch cannot survive in mono-culture
    1 Co cannot invade
    2 Only equilibrium of exclusion is stable
    3 Only equilibrium of coexistence is stable
    4 Two equilibria are UNstable
    5 two Equilibrium are stable (which may occur only when sCO vs rCh)
    """
    for tri in range(np.size(n_array)):
        counter=0
        n=n_array[tri]
        print(str("Hill coefficient is %d" %(n)))
        fate=np.zeros([class_number])#number of evolutionary fate should be reset
        if case==0 or case==1:
            fname=str('parameter-sweep-MC-n%d-case%d' %(n, case))
        else:
            print("Error in case")
            return 1
    
        for i in range(loop):
            if(i+1)%10000==0:
                print(i+1)
            Ks,cd,T0, alpha,=np.random.uniform(0,1,4)
            Kr,cr=np.random.uniform([Ks,0],[1,1],2)#Kr>Ks and cr.cd
            #check whether r is positive or not
            if case==0:
                r1=rmax*(1-cr-cd)#rCO
                r2=rmax#sCH
                W0Co=r1-dmax*T0**n/(T0**n+Kr**n)-alpha#initial growth of Cooperator
                W0Ch=r2-dmax*T0**n/(T0**n+Ks**n)-alpha#initial growth of Cheater
            elif case==1:
                r1=rmax*(1-cd)#sCo
                r2=rmax*(1-cr)#rCh
                W0Co=r1-dmax*T0**n/(T0**n+Ks**n)-alpha
                W0Ch=r2-dmax*T0**n/(T0**n+Kr**n)-alpha
                stab_e=0#initialize the falgs of stability
                stab_c=0
            if W0Co<0 or W0Ch<0:
                fate[0]+=1
                res=0
            else:
                #succeed in mono-culture             
                init=np.array([T0,10**(-6)])
                if case==0:                                        
                    solCo=odeint(DyCoop, init, time, args=(T0, r1, Kr, alpha, n))
                    Ts=solCo[-1,0]
                    #x1s=solCo[-1,1]
                    solCh=odeint(DyCheat, init, time, args=(T0, r2, Ks, alpha, n))
                    x2s=solCh[-1,1]
                else:
                    solCo=odeint(DyCoop, init, time, args=(T0, r1, Ks, alpha, n))
                    Ts=solCo[-1,0]
                    #x1s=solCo[-1,1]
                    solCh=odeint(DyCheat, init, time, args=(T0, r2, Kr, alpha, n))
                    x2s=solCh[-1,1]
                
                #Evolutionary dynamics 
                if case==0:
                    K=Kr
                else:
                    K=Ks
                if r1*(1-x2s)-dmax*T0**n/(T0**n+K**n)<alpha:
                    #Co cannot invade
                    fate[1]+=1
                    res=1
                else:
                    #Co can invade
                    #calculate Tdagger Td and check whether coexist or exclude
                    if case==0:
                        #rCo vs sCh
                        #in this case, at most one equilbrium is stable
                        tau=Quad(case,alpha,cr+cd,0,Kr, Ks, n)
                        Td=tau**(1/n)
                        if Td<Ts:
                            #Co exclude Ch
                            fate[2]+=1
                            res=2
                        else:
                            x1d=alpha*Kd*(T0-Td)/(fmax*Td-alpha*(T0-Td))
                            x2d=1-x1d-(dmax*Td**n/(Td**n+K**n)+alpha)/r1
                            #check the stability condition
                            stab=Stab_cond(alpha, T0, Td,x1d,x2d, r1,r2,n, K)
                            if stab==0:
                                #stable coexistence
                                fate[3]+=1
                                res=3
                            else:
                                #unstable coexistence nor exclusion
                                fate[4]+=1
                                res=4
                                print(Td, x1d, x2d)
                    else:
                        #sCo vs rCh
                        # in this case two equilibria can be stable at the same time
                        [tau_p,tau_m]=Quad(case,alpha,cd,cr,Ks, Kr, n)
                        if tau_m>Ts**n or tau_p<Ts**n:
                            # cexclusion is stable
                            stab_e=1
                        # stability in coexistence 
                        if tau_p<0:
                            stab_c=0
                        else:
                            Td=tau_p**(1/n)
                            x1d=alpha*Kd*(T0-Td)/(fmax*Td-alpha*(T0-Td))
                            x2d=1-x1d-(dmax*Td**n/(Td**n+K**n)+alpha)/r1
                            #check the stability condition
                            stab=Stab_cond(alpha, T0, Td,x1d,x2d, r1,r2,n, K)
                            if stab==0:
                                #stable coexistence
                                stab_c=1
                        #classify
                        if stab_e==1 and stab_c==1:
                            # two stable equilbria
                            fate[5]+=1
                            res=5
                        elif stab_e==1 and stab_c==0:
                            #only stable cexclusion
                            fate[2]+=1
                            res=2
                        elif stab_e==0 and stab_c==1:
                            #stable coexistence
                            fate[3]+=1
                            res=3
                        else:
                            #both unstable
                            fate[4]+=1
                            res=4
                        
            #save the results
            if counter==0:
                result=np.array([[Ks, Kr, cr, cd, alpha, T0,res]])
                #save the result with parameter values
                                            
            else:
                #add array of results
                R=np.array([[Ks, Kr, cr, cd, alpha, T0,res]])
                result=np.concatenate((result, R), axis=0)
            counter+=1
                
        #save csv file and graph
        np.savetxt(fname+'.csv',result, delimiter=',', header='Ks, Kr, cr, cd, alpha, T0, class', fmt='%.6f') 
        print(fate)
        fate_matrix[tri,:]=fate           
    if case==0:       
        np.savetxt('parameter_sweep_MC_total_case0.csv',fate_matrix, delimiter=',', header='cl0,l1,cl2,cl3,cl4', fmt='%d')
    else:
        np.savetxt('parameter_sweep_MC_total_case1.csv',fate_matrix, delimiter=',', header='cl0,l1,cl2,cl3,cl4,cl5', fmt='%d')
    Plot(case)    
        
def Plot(case):
    #plot the csv file and plot the results
    #rCo vs sCh case0
    #sCO vs rCh case1
    fname=str('parameter_sweep_MC_total_case%d'%(case))
    #read csv file
    data = np.loadtxt(fname+'.csv',delimiter=",",skiprows=1)   
    """
    data include
    n\class 0 1 2 3 4 5
    1 
    2
    3
    """
    left=np.array([1,2,3])#change here when analyzing more Hill parameters
    
    if case==0:
        class0=data[:,0].T#failure in mono culture
        class1=data[:,1].T#failure in invasion
        class2=data[:,2].T#only exclusion is stable
        class3=data[:,3].T#only coexistence is stable
        class4=data[:,4].T#both unstable
        
        #only exclusion or ceoxistence
        plt.figure(figsize=(8,6))
        p2=plt.bar(left, class2, color="blue")
        p3=plt.bar(left, class3, bottom=class2,  color="red")
        p4=plt.bar(left, class4, bottom=class2,  color="gray")
        plt.xlabel("Hill coefficient",fontsize=14)
        plt.xticks([1,2,3],["1","2","3"])
        plt.ylabel('frequency',fontsize=14)
        plt.legend((p2[0], p3[0], p4[0]), ("exclude", "coexist", "both unstable"))
        plt.savefig(fname+"_part.pdf")
        plt.show()
        
        #full bar plot
        plt.figure(figsize=(8,6))
        P0=plt.bar(left, class0, color="yellow")
        bt=class0
        P1=plt.bar(left, class1, bottom=bt, color="black")
        bt+=class1
        P2=plt.bar(left, class2, bottom=bt, color="blue")
        bt+=class2
        P3=plt.bar(left, class3, bottom=bt, color="red")
        bt+=class3
        P4=plt.bar(left, class4, bottom=bt, color="gray")
        plt.xlabel("Hill coefficient",fontsize=14)
        plt.xticks([1,2,3],["1","2","3"])
        plt.ylabel('frequency',fontsize=14)
        plt.legend((P0[0],P1[0], P2[0], P3[0], P4[0]), ("failure in mono-culture","failure in invasion","exclude", "coexist", "both unstable"))
        plt.savefig(fname+"_full.pdf")
        plt.show()
    elif case==1:
        # in this case two equilibrium may be stable at the same time
        class0=data[:,0].T#failure in mono culture
        class1=data[:,1].T#failure in invasion
        class5=data[:,5].T#both stable
        class2=data[:,2].T#only exclusion is stable
        class3=data[:,3].T#only coexistence is stable
        class4=data[:,4].T#both unstable
        
        #only exclusion or ceoxistence
        plt.figure(figsize=(8,6))
        p5=plt.bar(left, class5, color="purple")
        p2=plt.bar(left, class2, bottom=class5, color="blue")
        p3=plt.bar(left, class3, bottom=class5+class2,  color="red")
        p4=plt.bar(left, class4, bottom=class5+class2+class3,  color="gray")
        plt.xlabel("Hill coefficient",fontsize=12)
        plt.xticks([1,2,3],["1","2","3"])
        plt.ylabel('frequency',fontsize=12)
        plt.legend((p5[0], p2[0],p3[0], p4[0], ), ("both stable", "exclude", "coexist", "both unstable"))
        plt.savefig(fname+"_part.pdf")
        plt.show()
        
        #full bar plot
        plt.figure(figsize=(8,6))
        P0=plt.bar(left, class0, color="yellow")
        bt=class0
        P1=plt.bar(left, class1, bottom=bt, color="black")
        bt+=class1
        P5=plt.bar(left, class5, bottom=bt, color="purple")
        bt+=class5
        P2=plt.bar(left, class2, bottom=bt, color="blue")
        bt+=class2
        P3=plt.bar(left, class3, bottom=bt, color="red")
        bt+=class3
        P4=plt.bar(left, class4, bottom=bt, color="gray")
        plt.xlabel("Hill coefficient",fontsize=12)
        plt.xticks([1,2,3],["1","2","3"])
        plt.ylabel('frequency',fontsize=12)
        plt.legend((P0[0],P1[0], P5[0], P2[0], P3[0],P4[0]), ("failure in mono-culture","failure in invasion","both stable", "exclude","coexist", "both unstable"))
        plt.savefig(fname+"_full.pdf")
        plt.show()