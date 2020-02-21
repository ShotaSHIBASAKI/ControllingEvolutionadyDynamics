#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 15:23:31 2019
Figures for main text (Figs.1B-E, 2B-E, 4, and 5)
@author: shota
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib import colors
import sympy as sym


#fixed parameters
dmax=1 #maximum mortality rate caused by toxin
fmax=0.5 #maximum degradation efficiency
alpha=0.1#dilution rate
Kd=0.2 # median-effet density for the toxic degradation
Newton_size=10#step size for variable init condition in Newton method
rmax=1 #maximam intrisic growth rate
"""
Codes for Figs 1B-E
Run Fig1_ExampleODE(case)
"""

def DynamicsCoop(x, t, T0, cost, K, n, r):
    """
    ODE of mono-culture of cooperators to obtain the T^* and x_1^* 
    where mutation never occurs
    """
    
    #defune the degradation effect
    deg=fmax*x[1]/(x[1]+Kd)
    D=dmax*x[0]**n/(x[0]**n+K**n)
    #define ODEs
    y=np.zeros([np.size(x)])
    y[0]=alpha*T0-deg*x[0]-alpha*x[0] #dt/dt
    y[1]=x[1]*(r*(1-cost)*(1-x[1])-D-alpha)#d Co/dt
    
    return y

def DynamicsCheat(x, t, T0, cost, K, n, r):
    """
    ODE of mono-culture of cheaterss to obtain the T^* and x_1^* 
    where mutation never occurs
    """
    D=dmax*x[0]**n/(x[0]**n+K**n)
    #defune the degradation effect
    deg=0
    #define ODEs
    y=np.zeros([np.size(x)])
    y[0]=alpha*T0-deg*x[0]-alpha*x[0] #dt/dt
    y[1]=x[1]*(r*(1-cost)*(1-x[1])-D-alpha)#d Co/dt
    
    return y

def DynamicsCo(x, t, T0, cost_co, cost_ch, K_co, K_ch, n, r):
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
def DynamicsCC(x, t, T0, cost1, cost2, K1, K2, n, r):
    """
    Co-culture of the two coopertors 
    with the environmental feedback thorough changing the toxic concentration
    """
    y=np.zeros([np.size(x)])
    D=np.zeros([2])   
    #define fitnss
    D[0]=dmax*x[0]**n/(K1**n+x[0]**n) #cooperator
    D[1]=dmax*x[0]**n/(K2**n+x[0]**n) #cheater   
    #degradation
    deg=fmax*(x[1]+x[2])/(x[1]+x[2]+Kd)   
    #ODE of eco-evo dynamics
    y[0]=alpha*T0-deg*x[0]-alpha*x[0] #dt/dt
    y[1]=x[1]*(r*(1-cost1)*(1-x[1]-x[2])-D[0]-alpha)#d Co1/dt
    y[2]=x[2]*(r*(1-cost2)*(1-x[1]-x[2])-D[1]-alpha)#d Co2/dt
       
    return y

def Fig1_ExampleODE(case):
    T0=0.3
    time=np.linspace(0, 500, 50001)
    if case==0:
        """
        Mono-culture of Cooperators
        Same parameter values in Phase plain analysis
        """
        cd=0.15#cost of detofication
        Ks=0.3#median toxin concentration        
        n=3#Hill coeeficient
        init=np.array([T0, 10**(-6)])
        sol=odeint(DynamicsCoop, init, time, args=(T0, cd, Ks, n, rmax))
        fname='Mono-sCo.pdf'
        #plot the date
        plt.plot(time, sol[:,0], color='k', label='toxin')
        plt.plot(time, sol[:,1], color='b', linestyle='-', label='sCo')
        plt.xlim(0,time[-1])
    elif case==1:
        """
        sCh invade and exclude sCo
        same parameter values above
        begin from the equil of sCo
        """
        cd=0.15#cost of detofication
        Ks=0.3#median toxin concentration        
        n=3#Hill coeeficient
        fname='Co-sChsCo.pdf'
        time2=np.linspace(0, 2000, 200001)
        #mono- of sCo
        init=np.array([T0, 10**(-6)])
        sol=odeint(DynamicsCoop, init, time, args=(T0, cd, Ks, n, rmax))
        #invasion of sCh
        init_evo=np.array([sol[-1,0],sol[-1,1]-10**(-6), 10**(-6)])
        sol=odeint(DynamicsCo, init_evo, time2, args=(T0, cd, 0, Ks, Ks,  n, rmax))
        plt.plot(time2, sol[:,0], color='k', label='toxin')
        plt.plot(time2, sol[:,1], color='b', linestyle='-', label='sCo')
        plt.plot(time2, sol[:,2], color='r', linestyle='-', label='sCh')
        plt.xlim(0,time2[-1])
        
    elif case==2:
        """
        rCo invade and coexist with sCh
        """
        cd=0.15#cost of detofication
        cr=0.3
        Ks=0.3#median toxin concentration
        Kr=0.6        
        n=3#Hill coeeficient; if n is large, it is more likely to coexist (see Mon†´Carlo sim.)
        fname='Co-rCosCh-coexistence.pdf'
        time2=np.linspace(0, 2000, 200001)
        #mono- of sCh
        init=np.array([T0, 10**(-6)])
        sol=odeint(DynamicsCheat, init, time, args=(T0, 0, Ks, n, rmax))
        #invasion of sCh
        init_evo=np.array([sol[-1,0], 10**(-6), sol[-1,1]-10**(-6)])
        sol=odeint(DynamicsCo, init_evo, time2, args=(T0,cd+cr, 0, Kr, Ks,  n, rmax))
        plt.plot(time2, sol[:,0], color='k', label='toxin')
        plt.plot(time2, sol[:,2], color='r', linestyle='-', label='sCh')
        plt.plot(time2, sol[:,1], color='c', linestyle='--', label='rCo')
        plt.xlim(0,time2[-1])
    elif case==3:
        """
        sCo invade and exclude rCo
        from the equilibrium of the mono-clture of rCo
        """
        cd=0.15#cost of detofication
        cr=0.3
        Ks=0.3#median toxin concentration
        Kr=0.6        
        n=3#Hill coeeficient; if n is large, it is difficult to exclude sCh (see Mon†´Carlo sim.)
        
        fname='Co-sCorCo-exclusion.pdf'
        time2=np.linspace(0, 2000, 200001)
        #mono- of rCo
        init=np.array([T0, 10**(-6)])
        sol=odeint(DynamicsCoop, init, time, args=(T0, cr+cd, Kr, n, rmax))
        #invasion of sCh
        init_evo=np.array([sol[-1,0], 10**(-6), sol[-1,1]-10**(-6)])
        sol=odeint(DynamicsCC, init_evo, time2, args=(T0,cd, cd+cr, Ks, Kr,  n, rmax))
        plt.plot(time2, sol[:,0], color='k', label='toxin')
        plt.plot(time2, sol[:,1], color='b', linestyle='-', label='sCo')
        plt.plot(time2, sol[:,2], color='c', linestyle='--', label='rCo')
        plt.xlim(0,time2[-1])
    else:
        """
        rCo invade and exclude sCh
        small T0
        """
        T0=0.2
        cd=0.15#cost of detofication
        cr=0.2#if cr is larger, it is difficult to invade
        Ks=0.2#median toxin concentration
        Kr=0.6        
        n=1#Hill coeeficient; if n is large, it is more likely to coexist (see Mon†´Carlo sim.)
        fname='Co-rCosCh-exclusion.pdf'
        time2=np.linspace(0, 2000, 200001)
        #mono- of 
        init=np.array([T0, 10**(-6)])
        sol=odeint(DynamicsCheat, init, time, args=(T0, 0, Ks, n, rmax))
        #invasion of sCh
        init_evo=np.array([sol[-1,0], 10**(-6), sol[-1,1]-10**(-6)])
        sol=odeint(DynamicsCo, init_evo, time2, args=(T0,cd+cr, 0, Kr, Ks,  n, rmax))
        plt.plot(time2, sol[:,0], color='k', label='toxin')
        plt.plot(time2, sol[:,2], color='r', linestyle='-', label='sCh')
        plt.plot(time2, sol[:,1], color='c', linestyle='--', label='rCo')
        plt.xlim(0,time2[-1])
    plt.xlabel('time', fontsize=14)
    plt.ylabel('toxin or cell', fontsize=14)
    plt.ylim(0,1)
    plt.legend(loc='lower right', fontsize=12)
    if case==0:
        plt.plot(time, T0*np.ones(np.size(time)), color='k', linestyle='--')
    else:
        plt.plot(time2, T0*np.ones(np.size(time2)), color='k', linestyle='--')
    #plt.text(10, T0+0.05, r'$T_{in}$', fontsize=14)
    plt.savefig(fname)
    plt.show()
"""
Codes for Figs2B-E, G-J
Run Fig2_InvasionState(case,step_size)
"""
def Newton(x0, T0, r, K, n, alpha):
    """
    Newton method to find the equilibrium in mono-culture of cooperators
    """
    x=sym.Symbol('x')
    T=alpha*T0/(alpha+fmax*x/(x+Kd))
    F=r*(1-x)-(dmax*T**n/(T**n+K**n)+alpha)
    dF=sym.diff(F,x)
    #Newton method for obtaining equilibrium x
    k=50
    x_now=x0
    for i in range(k):
        x_next=x_now-F.subs(x,x_now)/dF.subs(x,x_now)
        if abs(x_next-x_now)<10**(-6):
            break
        else:
            x_now=x_next
    if abs(F.subs(x, x_next))<10**(-3):      
        return x_next
    else:
        # bad convergence
        return -1
def Stability(x_equil,T_equil, Tin, r, K, n, alpha):
    #stability in the case of mono-culture
    x=sym.Symbol('x')
    T=sym.Symbol('T')
    delta=dmax*T**n/(T**n+K**n)
    d_delta=sym.simplify(sym.diff(delta,T))
    #d_delta=dmax*n*T**(n-1)*K**n/(T**n+K**n)**2
    lhs=-r/d_delta.subs(T, T_equil)
    f=fmax*x/(x+Kd)
    df=sym.diff(f,x)
    rhs=-alpha*Tin*df.subs(x,x_equil)/((alpha+f.subs(x,x_equil))**2)
    if lhs<rhs:
        #stable
        return 0
    else:
        return 1
        
def Stability2(x1, x2, T_d, Tin, r1, r2, K1, n, alpha):
    #stability analysis when Co coexist Ch
    x=sym.Symbol('x')
    T=sym.Symbol('T')
    delta=dmax*T**n/(T**n+K1**n)
    d_delta=sym.simplify(sym.diff(delta,T))
    lhs=-(r1*x1+r2*x2)/(x1*d_delta.subs(T, T_d))
    f=fmax*x/(x+Kd)
    df=sym.diff(f,x)
    rhs=-alpha*Tin*df.subs(x,x1)/((alpha+f.subs(x,x1))**2)
    if lhs<rhs:
        #stable
        return 0
    else:
        return 1
    
def Quad(A, B, C, case):
    #Root of Quadratic function
    D=B**2-4*A*C
    if D>0:
        tau_p=(-B+D**0.5)/(2*A)
        tau_m=(-B-D**0.5)/(2*A)
        
        if case==0:
            # in this case smaller root can be stable
            return min(tau_p, tau_m)
        else:
            # in this case larger root can be stable
            return max(tau_p, tau_m)
    else:
        return -1
def Resident_Invasion(r_res, K_res, r_inv, K_inv, Tin_array, n, resident_coop):
    """
    boundaries where resident and invader persist, respectively
    r_res, K_res: r and K of resident
    r_inv, K_inv: r and K of invader
    Tin_array: given values of toxin flowing into the system
    n: Hill coefficient
    """
    
    #boundary 1: resident persists in absence of invader below this boundary
    boundary1=r_res-dmax*Tin_array**n/(Tin_array**n+K_res**n)
    
    #boundary 2: resident is not invaded by invader
    boundary2=np.ones(np.size(Tin_array))
    for i in range(np.size(Tin_array)):
        if resident_coop==0:
            # resident is not cooperator
            T_eq=Tin_array[i]
        else:
            stab=1
            while stab==1:
                # resident is cooperator
                # obtain toxin concentration when there exist only resident
                x0=np.random.random()
                x=Newton(x0, Tin_array[i], r_res, K_res, n, alpha)
                T_eq=alpha*Tin_array[i]/(alpha+fmax*x/(Kd+x))
                stab=Stability(x,T_eq, Tin_array[i], r_res, K_res, n, alpha)
        d_res=dmax*T_eq**n/(T_eq**n+K_res**n)
        d_inv=dmax*T_eq**n/(T_eq**n+K_inv**n)
        boundary2[i]=(r_inv*d_res-r_res*d_inv)/(r_res-r_inv)
        
        
        
    return [boundary1, boundary2]
def Fig2_InvasionState(case,step_size):
    """
    Case: which phenotype should be a resisdent?
    Step_size; step size of toxin concnetraion flowing into the system
    
    Classfy the results of invasion analysis according to alpha and T0
    case0: sCo vs rCo
    case1: sCo vs rCh
    case2: sCh vs rCo
    case3: sCh vs rCh
    
    Note: sCo vs sCh (or rCo vs rCh) is obvious
    """
    #set array of alpha and Tin, and x (init for Newton method)
    Tin_array=np.linspace(0,1,step_size+1)  
    #define parameter values
    r=rmax#maximum growth rate
    cd=0.15#cost of detofication
    cr=0.2#cost of resistance
    Ks=0.3#median toxin concentration for sensitive
    Kr=0.6#for resistant
    n=3#Hill coeeficient
    
    #sensitive strain
    if case==0 or case==1:
        #sCo
        r1=r*(1-cd)
        K1=Ks
        resident_coop1=1
    else:
        #sCh
        r1=r*1
        K1=Ks
        resident_coop1=0
    if case==0 or case==2:
       #rCo
        r2=r*(1-cd-cr)
        K2=Kr
        resident_coop2=1
    else:
        #rCh
        r2=r*(1-cr)
        K2=Kr
        resident_coop2=0
        
        
    if case==0:
        tname1='sCo → rCo'
        tname2='rCo → sCo'
        gname1='StateSpace_sCo-rCo.pdf'
        gname2='StateSpace_rCo-sCo.pdf'
    elif case==1:
        tname1='sCo → rCh'
        tname2='rCh → sCo'
        gname1='StateSpace_sCo-rCh.pdf'
        gname2='StateSpace_rCh-sCo.pdf'
    elif case==2:
        tname1='sCh → rCo'
        tname2='rCo → sCh'
        gname1='StateSpace_sCh-rCo.pdf'
        gname2='StateSpace_rCo-sCh.pdf'
    else:
        tname1='sCh → rCh'
        tname2='rCh → sCh'
        gname1='StateSpace_sCh-rCh.pdf'
        gname2='StateSpace_rCh-sCh.pdf'
         
    # simulation one: sensitive strain is resident and resistant strain is invader
    [boundary11, boundary12]=Resident_Invasion(r1, K1, r2, K2, 
    Tin_array, n, resident_coop1)
    
    # simulation two: resistnat strain is resident and sensitive strain is invader
    [boundary21, boundary22]=Resident_Invasion(r2, K2, r1, K1, 
    Tin_array, n, resident_coop2)
     
    plt.subplots(figsize=(6,6))
    plt.plot(Tin_array, boundary11, color='k', linestyle='-')
    plt.plot(Tin_array, boundary12, color='k', linestyle='--')
    Min=np.array([boundary11, boundary12])
    plt.fill_between(Tin_array, np.amin(Min, axis=0),0, color=(0.984375  , 0.55078125, 0.3828125), alpha=1.0) #invasion succed
    plt.fill_between(Tin_array, boundary11, boundary12, 
                     where=boundary11>boundary12, color=(0.99609375, 0.84765625, 0.18359375), alpha=1.0) #invasion fail
    plt.xlim(0, 1.0)
    plt.ylim(0, 1.0)
    ticks=np.array([0.0, 0.5, 1.0])
    locs=np.array([0.0, 0.5, 1.0])
    plt.yticks(locs,ticks, fontsize=18)
    plt.xticks(locs,ticks, fontsize=18)
    #plt.text(0.5, 0.9, tname1, fontsize=28)
    plt.savefig(gname1,bbox_inches='tight', pad_inches=0.1)
    plt.show()
    
    Min=np.array([boundary21, boundary22])
    plt.subplots(figsize=(6,6))
    plt.plot(Tin_array, boundary21, color='k', linestyle='-')
    plt.plot(Tin_array, boundary22, color='k', linestyle='--')
    plt.fill_between(Tin_array, boundary21, boundary22,
                     where=boundary21>boundary22, color=(0.984375  , 0.55078125, 0.3828125), alpha=1.0) #invasion succeed
    plt.fill_between(Tin_array, np.amin(Min, axis=0),0, color=(0.99609375, 0.84765625, 0.18359375), alpha=1.0) # invasion fails
    plt.xlim(0, 1.0)
    plt.ylim(0, 1.0)
    ticks=np.array([0.0, 0.5, 1.0])
    locs=np.array([0.0, 0.5, 1.0])
    plt.yticks(locs,ticks, fontsize=18)
    plt.xticks(locs,ticks, fontsize=18)
    #plt.text(0.5, 0.9, tname2, fontsize=28)
    plt.savefig(gname2,bbox_inches='tight', pad_inches=0.1)
    plt.show()
    

"""
Fig 4
Run Fig4_OptMultiPlot(cd, cr, Ks, Kr, ns, nr)
"""
def Productivity(alpha, T, Tin):
    return alpha*(Tin-T)

def OptPhi1(cost,K, n):
    """
    Optimize productivity in the mono-culture of cooperator
    cost: total cost cooperator has
    K; resistnace of cooperators to the toxin
    n: Hill coefficient
    """
    step_size=50#step size of T and alpha
    r=rmax*(1-cost)
    alpha_array=np.linspace(0,1,step_size+1)
    Tin_array=np.linspace(0,1,step_size+1)
    phi=np.zeros([step_size+1, step_size+1])
    x0_array=np.linspace(0,1,Newton_size)
    max_phi=0
    alpha_max=0#argmax phi(alpha, Tin)
    Tin_max=0
    for i in range(np.size(alpha_array)):
        alpha=alpha_array[i]         
        for j in range(np.size(Tin_array)):
            Tin=Tin_array[j]
            if alpha==0 or Tin==0:
                phi[i,j]=0#obviously productivity is zero and we should avoid devision  by zero
            else:
                T_equil=1
                #x_equil=0
                #check the condition
                W0=r-(dmax*Tin**n/(Tin**n+K**n)+alpha)
                if W0>0:
                    #print('success')
                    #calculate the euilibrium using Newton method
                    for k in range (Newton_size):
                        x0=x0_array[k]
                        x=Newton(x0, Tin, r, K, n, alpha)
                        if x>0:
                            #print(x)
                            T=alpha*Tin/(alpha+fmax*x/(Kd+x))
                            stab=Stability(x,T, Tin, r, K, n, alpha)
                            if stab==0:
                                #print('stable')
                                #linearly stable
                                if T_equil>T:
                                    #current T is lower than the current T_equil
                                    T_equil=T
                                    #x_equil=x
                    #calculate the productivity
                    #print(x_equil, T_equil)
                    phi[i,j]=Productivity(alpha, T_equil, Tin)
                    if phi[i,j]>max_phi:
                        max_phi=phi[i,j]
                        alpha_max=alpha
                        Tin_max=Tin
                else:
                    #out side of the constraint
                    phi[i,j]=0
    """
    #plot the results
    fig, ax = plt.subplots(figsize=(6,5))
    heatmap = ax.pcolor(phi, cmap=plt.cm.Blues)
    ax.xaxis.tick_bottom()
    xlocs=np.array([0.5,step_size/2+0.5, step_size+0.5])
    xticks=np.array([0.0, 0.5, 1.0])
    plt.xticks(xlocs,xticks)
    ylocs=np.array([0.5,step_size/2+0.5, step_size+0.5])
    yticks=np.array([0.0, 0.5,1.0])
    plt.yticks(ylocs,yticks)
    heatmap.set_clim(0,max_phi)
    plt.colorbar(heatmap)
    print(max_phi, Tin_max, alpha_max)
    
    
    #drow the boudary
    B=r-dmax*(Tin_array)**n/((Tin_array)**n+K**n)
    #B=r-dmax*(Tin_array*step_size+0.5)**n/((Tin_array*step_size+0.5)**n+K**n)
    plt.plot(Tin_array*step_size+0.5, B*step_size+0.5, color='gray', linestyle='--',linewidth=3)
    
    plt.xlabel(r'$T_{in}$',fontsize=14)
    plt.ylabel(r'$\alpha$',fontsize=14)
    plt.ylim(0, step_size+1)
    plt.scatter(Tin_max*step_size+0.5,alpha_max*step_size+0.5, color='r',s=300,marker='*', alpha=0.8)
    plt.savefig('Opt-productivity-Co.pdf')
    plt.show()
    """
    return [phi, max_phi, Tin_max, alpha_max]
def OptPhi2(case, cr, cd, Ks, Kr, n):
    """
    case 0 : rCo vs sCh
    case 1 : sCo vs rCh
    cr: cost of resistance
    cd: cost of detoxification
    Ks; medeian toxic concentration to sensitive
    Kr: to resistsnt
    n: Hill coefficient
    """
    step_size=50#step size of T and alpha
    if case==0:
        #rCo vs sCH
        c1=cr+cd
        K1=Kr
        c2=0
        K2=Ks
    elif case==1:
        c1=cd
        K1=Ks
        c2=cr
        K2=Kr
    else:
        print("Wrong case value. case should be 0 or 1")
        return 1
    r1=rmax*(1-c1)
    r2=rmax*(1-c2)
    if r1<0 or r2<0:
        print("Costs are too large")
        return 1
    #set the array of alpha, T0 etc...
    alpha_array=np.linspace(0,1,step_size+1)
    Tin_array=np.linspace(0,1,step_size+1)
    phi=np.zeros([step_size+1, step_size+1])
    max_phi=0
    alpha_max=0#argmax phi(alpha, Tin)
    Tin_max=0
    for i in range (step_size+1):
        alpha=alpha_array[i]
        for j in range(np.size(Tin_array)):
            Tin=Tin_array[j]
            if alpha==0 or Tin==0:
                phi[i,j]=0#obviously productivity is zero and we should avoid devision  by zero
            else:
                #check constrain 1: whether cheaters exist in the mono-culture
                Wch=r2-dmax*Tin**n/(Tin**n+K2**n)-alpha
                if Wch>0:
                   x_ch=1-(dmax*Tin**n/(Tin**n+K2**n)+alpha)/r2 #density of cheaters
                   #constrain 2 : invasion of cooperator succeed or not
                   Wco=r1*(1-x_ch)-dmax*Tin**n/(Tin**n+K1**n)-alpha
                   if Wco>0:
                       #calculate the equilibrium (at most one)
                       A=(dmax+alpha)*(c2-c1)
                       B=dmax*(1-c1)*K1**n-dmax*(1-c2)*K2**n+alpha*(c2-c1)*(K1**n+K2**n)
                       C=alpha*(c2-c1)*K1**n*K2**n
                       tau=Quad(A,B,C,case)
                       if tau>0:
                           T_d=tau**(1/n)#toxic concentration at the stable equilibrium
                           x_co_d=alpha*Kd*(Tin-T_d)/(fmax*T_d-alpha*(Tin-T_d))
                           x_ch_d=1-x_co_d-(dmax*T_d**n/(T_d**n+K2**n)+alpha)/r2
                           #check the stability of this equilibrium
                           stab=Stability2(x_co_d, x_ch_d, T_d, Tin, r1, r2, K1, n, alpha)
                           if stab==0:
                               phi[i,j]=Productivity(alpha, T_d, Tin)
                               if phi[i,j]>max_phi:
                                   max_phi=phi[i,j]
                                   alpha_max=alpha
                                   Tin_max=Tin
                           else:
                               phi[i,j]=0#unstable coexistence
                       else:
                            phi[i,j]=0 # unable to coexist
                   else:
                       phi[i,j]=0 # unable to invade
                    
                else:
                    phi[i,j]=0#Ch does not exist
    """
    #plot the results
    fig, ax = plt.subplots(figsize=(6,5))
    heatmap = ax.pcolor(phi, cmap=plt.cm.Blues)
    ax.xaxis.tick_bottom()
    xlocs=np.array([0.5,step_size/2+0.5, step_size+0.5])
    xticks=np.array([0.0, 0.5, 1.0])
    plt.xticks(xlocs,xticks)
    ylocs=np.array([0.5,step_size/2+0.5, step_size+0.5])
    yticks=np.array([0.0, 0.5,1.0])
    plt.yticks(ylocs,yticks)
    heatmap.set_clim(0,max_phi)
    plt.colorbar(heatmap)
    print(max_phi, Tin_max, alpha_max)
    #plot some boundaries
    """
    #boudary for existence of Ch
    B1=r2-dmax*Tin_array**n/(Tin_array**n+K2**n)
    plt.plot(Tin_array*step_size+0.5, B1*step_size+0.5, color='gray', linestyle='--',linewidth=3)
    #boundary for successful invasion
    B2=r1*(dmax*Tin_array**n/(Tin_array**n+K2**n)+alpha)/r2-dmax*Tin_array**n/(Tin_array**n+K1**n)
    plt.plot(Tin_array*step_size+0.5, B2*step_size+0.5, color='gray', linestyle='--',linewidth=3)
    #Other boudaries are very messy
    """
    #end of the plot the boudaries
    plt.scatter(Tin_max*step_size+0.5,alpha_max*step_size+0.5, color='r',s=300,marker='*', alpha=0.8)
    plt.ylim(0, step_size+1)
    plt.xlabel(r'$T_{in}$',fontsize=14)
    plt.ylabel(r'$\alpha$',fontsize=14)
    fname=str('Opt-productivity-CoCh-case%d' %(case))
    plt.savefig(fname+'.pdf')
    plt.show()
    """
    return [phi, Tin_max, alpha_max]
def Fig4_OptMultiPlot(cd, cr, Ks, Kr, ns, nr):
    #plot the Opt of four cases on the same panel
    step_size=50#step size of T and alpha
    
    #only sC0
    [sCo_phi, s_phi,sCo_T, sCo_alpha]=OptPhi1(cd,Ks,ns)
    
    #only rCo
    [rCo_phi, r_phi,rCo_T, rCo_alpha]=OptPhi1(cd+cr,Kr, nr)

    #note max phi is larger at the mono-culture
    #sCo vs rCh
    [sCorCh_phi,Co1_T, Co1_alpha]=OptPhi2(1, cr, cd, Ks, Kr, ns)
    #rCo vs sCh
    [rCosCh_phi,Co0_T, Co0_alpha]=OptPhi2(0, cr, cd, Ks, Kr, nr)
    
    #plot the results
    fig, ax= plt.subplots(figsize=(6,5))
    #LeftUpper sCo
    heatmap=ax.pcolor(sCo_phi, cmap=plt.cm.Blues)
    ax.xaxis.tick_bottom()
    xlocs=np.array([0.5,step_size/2+0.5, step_size+0.5])
    xticks=np.array([0.0, 0.5, 1.0])
    plt.xticks(xlocs,xticks)
    ylocs=np.array([0.5,step_size/2+0.5, step_size+0.5])
    yticks=np.array([0.0, 0.5,1.0])
    plt.yticks(ylocs,yticks)
    heatmap.set_clim(0,s_phi)
    clb=plt.colorbar(heatmap)
    clb.ax.set_title(r'$\phi$',fontsize=18)
    #drow the boudary
    r=rmax*(1-cd)
    Tin_array=np.linspace(0,1,step_size+1)
    B=r-dmax*(Tin_array)**ns/((Tin_array)**ns+Ks**ns)
    #B=r-dmax*(Tin_array*step_size+0.5)**n/((Tin_array*step_size+0.5)**n+K**n)
    plt.plot(Tin_array*step_size+0.5, B*step_size+0.5, color='gray', linestyle='--',linewidth=3)
    
    plt.xlabel(r'$T_{in}$',fontsize=18)
    plt.ylabel(r'$\alpha$',fontsize=18)
    plt.ylim(0, step_size+1)
    plt.scatter(sCo_T*step_size+0.5,sCo_alpha*step_size+0.5, color='r',s=300,marker='*', alpha=0.8)
    plt.savefig('Opt-productivity-sCo.pdf')
    plt.show()
    
    #RightUpper rCo
    fig, ax= plt.subplots(figsize=(6,5))
    heatmap=ax.pcolor(rCo_phi, cmap=plt.cm.Blues)
    ax.xaxis.tick_bottom()
    xlocs=np.array([0.5,step_size/2+0.5, step_size+0.5])
    xticks=np.array([0.0, 0.5, 1.0])
    plt.xticks(xlocs,xticks)
    ylocs=np.array([0.5,step_size/2+0.5, step_size+0.5])
    yticks=np.array([0.0, 0.5,1.0])
    plt.yticks(ylocs,yticks)
    heatmap.set_clim(0,r_phi)
    clb=plt.colorbar(heatmap) 
    clb.ax.set_title(r'$\phi$',fontsize=18)
    #drow the boudary
    r=rmax*(1-cd-cr)
    Tin_array=np.linspace(0,1,step_size+1)
    B=r-dmax*(Tin_array)**nr/((Tin_array)**nr+Kr**nr)
    #B=r-dmax*(Tin_array*step_size+0.5)**n/((Tin_array*step_size+0.5)**n+K**n)
    plt.plot(Tin_array*step_size+0.5, B*step_size+0.5, color='gray', linestyle='--',linewidth=3)
    
    plt.xlabel(r'$T_{in}$',fontsize=18)
    plt.ylabel(r'$\alpha$',fontsize=18)
    plt.ylim(0, step_size+1)
    plt.scatter(rCo_T*step_size+0.5,rCo_alpha*step_size+0.5, color='r',s=300,marker='*', alpha=0.8)
    plt.savefig('Opt-productivity-rCo.pdf')
    plt.show() 
    
    #Leftbottom sCo and rCh
    fig, ax= plt.subplots(figsize=(6,5))
    heatmap=ax.pcolor(sCorCh_phi,  cmap=plt.cm.Blues)
    ax.xaxis.tick_bottom()
    xlocs=np.array([0.5,step_size/2+0.5, step_size+0.5])
    xticks=np.array([0.0, 0.5, 1.0])
    plt.xticks(xlocs,xticks)
    ylocs=np.array([0.5,step_size/2+0.5, step_size+0.5])
    yticks=np.array([0.0, 0.5,1.0])
    plt.yticks(ylocs,yticks)
    heatmap.set_clim(0,r_phi)
    clb=plt.colorbar(heatmap)
    clb.ax.set_title(r'$\phi$',fontsize=18)
    plt.xlabel(r'$T_{in}$',fontsize=18)
    plt.ylabel(r'$\alpha$',fontsize=18)
    plt.ylim(0, step_size+1)
    plt.scatter(Co1_T*step_size+0.5,Co1_alpha*step_size+0.5, color='r',s=300,marker='*', alpha=0.8)
    plt.savefig('Opt-productivity-rCosCh.pdf')
    plt.show()
    #Rightbottom sCo and rCh
    fig, ax= plt.subplots(figsize=(6,5))
    heatmap = ax.pcolor(rCosCh_phi, cmap=plt.cm.Blues)
    ax.xaxis.tick_bottom()
    xlocs=np.array([0.5,step_size/2+0.5, step_size+0.5])
    xticks=np.array([0.0, 0.5, 1.0])
    plt.xticks(xlocs,xticks)
    ylocs=np.array([0.5,step_size/2+0.5, step_size+0.5])
    yticks=np.array([0.0, 0.5,1.0])
    plt.yticks(ylocs,yticks)
    heatmap.set_clim(0,s_phi)
    clb=plt.colorbar(heatmap)
    clb.ax.set_title(r'$\phi$',fontsize=18)     
    plt.xlabel(r'$T_{in}$',fontsize=18)
    plt.ylabel(r'$\alpha$',fontsize=18)
    plt.ylim(0, step_size+1)
    plt.scatter(Co0_T*step_size+0.5,Co0_alpha*step_size+0.5, color='r',s=300,marker='*', alpha=0.8)
    plt.savefig('Opt-productivity-sCorCh.pdf')
    plt.show()
    
"""
Fig.5
The code for Markov chain analysis with Dynamic programming 
(i.e., calculating the optimal intoruction rates of cooperators)
is performed by C language code named main.c in DynaProgMC
After running the code above, use this function
"""
def Fig5_DPMC(model):
    # read the csv file for plot the results
    if model==0:
        #eorgidc MC
        fname='optDPMC_ergodic.csv'
        gname='optDPMC_ergodic.pdf'
    else:
        #Non-eorgidc MC
        fname='optDPMC_Nonergodic.csv'
        gname='optDPMC_Nonergodic.pdf'
    data=np.loadtxt(fname, delimiter=",", skiprows=1)
    data_infty=data[-1,1:3]#subtract last row
    data=np.delete(data, -1,0)
    plt.plot(data[:,0], data[:,1], color='blue', label='$m_1$')
    plt.plot(data[:,0], data[:,2], color='cyan', label='$m_2$')
    plt.legend(loc='lower left')
    
    #data for t -> infinite
    infty1=np.ones(np.size(data[:,1]))*data_infty[0]
    infty2=np.ones(np.size(data[:,2]))*data_infty[1]
    
    plt.plot(data[:,0], infty1, color='blue', linestyle='--')
    plt.plot(data[:,0], infty2, color='cyan', linestyle='--')
    plt.xscale("log")
    plt.xlabel('time step (a.u.)',fontsize=14)
    plt.ylim(0,1.0)
    plt.ylabel('optimal '+'$m_i$',fontsize=14)
    plt.savefig(gname)
    plt.show()
