#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 14:37:07 2019
Code for Appendix 1
@author: shota
"""

import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
from sympy import init_printing

#fixed parameters
dmax=1.0 #maximum death rate caused by toxin
fmax=0.5#maximum degradation rate
alpha=0.1#dilution rate

def Newton(x0, T0, Kd, cd, Ks, n, alpha, r):
    #Newton method
    x=sym.Symbol('x')
    T=alpha*T0/(alpha+fmax*x/(x+Kd))
    D=dmax*T**n/(Ks**n+T**n)
    F=1-(alpha+D)/(r*(1-cd))-x
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

def Jacobian (T, x, Kd, cd, Ks, n, alpha, r):
    #Analyzing linear stability with Jacobian matrix
    D=dmax*T**n/(T**n+Ks**n)
    dD=dmax*n*T**(n-1)*Ks**n/(T**n+Ks**n)**2
    J=np.array([[-(fmax*x/(x+Kd)+alpha),-T*fmax*Kd/(x+Kd)**2],[-x*dD,r*(1-cd)*(1-2*x)-D-alpha]])
    
    if J[0,0]+J[1,1]<0 and J[0,0]*J[1,1]-J[0,1]*J[1,0]>0:
        return 0 #stable
    else:
        return 1
def Sys(x,T0, Kd, cd, Ks, n, alpha, r):
    #define ODE for quiver plot
    y=np.zeros([np.size(x)])
    deg=fmax*x[1]/(x[1]+Kd)
    
    D=dmax*x[0]**n/(x[0]**n+Ks**n)
    
    y[0]=alpha*T0-alpha*x[0]-deg*x[0]
    y[1]=x[1]*(r*(1-cd)*(1-x[1])-D-alpha)
    
    return y


def AnalyzeA1():
    init_printing()
    #some analytical calculations that might be complex in Appendix 1
    x=sym.Symbol('x_{i}')#density of microbes
    T=sym.Symbol('T')#toxin concentration
    T0=sym.Symbol('T_{in}')
    r=sym.Symbol('r_i')#intrinsic growth rate
    a=sym.Symbol('a')#dilution rate, instead of alpha
    f=sym.Symbol('f(x_{i})')#function of detoxification
    df=sym.Symbol("f'(x_{i})")#df/dx
    d=sym.Symbol('d(T)')#function of deah rate, instead of delta
    dd=sym.Symbol("d'(T)")#dd/dT 
    
    #Jacobian matrix
    J=sym.Matrix([[-(f+a), -T*df],[-dd*x, r*(1-2*x)-d-a]]) #gieven by Eq(A.5)
    print("In equation (A.5), the Jacobian matix, J is")
    display(J)
    
    #With the non-trivial equilibrium
    x_eq=sym.Symbol('x^{*}_{i}')
    T_eq=sym.Symbol('T^{*}')
    f_eq=sym.Symbol('f(x^{*}_{i})')
    df_eq=sym.Symbol("f'(x^{*}_i)")
    d_eq=sym.Symbol('d(T^{*})')
    dd_eq=sym.Symbol("d'(T^{*})") 
    T_eq=a*T0/(a+f_eq) #Eq (9a)
    x_eq=1-(d_eq+a)/r  #Eq (9b) <=> r(1-x_eq)= a +d(T)
    J[1,1]=-r*x
    print("\n Nontirivial equilibrium (T^*, x^*) from Eqs (9a) and (9b) are")
    display(T_eq, x_eq) 
    J=sym.simplify(J.subs([(x, x_eq), (T, T_eq), (f, f_eq), (df, df_eq), (d, d_eq), (dd, dd_eq)]))
    
    #det J >0 with non-trivial equilibrium
    detJ = sym.simplify(J.det())
    print('\n The diterminant of the Jacobian matrix (detJ) with the non-tirivial equilibrium is')
    display(detJ)
    print(r"Note that {a+d(T^*)-r}/[r{a+f(x^*)}] is negative.")
    print("Now, we obtain inequality (A.6) because detJ>0 <=>")
    display(sym.simplify(detJ/(a+d_eq-r)/(r*(a+f_eq))*r**2/dd_eq)<0)
    
    return
    
def FigA1():
    #Drawing fig A.1 (phase plain analysis)
    n=3 #Hill coefficient
    cd=0.15 # cost for detoxification
    r=1 #maximum growth rate
    Kd=0.2 # MM parameter for degradation
    Ks=0.3 #median-effect toxin conc. for sensitive
    T0=0.8 #toxin conc. flowing into chemostat
    
    x=np.linspace(0,1,1001)
    T=np.linspace(0,1,1001)
    F=alpha*T0/(alpha+fmax*x/(x+Kd))#nullcline for dT/dt
    G=1-(alpha+dmax*T**n/(T**n+Ks**n))/(r*(1-cd))
    Z=np.zeros(np.size(T))
    #plot the nullcline
    plt.plot(x, F, color='red')
    plt.plot(G, T, color='blue')
    plt.plot(Z,T,color='blue',linewidth=5)
    
    #plot the equilibrium point
    #search the equilibria with Newton method
    X=np.linspace(0,1,500)
    for i in range(np.size(X)): 
        x_eq=Newton(X[i],T0, Kd, cd, Ks, n, alpha, r)
        flag=0
        if i==0:
            Xeq_array=np.array([x_eq])
        else:
            for j in range (np.size(Xeq_array)):
                if abs(x_eq-Xeq_array[j])<0.03:
                    flag=1
                    break
            if flag==0:
                    Xeq_array=np.concatenate([Xeq_array,np.array([x_eq])],axis=0)
                   
    #plot and the stability
    for j in range (np.size(Xeq_array)):
        x_eq=Xeq_array[j]
        if x_eq>0:
            T_eq=alpha*T0/(alpha+fmax*x_eq/(x_eq+Kd))
            stab=Jacobian(T_eq, x_eq, Kd, cd, Ks, n, alpha, r)
            if stab==0:
                plt.scatter(x_eq,T_eq, c='k',edgecolors ='k',s=150)
            else:
                plt.scatter(x_eq,T_eq, c='w',edgecolors ='k',s=150)
            #plot the tangents at the equilibrium point
            tan1=-alpha*T0*fmax*Kd/(x_eq+Kd)**2/(alpha+fmax*x_eq/(x_eq+Kd))**2*(x-x_eq)+T_eq
            plt.plot(x, tan1, color='r', linestyle='--')
            tan2=-r/dmax*(Ks**n+T_eq**n)**2/(n*T_eq**(n-1)*Ks**n)*(x-x_eq)+T_eq
            plt.plot(x, tan2, color='b', linestyle='--')
            
    #plot the trivial equilibrium
    x_eq=0        
    T_eq=T0
    stab=Jacobian(T_eq, x_eq, Kd, cd, Ks, n, alpha, r)
    if stab==0:
        plt.scatter(x_eq,T_eq, c='k',edgecolors ='k',s=150)
    else:
        plt.scatter(x_eq,T_eq, c='w',edgecolors ='k',s=150)
    #plot the quiverplot
    a = np.linspace(0.02, 1, 15)#T
    b = np.linspace(0.02, 1, 15)#x

    A , B  = np.meshgrid(a, b)  
    DA=np.zeros([np.size(a), np.size(b)])
    DB=np.zeros([np.size(a), np.size(b)])
    for i in range(np.size(a)):
        for j in range(np.size(b)):
            DA[i,j], DB[i,j] = Sys([A[i,j], B[i,j]],T0, Kd, cd, Ks, n, alpha, r)#growth rate at each point
    M = (np.hypot(DA, DB))                        # norm growth rate 
    M[ M == 0] = 1.                                 # avoid zero division errors 
    DA /= M                                        # normalize each arrows
    DB /= M
    plt.quiver(b, a, DB.T, DA.T, M, pivot='mid')
    
    plt.xlabel('density of cooperators',fontsize=16)
    plt.ylabel('toxin concentration',fontsize=16)
    plt.xlim(0,1)
    plt.ylim(0,1)
    #plt.savefig('Popdynamics_PP2.pdf')
    plt.show()