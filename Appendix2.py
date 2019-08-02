#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 16:55:22 2019
Code for Appendix 2
@author: shota
"""

import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
from sympy import init_printing

#fixed parameters
dmax=0.5 #maximum mortality rate caused by toxin
fmax=0.5 #maximum degradation
r=1#maximum intrinsic growth rate
alpha=0.1#dilution rate

def Quadratic(a,b,c):
    #root of quadratic function
    x=sym.Symbol('x')
    func=a*x**2+b*x+c
    xl=(-b+(b**2-4*a*c)**(1/2))/(2*a)
    xs=(-b-(b**2-4*a*c)**(1/2))/(2*a)
    return [xs, xl]

def AnalyzeA2():
    init_printing()
    #Analytical calculation for A2
    n=sym.Symbol('n')#Hill coefficent
    t=sym.Symbol('t^{*}') #T^n where T is toxin concentration and n is Hill coefficient
    #note use t instead of tau here
    r=sym.Symbol('r')#maximum intrinsic growth rate
    ci=sym.Symbol('c_{i}')#total cost for i
    cj=sym.Symbol('c_{j}')#total cost for j
    ri=r*(1-ci)#intrisic growth rate of i
    rj=r*(1-cj)#intrisic growth rate of j
    a=sym.Symbol('a')#dilution rate, instead of alpha
    dmax=sym.Symbol('d_{max}')#maximum death rate
    Ki=sym.Symbol('K_{i}')#meadian-toxin-concentration for i
    Kj=sym.Symbol('K_{j}')#meadian-toxin-concentration for j    
    Wi=ri/(a+dmax*t/(t+Ki**n))#realtive fitness of i
    Wj=rj/(a+dmax*t/(t+Kj**n))#realtive fitness of i
    
    #define the condition for invasion
    F=Wi-Wj
    print("Strategy i can invade the population of j iff")
    display(F>0)
    #simplified the equation
    F=sym.simplify(F*(a+dmax*t/(t+Ki**n))*(a+dmax*t/(t+Kj**n))/r*(Ki**n+t)*(Kj**n+t))
    F=sym.expand(F)
    F=sym.collect(F,t)
    #simplify each coefficient
    F_coeff0=sym.simplify(F.coeff(t,0))
    F_coeff1=sym.simplify(F.coeff(t,1))
    F_coeff1=sym.factor(F_coeff1)
    F_coeff2=sym.simplify(F.coeff(t,2))
    F_coeff2=sym.factor(F_coeff2)
    F_new=F_coeff2*t**2+F_coeff1*t+F_coeff0
    print("Now inequality (A.13) is obtained. Each coefficient of quadratic function F is as below:")
    display(F_coeff2,+F_coeff1, +F_coeff0)#coefficients for inequality A.13
    
    
    #Example 1 When i=rCO invades j=sCh
    Ks=sym.Symbol("K_{s}")#meadian-toxin-concentration for sensitive
    Kr=sym.Symbol("K_{r}")#meadian-toxin-concentration for resistant
    c=sym.Symbol("c")#total cost c=cr+cd
    F1=F_new.subs([(ci,c),(cj,0), (Ki,Kr), (Kj,Ks)])
    print("\n As shown in inequality (A.16), rCo can invade sCh iff F1>0 whose coefficeints are")
    display(F1.coeff(t,2),F1.coeff(t,1),F1.coeff(t,0))
    #Condition of c where rCo can invade sCh
    #inequalities(A.20) and (A.21a) are not so difficult althogh they are messy...
    #Inequality (A.21b)
    A21b=(dmax*(1-c)*Kr**n-dmax*Ks**n-a*c*(Kr**n+Ks**n))#ignoring denominator (2*(d+a)*c)
    A21b=sym.expand(A21b)
    A21b=sym.collect(A21b,c)
    A21b_coeff1=sym.simplify(A21b.coeff(c,1))
    A21b_coeff0=sym.simplify(A21b.coeff(c,0))
    A21b_left=A21b_coeff0/-A21b_coeff1
    print("\n The left hand side of inequality (A.21b) is")
    display(A21b_left>c)
    A21b=(dmax*(1-c)*Kr**n-dmax*Ks**n-a*c*(Kr**n+Ks**n))-t*(2*(dmax+a)*c)
    A21b=sym.expand(A21b)
    A21b=sym.collect(A21b,c)
    A21b_coeff1=sym.simplify(A21b.coeff(c,1))
    A21b_coeff0=sym.simplify(A21b.coeff(c,0))
    A21b_right=A21b_coeff0/-A21b_coeff1
    print("\n The right hand side of inequality (A.21b) is")
    display(c>A21b_right)
   
    #Inequality (A.21c)
    """
    D=A^2-B^2>0 (B>0)
    <=> A>B>0 or A<-B<0
    but it is enough to analyze only A>B>0 case
    """
    sqrtD=dmax*(1-c)*Kr**n-dmax*Ks**n-a*c*(Kr**n+Ks**n)-2*c*(Kr*Ks)**(n/2)*(a*(a+dmax))**(1/2)
    sqrtD=sym.expand(sqrtD)
    sqrtD=sym.collect(sqrtD,c)
    sqrtD_coeff0=sym.simplify(sqrtD.coeff(c,0))
    sqrtD_coeff1=sym.simplify(sqrtD.coeff(c,1))  
    A21c_right=sqrtD_coeff0/-sqrtD_coeff1
    display(c<A21c_right)
    
    
    #Inequality (A.23)
    tin=sym.Symbol('t_{in}')
    A=dmax*Kr**n+a*(Kr**n+Ks**n)
    LHS=dmax*(Kr**n-Ks**n)/((dmax+a)*t+A+a*Kr**n*Ks**n/t)
    RHS=dmax*(Kr**n-Ks**n)/(2*(dmax+a)*t+A)
    A23=sym.expand(sym.simplify((LHS-RHS)*((dmax+a)*t+A+a*Kr**n*Ks**n/t)*(2*(dmax+a)*t+A))/(dmax*(Kr**n-Ks**n)))
    A23=sym.collect(A23, t)
    print("\n Inequality (A.23) is easily checked as below")
    display((sym.simplify(A23*t)).subs(t, tin)>0)
    
    #Example 2 When i=sCo invades j=rCh
    cd=sym.Symbol('c_{d}')
    cr=sym.Symbol('c_{r}')
    t=sym.Symbol('t^*')
    F_new=F_coeff2*t**2+F_coeff1*t+F_coeff0
    F2=F_new.subs([(ci,cd),(cj,cr), (Ki,Ks), (Kj,Kr)])
    F2=sym.collect(sym.expand(F2),t)
    print("\n Coefficients of F2 shown inequality(A.25) are")
    display(F2.coeff(t,2), F2.coeff(t,1), F2.coeff(t,0))
    
    
    return

def FigA2(cd, cr, Kr, Ks, n):
    """
    Draing Figure A2
    Given the values of 
    cd (cost of degradation),
    cr (cost of resistnace),
    Kr (median-effect concentration of resistant),
    Ks (median-effect conncentration of sensitive),
    and n (Hill"s coefficient),
    one can find the threshold of toxin concentration  T_th
    where cooperators' inasion succeed
    """
    #set the important parameters
    Sum_c=(cr+cd)
    Delta_c=(cr-cd)
    if Delta_c<0:
        print("Error in costs; cr>cd")
        return 1

    tau=np.linspace(0,1,1001) # tau = T^n
    Zero=np.zeros([np.size(tau)])
    tname=str(r'$c_d=%.1f, c_r=%.1f, K_r=%.1f, K_s=%.1f, n=%d$' %(cd, cr, Kr, Ks, n))
    """
    Quadratic functions
    F: rCo invades sCh
    G: sCo invades rCh
    """
    F=(dmax+alpha)*Sum_c*tau**2-(dmax*(1-Sum_c)*Kr**n-dmax*Ks**n-alpha*Sum_c*(Kr**n+Ks**n))*tau+alpha*Sum_c*Kr**n*Ks**n
    G=(dmax+alpha)*Delta_c*tau**2+(dmax*(1-cd)*Ks**n-dmax*(1-cr)*Kr**n+alpha*Delta_c*(Kr**n+Ks**n))*tau+alpha*Delta_c*Kr**n*Ks**n
    
    #find the threshold
    tau_F_small, tau_F_large= Quadratic((dmax+alpha)*Sum_c,-(dmax*(1-Sum_c)*Kr**n-dmax*Ks**n-alpha*Sum_c*(Kr**n+Ks**n)),alpha*Sum_c*Kr**n*Ks**n)
    print(tau_F_small, tau_F_large)
    tau_G_small, tau_G_large= Quadratic((dmax+alpha)*Delta_c,dmax*(1-cd)*Ks**n-dmax*(1-cr)*Kr**n+alpha*Delta_c*(Kr**n+Ks**n),alpha*Delta_c*Kr**n*Ks**n)
    print(tau_G_small, tau_G_large)
    #plot the range of tau
    plt.plot(tau,-F, color='b',label=r'$F_{\rm{rCo}, \rm{sCh}}(\tau_{\rm{in}})$')
    plt.plot(tau, G, color='r', label=r'$F_{\rm{sCo},\rm{rCh}}(\tau_{\rm{in}})$')
    plt.plot(tau,Zero, color='k',linestyle="--")
    plt.fill_between(tau, -F, Zero, where=-F>Zero,facecolor='b', alpha=0.5)
    plt.fill_between(tau, Zero, G, where=G>Zero,facecolor='r', alpha=0.5)
    plt.scatter(tau_F_small,0,color='w',edgecolors='b',s=80)
    plt.scatter(tau_F_large,0,color='w',edgecolors='b',s=80)
    plt.scatter(tau_G_small,0,color='w',edgecolors='r',s=80)
    plt.scatter(tau_G_large,0,color='w',edgecolors='r',s=80)
    #plt.text(tau_F_small,0.025,str(r'$\tau=%.4f$' %(tau_F_small)),fontsize=14)
    #plt.text(tau_G_small,-0.02,str(r'$\tau=%.4f$' %(tau_G_small)),fontsize=14)
    #plt.text(tau_F_large-0.05,0.015,str(r'$\tau=%.4f$' %(tau_F_large)),fontsize=14)
    #plt.text(tau_G_large,-0.015,str(r'$\tau=%.4f$' %(tau_G_large)),fontsize=14)
    plt.xlim(0,1)
    plt.ylim(-0.02,0.02)
    plt.legend(loc='upper left',fontsize=12)
    plt.xlabel(r'$\tau_{\rm{in}}=T^{n}_{\rm{in}}$', fontsize=14)
    plt.yticks([0,0.02])
    plt.ylabel(r'$F(\tau_{\rm{in}})$', fontsize=14)
    plt.title(tname,fontsize=12)
    plt.savefig('range-tau.pdf')
    plt.show()
    
    #plot the range of T
    rCo=np.array([tau_F_small**(1/n), tau_F_large**(1/n)])
    sCo1=np.array([0, tau_G_small**(1/n)])
    sCo2=np.array([tau_G_large**(1/n),1])
    plt.plot(rCo, np.ones([np.size(rCo)])*0.1, color='b', marker='o', markerfacecolor='w')
    plt.plot(sCo1, np.ones([np.size(rCo)])*0.15, color='r', marker='o', markerfacecolor='w')
    plt.plot(sCo2, np.ones([np.size(rCo)])*0.15, color='r', marker='o', markerfacecolor='w')
    plt.xlim(0,1)
    plt.ylim(0.09,0.18)
    plt.text(tau_F_small**(1/n),0.11,'rCo invades sCh', fontsize=14, color='b')
    plt.text(tau_G_large**(1/n)-0.15,0.16,'sCo invades rCh', fontsize=14, color='r')
    plt.xlabel(r'$T_{\rm{in}}$',fontsize=16)
    plt.yticks([])
    plt.title(tname,fontsize=12)
    plt.savefig('range-T.pdf')
    plt.show()  
    
def FigA3(n,Ks,Kr):
    """
    Codes for Fig A3 (given values of Ks and Kr)
    Plot the range of cr and cr where two types of invasion succeed
    """
    ntext=str('_n%d_Ks%d_Kr_%d' %(n,10*Ks, 10*Kr))
    cr=np.linspace(0,1-alpha,1001)
    #rCo invade sCh
    #determine the uper boundary of tau
    a1=r-alpha
    if dmax-a1<a1*Ks**n:
        tau_h=1
    else:
        tau_h=a1*Ks**n/(dmax-a1)
    print(tau_h)
    #condition for the invasion is affected by tau_h
    Th=(alpha*Ks**n*Kr**n/(dmax+alpha))**(1/2)#threshold
    print(Th)
    A=dmax*Kr**n+alpha*(Kr**n+Ks**n)
    cd1a=-cr+dmax*(Kr**n-Ks**n)/(A+2*(alpha*(alpha+1)*Kr**n*Ks**n))  
    plt.plot(cr, cd1a, color='black')
    if tau_h<=Th:
        flag=1
        cd1b=-cr+dmax*(Kr**n-Ks**n)/(2*(dmax+alpha)+A)
        plt.plot(cr, cd1b, color='blue')
        cd1c=-cr+dmax*(Kr**n-Ks**n)/((dmax+alpha)*tau_h+A+alpha*Kr**n*Ks**n/tau_h)
        plt.plot(cr, cd1c, color='black')
        plt.fill_between(cr,cd1a,cd1c,where=cd1a>cd1c,facecolor='k',alpha=0.5)
        plt.fill_between(cr,cd1b,0,where=cd1b>0,facecolor='b',alpha=0.5)
    else:
        flag=0
    
    
    #sCO invade sCh
    tau_h=np.zeros(np.size(cr))
    a2=np.zeros(np.size(cr))
    for i in range (np.size(cr)):
        a2[i]=r*(1-cr[i])-alpha
        if dmax-a2[i]<a2[i]*Kr**n:
            tau_h[i]=1
        else:
            tau_h[i]=a2[i]*Kr**n/(dmax-a2[i])
    B=(dmax+alpha)*tau_h+alpha*(Kr**n+Ks**n)+alpha*Kr**n*Ks**n/tau_h
    cd2=((B+dmax*Kr**n)*cr+dmax*(Ks**n-Kr**n))/(B+dmax*Ks**n)
    plt.plot(cr, cd2, color='black')
        
        
    if flag==0:
        plt.fill_between(cr,cd1a,0,where=cd1a>0,facecolor='b',alpha=0.5)
        plt.fill_between(cr,cd1a,cd2,where=cd1a>cd2,facecolor='w',alpha=1)
        plt.xlim(0,1-alpha/r)
        plt.ylim(0,1)
        plt.title(str(r'$K_r= %.1f, K_s=%.1f, n= %d$' %(Kr, Ks, n)), fontsize=14)
        plt.ylabel('cost of detoxification '+r'$c_d$', fontsize=16)
        plt.xlabel('cost of resistance 'r'$c_r$', fontsize=16)
        plt.savefig('range-cd-cr'+ntext+'.pdf')
    plt.show()