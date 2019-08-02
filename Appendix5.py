#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 15:10:27 2019
Code for Appendix 5
@author: shota
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
#fixed parameter
fmax=0.5#maximum detoxification efficiency
dmax=1.0#maximum death rate
Kd=0.2#medain effect cooperator density


def death_rate(T, n, K):
    return dmax*T**n/(T**n + K**n)
def DynamicsMut(x, t, Tin, r, K, n, alpha, mut1, mut2):
    
    """
    Mutation matrix Q is intrduced
    mut1 changing cooperator <-> cheater
    mut2 changing sensitive <-> resistance
    r and K should be array
    """
    
    Q=np.zeros([4,4])#mutation matrix
    Q[0,0]=1-mut1-mut2+mut1*mut2
    Q[0,1]=mut1*(1-mut2)
    Q[0,2]=(1-mut1)*mut2
    Q[0,3]=mut1*mut2
    Q[1,0]=mut1*(1-mut2)
    Q[1,1]=1-mut1-mut2+mut1*mut2
    Q[1,2]=mut1*mut2
    Q[1,3]=(1-mut1)*mut2
    Q[2,0]=(1-mut1)*mut2
    Q[2,1]=mut1*mut2
    Q[2,2]=1-mut1-mut2+mut1*mut2
    Q[2,3]=mut1*(1-mut2)
    Q[3,0]=mut1*mut2
    Q[3,1]=(1-mut1)*mut2
    Q[3,2]=mut1*(1-mut2)
    Q[3,3]=1-mut1-mut2+mut1*mut2
    
    y=np.zeros([np.size(x)])
    D=np.zeros([4])   
    #define fitnss
    for i in range (4):
        D[i]=death_rate(x[0], n, K[i])
    #degradation of toxin
    deg=fmax*(x[1]+x[3])/(x[1]+x[3]+Kd)   
    #ODE of eco-evo dynamics
    y[0]=alpha*Tin-deg*x[0]-alpha*x[0] #dt/dt
    y[1]=(1-sum(x[1:5]))*np.dot(r*x[1:5],Q[:,0])-(D[0]+alpha)*x[1]#d sCo/dt
    y[2]=(1-sum(x[1:5]))*np.dot(r*x[1:5],Q[:,1])-(D[1]+alpha)*x[2]#d sCh/dt
    y[3]=(1-sum(x[1:5]))*np.dot(r*x[1:5],Q[:,2])-(D[2]+alpha)*x[3]#d rCo/dt
    y[4]=(1-sum(x[1:5]))*np.dot(r*x[1:5],Q[:,3])-(D[3]+alpha)*x[4]#d rCoh/dt
    return y


def ODE_Mutation(mu1, mu2):
    """
    Simulation for ODE given mutation rates
    
    """
    #Given Parameters
    r=1#maximum growth rate
    cd=0.15#cost of producing degradation
    cr=0.3 #cost of resistance
    K_s=0.3#median effect toxic for sensitive
    K_r=0.5#median effect toxic for resistant
    alpha=0.1#delution rate
    Tin=0.3#toxin flowing into the system
    n=3 #hill coefficient
    fname=str('replicator-mutator-rebuild_mu1-%d_mu2-%d' %(mu1, mu2))
    #assuming the equilibrlium state of the mono-cukture of sCo
    init=np.array([Tin,0.01,0.0,0.0,0.0])
    if mu1==0:
        mut1=0
    else:    
        mut1=10**(-mu1)
    if mu2==0:
        mut2=0
    else:
        mut2=10**(-mu2)
    r_array=r*np.array([1-cd, 1, 1-cr-cd, 1-cr])#growth rate of sCo, sCh, rCo, rCh
    K_array=np.array([K_s, K_s, K_r, K_r])#median toxin concentration of sCo, sCh, rCo, rCh
    time=np.linspace(0,500,50001)
    
    #no  mutation test
    sol=odeint(DynamicsMut, init, time, args=(Tin, r_array, K_array, n, alpha,  0, 0))
    """
    #plot the results
    plt.plot(time, sol[:,0])
    plt.xlabel('time', fontsize=14)
    plt.ylabel('toxic concentration', fontsize=14)
    plt.ylim(0,Tin)
    plt.show()
    
    
    plt.plot(time, sol[:,1], color='b', linestyle='-', label='sCo')
    plt.plot(time, sol[:,2], color='r', linestyle='-', label='sCh')
    plt.plot(time, sol[:,3], color='c', linestyle='--', label='rCo')
    plt.plot(time, sol[:,4], color='m', linestyle='--', label='rCh')
    plt.xlabel('time', fontsize=14)
    plt.ylabel('density', fontsize=14)
    plt.ylim(0,1)
    plt.legend(loc='lower right')
    plt.show()
    """
    #now let start the simulation with mutation
    #init=sol[-1,:]
    #print(init)
    time=np.linspace(0,1200,120001)    
    sol=odeint(DynamicsMut, init, time, args=(Tin, r_array, K_array, n, alpha, mut1, mut2))
    """
    #plot the results
    plt.plot(time, sol[:,0])
    plt.xlabel('time', fontsize=14)
    plt.ylabel('toxic concentration', fontsize=14)
    plt.ylim(0,Tin)
    #plt.title(r'$\mu_1=$'+str('%.4f and '%(mut1))+r'$\mu_2$='+str('%.4f' %(mut2)))
    plt.title(str(r'$\mu_1=10^{-%d},\quad \mu_2=10^{-%d}$' %(mu1, mu2)))
   # plt.savefig(fname+'_toxin.pdf')
    plt.show()
    
    
    plt.plot(time, sol[:,1], color='b', linestyle='-', label='sCo')
    plt.plot(time, sol[:,2], color='r', linestyle='-', label='sCh')
    plt.plot(time, sol[:,3], color='c', linestyle='--', label='rCo')
    plt.plot(time, sol[:,4], color='m', linestyle='--', label='rCh')
    plt.xlabel('time', fontsize=14)
    plt.ylabel('density', fontsize=14)
    plt.ylim(0,1)
   # plt.title(r'$\mu_1=$'+str('%.4f and '%(mut1))+r'$\mu_2$='+str('%.4f' %(mut2)))
    plt.title(str(r'$\mu_1=10^{-%d},\quad \mu_2=10^{-%d}$' %(mu1, mu2)))
    plt.legend(loc='lower right')
    #plt.savefig(fname+'_cell.pdf')
    plt.show()
    """
    return [time,sol]

def FigA5():
    #Draw Fig. A5  
    Tin=0.3
    Mu1=np.array([1.0,4.0])#array for the power of mu1
    Mu2=np.array([1.0,4.0])#array for the power of mu2
    count=0
    plt.figure(figsize=(12, 9))
    for i in range (np.size(Mu1)):
        mu1=Mu1[i]
        for j in range (np.size(Mu2)):
            mu2=Mu2[j]
            time, sol=ODE_Mutation(mu1, mu2)
            count+=1
            plt.subplot(np.size(Mu1), np.size(Mu2), count)
            plt.xlim(0, time[-1])
            plt.ylim(0,1)
            plt.plot(time, sol[:,0],color='k', label='toxin')
            plt.plot(time, sol[:,1], color='b', linestyle='-', label='sCo')
            plt.plot(time, sol[:,2], color='r', linestyle='-', label='sCh')
            plt.plot(time, sol[:,3], color='c', linestyle='--', label='rCo')
            plt.plot(time, sol[:,4], color='m', linestyle='--', label='rCh')
            plt.legend(loc='upper right',fontsize=12)
            plt.plot(time, np.ones([np.size(time)])*Tin,color='k', linestyle='--',label='toxin')
            plt.title(str(r'$\mu_1=10^{-%d},\quad \mu_2=10^{-%d}$' %(mu1, mu2)), fontsize=16)
    plt.savefig('replicator-mutator.pdf')
    plt.show()
    