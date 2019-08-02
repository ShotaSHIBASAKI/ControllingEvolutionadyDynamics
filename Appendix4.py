#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 14:21:29 2019
Code for Appendix 4
@author: shota
"""

import sympy as sym
from sympy import init_printing

def AnalyzeA4():
    init_printing()
    """
    From Eq (A.52)to Eq (A.54)
    Proove that the eigenvaues of the 4x4 Jacobian matrix are
    diagonal elements of the Jacobian when there exist 3 strategies i,j, and k
    (e.g. rCo, sCh, and rCh)
    """
    
    a=sym.Symbol('a')#dilution rate instead of alpha
    T=sym.Symbol('T')#toxin concentration
    df=sym.Symbol("f'")#df/dxi
    xk=sym.Symbol('x_{k}')#density of strategy k
    ri=sym.Symbol('r_{i}')#intrinsic growth rate i
    rj=sym.Symbol('r_{j}')#intrinsic growth rate i
    rk=sym.Symbol('r_{k}')#intrinsic growth rate i
    di=sym.Symbol('d_{i}')#death rate of strategy i
    dj=sym.Symbol('d_{j}')#death rate of strategy i
    dk=sym.Symbol('d_{k}')#death rate of strategy i
    ddk=sym.Symbol("d_{k}'")# d d_k/dT
    lam=sym.Symbol('lambda') #eigenvalue of Jacobian
    
    #define Jacobian matirx
    J=sym.Matrix([[-a, -T*df, 0, 0],
                  [0, ri*(1-xk)-di-a,0,0],
                  [0,0,rj*(1-xk)-dj-a,0],
                  [-xk*ddk, -rk*xk, -rk*xk, -rk*xk]])
    
    print('Jacobian matrix')
    display(J)
    I=lam*sym.eye(4)
    
    #Charatcteristic equation
    C=sym.det(I-J)
    #display(C)
    C=sym.expand(C)
    C=sym.collect(C, lam)
    C=sym.factor(C)
    print("As we can see, the eivengaleus are the diagonal elements of the Jacobian matrix above")
    display(C)
    return
