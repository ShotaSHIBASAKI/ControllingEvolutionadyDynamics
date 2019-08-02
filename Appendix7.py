#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 15:36:25 2019
Codes for Apendix 7
Analyzing Markov chian
For drawing Fig A.10,use FigMaintext.py
@author: shota
"""

import sympy as sym
from sympy import init_printing

def AnalyzeA7():
    init_printing()
    """
    This code calculates the stationary distribution  Eq (A.63)
    when a Markov chain is given by Fig, 5 in the main text.
    """
    #Set probability distribution pi
    pi1=sym.Symbol('pi^*_{1}')#mono-culture of sCo
    pi2=sym.Symbol('pi^*_{2}')#sCo->sCh
    pi3=sym.Symbol('pi^*_{3}')#mono-culture of sCh
    pi4=sym.Symbol('pi^*_{4}')#sCh -> rCo
    pi5=sym.Symbol('pi^*_{5}')#mono-culture of rCo
    pi6=sym.Symbol('pi^*_{6}')#rCo -> rCh
    pi7=sym.Symbol('pi^*_{7}')#mono-culture of rCh
    pi8=sym.Symbol('pi^*_{8}')#rCh->sCo
    pi9=sym.Symbol('pi^*_{9}')#sCo ->rCo
    pi10=sym.Symbol('pi^*_{10}')#rCo-> sCo
    pi11=sym.Symbol('pi^*_{11}')#sCo -> sCo
    pi12=sym.Symbol('pi^*_{12}')#sCh -> sCh (failure of introduction of sCo)
    pi13=sym.Symbol('pi^*_{13}')#rCO -> rCo
    pi14=sym.Symbol('pi^*_{14}')#rCh -> rCh (failure of introduction of rCo)
    Pi=sym.Matrix([pi1, pi2, pi3, pi4, pi5, pi6, pi7,
                   pi8, pi9, pi10, pi11, pi12, pi13, pi14]).T
    
    display(Pi)
    #set transition matrix P
    mu1=sym.Symbol('mu1')#mutation rate from Co to Ch
    m1=sym.Symbol('m1')#introduction rate of sCo
    m2=sym.Symbol('m2')#introduction rate of rCo
    
    P=sym.Matrix([
    [1-mu1-m1-m2,mu1,0,0,0,0,0,0,m2,0,m1,0,0,0] ,       
    [0,0,1,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,1-m1-m2,m2,0,0,0,0,0,0,0,m1,0,0],
    [0,0,0,0,1,0,0,0,0,0,0,0,0,0],        
    [0,0,0,0,1-mu1-m1-m2,mu1,0,0,0,m1,0,0,m2,0],
    [0,0,0,0,0,0,1,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,1-m1-m2,m1,0,0,0,0,0,m2],
    [1,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,1,0,0,0,0,0,0,0,0,0],
    [1,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [1,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,1,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,1,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,1,0,0,0,0,0,0,0]])
    print('Transition matrix P of Fig.5 is given by as below')
    display(P)
    print('Sow the right-hand side in Eqs (A.62b) -(A.62o) which is optained by Pi*P')
    RHS=Pi*P
    display(RHS.T)
    """
    Note it would be time consuming to directly solve the simultaneous linear equations.
    to avoind this problem, we show that Pi_st-Pi_st*P=0 
    where Pi_st is the stationary probability distribution given by Eq(A.63)
    """
    #convert pi_i to a function of pi_1 where i=2,3,...,14
    Pi_st=Pi.subs([(pi2, mu1*pi1), (pi3,mu1/m2*pi1), (pi4, mu1*pi1), 
                   (pi5,(mu1+m2)/(mu1+m1)*pi1),(pi6,mu1*(mu1+m2)/(mu1+m1)*pi1),
                   (pi7,mu1*(mu1+m2)/m1/(mu1+m1)*pi1),(pi8,mu1*(mu1+m2)/(mu1+m1)*pi1),
                   (pi9,m2*pi1), (pi10,m1*(mu1+m2)/(mu1+m1)*pi1),(pi11,m1*pi1), 
                   (pi12,mu1*m1/m2*pi1), (pi13,m2*(mu1+m2)/(mu1+m1)*pi1),
                   (pi14,mu1*m2*(mu1+m2)/m1/(mu1+m1)*pi1)])
    print('The stationary distribution Pi_st is as shown in Eq (A.63)')
    display(Pi_st.T)
    #sum pi =1
    print('pi^*_1 is asshown in Eq (A.64)')
    pi1_st=sym.simplify(1/sum(Pi_st)*pi1)
    display(pi1_st)
    """
    A=1/(pi1_st/(m1*m2*(m1+mu1)))#denominator of Eq (A.63)
    display(sym.expand(A))
    print(sym.latex(A))
    """
    print('If this distribution is really a stationary distribution, Pi_st-Pi_st*P should be zero.')
    print('Indeed, Pi_st - Pi_st*P equals to')
    V=Pi_st-Pi_st*P
    V=V.subs([(pi1, pi1_st)])
    display(sym.simplify(V.T))
    print('Therefore, Eq (A.63) is a stationary distibution')
    
    
    