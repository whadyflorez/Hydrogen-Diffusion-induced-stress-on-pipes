#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 07:06:37 2024

@author: whadyimac
Calculo dependencia de la temperatura del coeficiente de difusion
"""
import numpy as np
from scipy.optimize import minimize

T=np.array([313.15,303.15])
D=np.array([1.68e-10,1.40e-10])
R=8.31446261815324
n=np.shape(T)[0]

def m(T,p):
    y=p[0]-p[1]/(R*T)
    return y

def f(p):
    E=0.0
    for i in range(n):
      E+=(np.log(D[i])-m(T[i],p))**2
    return E  

p_initialize=np.array([0.5,0.5])

res=minimize(f,p_initialize,method='CG')  

p=res.x
p[0]=np.exp(p[0])


    
