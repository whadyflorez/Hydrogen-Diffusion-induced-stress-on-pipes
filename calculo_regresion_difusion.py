#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 07:06:37 2024

@author: whadyimac
Calculo dependencia de la temperatura del coeficiente de difusion
"""
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

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

p_initialize=np.array([-1.0,1.0])

res=minimize(f,p_initialize,method='Powell',tol=1.0e-8)  

p=res.x
p[0]=np.exp(p[0])

print('Ecuacion de Arrhenius:')
print(f"{p[0]:.3E} exp(-{p[1]:.3E}/(R T))")



    
