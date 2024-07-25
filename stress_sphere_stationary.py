#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 08:19:53 2024
esfuerzo termico solo en esferas
@author: whadymacbook2016
"""
import numpy as np
import matplotlib.pyplot as plt

b=2.0e-2
a=b-3.0e-3
Ta=0.1
Tb=0.0
E=10e9 #Pa
nu=0.3
alfa=2.0
n=100 #nodes
dr=(b-a)/(n-1)


T=np.zeros(n)
sigma_r=np.zeros(n)
sigma_t=np.zeros(n)
epsi_r=np.zeros(n)
epsi_t=np.zeros(n)
disp=np.zeros(n)
r=np.linspace(a,b,n)

def temp(r):
    y=(Tb*b-Ta*a)*(b-a)+(Ta-Tb)*a*b/((b-a)*r)
    return y

def u(r):
    disp=(alfa/(1-nu))*(Ta-Tb)*a/(2*(b**3-a**3))*\
    ((1+nu)*b*(a**2+a*b+b**2)+2*(nu*a**2-a**2-nu*a*b-nu*b**2)*r-(1+nu)*a**2*b**3/r**2)+Tb*alfa*r
    return disp

def e_r(r):
    y=(alfa/(1-nu))*(Ta-Tb)*a/(b**3-a**3)*\
    ((nu*a**2-a**2-nu*a*b-nu*b**2)+(1+nu)*a**2*b**3/r**3)+Tb*alfa
    return y

def e_t(r):
    y=(alfa/(1-nu))*(Ta-Tb)*a/(2*(b**3-a**3)*r)*\
    (((1+nu)*b*(a**2+a*b+b**2)+2*(nu*a**2-a**2-nu*a*b-nu*b**2)*r-(1+nu)*a**2*b**3/r**2))+Tb*alfa
    return y

def s_r(r):
    y=(E*alfa*nu/(1-nu))*(Ta-Tb)*a*b/(b**3-a**3)*(r-a)*(r-b)*(r*a+r*b+a*b)
    return y

def s_t(r):
    y=(E*alfa/(2*(1-nu)))*(Ta-Tb)*a*b/(b**3-a**3)*\
    (2*(a+b)-(a**2+a*b+b**2)/r-a**2*b**2/r**3)
    return y

for i in range(n):
    T[i]=temp(r[i])
    disp[i]=u(r[i])
    sigma_r[i]=s_r(r[i])
    sigma_t[i]=s_t(r[i])
    epsi_r[i]=e_r(r[i])
    epsi_t[i]=e_t(r[i])
    
#grafica 
plt.figure()
plt.plot(r*1e3,T) 
plt.xlabel('r [mm]')
plt.ylabel('T') 
      
plt.figure()
plt.plot(r*1e3,disp*1e3) 
plt.xlabel('r [mm]')
plt.ylabel('u mm') 
    

plt.figure()
plt.plot(r*1e3,sigma_r) 
plt.xlabel('r [mm]')
plt.ylabel('$\sigma_r$ Pa') 

plt.figure()
plt.plot(r*1e3,sigma_t) 
plt.xlabel('r [mm]')
plt.ylabel('$\sigma_t$ Pa') 

plt.figure()
plt.plot(r*1e3,epsi_r) 
plt.xlabel('r [mm]')
plt.ylabel('$\epsilon_r$') 

plt.figure()
plt.plot(r*1e3,epsi_t) 
plt.xlabel('r [mm]')
plt.ylabel('$\epsilon_t$') 
