#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 23 10:15:49 2025

@author: whadymacbook2016
"""
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt

Ci=400.0
Co=0.0
Dif=1.7e-11
L=1.0
ri=10.2e-3
ro=12.3e-3
Ru=8.314  #m3⋅Pa⋅K−1⋅mol−1
T=300.0 
r_pos=11.0e-3
ro_pore=1.0e-6
a=0.2476*1e-1 #H2 m**6*Pa/mol**2
b=0.02661*1e-3 #H2 m**3/mol
E=0.565e9 #elastic modulus HDPE Pa
nu=0.4  #poisson ration HDPE
Sc_t=24*3600 #time scale days
t_end=360*10 #days
dt=0.025 #days
n_step=int(t_end/dt)

n=np.zeros(n_step+1) #moles inside 
P=np.zeros(n_step+1) #Pore pressure
r=np.zeros(n_step+1) # pore radius
V=np.zeros(n_step+1) # pore volume
c_pore=np.zeros(n_step+1) # concetration in pore

t=np.linspace(0,t_end,n_step+1)

#flux through pipe wall
def flux_pipe(r):
    n=(Ci-Co)/(np.log(r/ri)/(2*pi*Dif*L))
    return n

#shape factor for the defect
def SF(D,z):
    s=2*pi*D/(1-D/(4*z))
    return s

#staate equation for pressure @V,t=constant
def P_se(n,V,T):
#    P=n*Ru*T/V #ideal gas
    P=n*Ru*T/(V-n*b)-a*n**2/V**2 #Van der walls
    return P

#time steps
r[0]=ro_pore
n_in=flux_pipe(r_pos)
V[0]=4.0/3.0*pi*r[0]**2
S1=SF(2*r[0],r_pos-ri)*0
S2=SF(2*r[0],ro-r_pos)
for j in range(n_step):
    denom=1.0+dt*Sc_t*(0.5*S1*Dif*1.0/V[j]+0.5*S2*Dif*1.0/V[j])
    nume=n[j]+dt*Sc_t*(4*pi*r[j]**2*n_in+0.5*S1*Dif*Ci+0.5*S2*Dif*Co\
    -0.5*S1*Dif*(n[j]/V[j]-Ci)-0.5*S2*Dif*(n[j]/V[j]-Co))                   
    n[j+1]=nume/denom
    dr=(1+nu)/E*r[0]*P[j]
    r[j+1]=r[0]+dr
    V[j+1]=4.0/3.0*pi*r[j+1]**2
    P[j+1]=P_se(n[j+1],V[j+1],T) 
    c_pore[j+1]=n[j+1]/V[j+1]
    S1=SF(2*r[j+1],r_pos-ri)*0
    S2=SF(2*r[j+1],ro-r_pos)       

#moles
fig, ax = plt.subplots()
ax.plot(t,n)   
ax.set(xlabel='t (days)')
ax.set(ylabel='n (mol)')   

#Pressure   
fig, ax = plt.subplots()
ax.plot(t,P)   
ax.set(xlabel='t (days)')
ax.set(ylabel='P (Pa)')   
 
#Radius     
fig, ax = plt.subplots()
ax.plot(t,r)   
ax.set(xlabel='t (days)')
ax.set(ylabel='r (mm)')      
 
#concentration  
fig, ax = plt.subplots()
ax.plot(t,c_pore)   
ax.set(xlabel='t (days)')
ax.set(ylabel=r'$c\quad(mol/m3)$')      
      
   
   
   
   
   
   


