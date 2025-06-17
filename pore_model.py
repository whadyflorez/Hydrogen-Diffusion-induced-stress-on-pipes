#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 23 10:15:49 2025

@author: whadymacbook2016
"""
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

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
k_trap=1.0e-5
Sc_t=3600 #time scale days
t_end=5 #days
dt=0.01 #days
n_step=int(t_end/dt)

n=np.zeros(n_step+1) #moles inside 
P=np.zeros(n_step+1) #Pore pressure
r=np.zeros(n_step+1) # pore radius
V=np.zeros(n_step+1) # pore volume
c_pore=np.zeros(n_step+1) # concetration in pore

t=np.linspace(0,t_end,n_step+1)

#flux through pipe wall
def flux_pipe(r):
    n=(Ci-Co)/((r/Dif)*np.log(ro/ri))
    return n

#shape factor for the defect
def SF(D,z):
    s=2*pi*D/(1-D/(4*z))*k_trap  #k_trap simula el atrapamiento de H2
    return s

#staate equation for pressure @V,t=constant
def P_se(n,V,T):
    P=n*Ru*T/V #ideal gas
#    P=n*Ru*T/(V-n*b)-a*n**2/V**2 #Van der walls
    return P

#nonlinear implicit equations
def f(x,n_old,dt,n_in,S1,S2,Ci,Co,Dif,Ru,T,ro_pore,nu,E):
    n,p,r,v,c=x
    z=np.zeros(5)
    As=4*pi*r**2
    z[0]=n-n_old-dt*Sc_t*(As*n_in-S1*Dif*(c-Ci)-S2*Dif*(c-Co))
    z[1]=p*v-n*Ru*T
    z[2]=r-ro_pore-(1+nu)/E*ro_pore*p
    z[3]=c*v-n
    z[4]=v-(4/3)*pi*r**3
    return z

#time steps
r[0]=ro_pore
n_in=flux_pipe(r_pos)
V[0]=4.0/3.0*pi*r[0]**3
S1=SF(2*r[0],r_pos-ri)
S2=SF(2*r[0],ro-r_pos)
n_old=0
p_old=0
c_old=0
x0=np.array([n_old,p_old,r[0],V[0],c_old])
for j in range(n_step):
    x=fsolve(f,x0,args=(n_old,dt,n_in,S1,S2,Ci,Co,Dif,Ru,T,ro_pore,nu,E))
    n[j+1]=x[0]
    P[j+1]=x[1]
    r[j+1]=x[2]
    V[j+1]=x[3]
    c_pore[j+1]=x[4]
    S1=SF(2*r[j+1],r_pos-ri)
    S2=SF(2*r[j+1],ro-r_pos) 
    n_old=n[j+1]
    p_old=P[j+1]
    X0=np.array([n_old,p_old,r[j+1],V[j+1],c_pore[j+1]])     

#moles
fig, ax = plt.subplots()
ax.plot(t,n)   
ax.set(xlabel='t (days)')
ax.set(ylabel='n (mol)')   
fig.savefig("moles.png")
            
#Pressure   
fig, ax = plt.subplots()
ax.plot(t,P)   
ax.set(xlabel='t (days)')
ax.set(ylabel='P (Pa)')   
fig.savefig("presion.png")
 
#Radius     
fig, ax = plt.subplots()
ax.plot(t,r)   
ax.set(xlabel='t (days)')
ax.set(ylabel='r (mm)')   
fig.savefig("radio.png")   
 
#concentration  
fig, ax = plt.subplots()
ax.plot(t,c_pore)   
ax.set(xlabel='t (days)')
ax.set(ylabel=r'$c\quad(mol/m3)$')  
fig.savefig("concentracion.png")    
      
   
   
   
   
   
   


