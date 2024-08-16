#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 08:19:53 2024
Diffusion cylinder implicit FD constant diffusivity
@author: whadymacbook2016
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded

Ro=2.0e-2
Ri=Ro-3.0e-3
Cin=0.1
Cout=0.0
D=2.2e-15
C0=0.0
n=100 #nodes
dr=(Ro-Ri)/(n-1)
#dt=0.083
S=31536000 # seconds in a year
t_end=100.0
#nt=int(t_end/dt)
nt=300
dt=t_end/(nt-1)


Cold=np.zeros(n)
C=np.zeros(n)
A=np.zeros((3,n))
rhs=np.zeros(n)
r=np.linspace(Ri,Ro,n)
t=np.linspace(0,nt*dt,nt+1)
H=np.zeros((n,nt+1))
Cflux=np.zeros(n)
Hflux=np.zeros((n,nt+1))
Table_tot_flux=np.zeros((nt,3))

Cold[:]=C0
H[:,0]=Cold
Hflux[:,0]=0

for i in range(1,n-1):
    A[1+i-i,i]=-2*S*D/dr**2+S*D/(r[i]*dr)-1/dt
    A[1+i-(i+1),i+1]=S*D/dr**2
    A[1+i-(i-1),i-1]=S*D/dr**2-S*D/(r[i]*dr)
A[1+0-0,0]=1 
A[1+n-1-(n-1),n-1]=1

for j in range(nt):
    for i in range(1,n-1):
        rhs[i]=-Cold[i]/dt
    rhs[0]=Cin
    rhs[n-1]=Cout
    C=solve_banded((1,1),A,rhs)
    for i in range(1,n-1):
        Cflux[i]=-2.0*np.pi*r[i]*D*(C[i]-C[i-1])/dr
    Cflux[0]=-2.0*np.pi*r[0]*D*(C[1]-C[0])/dr  
    Cflux[n-1]=-2.0*np.pi*r[n-1]*D*(C[n-1]-C[n-2])/dr
    H[:,j+1]=C
    Hflux[:,j+1]=Cflux
    Cold[:]=C
    
tot_flux_in=0
tot_flux_out=0

for i in range(nt):
    tot_flux_in+=0.5*(Hflux[0,i]+Hflux[0,i+1])*S*dt
    tot_flux_out+=0.5*(Hflux[n-1,i]+Hflux[n-1,i+1])*S*dt
    tot_acum=tot_flux_in-tot_flux_out
    Table_tot_flux[i,0:3]=[tot_flux_in,tot_flux_out,tot_acum]
    

#concentracions vs radio    
plt.figure()
for i in range(0,nt+1,20):
    plt.plot(r*1e3,H[:,i],label=f't={i*dt:.1f} years') 
plt.xlabel('r [mm]')
plt.ylabel('C') 
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    
#concentraciones vs tiempo
plt.figure()
for i in range(0,n,10):
    plt.plot(t[:],H[i,:],label=f'r={i*dr*1e3:.1f} mm') 
plt.xlabel('t [years]')
plt.ylabel('C')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))

#flujos vs radio
plt.figure()
for i in range(1,nt+1,1):
    plt.plot(r*1e3,Hflux[:,i],label=f't={i*dt:.1f} years') 
plt.xlabel('r [mm]')
plt.ylabel('Flux m^3/s') 
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))

#flujos interno y externo
plt.figure()
plt.plot(t,Hflux[0,:],label=f'inner radius flow') 
plt.plot(t,Hflux[n-1,:],label=f'outer radius flow') 
plt.xlabel('t [years]')
plt.ylabel('Flux m^3/s') 
plt.legend(loc='upper right')

#flujos acumulados
plt.figure()
plt.plot(t[1:],Table_tot_flux[:,0],label=f'inner radius total mass input') 
plt.plot(t[1:],Table_tot_flux[:,1],label=f'outer radius total mass output') 
plt.plot(t[1:],Table_tot_flux[:,2],label=f'Accumulated moles in pipe') 
plt.xlabel('t [years]')
plt.ylabel('total mass mol/m') 
plt.legend(loc='upper right')



