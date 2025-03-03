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
from scipy.optimize import minimize_scalar
from scipy.interpolate import CubicSpline


Ro=2.0e-2
Ri=Ro-5.0e-3
Cout=0.0
D=2.0e-10
C0=0.320034154  #concet iniciall mol/m3
n=100 #nodes
dr=(Ro-Ri)/(n-1)
#S=31536000 # seconds in a year
#S=3600*24*30 #seconds per month
#S=3600*24 #seconds per day
S=60 #seconds per minute
t_end=120.0
#nt=int(t_end/dt)
nt=300
dt=t_end/(nt-1)
ndata=16 #data


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
Hfluxout=np.zeros(nt)


desor_data = np.loadtxt("data.csv", delimiter=",") 
desor_data[:,1]*=np.pi*(Ro**2-Ri**2)
ndata=desor_data.shape[0]
tot_flux_interp=np.zeros(ndata)

def Cbnd_in(t): 
    y=0.0
    return y

def desor_model(D):

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
        time=j*dt
        Cin=Cbnd_in(time)
        for i in range(1,n-1):
            rhs[i]=-Cold[i]/dt
        rhs[0]=Cin
        rhs[n-1]=Cout
        C=solve_banded((1,1),A,rhs)
        for i in range(1,n-1):
            Cflux[i]=-2.0*np.pi*r[i]*D*(C[i]-C[i-1])/dr
        Cflux[0]=-2.0*np.pi*r[0]*D*(C[1]-C[0])/dr*(-1.0) #el -1.0 del vector normal
        Cflux[n-1]=-2.0*np.pi*r[n-1]*D*(C[n-1]-C[n-2])/dr
        H[:,j+1]=C
        Hflux[:,j+1]=Cflux
        Cold[:]=C
        
    tot_flux_in=0
    tot_flux_out=0
    Hfluxout=Hflux[0,:]+Hflux[n-1,:]
    tot_acum=C0*np.pi*(Ro**2-Ri**2)
    
    for i in range(nt):
        tot_flux_in+=0.5*(Hflux[0,i]+Hflux[0,i+1])*S*dt
        tot_flux_out+=0.5*(Hflux[n-1,i]+Hflux[n-1,i+1])*S*dt  
        tot_acum+=-0.5*(Hfluxout[i]+Hfluxout[i+1])*S*dt #el - indica salida o desorpcion
    #    tot_acum=-(tot_flux_in+tot_flux_out) #el - indica salida o desorpcion
        Table_tot_flux[i,0:3]=[tot_flux_in,tot_flux_out,tot_acum]
        
    for i in range(ndata):
        spl=CubicSpline(t[1:],Table_tot_flux[:,1])
        tot_flux_interp[i]=spl(desor_data[i,0])
    
        
    return [t[1:],Table_tot_flux[:,2],tot_flux_interp]


def E2loss(D):
    y=desor_model(D)[2]-desor_data[:,1]
    y=np.dot(y,y)
    return y

solmin=minimize_scalar(E2loss,bounds=(1e-12,1e-6))
D_fit=solmin.x
print('D fitted=',D_fit)
print('Optimization result=',solmin.message)


t_model=desor_model(D_fit)[0]
C_model=desor_model(D_fit)[1]
plt.figure()
plt.plot(t_model,C_model,'--')
plt.plot(desor_data[:,0],desor_data[:,1],'o',alpha=0.5)


    


