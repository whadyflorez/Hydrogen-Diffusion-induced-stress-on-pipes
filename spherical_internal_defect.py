#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 08:19:53 2024
stresses in internal spherical cavity
@author: whadymacbook2016
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded
import matplotlib.ticker as ticker

Ro=2.0e-2
Ri=Ro-3.0e-3
Cin=0.1
Cout=0.0
D=1.0e-15
E=10e9 #Pa
nu=0.3
Omega=5.0e-3
C0=0.0
pin=400.0e3 #Pa https://www.ugc.edu.co/pages/juridica/documentos/institucionales/NTC_2505_Instalaciones_Suministro_De_Gas.pdf
R_cavity=1.0e-4
r_cavity=0.5*(Ro+Ri) #position of cavity center
Rg=8.31446261815324 
T_ref=25.0+273.15 #reference temperature in K


n=2000 #nodes
dr=(Ro-Ri)/(n-1)
#dt=0.083
S=31536000 # seconds in a year
t_end=100.0
#nt=int(t_end/dt)
nt=10000
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
H_mol_cavity=np.zeros(nt+1)
H_pressure_cavity=np.zeros(nt+1)
H_r_cavity=np.zeros(nt+1)
H_stressr_cavity=np.zeros(nt+1)
H_strainr_cavity=np.zeros(nt+1)


Cold[:]=C0
H[:,0]=Cold
Hflux[:,0]=0

index_cavity=np.argmin(np.abs(r-r_cavity)) #the closest index to the value r_cavity in r

# modelo de difusion y matrices para ambos problemas
for i in range(1,n-1):
    A[1+i-i,i]=-2*S*D/dr**2-1/dt
    A[1+i-(i+1),i+1]=S*D/dr**2+S*D/(2*r[i]*dr)
    A[1+i-(i-1),i-1]=S*D/dr**2-S*D/(2*r[i]*dr)
A[1+0-0,0]=1 
A[1+n-1-(n-1),n-1]=1

mol_cavity=0.0
H_r_cavity[0]=R_cavity
for j in range(nt):
    for i in range(1,n-1):
        rhs[i]=-Cold[i]/dt
    rhs[0]=Cin
    rhs[n-1]=Cout
    C=solve_banded((1,1),A,rhs)
    
    for i in range(1,n-1): # post processing de flujos
        Cflux[i]=-D*(C[i+1]-C[i-1])/(2*dr)
    Cflux[0]=D*(-3*C[0]+4*C[1]-C[2])/(2*dr)  
    Cflux[n-1]=-D*(3*C[n-1]-4*C[n-2]+C[n-3])/(2*dr)
    H[:,j+1]=C
    Hflux[:,j+1]=Cflux
    Cold[:]=C
    mol_cavity+=Cflux[index_cavity]*4.0*np.pi*R_cavity**2*S*dt
    Vol_cavity=(4.0/3.0)*np.pi*R_cavity**3
    pressure_cavity=mol_cavity*Rg*T_ref/Vol_cavity
    sigmar_cavity=pressure_cavity
    deltaR_cavity=R_cavity*pressure_cavity*(1.0-nu)/(2*E)
    R_cavity+=deltaR_cavity
    H_mol_cavity[j+1]=mol_cavity
    H_pressure_cavity[j+1]=pressure_cavity
    H_r_cavity[j+1]=R_cavity
    H_stressr_cavity[j+1]=sigmar_cavity
    
#post processing de flujos totales    
# tot_flux_in=0
# tot_flux_out=0
# for i in range(nt):
#     tot_flux_in+=2*0.5*(Hflux[0,i]+Hflux[0,i+1])*S*dt
#     tot_flux_out+=0.5*(Hflux[n-1,i]+Hflux[n-1,i+1])*S*dt
   
#grafica Concentracion vs r para cada tiempo    
plt.figure()
for i in range(0,nt+1,1):
    plt.plot(r*1e3,H[:,i],label=f't={i*dt:.1f} years') 
plt.xlabel('r [mm]')
plt.ylabel('C') 
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    
#grafica concentracion vs tiempo para cada radio
plt.figure()
for i in range(0,n,10):
    plt.plot(t[:],H[i,:],label=f'r={i*dr*1e3:.1f} mm') 
plt.xlabel('t [years]')
plt.ylabel('C')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))

#grafica flujo vs r para cada tiempo
plt.figure()
for i in range(0,nt+1,1):
    plt.plot(r*1e3,Hflux[:,i],label=f't={i*dt:.1f} years') 
plt.xlabel('r [mm]')
plt.ylabel('Flux $m^3/s$') 
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))

#grafica flujos en radio interno y externo vs tiempo
plt.figure()
plt.plot(t,Hflux[0,:],label=f'inner radius flow') 
plt.plot(t,Hflux[n-1,:],label=f'outer radius flow') 
plt.xlabel('t [years]')
plt.ylabel('Flux $m^3/s$') 
plt.legend(loc='upper right')



#grafica de la presion en la cavidad vs tiempo
plt.figure()
plt.plot(t,H_pressure_cavity) 
plt.xlabel('t ')
plt.ylabel('p (Pa)') 
#grafica del radio en la cavidad vs tiempo
plt.figure()
plt.plot(t,H_r_cavity*1e3) 
plt.xlabel('t ')
plt.ylabel('$r_{cavity}$ (mm)') 
# Formatear el eje Y en notación científica
ax = plt.gca()  # Obtener el eje actual
formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-2, 2))  # Controla cuándo activar la notación científica
ax.yaxis.set_major_formatter(formatter)

plt.show()


