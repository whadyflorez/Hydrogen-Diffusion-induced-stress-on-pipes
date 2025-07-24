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
from scipy.stats import gaussian_kde
import pandas as pd

Ro=2.0e-2
Ri=Ro-3.0e-3
Cin=80.0
Cout=0.0
D=1.4e-10
E=10e9 #Pa
pin=3e5 #Pa https://www.ugc.edu.co/pages/juridica/documentos/institucionales/NTC_2505_Instalaciones_Suministro_De_Gas.pdf

def funcion_sigmat_max(D, Cin, E, pin, Ri, Ro):
    Cout=0.0
    nu=0.3
    Omega=44.5e-6 #2.947e-6
    C0=0.0
    n=250 #nodes
    dr=(Ro-Ri)/(n-1)
    #dt=0.083
    #S=31536000 # seconds in a year
    S=60*60 # seconds per hour
    t_end=5.0
    #nt=int(t_end/dt)
    nt=1000
    dt=t_end/(nt-1)
    
    
    Cold=np.zeros(n)
    C=np.zeros(n)
    Disp=np.zeros(n)
    sigma_r=np.zeros(n)
    sigma_t=np.zeros(n)
    epsi_r=np.zeros(n)
    epsi_t=np.zeros(n)
    A=np.zeros((3,n))
    ADisp=np.zeros((5,n))
    rhs=np.zeros(n)
    rhsDisp=np.zeros(n)
    r=np.linspace(Ri,Ro,n)
    t=np.linspace(0,nt*dt,nt+1)
    H=np.zeros((n,nt+1))
    Cflux=np.zeros(n)
    Hflux=np.zeros((n,nt+1))
    HDisp=np.zeros((n,nt+1))
    HStress_r=np.zeros((n,nt+1))
    HStress_t=np.zeros((n,nt+1))
    HStrain_r=np.zeros((n,nt+1))
    HStrain_t=np.zeros((n,nt+1))
    
    Cold[:]=C0
    H[:,0]=Cold
    Hflux[:,0]=0
    
    # modelo de difusion y matrices para ambos problemas
    for i in range(1,n-1):
        A[1+i-i,i]=-2*S*D/dr**2-1/dt
        A[1+i-(i+1),i+1]=S*D/dr**2+S*D/(2*r[i]*dr)
        A[1+i-(i-1),i-1]=S*D/dr**2-S*D/(2*r[i]*dr)
        ADisp[2+i-(i+1),i+1]=1/dr**2+1/(2*r[i]*dr)
        ADisp[2+i-(i-1),i-1]=1/dr**2-1/(2*r[i]*dr)
        ADisp[2+i-i,i]=-2/dr**2-1./r[i]**2
    A[1+0-0,0]=1 
    A[1+n-1-(n-1),n-1]=1
    ADisp[2+0-0,0]=-3*(nu-1.)/(2*dr)-nu/r[0]
    ADisp[2+0-1,1]=4*(nu-1.)/(2*dr)
    ADisp[2+0-2,2]=-(nu-1.)/(2*dr)
    ADisp[2+n-1-(n-3),n-3]=(nu-1.)/(2*dr)
    ADisp[2+n-1-(n-2),n-2]=-4*(nu-1.)/(2*dr)
    ADisp[2+n-1-(n-1),n-1]=3*(nu-1.)/(2*dr)-nu/r[n-1]
    
    
    for j in range(nt):
        for i in range(1,n-1):
            rhs[i]=-Cold[i]/dt
        rhs[0]=Cin
        rhs[n-1]=Cout
        C=solve_banded((1,1),A,rhs)
        for i in range(1,n-1): # post processing de flujos
            Cflux[i]=-2.0*np.pi*r[i]*D*(C[i+1]-C[i-1])/(2*dr)
        Cflux[0]=-2.0*np.pi*r[0]*D*(-3*C[0]+4*C[1]-C[2])/(2*dr)  
        Cflux[n-1]=-2.0*np.pi*r[n-1]*D*(3*C[n-1]-4*C[n-2]+C[n-3])/(2*dr)
        H[:,j+1]=C
        Hflux[:,j+1]=Cflux
        Cold[:]=C
    #post processing de flujos totales    
    tot_flux_in=0
    tot_flux_out=0
    for i in range(nt):
        tot_flux_in+=0.5*(Hflux[0,i]+Hflux[0,i+1])*S*dt
        tot_flux_out+=0.5*(Hflux[n-1,i]+Hflux[n-1,i+1])*S*dt
    
    # modelo de esfuerzos    
    for j in range(nt+1):
        C[:]=H[:,j]
        for i in range(1,n-1):
            rhsDisp[i]=-(1.0/3.0)*(Omega/(nu-1.0))*(C[i+1]-C[i-1])/(2.*dr)
        rhsDisp[0]=-pin*(1.+nu)*(2*nu-1.)/E-(1./3.)*Omega*C[0]
        rhsDisp[n-1]=-(1./3.)*Omega*C[n-1]
        Disp=solve_banded((2,2),ADisp,rhsDisp)
        HDisp[:,j]=Disp
    #post processing de stress and strain    
        epsi_t=Disp/r
        for i in range(n):
            if i>0 and i<n-1:
                epsi_r[i]=(Disp[i+1]-Disp[i-1])/(2*dr)
            if i==0:
                epsi_r[i]=(-3*Disp[i]+4*Disp[i+1]-Disp[i+2] )/(2*dr)
            if i==n-1:    
                epsi_r[i]=(3*Disp[i]-4*Disp[i-1]+Disp[i-2] )/(2*dr)
            sigma_r[i]=(E/((1.+nu)*(2.*nu-1.)))*\
            ((nu-1.)*epsi_r[i]-nu*epsi_t[i]+(1./3.)*Omega*C[i])
            sigma_t[i]=(E/((1.+nu)*(2.*nu-1.)))*\
            ((nu-1.)*epsi_t[i]-nu*epsi_r[i]+(1./3.)*Omega*C[i])
        HStress_r[:,j]=sigma_r
        HStress_t[:,j]=sigma_t
        HStrain_r[:,j]=epsi_r
        HStrain_t[:,j]=epsi_t
    
    #calculo del esfuerzo maximo tangencial durante toda la historia del proceso
    indice = np.argmax(np.abs(HStress_t))
    posicion = np.unravel_index(indice, HStress_t.shape)
    sigma_t_max=np.abs(HStress_t[posicion])
        
    return sigma_t_max

sigmat_max=funcion_sigmat_max(D, Cin, E, pin, Ri, Ro)

# Número de simulaciones
N = 1000
# Distribuciones aleatorias
D_med = 1.4e-10
g_sd_D = 2.0
ln_sigma_D = np.log(g_sd_D)
D_samples = np.random.lognormal(mean=np.log(D_med), sigma=ln_sigma_D, size=N)

E_med = 10e9
E_std = 1e9
E_samples = np.random.normal(loc=E_med, scale=E_std, size=N)
E_samples = E_samples[E_samples > 0]
while len(E_samples) < N:
    extra = np.random.normal(loc=E_med, scale=E_std, size=N - len(E_samples))
    E_samples = np.concatenate([E_samples, extra[extra > 0]])
E_samples = E_samples[:N]

# Ejecutar Monte Carlo
sigma_max_list = []
for D_i, E_i in zip(D_samples, E_samples):
    sigma_max = funcion_sigmat_max(D_i, Cin, E_i, pin, Ri, Ro)
    sigma_max_list.append(sigma_max)

# Mostrar resultados como tabla
resultados = pd.DataFrame({
    "D (m²/s)": D_samples,
    "E (Pa)": E_samples,
    "σ_t max (Pa)": sigma_max_list
})
print("\nResultados Monte Carlo:\n")
print(resultados)

# Histograma como PDF (área total = 1)
counts, bins, _ = plt.hist(sigma_max_list, bins=10, alpha=0.6, density=True, label="Histograma (PDF)")

# KDE ya es una PDF
kde = gaussian_kde(sigma_max_list)
x_vals = np.linspace(min(sigma_max_list), max(sigma_max_list), 200)
plt.plot(x_vals, kde(x_vals), label="KDE", color="red")

plt.xlabel(r"$\sigma_{\theta}^{\mathrm{max}}$ [Pa]")
plt.ylabel("Densidad")
plt.title("Distribución Monte Carlo normalizada (PDF)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
