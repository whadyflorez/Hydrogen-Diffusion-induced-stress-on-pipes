#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 14:18:18 2024

@author: whadymacbook2016
Modelo difusiond e hidrogeno cilindros analitico
"""
import numpy as np
from scipy.special import j0,j1,y0,y1
from scipy.optimize import root_scalar
from scipy.integrate import quad,dblquad
import matplotlib.pyplot as plt

Ro=2.0e-2
Ri=Ro-3.0e-3
g0=0.1
h0=0.0
D=2.2e-15
V0=0.0

# Intervalo de busqueda raices
a = 800.0
b = 10000.0
num_points=2000

limit=200 #integration intervals

t_end=3600*24*365*100
secs_year=3600*24*365

def find_sign_changes(f, a, b, num_points):
    """
    Encuentra los subintervalos donde la función f(x) cambia de signo en el intervalo [a, b].

    Parámetros:
    f          : función a evaluar
    a, b       : límites del intervalo [a, b]
    num_points : número de puntos para dividir el intervalo

    Retorna:
    Una lista de subintervalos donde la función cambia de signo.
    """
    # Generar puntos equidistantes en el intervalo [a, b]
    x = np.linspace(a, b, num_points)
    y = f(x)

    # Lista para almacenar los subintervalos donde la función cambia de signo
    sign_change_intervals = []
    roots=[]

    # Revisar cambios de signo entre puntos consecutivos
    for i in range(len(x) - 1):
        if y[i] * y[i + 1] < 0:
            sign_change_intervals.append((x[i], x[i + 1]))
            root=root_scalar(f,x0=x[i],x1=x[i + 1],method='newton')
            roots.append(root.root)
    return sign_change_intervals,roots

# Ejemplo de uso:
def funcion_m(m):
    r=j0(m*Ro)*y0(m*Ri)-j0(m*Ri)*y0(m*Ro)
    return r


# Encontrar los subintervalos donde cambia de signo la función ejemplo_funcion
intervalos = find_sign_changes(funcion_m, a, b, num_points)
mn=intervalos[1]
Bn=[]
for i in mn:
    Bn.append(-j0(i*Ri)/y0(i*Ri))    

#args[0]=mn  args[1]=Bn
n=len(mn)   
Ln=np.zeros(n)
An=np.zeros(n)
Sn=np.zeros(n)

for i in range(n):
   def Qn_fun(r,*args):
      mn=args[0]
      Bn=args[1]
      y=j0(mn*r)+Bn*y0(mn*r)
      return y
   def flux_Qn_fun(r,*args):
      mn=args[0]
      Bn=args[1]
      y=-mn*j1(mn*r)-Bn*mn*y1(mn*r)
      y=-D*y
      return y
   def Ln_fun(r,*args):
      y=r*Qn_fun(r,*args)**2
      return y
   result_int=quad(Ln_fun,Ri,Ro,limit=limit,args=(mn[i],Bn[i]))
   Ln[i]=result_int[0]
   def V_fun(r):
       y=V0-np.log(r/Ro)/np.log(Ri/Ro)*g0-np.log(Ri/r)/np.log(Ri/Ro)*h0
       return y
   def An_fun(r,*args):
       y=r*Qn_fun(r,*args)*V_fun(r)
       return y
   result_int=quad(An_fun,Ri,Ro,limit=limit,args=(mn[i],Bn[i]))
   An[i]=(1.0/Ln[i])*result_int[0]
   

def C_series(r,t):
    for i in range(n):
       Sn[i]=An[i]*Qn_fun(r,mn[i],Bn[i])*np.exp(-mn[i]**2*D*t)    
    S=np.sum(Sn)+np.log(r/Ro)/np.log(Ri/Ro)*g0+np.log(Ri/r)/np.log(Ri/Ro)*h0    
    return S   

def flux_C_series(r,t):
    for i in range(n):
       Sn[i]=An[i]*flux_Qn_fun(r,mn[i],Bn[i])*np.exp(-mn[i]**2*D*t)    
    S=np.sum(Sn)-D*(1.0/r)/np.log(Ri/Ro)*g0+D*(1.0/r)/np.log(Ri/Ro)*h0    
    S=2.0*np.pi*r*S
    return S  

def flux_Ci_series(t):
  return flux_C_series(Ri,t)   

def flux_Co_series(t):
  return flux_C_series(Ro,t)  
    

nr,nt=1000,300
r = np.linspace(Ri, Ro, nr)
t = np.linspace(0, t_end, nt)
dt=t_end/(nt-1)
dr=(Ro-Ri)/(nr-1)
C=np.zeros((nr,nt))
Ci=np.zeros(nt)
Co=np.zeros(nt)
flux_C=np.zeros((nr,nt))
for i in range(nt):
    for j in range(nr):
        C[j,i]=C_series(r[j],t[i])
        flux_C[j,i]=flux_C_series(r[j],t[i])
Ci=flux_C[0,:] 
Co=flux_C[nr-1,:]       

#plt.rcParams['text.usetex'] = True        
plt.figure()
for i in range(0,nt,20):
    plt.plot(r*1e3,C[:,i],label=f't={i*dt/secs_year:.1f} years') 
plt.xlabel('r [mm]')
plt.ylabel('C') 
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))


plt.figure()
for i in range(0,nr,100):
    plt.plot(t[1:]/secs_year,C[i,1:],label=f't={i*dr*1e3:.1f} mm') 
plt.xlabel('t [years]')
plt.ylabel('C')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))

plt.figure()
for i in range(1,nt,1):
    plt.plot(r*1e3,flux_C[:,i],label=f't={i*dt/secs_year:.1f} years') 
plt.xlabel('r [mm]')
plt.ylabel('Flux m^3/s') 
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))

plt.figure()
plt.plot(t/secs_year,Ci,label=f'inner radius flow') 
plt.plot(t/secs_year,Co,label=f'outer radius flow') 
plt.xlabel('t [years]')
plt.ylabel('Flux m^3/s') 
plt.legend(loc='upper right')

total_C_in=quad(flux_Ci_series,0,t_end)[0]
total_C_out=quad(flux_Co_series,0,t_end)[0]
total_C_variation=total_C_in-total_C_out
print('H2 acumulado',total_C_variation)


# plt.plot(r,C[:,1:nt:10])  
# plt.xlabel('r [m]')
# plt.ylabel('C')  

# plt.figure()
# plt.plot(t[1:]/secs_year,np.transpose(C[0:nr:50,1:]))   
# plt.xlabel('t [years]')
# plt.ylabel('C')


