import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

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
k_trap=1.0e-4 # papers bathia farid porosity epsilo**2 epsilon=0.01 hdpe
#Sc_t=3600*24 #time scale 
#t_end=30*12*10
Sc_t=60.0
t_end=10
n_step=1000

#flux through pipe wall
def flux_pipe(r):
    n=(Ci-Co)/((r/Dif)*np.log(ro/ri))*k_trap
    return n

#shape factor for the defect
def SF(D,z):
    s=2*pi*D/(1-D/(4*z))*k_trap  #k_trap simula el atrapamiento de H2
    return s


def df(t,y):
    n,p,r=y
    df=np.zeros(3)
    n_in=flux_pipe(r_pos)
    S1=SF(2*r,r_pos-ri)
    S2=SF(2*r,ro-r_pos)
    V=(4.0/3.0)*pi*r**3
    As=4*pi*r**2
    c=n/V
    df[0]=(As*n_in-S1*Dif*(c-Ci)-S2*Dif*(c-Co))*Sc_t
    df[1]=Ru*T*df[0]/\
        (4*pi*r**2*p*ro_pore*((1+nu)/E)+(4/3)*pi*r**3)           
    df[2]=ro_pore*((1+nu)/E)*df[1]   
    return df

n_ini=0.0
p_ini=0.0
y0=np.array([n_ini,p_ini,ro_pore])
t_eval=np.linspace(0,t_end,n_step+1)
t_span=[0.0,t_end]
    
sol=solve_ivp(df, t_span, y0, method='Radau', t_eval=t_eval) 
n=sol.y[0,:]
p=sol.y[1,:]
r=sol.y[2,:]
V=(4/3)*pi*r**3
c=n/V
t=sol.t

r_cr=min(r_pos-ri,ro-r_pos)

   
#figures
#moles
fig, ax = plt.subplots()
ax.plot(t,n)   
ax.set(xlabel='t (days)')
ax.set(ylabel='n (mol)') 

#Pressure   
fig, ax = plt.subplots()
ax.plot(t,p)   
ax.set(xlabel='t (days)')
ax.set(ylabel='P (Pa)')    

#Radius     
fig, ax = plt.subplots()
ax.plot(t,r)   
ax.set(xlabel='t (days)')
ax.set(ylabel='r (mm)')      
 
#concentration  
fig, ax = plt.subplots()
ax.plot(t,c)   
ax.set(xlabel='t (days)')
ax.set(ylabel=r'$c\quad(mol/m3)$')      
      
      
    
    
    
    
    
    
    
    
    
    
    
    
    
    
            