import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from assimulo.problem import Implicit_Problem
from assimulo.solvers import IDA

Ci=400.0
Co=0.0
Dif=1.7e-11
L=1.0
ri=10.2e-3
ro=12.3e-3
Ru=8.314  #m3⋅Pa⋅K−1⋅mol−1
T=300.0 
r_pos=11.0e-3
ro_pore=1.0e-3
a=0.2476*1e-1 #H2 m**6*Pa/mol**2
b=0.02661*1e-3 #H2 m**3/mol
E=0.565e9 #elastic modulus HDPE Pa
nu=0.4  #poisson ration HDPE
k_trap=1.0e-4 # papers bathia farid porosity epsilo**2 epsilon=0.01 hdpe
Deff=Dif*k_trap #modelos de difusividad efectiva
Sc_t=3600.0*24 #time scale 
t_end=30*12*10
n_step=1000

#flux through pipe wall
def flux_pipe(r):
    n=(Ci-Co)/((r/Deff)*np.log(ro/ri))
    return n

#shape factor for the defect
def SF(D,z):
    s=2*pi*D/(1-D/(4*z))  #k_trap simula el atrapamiento de H2
    return s


def residual(t, y, ydot):
    n,p,r=y
    res=np.zeros(3)
    n_in=flux_pipe(r_pos)
    S1=SF(2*r,r_pos-ri)
    S2=SF(2*r,ro-r_pos)
    V=(4.0/3.0)*pi*r**3
    As=4*pi*r**2
    c=n/V
    res[0]=ydot[0]-(As*n_in-S1*Deff*(c-Ci)-S2*Deff*(c-Co))*Sc_t
    res[1]=p*V-n*Ru*T           
    res[2]=r-ro_pore*(1.0+(1.0+nu)*p/E)  
    return res

n_ini=0.0
p_ini=0.0

V0 = (4.0/3.0)*pi*ro_pore**3
c0 = n_ini / V0
n_in = flux_pipe(r_pos)
S1 = SF(2*ro_pore, r_pos - ri)
S2 = SF(2*ro_pore, ro - r_pos)
As = 4*pi*ro_pore**2

y0=np.array([n_ini,p_ini,ro_pore])

dn_dt_0 = (As*n_in - S1*Deff*(c0 - Ci) - S2*Deff*(c0 - Co)) * Sc_t
yd0 = np.array([dn_dt_0, 0.0, 0.0])

problem = Implicit_Problem(residual, y0, yd0, t0=0.0)
solver = IDA(problem)
solver.algvar = [1, 0, 0]
solver.make_consistent('IDA_YA_YDP_INIT')
solver.inith = 1e-10  # crítico para DAE rígidas con y0=0

solver.atol = [1e-6]*3
solver.rtol = 1e-6

t, y, yd = solver.simulate(t_end)

n=y[:,0]
p=y[:,1]
r=y[:,2]
V=(4.0/3.0)*pi*r**3
c=n/V


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
      
      
    
    
    
    
    
    
    
    
    
    
    
    
    
    
            