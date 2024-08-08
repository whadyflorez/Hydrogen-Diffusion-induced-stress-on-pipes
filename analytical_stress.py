#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 28 17:40:41 2024
analytical solution for stresses in cilinders
@author: whadymacbook2016
"""
import numpy as np
Ro=2.0e-2
Ri=Ro-3.0e-3
delta=Ro-Ri
Cin=0.1
Cout=0.0
D=2.2e-15
E=10e9 #Pa
nu=0.3
Omega=0.1
C0=0.0
pin=400.0e3 #Pa https://www.ugc.edu.co/pages/juridica/documentos/institucionales/NTC_2505_Instalaciones_Suministro_De_Gas.pdf
n=2000 #nodes
dr=(Ro-Ri)/(n-1)

#thick wall cylinder https://en.wikipedia.org/wiki/Cylinder_stress
A=pin*Ri**2/(Ro**2-Ri**2)
B=Ri**2*Ro**2*pin/(Ro**2-Ri**2)
sigma_theta=A+B/Ro**2
sigmamax=pin*(Ri**2+Ro**2)/(delta*(Ri+Ro))
