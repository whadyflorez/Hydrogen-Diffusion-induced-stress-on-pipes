#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 21 10:01:02 2024

@author: whadyimac
Ejemplo calculo permeabilidad y fugas

"""
import numpy as np

R=8.2057338e-5 #m3.atm.K-1.mol-1 https://byjus.com/physics/value-of-r-in-atm/#:~:text=The%20gas%20constant%20value%20is,1⋅K%5E%E2%88%921.
P=1.1 #atm
T=293.15 #https://www.ferrovial.com/es/recursos/condiciones-normalizadas-de-presion-y-temperatura/#:~:text=Las%20Condiciones%20Normalizadas%20de%20Presión%20y%20Temperatura%20(CNPT)%2C%20comúnmente,y%20consistentes%20en%20diferentes%20aplicaciones.
c=P/(R*T) #mol/m3
r_in=18e-3 
r_out=20e-3
Perm=8.21e-16
j_perm_H2=P*1.01e5/(np.log(r_out/r_in)/(2*np.pi*Perm)) #mol/s
S=3600*24*365  #segundos en un año
J_perm_H2=j_perm_H2*S*22.4 #litros/(año.metro)





