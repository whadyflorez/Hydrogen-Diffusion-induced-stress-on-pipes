#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 21 10:01:02 2024

@author: whadyimac
Ejemplo calculo permeabilidad y fugas

"""
import numpy as np

R= 8.31446261815324 #joules/K mol https://byjus.com/physics/value-of-r-in-atm/#:~:text=The%20gas%20constant%20value%20is,1⋅K%5E%E2%88%921.
P=1.1*1.01e5 #Pa
T=293.15 #https://www.ferrovial.com/es/recursos/condiciones-normalizadas-de-presion-y-temperatura/#:~:text=Las%20Condiciones%20Normalizadas%20de%20Presión%20y%20Temperatura%20(CNPT)%2C%20comúnmente,y%20consistentes%20en%20diferentes%20aplicaciones.
c=P/(R*T) #mol/m3
S=3600*24*365  #segundos en un año
r_out=20e-3
r_in=r_out-3e-3 

Perm=1.15e-16
Dif=6.6e-12

#perdidas por permebilidad
Res_perm=np.log(r_out/r_in)/(2*np.pi*Perm)
j_perm=P/Res_perm #ppm mol/s 
j_perm_tot=j_perm*S*22.4 #litros/(año.metro)

#perdidas pro difusion
Res_dif=np.log(r_out/r_in)/(2*np.pi*Dif)
j_dif=c/Res_dif #ppm mol/s 
j_dif_tot=j_dif*S*22.4 #litros/(año.metro)







