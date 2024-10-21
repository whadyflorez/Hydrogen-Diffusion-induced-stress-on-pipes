#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 21 10:01:02 2024

@author: whadyimac
Ejemplo calculo permeabilidad y fugas

"""
import numpy as np

R= 8.31446261815324 #joules/K mol https://byjus.com/physics/value-of-r-in-atm/#:~:text=The%20gas%20constant%20value%20is,1⋅K%5E%E2%88%921.
MWH2=2.016e-3 #peso molecular h2 kg/mol
MWCH4=16.04e-3 #peso molecular metano kg/mol
MWblend=0.2*MWH2+0.8*MWCH4
unit_cost=12 #costo unitario del hidrogeno USD/kg
P=3*1.01e5 #Pa
T=273.15 #https://www.ferrovial.com/es/recursos/condiciones-normalizadas-de-presion-y-temperatura/#:~:text=Las%20Condiciones%20Normalizadas%20de%20Presión%20y%20Temperatura%20(CNPT)%2C%20comúnmente,y%20consistentes%20en%20diferentes%20aplicaciones.
c=P/(R*T) #mol/m3
S=3600*24*365  #segundos en un año
r_out=20e-3
r_in=r_out-3e-3 

Perm=3.2e-16
Dif=3.5e-11
L=100 #mts

#perdidas por permebilidad
Res_perm=np.log(r_out/r_in)/(2*np.pi*Perm)
j_perm=P/Res_perm*L # mol/s 
j_perm_hour=j_perm*3600 #mol/hour
j_perm_hour_lt=j_perm_hour*22.4 #lt/h
j_perm_year=j_perm*S*MWblend #kg/(año.metro)


#perdidas pro difusion
Res_dif=np.log(r_out/r_in)/(2*np.pi*Dif)
j_dif=c/Res_dif*L # mol/s 
j_dif_hour=j_dif*3600 #mol/h
j_dif_hour_lt=j_dif_hour*22.4
j_dif_year=j_dif*S*MWblend #kg/(año.metro)

#ejemplo calculo de cost de gas a condiciones normales
costo_perdidas_h2=(j_perm_year+j_dif_year)*unit_cost





