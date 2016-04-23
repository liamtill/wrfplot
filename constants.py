# -*- coding: utf-8 -*-
"""
Liam Till
University of Oklahoma / Unversity of Reading
constants.py

Constants for wrfplot package
Last modified: 22/04/16
"""

#define constants
g = 9.81 # gravity
R = 287. # gas constant for dry air
Rg = 8.31447 # ideal gas constant
M = 0.0289644 # molar mass of air
N0rain = 8.0e6 # number conc. for rain
#N0snow = 2e7 # number conc. for snow
N0gra = 4.0e6 # number conc. for graupel
#gamma7 = 720. # gamma 7 function
rhol = 1000. # density of liquid
rhos = 100. # density of snow
rhog = 400. # density of graupel
rhoi = 917. # density of ice
alpha = 0.224 # alpha (rhol/rhoi)^2*(mod(Ki)^2/mod(Kl)^2) for composite refl.
Cpd = 1005.7 # J/kg*K specific heat capacity of dry air at constant pressure
Lv = 2264760. # J/kg # latent heat of vaporisation