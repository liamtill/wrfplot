# -*- coding: utf-8 -*-
"""
Liam Till
University of Oklahoma / Unversity of Reading
conversions.py

Conversions of variables for wrfplot package
Last modified: 23/04/16
"""
import numpy as np
import constants as const

#convert geopotential height to height in m
def gphgt_to_hgt(totalgp):
    """
    Parameters
    Geopotential (m^2/s^2)
    
    Returns
    Geopotential Height (m)
    """
    hgt = (totalgp / const.g)
    return hgt

#convert pa to hPa    
def pa_to_hpa(pres):
    """
    Parameters
    Pressure (Pa)
    
    Returns
    Pressure (hPa)
    """
    pres_hpa = pres * 0.01
    return pres_hpa

#convert kelvin to celcius    
def k_to_c(tempk):
    """
    Parameters
    Temp (K)
    
    Returns
    Temp(C)
    """
    tempc = tempk - 273.15
    return tempc

#convert kelvin to fahreinheit    
def k_to_f(tempk):
    """
    Parameters
    Temp (K)
    
    Returns
    Temp (F)
    """
    tempc = k_to_c(tempk)
    tempf = tempc * (9./5.) + 32.
    return tempf

# convert c to f
def c_to_f(tempc):
    """
    Parameters
    Temp (C)
    
    Returns
    Temp (F)
    """
    tempf = tempc * (9./5.) + 32.
    return tempf
    
# ms^-1 to knots    
def ms_to_kts(spd):
    """
    Parameters
    Speed (ms^-1)
    
    Returns
    Speed (knots)
    """
    spdkts = spd * 1.94384449
    return spdkts
    
def mm_to_in(ppn):
    """
    Parameters
    Precipitation (mm)
    
    Returns
    Precipitation (inches)
    """
    ppnin = ppn * 0.0393701
    return ppnin

