# -*- coding: utf-8 -*-
"""
Liam Till
University of Oklahoma / Unversity of Reading
calc_vars.py

Calculate variables for wrfplot package
Last modified: 23/04/16
"""

import numpy as np
import constants as const

# mean sea level pressure
def calc_mslp(psfchpa, thgt, t2):
    """
    Parameters
    SFC Pres (hPa), Terrain Height (m), 2m Temps (K)
    
    Returns
    MSLP (hPa)
    """
    #mslp = psfchpa*np.exp((g*M*thgt[0])/(Rg*t2[time])) 
    mslp = psfchpa*(np.exp((const.g*thgt)/(const.R*t2))) # calc MSLP
    return mslp
    
# saturation vapour pressure    
def calc_es(tempc):
    """
    Parameters
    Temps (C)
    
    Returns
    Saturation Vapour Pressure (hPa)
    """
    #using equations from NOAA NWS page at http://www.srh.noaa.gov/epz/?n=wxcalc
    #es = 6.11*10**((7.5*tempc)/(237.3+tempc)) #calculate saturation vapour pressure 
    #es = 6.112*np.exp(17.67*(temp-273.15)/(temp-29.65))  # temp in K      
    # from Rogers and Yau : A Short Course in Cloud Physics 
    es = 6.112*np.exp(17.67*tempc/(tempc+243.5))
    return es
    
def calc_vappres(w, pres):
    """
    Parameters
    Saturation mixing ratio (kg/kg)
    Pressure (Pa)
    
    Returns
    Vapour Pressure (hPa)
    """
    e = (w*pres/(.622+w))/100
    #e = es*(rh/100)
    return e
    
# mixing ratio
def calc_ws(es, preshpa):
    """
    Parameters
    Sat. Vapour Pressure (hPa)
    Pressure (hPa)
    
    Returns
    Mixing Ratio (kg/kg)
    """
    ws = 0.622*(es/(preshpa-es))
    return ws

# saturation mixing ratio
def calc_w(q):  
    """
    Parameters
    Water vapour mixing ratio (kg/kg)
    
    Returns
    Saturation mixing ratio (kg/kg)
    """
    #from Wallace and Hobbs
    w = q/(1-q)
    return w
    
# relative humidity    
def calc_rh(q, ws):
    """
    Parameters
    Water vapour mixing ratio (kg/kg)
    Mixing ratio (kg/kg)
    
    Returns
    Relative Humidity (%)
    """
    # q is the water vapour mixing ratio
    #rh = (e/es)*100
    #rh = 100*(w/ws) #relative humidty
    rh = (q/ws)*100 #calc relative humidity
    rh = np.where(rh > 100,100, rh) #remove values > 100
    rh = np.where(rh < 0,0, rh) #remove values < 0 (just in case)
    return rh

# dewpoint    
def calc_dewpoint(es, rh):
    """
    Parameters
    Saturation vapour pressure (hPa)
    Relative Humidity (%)

    Returns
    Dewpoint (C)
    """
    #td = (243.5*np.log(e/6.112))/(17.67-np.log(e/6.112))
    #td = np.log10(e/6.112) * (243.5/(17.67 - np.log10(e/6.112))) # calc Td another method using e **UNDERESTIMATES**
    td1 = 237.3*np.log10( (es*rh) / 611) #calc numerator
    td2 = 7.5*np.log10(10)-np.log10( (es*rh) / 611) #calc denominator
    td = td1/td2 #calc Td . GOOD METHOD in degC
    return td

# wind direction    
def calc_wdir(u, v):
    """
    Parameters
    U (ms^-1)
    V (ms^-1)
    
    Returns
    Wind direction (Meteorological Degrees)
    """
    r2d = 45.0/np.arctan(1.0)
    wdir = np.arctan2(u, v) * r2d + 180
    return wdir

# wind speed (magnitude)    
def calc_wspeed(u, v):
    """
    Parameters
    U (ms^-1)
    V (ms^-1)
    
    Returns
    Wind Speed (Magnitude) (ms^-1)
    """
    wspd = np.sqrt(u*u + v*v)
    return wspd
    
# potential temp to temp    
def theta_to_temp(theta, totalp):
    """
    Parameters
    Potential temperature (K)
    Total Pressure (Pa)
    
    Returns
    Temperature (K)
    """
    tempsfac = ( (totalp*0.01) / 1000. )**(const.R/1004.) # factor to multiply theta by
    temps = (theta*tempsfac)
    return temps

# density of air    
def calc_rhoa(pres, temps):
    """
    Parameteres
    Pressure (Pa)
    Temps (K)
    
    Returns
    Density (Kg/m^3)
    """
    rhoa = ( pres / (const.R*temps) )
    return rhoa
    
####
    # Reference for simulated and composite reflectivity:   
    # Koch, S., Ferrier, B., Stoelinga, M., Szoke, E., Weiss, S., Kain, J., 2005:
    # THE USE OF SIMULATED RADAR REFLECTIVITY FIELDS IN THE DIAGNOSIS OF
    # MESOSCALE PHENOMENA FROM HIGH-RESOLUTION WRF MODEL FORECASTS
    # J4J.7
####
# calculate dBZ (simulated reflectivity)    
def calc_dbz(t2c, rhoa, Qrain, Qsnow):
    """
    Parameters
    Temps (C)
    Density (kg/m^3)
    Rain mixing ratio (kg/kg)
    Snow mixing ratio (kg/kg)
    
    Returns
    Relfectivity (dBZ)    
    """
    N0snow = 2.0e6*np.exp(-0.12*t2c) # N0 for snow (intercept parameter)
    lambr = np.divide((np.pi*const.N0rain*const.rhol), np.multiply(rhoa,Qrain))**0.25 # rain lambda slope factor
    lambs = np.exp(-0.0536*t2c) # snow lambda slope factor
    Zer = ( (720.0*const.N0rain)*(lambr**-7.0) )*1e18 # rain equivalent reflectivity factor
    Zes_int = np.divide((lambs*Qsnow*rhoa), N0snow)
    Zes = ( (0.224*720.*1e18) / (np.pi*const.rhol)**2 )*(Zes_int**2) # snow equivalent reflectivity factor
    Ze = np.add(Zer,Zes) # total reflectivity factor
    dBZ = 10*np.log10(Ze) # compare apples to apples / convert Ze to dBZ
    dBZ = np.nan_to_num(dBZ) # get rid of NaN's to zeros
    return dBZ
    
# calc lifted condensation level
def calc_lcl(t2c, td):
    """
    Parameters
    Temps (C)
    Dewpoint Temps (C)
    
    Returns
    LCL Height (m)
    """
    lcl = 125.*(t2c-td) # calc lcl
    lcl = np.where(lcl < 0, 0, lcl)    
    lcl = np.nan_to_num(lcl) # get rid of NaNs
    return lcl
    
# equivalent potential temperature
def calc_thetae(theta, temp, ws):
    """
    Parameters
    Potential temp (K)
    Temp (K)
    Mixing ratio (kg/kg)
    
    Returns
    Equivalent Potential Temperature (K)
    """
    # Version from Wallace and Hobbs
    thetae = theta*np.exp( (const.Lv*ws) / (1004.*temp) ) # calc equivalent potential temperature
    #Simplified version from Stull
    #thetae = (t2[time]+(Lv/Cpd)*q2[time])*(1000/psfchpa)**(R/Cpd)
    return thetae

# coriolis parameter    
def calc_fcoriolis(lats):
    """
    Parameters
    Latitudes
    
    Returns
    f (s^-1)
    """
    fcoriolis = 2 * ( 2*np.pi/86400 ) * np.sin( lats*(np.pi/180.) ) # calc coriolis term
    return fcoriolis
    
# vertical vorticity
def calc_vertvort(u, v, dx):
    """
    Parameters 
    u (ms^-1)
    v (ms^-1)
    dx (m)
    
    Returns
    Vertical vorticity (s^-1) 
    """
    dvdx = np.gradient(v,dx,dx)[1] # calc dvdx
    dudy = np.gradient(u,dx,dx)[0] # calc dudy
    vertvort = dvdx - dudy
    return vertvort