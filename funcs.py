# -*- coding: utf-8 -*-
"""
Liam Till
University of Oklahoma / Unversity of Reading
funcs.py

Functions for interpolating and other stuff for wrfplot package
Last modified: 05/05/16
"""

import numpy as np
import sys
from scipy.ndimage.filters import minimum_filter, maximum_filter

# write skew-t data in SPC file format
def skewt_data(timestamp, level, height, tempc, dwpt, winddir, wspd, inlat, inlon):
    """
    Parameters
    Timestamp, level (hPa), height (m), Temps (C), Dewpoints (C), Wind dir (degrees), Wind Speed (knots), Lat/Lon
    
    Returns
    Writes a file of sounding data in SPC format
    """
    top='%TITLE% \n '+str(inlat)+' '+str(inlon)+'   '+timestamp+'\n\n'
    header='   LEVEL       HGHT       TEMP       DWPT       WDIR       WSPD\n-------------------------------------------------------------------\n%RAW%\n'
    bottom='%END%'
    with open('skewt_data', 'w') as f:
        f.write(top+header)
        for line in range(len(level)):
            #kinda get the UWYO format
            eachline = '{:>8.2f},  {:>8.2f},  {:>8.2f},  {:>8.2f},  {:>8.2f},  {:>8.2f} \n'.format(level[line], height[line], tempc[line], dwpt[line], winddir[line], wspd[line])
            f.write(eachline)
        f.write(bottom)
    f.close()
    
# find the nearest lat/lon to entered lat/lon
def latlon_ij(inlat, inlon, xlat, xlong): 
    """
    Parameters
    inlat, inlon (Lat/Lon of interest)
    xlat, xlong (Array of Lats/Lons)
    
    Returns
    i, j index of nearest Lat/Lon to entered Lat/Lon
    """
    #inlat = opt.slat
    #inlon = opt.slon
    dlat = dlon = 2.0
    dfactor = 2.0
    if (inlat < np.min(xlat[0])) or (inlat > np.max(xlat[0])) or (inlon < np.min(xlong[0])) or (inlon > np.max(xlong[0])):
        print "Lat/Lon entered out of range of domain"
        sys.exit()
    else:
        loop = 1
        loopcount = 0
        i = j = [0,0]
        while loop == 1:
            while (len(i) > 1) or (len(j) > 1):
                latij = np.where((xlat[0] < (inlat+dlat)) & (xlat[0] >= (inlat-dlat)),1,0)
                lonij = np.where((xlong[0] <= (inlon+dlon)) & (xlong[0] > (inlon-dlon)),1,0)
                i,j = np.where((latij == 1) & (lonij == 1))
                dlat = dlat / dfactor
                dlon = dlon / dfactor
                if loopcount > 2000:
                    print "loopcount= ", loopcount
                    print "TOO MANY LOOPS. QUITTING"                   
                    sys.exit()
                loopcount += 1
            if (np.shape(i)[0] == 0) or (np.shape(j)[0] == 0):
                i=[0,0]
                j=[0,0]
            if (np.shape(i)[0] == 1) & (np.shape(j)[0] == 1): 
                loop = 0
        #print "i,j found is: ", i,j
        print "Nearest Lat/Lon is: ", xlat[0,i[0],j[0]], xlong[0,i[0],j[0]]
        
        return i[0], j[0]

# return i, j of input lat and lon SLIGHTLY DIFFERENT FROM ABOVE
def latlon_ij2(inlat, inlon, xlat, xlong): 
    """
    THIS METHOD LESS ACCURATE
    Parameters
    inlat, inlon (Lat/Lon of interest)
    xlat, xlong (Array of Lats/Lons)
    
    Returns
    i, j index of nearest Lat/Lon to entered Lat/Lon
    """
    if (inlat < np.min(xlat[0])) or (inlat > np.max(xlat[0])) or (inlon < np.min(xlong[0])) or (inlon > np.max(xlong[0])):
        print "Lat/Lon entered out of range of domain"
        sys.exit()
    i = np.argmax(xlat[0] > inlat)
    j = np.argmax(xlong[0] > inlon)
    #print len(xlat[0].flatten())
    #print "i,j of flattened array: ", i,j
    newlat = xlat[0].flatten()[i]
    newlon = xlong[0].flatten()[j]
    #print "Lat/Lon from array at this i,j: ", newlat, newlon
    newi = np.unravel_index(i, xlat[0].shape)
    newj = np.unravel_index(j, xlong[0].shape)
    #print "i,j of original array: ", newi, newj
    print "Nearest Lat/Lon is: ", xlat[0,newi[0],newi[1]], xlong[0,newj[0],newj[1]]

    return newi,newj

# weighted linear interpolation of 3d array of data to pressure level
def linear_interp(data, totalp, plev):
    """
    Interpolates data to a specified pressure level
    
    Parameters
    data (3d array)
    total pressure (3d array)
    Pressure level (Pa)
    
    Returns
    Data interpolated to specified pressure level
    """
    """
    # Slow iterative way of doing it
    weight = np.zeros(np.shape(data[0]))
    outVal = np.zeros(np.shape(data[0]))
    #print below.shape, above.shape, weight.shape, value.shape
    for i in range(totalp.shape[0]):
        for j in range(totalp.shape[1]):
            for k in range(totalp.shape[2]):
                currP = totalp[i,j,k]
                prevP = totalp[i-1,j,k]
                if plev*100. > currP and plev*100. < prevP:
                    upper = i
                    lower = upper-1
                    #print upper,lower
                    upperP = totalp[upper,j,k]
                    lowerP = totalp[lower,j,k]
                    upperVal = data[upper,j,k]
                    lowerVal = data[lower,j,k]
                    totaldist = upperP - lowerP
                    weight[j,k] = np.abs( ( (plev*100.) - lowerP) / (totaldist) )
                    currweight = weight[j,k]
                    #print currweight
                    outVal[j,k] = ( lowerVal * (1 - currweight) ) + ( upperVal * currweight )
    #print weight.shape, outVal.shape
    """         
    # more efficient way of doing it
    #find index of level above plev        
    above = np.argmax(totalp < plev*100., axis=0)
    below = above - 1 # index of pressure below
    # pressure at level above plev 
    nz,ny,nx = totalp.shape
    upperP = totalp.reshape(nz,ny*nx)[above.flatten(),range(ny*nx)].reshape(ny,nx)
    #pressure at level below plev
    lowerP = totalp.reshape(nz,ny*nx)[below.flatten(),range(ny*nx)].reshape(ny,nx)
    # value above plev
    nz,ny,nx = data.shape
    aboveVal = data.reshape(nz, ny*nx)[above.flatten(),range(ny*nx)].reshape(ny,nx)
    #value below plev
    belowVal = data.reshape(nz, ny*nx)[below.flatten(),range(ny*nx)].reshape(ny,nx)
    # calc total dist betweek upper and lower
    totaldist = upperP - lowerP
    #calc weighting
    weight = np.abs( ( (plev*100.) - lowerP) / (totaldist) )
    #calc interpolated value    
    outVal = ( belowVal * (1 - weight) ) + ( aboveVal * weight )
    
    return outVal
    
#function to unstagger data from staggered grid given dimension
def unstagger(data,dim): 
    """
    Unstagger data from staggered grid
    
    Parameter
    data (Nd array)
    dim
    
    Returns
    Unstaggered data array of Nd shape
    """
    nd = len(data.shape)
    if dim == 'X' or dim == 'U':
        if nd == 5:
            dataout = ( data[:,:,:,:,:-1] + data[:,:,:,:,1:] ) / 2.
        elif nd == 4:
            dataout = ( data[:,:,:,:-1] + data[:,:,:,1:] ) / 2.
        elif nd == 3:
            dataout = ( data[:,:,:-1] + data[:,:,1:] ) / 2.
        elif nd == 2:
            dataout = ( data[:,:-1] + data[:,1:] ) / 2.
        elif nd == 1:
            dataout = ( data[:-1] + data[1:] ) / 2.
        else: pass
    if dim == 'Y' or dim == 'V':
        if nd == 5:
            dataout = ( data[:,:,:,:-1,:] + data[:,:,:,1:,:] ) / 2.
        elif nd == 4:
            dataout = ( data[:,:,:-1,:] + data[:,:,1:,:] ) / 2.
        elif nd == 3:
            dataout = ( data[:,:-1,:] + data[:,1:,:] ) / 2.
        elif nd == 2:
            dataout = ( data[:-1,:] + data[1:,:] ) / 2.
        elif nd == 1:
            dataout = ( data[:-1] + data[1:] ) / 2.
        else: pass
    if dim == 'Z' or dim == 'W':
        if nd == 5:
            dataout = ( data[:,:,:-1,:,:] + data[:,:,1:,:,:] ) / 2.
        elif nd == 4:
            dataout = ( data[:,:-1,:,:] + data[:,1:,:,:] ) / 2.
        elif nd == 3:
            dataout = ( data[:-1,:,:] + data[1:,:,:] ) / 2.
        elif nd == 1:
            dataout = ( data[:-1] + data[1:] ) / 2.
        else: pass
    return dataout

def extrema(mat,mode='wrap',window=10): # function to find the pressure extrema
    """
    Find the indices of local extrema (min and max) in the input array.
    
    Parameters 
    mat (input array)
    mode
    window (sensitivity)
    
    Returns 
    Indices of extrema
    
    """
    mn = minimum_filter(mat, size=window, mode=mode)
    mx = maximum_filter(mat, size=window, mode=mode)
    # (mat == mx) true if pixel is equal to the local max
    # (mat == mn) true if pixel is equal to the local in
    # Return the indices of the maxima, minima
    return np.nonzero(mat == mn), np.nonzero(mat == mx)
   
#function to interpolate data to given level(s) using np.interp
def interp_to_level(level, coords, data):
    """
    Parameters
    level (int/array of level(s) in INCREASING order)
    data (3d array)
    coords (3d array of levels in INCREASING order)
    
    Returns
    Data interpolated onto given level
    """
    #set up dimensions of interpolated output data
    #height_dim = np.shape(data)[0]
    we_dim = np.shape(data)[1]
    sn_dim = np.shape(data)[2]
    
    interpdata = np.zeros((we_dim,sn_dim)) #array for new data
    coords = np.sort(coords, axis=0) # sort coords into increasing order
        
    for i in range(we_dim): #loop to do interpolation
        for j in range(sn_dim):
            interpdata[i,j] = np.interp(level, coords[:,i,j], data[:,i,j], left=np.nan, right=np.nan)
        
    return interpdata

#generic linear interpolation function using scipy.interpolate.interp1d
def interp_generic(level, coords, data):
    """
    Parameters
    level - int/float of level to interpolate too
    coords - 3d array of co-ordinates (height/pres, lat, lon)
    data - 3d array of data on coords
    
    Returns
    Value at interpolated level in 2d array
    """
    from scipy.interpolate import interp1d
    out = np.zeros((np.shape(data)[1],np.shape(data)[2])) # create array to hold interpolated data
    for i in range(np.shape(data)[1]):
        for j in range(np.shape(data)[2]):
            f = interp1d(coords[:,i,j], data[:,i,j], kind='linear', fill_value=np.nan, bounds_error=False, assume_sorted=False)
            out[i,j] = f(level) # do interpolation
    return out
    
def nine_point_smooth(data, p=0.50, q=-0.25):
    """
    
    """
    smoothed_data = data[1:-1,1:-1] + ( p / 4. ) * ( data[:-1,1:-1] + data[1:-1,1:] + data[1:,1:-1] + data[1:-1,:-1] - ( 4 * data[1:-1,1:-1] ) \
                  + ( q / 4. ) * ( data[:-1,1:-1] + data[:-1,1:-1] + data[1:,1:-1] + data[1:,1:-1] - (4 * data[1:-1,1:-1]) ) )
    return smoothed_data