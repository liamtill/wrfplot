#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Liam Till
University of Oklahoma 2016
METR4330 and personal project
wrfplot.py
Python script to plot various WRF model output. Plots are saved as PNG.
example usage: wrfplot.py --infile filename.nc --sfc --tunit C --ppn -punit mm --td
Will plot surface chart and dewpoint in Celcius and precipitation in mm.
Use wrfplot.py --help to list all options
Last modified: 03/04/16
"""

import matplotlib
#matplotlib.use('Agg') # UNCOMMENT THIS ONLY WHEN INVOKING FROM CRON SCRIPT 
from scipy.io import netcdf # USE SCIPY MODULE
#from netCDF4 import Dataset # UNCOMMENT TO USE NETCDF 4 MODULE
from mpl_toolkits.basemap import Basemap, cm
from matplotlib import cm as cm2
import matplotlib.pyplot as plt
import numpy as np
import datetime
from optparse import OptionParser
import os.path
import sys
from scipy.ndimage.filters import gaussian_filter

# option parser
usage="usage: %prog [options] \n example usage: wrfplot.py --infile filename.nc --sfc --tunit C --td --ppn --punit mm"
parser = OptionParser(usage=usage, version="%prog 5.2 by Liam Till")
parser.add_option("--sfc", dest="sfc",action="store_true",help="Plot surface chart with 2m temp, wind barbs and MSLP")
parser.add_option("--t2", dest="t2", action="store_true", help="Plot 2m temp and wind barbs only")
parser.add_option("--mslp", dest="mslp", action="store_true", help="Plot MSLP only")
parser.add_option("--ppnaccum", dest="ppnaccum", action="store_true", help="Plot total accumulated precipitation")
parser.add_option("--ppn", dest="ppn", action="store_true", help="Plot total precipitation")
parser.add_option("--convppn", dest="convppn", action="store_true", help="Plot convective precipitation")
parser.add_option("--td", dest="td", action="store_true", help="Plot 2m dew point temperature")
parser.add_option("--rh", dest="rh", action="store_true", help="Plot relative humidity")
parser.add_option("--snow", dest="snow", action="store_true", help="Plot snow accumulation")
parser.add_option("--hail", dest="hail", action="store_true", help="Plot hail accumulaton")
parser.add_option("--simdbz", dest="simdbz", action="store_true", help="Plot simulated reflectivity")
parser.add_option("--compdbz", dest="compdbz", action="store_true", help="Plot composite reflectivity")
parser.add_option("--lcl", dest="lcl", action="store_true", help="Plot LCL (lifted condensation level)")
parser.add_option("--thetae", dest="thetae", action="store_true", help="Plot Theta-e (equivalent potential temperature)")
parser.add_option("--ua", dest="ua", action="store_true", help="Plot chart at given pressure levels (hPa)")
parser.add_option("--lvl", dest="lvl", help="Specify levels for upper level charts option --ua. Comma seperated e.g 250,500", default="500")
parser.add_option("--run", dest="run", type="string", help="Model initialisation time", default="00")
parser.add_option("--indir", dest="indir", type="string", help="Directory of the NetCDF file", default="")
parser.add_option("--outdir", dest="outdir", type="string", help="Directory to save plots too", default="")
parser.add_option("--infile", dest="infile", type="string", help="NetCDF filename", default="")
parser.add_option("--thin", dest="thin", type="int", help="Thinning factor for wind barbs", default=5)
parser.add_option("--tunit", dest="tunit", type="string", help="Unit of temperature (C or F)", default="C")
parser.add_option("--punit", dest="punit", type="string", help="Unit of precipitation (mm or inches)", default="mm")
parser.add_option("--save", dest="save", action="store_true", help="Save plots as png files")
parser.add_option("--v", dest="verbose", action="store_true", help="Enable verbose")
parser.add_option("--auto", dest="auto", action="store_true", help="Enable auto file input for daily WRF runs")
parser.add_option("--barbsize", dest="barbsize", type="int", help="Set the length of the wind barbs", default=7)
parser.add_option("--75lr", dest="lr75", action="store_true", help="Plot the H7-H5 lapse rates")
parser.add_option("--vort500", dest="vort500", action="store_true", help="Plot the 500mb absolute vorticity")
(opt, arg) = parser.parse_args()

indir = opt.indir # dir of input file
filein = opt.infile
if opt.auto: # for auto file input for daily runs
    run = opt.run # model init time
    filein = indir+'wrfout_d01_'+datetime.datetime.utcnow().strftime('%Y-%m-%d')+'_'+run+':00:00' # auto filename for current days run
while os.path.isfile(filein) is False and not opt.auto: #if file doesnt exist get filename
    print "File", filein, "not found! in directory:", indir
    indir = raw_input("Please enter a directory (blank for current dir): ")
    filein = raw_input("Please enter a filename: ")
try: #check if file exists and read in
    print "Reading in file: ", filein
    #nc = Dataset(filein) # for netcdf 4
    nc = netcdf.netcdf_file(filein,'r') # for scipy
except: # quit if cant read file
    print "Something went wrong reading in the file"
    print "QUITTING"
    sys.exit()

outdir = opt.outdir # output image dir

## BASEMAP STUFF
# dimensions of domain in staggered unit / grid points
x_dim = nc.dimensions['west_east']
y_dim = nc.dimensions['south_north']

#thin factor for wind barbs
thin = opt.thin

#get lats and lons for map projection
cen_lat  = float(nc.CEN_LAT)
cen_lon  = float(nc.CEN_LON)
truelat1 = float(nc.TRUELAT1)
truelat2 = float(nc.TRUELAT2)
standlon = float(nc.STAND_LON)
xlat     = nc.variables['XLAT']
xlong    = nc.variables['XLONG']
map_proj = int(nc.MAP_PROJ)

# Get dx and dy. Grid size
dx = float(nc.DX)
dy = float(nc.DY)

#calculate plot width and height from grid size * dimension. Domain size
width_meters  = dx * (x_dim - 1)
height_meters = dy * (y_dim - 1) 

# Define gridlines
parallels = np.arange(-90,90,10)
meridians = np.arange(0,360,10)

# find projection and create map. Only LCC tested.
if map_proj == 1: #lambert conformal. 
    proj = 'lcc'
    m = Basemap(resolution='i',projection=proj,width=width_meters,height=height_meters,lat_0=cen_lat,lon_0=cen_lon,lat_1=truelat1,lat_2=truelat2)
    # get lat/lons of ny by nx evenly space grid              
    lons, lats = m.makegrid(x_dim, y_dim)
    x, y = m(lons, lats) # compute map proj coordinates.
    proj = 'Lambert Conformal'
elif map_proj == 2: # polar stereographic
    proj = 'npstere'
    m = Basemap(resolution='i',projection=proj,llcrnrlon=xlong[0,0,0],llcrnrlat=xlat[0,0,0],urcrnrlon=xlong[0,-1,-1],urcrnrlat=xlat[0,-1,-1],lat_0=cen_lat,lon_0=cen_lon)    
    x, y = m(xlong[0,:,:],xlat[0,:,:])
    proj = 'Polar Stereographic'
elif map_proj == 3: # mercator
    proj = 'merc'
    m = Basemap(resolution='i',projection=proj,llcrnrlon=xlong[0,0,0],llcrnrlat=xlat[0,0,0],urcrnrlon=xlong[0,-1,-1],urcrnrlat=xlat[0,-1,-1],lat_0=cen_lat,lon_0=cen_lon)    
    x, y = m(xlong[0,:,:],xlat[0,:,:])
    proj = 'Mercator'
else: # not supported and quit
    print "Projection ", map_proj, "unknown"
    print "QUITTING"
    sys.exit()
print "Using map projection: ", proj

## GET DATA
times = nc.variables['Times'] #each time output in wrf nc file
t2 = nc.variables['T2'] #temp at 2m / Kelvin
u10 = nc.variables['U10'] #u10 wind / ms/s
v10 = nc.variables['V10'] #v10 wind / ms/s
psfc = nc.variables['PSFC'] #surface pressure / Pascals
rainc = nc.variables['RAINC'] # accumulated total cumulus precip
rainnc = nc.variables['RAINNC'] # accumulated total grid scale precip
#qvapor = nc.variables['QVAPOR'] #water vapor mixing ratio kg/kg
thgt = nc.variables['HGT'] #terrain height

#define constants
g = 9.81 #gravity
R = 287. #gas constant for dry air
Rg = 8.31447 #ideal gas constant
M = 0.0289644 #molar mass of air
N0rain = 8.0e6 # number conc. for rain
#N0snow = 2e7 #constant number conc. for snow
N0gra = 4.0e6 # number conc. for graupel
#gamma7 = 720. # gamma 7 function
rhol = 1000. # density of liquid
rhos = 100. # density of snow
rhog = 400. # density of graupel
rhoi = 917. # density of ice
alpha = 0.224 # alpha (rhol/rhoi)^2*(mod(Ki)^2/mod(Kl)^2)
Cpd = 1005.7 # J/kg*K
Lv = 2264760. # J/kg

#general info
init = str(''.join(times[0])).replace('_',' ') # model init time
alltimes = [] #list to hold all times

## COLOR MAPS
# define color maps for when the matplotlib ones arent very good
# color maps created with this amazing online hex gradient color generator
# http://www.strangeplanet.fr/work/gradient-generator/index.php
nws_precip_colors = [
    "#04e9e7",  # 0.01 - 0.10 inches
    "#019ff4",  # 0.10 - 0.25 inches
    "#0300f4",  # 0.25 - 0.50 inches
    "#02fd02",  # 0.50 - 0.75 inches
    "#01c501",  # 0.75 - 1.00 inches
    "#008e00",  # 1.00 - 1.50 inches
    "#fdf802",  # 1.50 - 2.00 inches
    "#e5bc00",  # 2.00 - 2.50 inches
    "#fd9500",  # 2.50 - 3.00 inches
    "#fd0000",  # 3.00 - 4.00 inches
    "#d40000",  # 4.00 - 5.00 inches
    "#bc0000",  # 5.00 - 6.00 inches
    "#f800fd",  # 6.00 - 8.00 inches
    "#9854c6",  # 8.00 - 10.00 inches
    "#fdfdfd"]   # 10.00+

precip_colormap = matplotlib.colors.ListedColormap(nws_precip_colors)

nws_dbz_colors = [
    "#FFFFFF", # 0
    "#808080", # 5
    "#00FFFF", # 5-10
    "#00BFFF", # 10-15
    "#0000FF", # 15-20
    "#00FF00", # 20-25
    "#32CD32", # 25-30
    "#008000", # 30-35
    "#FFFF00", # 35-40
    "#DAA520", # 40-45
    "#FFA500", # 45-50
    "#FF0000", # 50-55
    "#8B0000", # 55-60
    "#800000", # 60-65
    "#FF00FF", # 65-70
    "#8A2BE2", # 70-75
    "#FFFFFF"] # 75+ (usually white but map is white)

dbz_colormap=matplotlib.colors.ListedColormap(nws_dbz_colors)

lcl_colors = ["#FFFFFF", #0 m
              "#867FFF", #500
              "#0D00FF", #1000
              "#0B7F7F", #1500
              "#09FF00", #2000
              "#84D200", #2500
              "#FFA600", #3000
              "#FF5300", #3500
              "#FF0000", #4000
              "#FF0075", #4500
              "#FF00EA"] #5000
              
lcl_colormap=matplotlib.colors.ListedColormap(lcl_colors)

snow_colors = ["#FFFFFF", "#C2BFFF", "#867FFF", "#493FFF", "#0D00FF", "#0C3FBF", "#0B7F7F", "#0ABF3F", "#09FF00", "#46E800", "#84D200", "#C1BC00", "#FFA600", "#FF7C00", "#FF5300", "#FF2900", "#FF0000", "#FF004E", "#FF009C", "#FF00EA"]
snow_colormap=matplotlib.colors.ListedColormap(snow_colors)

hail_colors = ["#FFFFFF", "#AEAAFF", "#5D55FF", "#0D00FF", "#0B55AA", "#0AAA55", "#09FF00", "#5BE100", "#ADC300", "#FFA600", "#FF6E00", "#FF3700", "#FF0000", "#FF0075", "#FF00EA"]
hail_colormap=matplotlib.colors.ListedColormap(hail_colors)

### BEGIN FUNCTIONS ###

def savefile(filename): #save plot image as png
    print "Saving file: ", filename           
    #print filename
    plt.savefig(outdir+filename)
    
def makeplot(data,title,cblabel,clevs,cbticks,ftitle): # function to make plots
    fig = plt.gcf() #get current fig
    ax = plt.gca()  #get current axis
    #ax  = fig.add_axes([0.1,0.1,0.8,0.8])
    # draw parallels and meridians
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
    # draw coastlines, state and country boundaries
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()

    # set plot title
    ax.set_title(title+currtime)
    fig.suptitle('Init: '+init+'\n', fontsize=12) #init title
    if clevs is False:
        # No color bar
        pass
    else: #create color bar
        cbar = m.colorbar(data,location='bottom',pad="5%")
        cbar.set_label(cblabel)
        if cbticks:
            cbar.set_ticks(clevs)
            cbar.ax.tick_params(labelsize=8)
    if opt.save:        
        #create filename for image and save file
        filename = ftitle+filetime+'.png'
        #filename = ftitle+str(time)+'.png' #save file with number instead of date and time
        savefile(filename) #save image file
    else:
        plt.show()

def unstagger(data,dim): #function to unstagger data from staggered grid given dimension
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

# weighted linear interpolation of 3d array of data to pressure level
def linear_interp(data, totalp, plev):
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
    outVal = ( belowVal * (1 - weight) ) + ( aboveVal * weight)

    return outVal

## THIS FUNCTION NOT WORKING CORRECTLY ##        
def interp_pres_lvls(preslvls,data,modelpreslvls): #function to interpolate data to given pressure levels
    #set up dimensions of interpolated output data
    #height_dim = np.shape(data)[0]
    we_dim = np.shape(data)[1]
    sn_dim = np.shape(data)[2]
    
    interpdata = np.zeros((len(preslvls),we_dim,sn_dim)) #array for new data

    for i in range(we_dim): #loop to do interpolation
        for j in range(sn_dim):
            interpdata[:,i,j] = np.interp(preslvls,data[:,i,j],modelpreslvls[:,i,j])
        
    return interpdata
##    
    
def t2wind(): # plot t2 and wind barbs
    # create figure
    plt.figure(figsize=(8,8))
    t2c = t2[time]-273.15 #convert temp to celcius
    if opt.tunit == 'F':        
        t2f = (9./5. * t2c)+32 #convert celcius to fahrenheit
        clevs = np.arange(-30,115,5) # levels / degF
        cs = m.contourf(x,y,t2f,clevs,cmap=cm2.get_cmap('gist_ncar'))
    elif opt.tunit == 'C':
        clevs = np.arange(-40,55,5) # levels / degC
        cs = m.contourf(x,y,t2c,clevs,cmap=cm2.get_cmap('gist_ncar'))
        
    #make x and y grid points for barbs
    #yy = np.arange(0, len(y), 8)
    #xx = np.arange(0, len(x), 8)
    #gp = np.meshgrid(yy, xx)

    #print x[::thin,::thin].shape #check x co-ord thinning
    #print u[time,::thin,::thin].shape #check u10 thinning
    #x_th,y_th = m(xlong[0,::thin,::thin],xlat[0,::thin,::thin]) #another method to thin barbs
    
    #convert wind to kts        
    u10kts = u10[time]*1.94384449
    v10kts = v10[time]*1.94384449    
    m.barbs(x[::thin,::thin], y[::thin,::thin], u10kts[::thin,::thin], v10kts[::thin,::thin],length=opt.barbsize) #plot barbs
    title = "2m Temperature and Wind Barbs (kts) \n Valid: " 
    ftitle = "t2-wind-"   
    if opt.tunit == 'C':        
        cblabel = r'$\degree$C'
    elif opt.tunit == 'F':
        cblabel = r'$\degree$F'
    cbticks = True
    makeplot(cs,title,cblabel,clevs,cbticks,ftitle)
    
def mslponly(): # plot MSLP only
    #create figure
    plt.figure(figsize=(8,8))
    psfchpa = psfc[time]*0.01 #convert Pa to hPa
    #t2c = t2[time]-273.15 #convert temp to celcius      
    mslp = psfchpa*(np.exp((g*thgt[0])/(R*t2[time]))) # calc MSLP
    mslp = gaussian_filter(mslp, sigma=3) #smooth wiggles
    #mslp = psfchpa*np.exp((g*M*thgt[0])/(Rg*t2[time]))        
    clevs = np.arange(900,1055,2.)        
    cs = m.contour(x,y,mslp,clevs,colors='k',linewidths=2.)
    plt.clabel(cs, inline=True, fmt='%1.0f', fontsize=12, colors='k')
    title = "MSLP (hPa) \n Valid: "
    ftitle = 'mslp-'
    cblabel = ''
    clevs = False # no color bar levels
    cbticks = False
    makeplot(cs,title,cblabel,clevs,cbticks,ftitle)
    
def precipaccum(): # plot total precip accumulation
    # create figure
    plt.figure(figsize=(8,8))
    ppn = rainc[time]+rainnc[time] #ppn / mm
    if opt.punit == 'mm':
        clevs = [0.1,0.5,1,2,5,10,15,20,30,40,50,80,100,200,300,500] #levels / mm
    elif opt.punit == 'in':
        clevs = [0.01, 0.1, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, \
                6.0, 8.0, 10., 20.0] # levels / in
        ppn *= 0.0393701 # convert ppn to inches . comment out to have mm
    norm = matplotlib.colors.BoundaryNorm(clevs, 15) # set boundary of data by normalizing (0,1)
    cs = m.contourf(x,y,ppn,clevs,norm=norm,cmap=precip_colormap) #plot total
    title = "Precipitation Accumulation \n Valid:  "
    ftitle = 'ppnaccum-'
    if opt.punit == 'mm':        
        cblabel = 'mm'
    elif opt.punit == 'in':
        cblabel = 'inches'
    cbticks = True    
    makeplot(cs,title,cblabel,clevs,cbticks,ftitle)
    
def precip(): # plot current precip at each time
    # create figure
    plt.figure(figsize=(8,8))
    ppn = rainc[time]+rainnc[time] # total ppn / mm
    currppn = np.array(ppn.shape)
    if time == 0: # initial amount
        currppn = ppn
    else: # current amount
        prev = rainc[time-1]+rainnc[time-1]
        currppn = ppn-prev
    if opt.punit == 'mm':
        clevs = [0.1,0.5,1,2,5,10,15,20,30,40,50,80,100,200,300,500] #levels / mm
    elif opt.punit == 'in':
        clevs = [0.01, 0.1, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, \
                6.0, 8.0, 10., 20.0] # levels / in
        currppn *= 0.0393701 # convert ppn to inches . comment out to have mm
    norm = matplotlib.colors.BoundaryNorm(clevs, 15) # set boundary of data by normalizing (0,1)
    cs = m.contourf(x,y,currppn,clevs,norm=norm,cmap=precip_colormap) #plot total
    title = "Precipitation \n Valid:  "
    ftitle = 'ppn-'
    if opt.punit == 'mm':        
        cblabel = 'mm'
    elif opt.punit == 'in':
        cblabel = 'inches'
    cbticks = True    
    makeplot(cs,title,cblabel,clevs,cbticks,ftitle)
    
def convprecip(): # plot current convective precip at each time
    # create figure
    plt.figure(figsize=(8,8))
    convppn = rainc[time] #ppn / mm
    currppn = np.array(convppn.shape)
    if time == 0:
        currppn = convppn
    else:
        prev = rainc[time-1]
        currppn = convppn-prev
    if opt.punit == 'mm':
        clevs = [0.1,0.5,1,2,5,10,15,20,30,40,50,80,100,200,300,500] #levels / mm
    elif opt.punit == 'in':
        clevs = [0.01, 0.1, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, \
                6.0, 8.0, 10., 20.0] # levels / in
        currppn *= 0.0393701 # convert ppn to inches . comment out to have mm
    norm = matplotlib.colors.BoundaryNorm(clevs, 15) # set boundary of data by normalizing (0,1)
    cs = m.contourf(x,y,currppn,clevs,norm=norm,cmap=precip_colormap) #plot total
    title = "Convective Precipitation \n Valid:  "
    ftitle = 'convppn-'
    if opt.punit == 'mm':        
        cblabel = 'mm'
    elif opt.punit == 'in':
        cblabel = 'inches'
    cbticks = True    
    makeplot(cs,title,cblabel,clevs,cbticks,ftitle)
    
def tdrh(): # plot td and rh
    # create figure
    plt.figure(figsize=(8,8))
    q2 = nc.variables['Q2'] # water vapour mixing ratio at 2m
    t2c = t2[time]-273.15 #convert temp to celcius
    psfchpa = psfc[time]*0.01
    #using equations from NOAA NWS page at http://www.srh.noaa.gov/epz/?n=wxcalc
    #es = 6.11*10**((7.5*t2c)/(237.3+t2c)) #calculate saturation vapour pressure        
    # from Rogers and Yau : A Short Course in Cloud Physics        
    es = 6.112*np.exp(17.67*t2c/(t2c+243.5)) #saturation vapour pressure in hPa       
    #es = 6.112*np.exp(17.67*(t2[time]-273.15)/(t2[time]-29.65))    
    #from Wallace and Hobbs
    #w = q2[time]/(1-q2[time]) # saturation mixing ratio
    #print w
    ws = 0.622*(es/(psfchpa-es)) #calc mixing ratio
    #print ws
    #print es
    #rh = 100*(w/ws) #relative humidty
    #e = (w*psfc[time]/(.622+w))/100
    #e = es*(rh/100) #calc vapour pressure / Pa
    #rh = (e/es)*100
    #td = (243.5*np.log(e/6.112))/(17.67-np.log(e/6.112))
    #td = np.log10(e/6.112) * (243.5/(17.67 - np.log10(e/6.112))) # calc Td another method using e **UNDERESTIMATES**
    if opt.rh:
        rh = (q2[time]/ws)*100 #calc relative humidity
        rh = np.where(rh > 100,100, rh) #remove values > 100
        rh = np.where(rh < 0,0, rh) #remove values < 0 (just in case)        
        clevs = np.arange(0,105,5)     
        cs = m.contourf(x,y,rh,clevs,cmap=cm2.get_cmap('jet')) #plot RH
        #print np.min(rh), np.max(rh)
        #print np.max(ws)
        cblabel='RH \ %'
        title = "Relative Humidity \n Valid: "
        ftitle = 'rh-'
        cbticks = True
    elif opt.td: #more efficient to have it here instead of above (same as two funcs)
        rh = (q2[time]/ws)*100 #calc relative humidity
        rh = np.where(rh > 100,100, rh) #remove values > 100
        rh = np.where(rh < 0,0, rh) #remove values < 0 (just in case) 
        td1 = 237.3*np.log10( (es*rh) / 611) #calc numerator
        td2 = 7.5*np.log10(10)-np.log10( (es*rh) / 611) #calc denominator
        td = td1/td2 #calc Td . GOOD METHOD
        title = "2m Dew Point \n Valid: "
        ftitle = 'td-'
        if opt.tunit == 'C':
            clevs = np.arange(-30,65,5) # levels / degC
            cblabel = r'$\degree$C'
        elif opt.tunit == 'F':
            clevs = np.arange(-20,125,5) # levels / degF
            td = (9./5. * td)+32 #convert celcius to fahrenheit
            cblabel = r'$\degree$F'
        cs = m.contourf(x,y,td,clevs,cmap=cm2.get_cmap('gist_ncar')) #plot Td
        cbticks=True
    
    makeplot(cs,title,cblabel,clevs,cbticks,ftitle)
    
def upperair(): # plot upper air chart for given level. geopotential height, wind bards and temp
    pb = nc.variables['PB'] #base state pressure, Pa
    p = nc.variables['P'] # perturbation pressure, Pa
    totalp = pb[time,:,:,:]+p[time,:,:,:] # total pressure in Pa
    U = nc.variables['U'][time] # U wind component
    V = nc.variables['V'][time] # V wind component
    Unew = unstagger(U,'U') # unstagger u
    Vnew = unstagger(V,'V') # unstagger v
    ph = nc.variables['PH'] #perturbation geopotential
    phb = nc.variables['PHB'] #base state geopotential
    totalgp = phb[time,:,:,:]+ph[time,:,:,:] # total geopotential
    totalgp = unstagger(totalgp,'Z') #total geopotential unstaggered
    t = nc.variables['T'] #perturbation potential temperature (theta-t0)
    t00 = nc.variables['T00'] #base state theta
    #print t.shape, t00.shape
    totalTheta = t[time,:,:,:]+t00[0] # total potential temp
    totalTfac = ( (totalp*0.01) / 1000. )**(R/1004.) # factor to multiply theta by
    totalT = (totalTheta*totalTfac)-273.15 # calc temp in deg C
    #print totalT.shape
    #print pb.shape, p.shape, totalp.shape
    #print U.shape, V.shape, Unew.shape, Vnew.shape
    #print totalgp.shape
    
    #pressure levels to interpolate to in Pa for interp_pres_lvls
    #preslvls = np.array([100000,95000,92500,90000,85000,80000,75000,70000,65000,60000,55000,50000,45000,40000,35000,30000,25000,20000,15000,10000,5000,1000])
    #preslvls = np.array([1000,5000,10000,15000,20000,25000,30000,35000,40000,45000,50000,55000,60000,65000,70000,75000,80000,85000,90000,92500,95000,100000])
    
    # interpolate using np.interp 
    #gphgt = interp_pres_lvls(preslvls,totalgp,totalp) 
    #totalTfinal = interp_pres_lvls(preslvls,totalp,totalT)
    
    levels = opt.lvl.split(',') # get list of levels
    for level in levels: 
        plt.figure(figsize=(8,8)) #create fig for each plot
        level = int(level) # make it int
        #interp data for level
        gphgt = linear_interp(totalgp,totalp,level)
        totalTfinal = linear_interp(totalT,totalp,level)
        uinterp = linear_interp(Unew,totalp,level)
        vinterp = linear_interp(Vnew,totalp,level) 
        Ufinal = uinterp*1.94384449  #convert to kts
        Vfinal = vinterp*1.94384449
        gphgt = (gphgt/9.81) # convert to height (m)
        gphgt = gaussian_filter(gphgt, sigma=3) # smooth wiggles  
        totalTfinal = gaussian_filter(totalTfinal, sigma=2)
        # set gpheight levels for common pressure levels        
        if level == 250:
            gpclevs = np.arange(np.min(gphgt),np.max(gphgt),60)
        elif level == 500:
            gpclevs = np.arange(np.min(gphgt),np.max(gphgt),60)
        elif level == 700:
            gpclevs = np.arange(np.min(gphgt),np.max(gphgt),30)
        elif level == 850:
            gpclevs = np.arange(np.min(gphgt),np.max(gphgt),30)
        elif level == 925:
            gpclevs = np.arange(np.min(gphgt),np.max(gphgt),30)
        else: # use generic 30m spacing
            gpclevs = np.arange(np.min(gphgt),np.max(gphgt),30)
        
        #plot all this up    
        cs = m.contour(x,y,gphgt,gpclevs,colors='k',linewidths=2.)
        plt.clabel(cs, inline=True, fmt='%1.0f', fontsize=12, colors='k')
        tclevs = np.arange(np.min(totalTfinal),np.max(totalTfinal),4)
        cs2 = m.contour(x,y,totalTfinal,tclevs,colors='r',linestyles='-',linewidths=2.)
        plt.clabel(cs2,inline=True,fmt='%1.0f',fontsize=12,colors='r')
        m.barbs(x[::thin,::thin], y[::thin,::thin], Ufinal[::thin,::thin], Vfinal[::thin,::thin],length=opt.barbsize) #plot barbs
        level = str(level)
        title = level+'mb Height (m), Temp (C), Wind Barbs (kts) \n Valid: '
        ftitle = level+'mb-' 
        cblabel = 'kts'
        clevs = False
        cbticks = False
        makeplot(cs,title,cblabel,clevs,cbticks,ftitle)
        
def surface(): # plot surface chart. t2, wind barbs and mslp
    # create figure
    plt.figure(figsize=(8,8))
    t2c = t2[time]-273.15 #convert temp to celcius
    if opt.tunit == 'F':
        t2f = (9./5. * t2c)+32 #convert celcius to fahrenheit
        clevs = np.arange(-30,115,5) # levels / degF
        cs = m.contourf(x,y,t2f,clevs,cmap=cm2.get_cmap('gist_ncar'))
        cblabel = r'$\degree$F'
    elif opt.tunit == 'C':
        clevs = np.arange(-40,55,5) # levels / degC
        cs = m.contourf(x,y,t2c,clevs,cmap=cm2.get_cmap('gist_ncar'))
        cblabel = r'$\degree$C'
    cbticks = True    
    
    psfchpa = psfc[time]*0.01 #convert Pa to hPa
    mslp = psfchpa*(np.exp((g*thgt[0])/(R*t2[time]))) # calc MSLP
    mslp = gaussian_filter(mslp, sigma=3) # smooth wiggles
    #make x and y grid points for barbs
    #yy = np.arange(0, len(y), 8)
    #xx = np.arange(0, len(x), 8)
    #gp = np.meshgrid(yy, xx)

    #print x[::thin,::thin].shape #check x co-ord thinning
    #print u[time,::thin,::thin].shape #check u10 thinning
    #x_th,y_th = m(xlong[0,::thin,::thin],xlat[0,::thin,::thin]) #another method to thin barbs
        
    #convert wind to kts        
    u10kts = u10[time]*1.94384449
    v10kts = v10[time]*1.94384449      
    m.barbs(x[::thin,::thin], y[::thin,::thin], u10kts[::thin,::thin], v10kts[::thin,::thin],length=opt.barbsize) #plot barbs
    title = "2m Temp, Wind Barbs (kts), MSLP (hPa) \n Valid: "
    ftitle = 'sfc-'       
    pclevs = np.arange(900,1055,2.)        
    pcs = m.contour(x,y,mslp,pclevs,colors='k',linewidths=2.)
    plt.clabel(pcs, inline=True, fmt='%1.0f', fontsize=12, colors='k')
   
    makeplot(cs,title,cblabel,clevs,cbticks,ftitle)
    
def snowaccum(): # plot snow accumulation
    # create figure
    plt.figure(figsize=(8,8))
    snowin = nc.variables['SNOWNC'] # total accumulated grid scale snow and ice / mm
    snow = snowin[time] # snow accum in mm for each time
    if opt.punit == 'mm':
        clevs = [0,0.5,1,2.5,3,4,5,8,10,15,20,30,40,50,80,100,150,200,250,500]
        cblabel = 'mm'
    elif opt.punit == 'in':
        snow = snow*0.0393701 # convert to inches
        clevs = [0.25,0.5,0.75,1,1.5,2,2.5,3,4,5,6,8,10,12,14,16,18,20,22,24]
        cblabel = 'inches'
    cbticks = True
    norm = matplotlib.colors.BoundaryNorm(clevs, 19) # set boundary of data by normalizing (0,1)    
    cs = m.contourf(x,y,snow,clevs,norm=norm,cmap=snow_colormap)
    title = "Snow Accumulation \n Valid: "
    ftitle = 'snow-'

    makeplot(cs,title,cblabel,clevs,cbticks,ftitle)
    
def hailaccum(): # plot hail accumulation
    # create figure
    plt.figure(figsize=(8,8))
    hailin = nc.variables['HAILNC'] # accimulated total grid scale hail / mm
    hail = hailin[time] #hail accum at each time
    if opt.punit == 'mm':
        clevs = [0.5,1.,1.5,2.,2.5,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.]
        cblabel = 'mm'
    elif opt.punit == 'in':
        hail = hail*0.0393701 # convert to inches
        clevs = [0.1,0.2,0.3,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.]
        cblabel = 'inches'
    cbticks = True
    norm = matplotlib.colors.BoundaryNorm(clevs, 14) # set boundary of data by normalizing (0,1)    
    cs = m.contourf(x,y,hail,clevs,norm=norm,cmap=hail_colormap)
    title = "Hail Accumulation \n Valid: "
    ftitle = 'hail-'

    makeplot(cs,title,cblabel,clevs,cbticks,ftitle)
    
    ####
    # Reference for simulated and composite reflectivity:   
    # Koch, S., Ferrier, B., Stoelinga, M., Szoke, E., Weiss, S., Kain, J., 2005:
    # THE USE OF SIMULATED RADAR REFLECTIVITY FIELDS IN THE DIAGNOSIS OF
    # MESOSCALE PHENOMENA FROM HIGH-RESOLUTION WRF MODEL FORECASTS
    # J4J.7
    ####
def simudbz(): # plot simulated reflectivity, mp_physics dependent
    # create figure
    plt.figure(figsize=(8,8))
    qrain = nc.variables['QRAIN'] # rain water mixing ratio
    t2c = t2[time]-273.15 #convert temp to celcius
    rhoa = ( psfc[time] / (R*t2[time]) ) # density of air
    N0snow = 2.0e6*np.exp(-0.12*t2c) # N0 for snow (intercept parameter)
    Qrain = qrain[time,1] # rain mixing ratio
    Qrain = np.nan_to_num(Qrain) # change NaN to zeroes, changge infs to nums
    try: #depends on MP scheme
        Qsn = nc.variables['QSNOW'] # try to get snow mixing ratio
    except:
        Qsn = np.zeros(np.shape(qrain)) # else create zeros array same shape as qrain
    Qsnow = Qsn[time,1] # snow mixing ratio
    Qsnow = np.nan_to_num(Qsnow) # change NaN to zeros
    # using these seems tidier for these calcs
    lambr = np.divide((np.pi*N0rain*rhol), np.multiply(rhoa,Qrain))**0.25 # rain lambda slope factor
    lambs = np.exp(-0.0536*t2c) # snow lambda slope factor
    Zer = ( (720.0*N0rain)*(lambr**-7.0) )*1e18 # rain equivalent reflectivity factor
    Zes_int = np.divide((lambs*Qsnow*rhoa), N0snow)
    Zes = ( (0.224*720.*1e18) / (np.pi*rhol)**2 )*(Zes_int**2) # snow equivalent reflectivity factor
    Ze = np.add(Zer,Zes) # total reflectivity factor
    dBZ = 10*np.log10(Ze) # compare apples to apples / convert Ze to dBZ
    dBZ = np.nan_to_num(dBZ) # get rid of NaN's to zeros
    #Zeraindbz = 10*np.log10(Zerain)
    #Zeraindbz = np.nan_to_num(Zeraindbz)
    clevs = np.arange(0,85,5)
    norm = matplotlib.colors.BoundaryNorm(clevs, 17) # normalize levels
    cs = m.contourf(x,y,dBZ,clevs,norm=norm,cmap=dbz_colormap)
    title = "Simulated Reflectivity \n Valid: "
    ftitle = 'simdbz-'
    cblabel = 'dBZ'
    cbticks = True
    makeplot(cs,title,cblabel,clevs,cbticks,ftitle)

def compodbz(): # plot composite reflectivity, mp_physics dependent
    # create figure
    plt.figure(figsize=(8,8))
    try: #get refl from do_radar_ref=1
        refl = nc.variables['REFL_10CM'][time]
        dBZ = np.zeros(refl[0,0].shape)
        dBZ = np.max(refl, axis=0)
        #for i in range(len(refl[1,:,1])):
         #   for j in range(len(refl[1,1,:])):
          #      dBZ[i,j]=np.max(refl[:,i,j])
    except: # calculate reflectivity
        qrain = nc.variables['QRAIN'] # rain water mixing ratio
        Qrainall = qrain[time] # get all Qrain values at all levels for each time
        t2c = t2[time]-273.15 #convert temp to celcius
        rhoa = ( psfc[time] / (R*t2[time]) ) # density of air
        N0snow = 2.0e6*np.exp(-0.12*t2c) # N0 for snow (intercept parameter)
        try: # depends on MP scheme
            Qsn = nc.variables['QSNOW'] # try to get snow mixing ratio
        except:
            Qsn = np.zeros(np.shape(qrain)) # else create zeros array same shape as qrain
        Qsnowall = Qsn[time] # get all Qsnow values at all levels for each time
        # using these seems tidier for these calcs
        Qrainmax = np.max(Qrainall, axis=0) #max rain QV
        Qsnowmax = np.max(Qsnowall, axis=0) #max snow QV
        lambr = np.divide((np.pi*N0rain*rhol), np.multiply(rhoa,Qrainmax))** 0.25 # rain lambda slope factor
        lambs = np.exp(-0.0536 * t2c) # snow lambda slope factor
        Zer = ( (720.0*N0rain)*(lambr ** -7.0) ) * 1e18 # rain equivalent reflectivity factor
        Zes_int = np.divide((lambs * Qsnowmax * rhoa), N0snow)
        Zes = ( (0.224 * 720. * 1e18) / (np.pi * rhol) ** 2 ) * (Zes_int **2) # snow equivalent reflectivity factor
        Ze = np.add(Zer,Zes) # total reflectivity factor
        dBZ = 10.0*np.log10(Ze) # compare apples to apples / convert Ze to dBZ
        dBZ = np.nan_to_num(dBZ) # get rid of NaN's to zeros
    clevs = np.arange(0,85,5)
    norm = matplotlib.colors.BoundaryNorm(clevs, 17) # normalize levels
    cs = m.contourf(x,y,dBZ,clevs,norm=norm,cmap=dbz_colormap)
    title = "Composite Reflectivity \n Valid: "
    ftitle = 'compdbz-'
    cblabel = 'dBZ'
    cbticks = True
    makeplot(cs,title,cblabel,clevs,cbticks,ftitle)
    
def lclhgt(): # plot lcl height
    # create figure
    plt.figure(figsize=(8,8))
    q2 = nc.variables['Q2'] # water vapour mixing ratio at 2m
    t2c = t2[time]-273.15 #convert temp to celcius
    psfchpa = psfc[time]*0.01
    #using equations from NOAA NWS page at http://www.srh.noaa.gov/epz/?n=wxcalc
    #es = 6.11*10**((7.5*t2c)/(237.3+t2c)) #calculate saturation vapour pressure        
    # from Rogers and Yau : A Short Course in Cloud Physics        
    es = 6.112*np.exp(17.67*t2c/(t2c+243.5)) #saturation vapour pressure in hPa       
    #from Wallace and Hobbs
    #w = q2[time]/(1-q2[time]) # saturation mixing ratio
    #print w
    ws = 0.622*(es/(psfchpa-es)) #calc mixing ratio
    #print ws
    #print es
    #rh = 100*(w/ws) #relative humidty
    #e = es*(rh/100) #calc vapour pressure / Pa
    rh = (q2[time]/ws)*100. #calc relative humidity
    rh = np.where(rh > 100.,100., rh) #remove values > 100
    rh = np.where(rh < 0,0, rh) #remove values < 0 (just in case)
    #td = np.log10(e/6.112)*(243.5/(17.67-np.log10(e/6.112))) # calc Td another method using e UNDERESTIMATES
    td1 = 237.3*np.log10((es*rh)/611.) #calc numerator
    td2 = 7.5*np.log10(10.)-np.log10( (es*rh) / 611. ) #calc denominator
    td = td1/td2 #calc Td . GOOD METHOD
    lcl = 125.*(t2c-td) # calc lcl
    lcl = np.where(lcl < 0, 0, lcl)    
    #print np.min(lcl), np.max(lcl)
    lcl = np.nan_to_num(lcl) # get rid of NaNs
    clevs = np.arange(0,6000,500)
    cs = m.contourf(x,y,lcl,clevs,cmap=lcl_colormap)
    title = "LCL Height \n Valid: "
    ftitle = 'lcl-'
    cblabel = 'm'
    cbticks = True
    makeplot(cs,title,cblabel,clevs,cbticks,ftitle)
    
def thetaE(): # plot theta-e
    # create figure
    plt.figure(figsize=(8,8))
    t = nc.variables['T'] #perturbation potential temperature (theta-t0)
    t00 = nc.variables['T00'] #base state theta
    psfchpa = psfc[time]*0.01
    t2c = t2[time]-273.15 #convert temp to celcius
    theta = t[time,0]+t00[0] # total theta
    #es = 6.112*np.exp((17.27*(t2[time]-273.16))/(t2[time]-35.86))
    es = 6.112*np.exp( 17.67*t2c / (t2c+243.5) ) #saturation vapour pressure in hPa       
    ws = 0.622*(es/(psfchpa-es)) # water vapour mixing ratio
    # Version from Wallace and Hobbs
    thetae = theta*np.exp( (Lv*ws) / (1004.*t2[time]) ) # calc equivalent potential temperature
    #Simplified version from Stull
    #thetae = (t2[time]+(Lv/Cpd)*q2[time])*(1000/psfchpa)**(R/Cpd)
    clevs = np.arange(260,372,4) # set by max and min of data    
    cs = m.contourf(x,y,thetae,clevs,cmap=cm2.get_cmap('gist_ncar'))
    title = "Theta-e \n Valid: "
    ftitle = 'thetae-'
    cblabel = 'K'
    cbticks = True
    
    makeplot(cs,title,cblabel,clevs,cbticks,ftitle)
    
def h75lr(): # 700-500mb lapse rates
    # create figure
    plt.figure(figsize=(8,8))
    pb = nc.variables['PB'] #base state pressure, Pa
    p = nc.variables['P'] # perturbation pressure, Pa
    totalp = pb[time,:,:,:]+p[time,:,:,:] # total pressure in Pa
    t = nc.variables['T'] #perturbation potential temperature (theta-t0)
    t00 = nc.variables['T00'] #base state theta
    #print t.shape, t00.shape
    totalTheta = t[time,:,:,:]+t00[0] # total potential temp
    totalTfac = ( (totalp*0.01) / 1000. )**(R/1004.) # factor to multiply theta by
    totalT=(totalTheta*totalTfac)-273.15 # calc temp in deg C
    # interp temps to levels
    totalT700 = linear_interp(totalT,totalp,700)
    totalT500 = linear_interp(totalT,totalp,500)
    # calc h7-h5 lapse rates
    lr = totalT700 - totalT500
    #print np.min(lr), np.max(lr)
    lr = np.round(lr, decimals=2) #round off
    
    clevs = np.arange(5,10.5,.5) # these levels seem to be of importance
    #clevs=False
    cs = m.contourf(x,y,lr,clevs,cmap=cm2.get_cmap('gist_ncar'))
    title = "H7-H5 Lapse Rates \n Valid: "
    ftitle = 'h75lr-'
    cblabel = r'$\degree$C'
    cbticks = True
    makeplot(cs,title,cblabel,clevs,cbticks,ftitle)
    
def absvort500(): # plot 500mb absolute vorticity
    # create figure
    plt.figure(figsize=(8,8))
    pb = nc.variables['PB'] #base state pressure, Pa
    p = nc.variables['P'] # perturbation pressure, Pa
    totalp = pb[time,:,:,:]+p[time,:,:,:] # total pressure in Pa
    U = nc.variables['U'][time] # U wind component
    V = nc.variables['V'][time] # V wind component
    Unew = unstagger(U,'U') # unstagger u
    Vnew = unstagger(V,'V') # unstagger v
    fcoriolis = 2 * ( 2*np.pi/86400 ) * np.sin( xlat[0]*(np.pi/180.) ) # calc coriolis term
    uinterp = linear_interp(Unew,totalp,500) #interp to 500mb
    vinterp = linear_interp(Vnew,totalp,500) 
    dvdx = np.gradient(vinterp,dx,dx)[1] # calc dvdx
    dudy = np.gradient(uinterp,dx,dx)[0] # calc dudy
    avort = dvdx - dudy + fcoriolis # absolute vorticity
    clevs = np.arange(-4,50,2) # levels used on COD
    cs = m.contourf(x,y,avort,clevs,cmap=cm2.get_cmap('gist_ncar'))
    title = '500mb Absolute Vorticity \n Valid: '
    ftitle = '500absvort-' 
    cblabel = r'$s^{-1}$'
    cbticks = True
    makeplot(cs,title,cblabel,clevs,cbticks,ftitle)
    
#### END FUNCTIONS ####
flag = False # to check for plotting options
#### BEGIN TIME LOOP ####    
for time in range(times.shape[0]):
    currtime = str(''.join(times[time])).replace('_', ' ') #get current model time    
    filetime = currtime.translate(None, ':').replace(' ', '_') # time for filename
    alltimes.append(currtime) # all times in output    
    
    if opt.t2: #plot 2m temp and wind barbs
        print "Plotting Temperature and Wind Barbs for time: ", currtime       
        t2wind()
        flag = True
            
    if opt.mslp: #plot surface pressure only
        print "Plotting MSLP for time: ", currtime
        mslponly()
        flag = True

    if opt.ppnaccum: #plot total precipitation
        print "Plotting Precipitation Accumulation for time: ", currtime
        precipaccum()
        flag = True
        
    if opt.ppn: # plot current ppn
        print "Plotting Precipitation for time: ", currtime
        precip()
        flag = True

    if opt.convppn: # plot convective ppn
        print "Plotting Convective Precipitation for time: ", currtime
        convprecip()
        flag = True        
        
    if opt.td or opt.rh: #plot dew point or RH
        flag = True        
        if opt.td:        
            print "Plotting Dew Point for time: ", currtime
        elif opt.rh:
            print "Plotting RH for time: ", currtime
        tdrh()
        
    if opt.ua: #plot upper air charts
        print "Plotting upper level chart for time: ", currtime
        upperair()
        flag = True
        
    if opt.sfc: #plot surface chart. t2, wind and mslp
        print "Plotting Surface Chart for time: ", currtime
        surface()
        flag = True
        
    if opt.snow: #plot snow accumulation
        print "Plotting Snow Accumulation for time: ", currtime
        snowaccum()
        flag = True        
    if opt.hail: #plot hail accumulation
        print "Plotting Hail Accumulation for time: ", currtime
        hailaccum()
        flag = True
        
    if opt.simdbz: #simulated reflectivity
        print "Plotting Simulated Reflectivity for time: ", currtime
        simudbz()
        flag = True
        
    if opt.compdbz: #composite reflectivity
        print "Plotting Composite Reflectivity for time: ", currtime
        compodbz()
        flag = True 
        
    if opt.lcl: #plot LCL
        print "Plotting LCL for time: ", currtime
        lclhgt()
        flag = True
        
    if opt.thetae: #plot theta-e
        print "Plotting Theta-e for time: ", currtime
        thetaE()  
        flag= True
        
    if opt.lr75: #plot h7-h5 lapse rates
        print "Plotting H7-H5 lapse rates for time: ", currtime
        h75lr()
        flag = True
        
    if opt.vort500: # plot 500mb absolute vorticity
        print "Plotting 500mb absolute vorticity for time: ", currtime
        absvort500()
        flag = True
        
    if flag is False: # do this when no options given
        print "Please provide options to plot. Use wrfplot.py --help"
        print "QUITTING"
        sys.exit()
        #pass
#### END TIME LOOP ####
    
if opt.verbose: #verbose output
    print "\n*VERBOSE OUTPUT*"
    print "\nindir= ", indir
    print "infile= ", filein
    print "outdir=", outdir
    print "Model initialisation time: ", init
    print "Timestep: ", nc.variables['ITIMESTEP'][1]
    print "Times in file: ", alltimes
    print "west_east: ", x_dim
    print "south_north: ", y_dim
    print "Model dimentions (metres): ", width_meters, height_meters
    print "dx, dy: ", dx, dy
    print "Center lat: ", cen_lat
    print "Center lon: ", cen_lon
    print "Model top: ", nc.variables['P_TOP'][0]
    print "Map projection: ", proj

nc.close() # close netcdf file