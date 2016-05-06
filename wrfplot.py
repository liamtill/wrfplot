#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Liam Till
University of Oklahoma / Unversity of Reading
wrfplot.py
Python 2.x
Python script to plot various WRF model output. Plots are saved as PNG.
example usage: wrfplot.py --infile filename.nc --sfc --tunit C --ppn -punit mm --td
Will plot surface chart and dewpoint in Celcius and precipitation in mm.
Use wrfplot.py --help to list all options
Last modified: 05/05/16

Skew-T plotting with the pyMeteo package available at: https://github.com/cwebster2/pyMeteo
Credit to Casey Webster

Skew-t plotting with SHARPpy package available at: https://github.com/sharppy/SHARPpy
Credit to: Patrick Marsh (SPC), Kelton Halbert (OU School of Meteorology), Greg Blumberg (OU/CIMMS), Tim Supinie (OU School of Meteorology)
"""

import matplotlib
#matplotlib.use('Agg') # UNCOMMENT THIS ONLY WHEN INVOKING FROM CRON SCRIPT 
from scipy.io import netcdf # USE SCIPY MODULE
#from netCDF4 import Dataset # UNCOMMENT TO USE NETCDF 4 MODULE
from mpl_toolkits.basemap import Basemap
from matplotlib import cm
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter
import numpy as np
import datetime
from optparse import OptionParser
import os.path
import sys
import conversions as conv
import calc_vars as calc
import plot_funcs as pltfuncs
import funcs
import colormaps as cmap

# option parser
usage="usage: %prog [options] \n example usage: wrfplot.py --infile filename.nc --sfc --tunit C --td --ppn --punit mm"
parser = OptionParser(usage=usage, version="%prog 6.0 by Liam Till")
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
parser.add_option("--ua", dest="ua", action="store_true", help="Plot geopotential height, temperature and wind barbs at given pressure levels (hPa), --lvl")
parser.add_option("--lvl", dest="lvl", help="Pressure levels to interpolate to for upper level charts option --ua, --vv. Comma seperated e.g 250,500", default="500")
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
parser.add_option("--shear06", dest="shear06", action="store_true", help="Plot the 0-6km shear")
parser.add_option("--vv", dest="vv", action="store_true", help="Plot vertical velocity at specified levels --lvl")
parser.add_option("--irtemp", dest="irtemp", action="store_true", help="Plot IR Brightness Temperature")
parser.add_option("--skewt", dest="skewt", action="store_true", help="Plot Skew-t for a location. Uses pyMeteo package.")
parser.add_option("--slat", dest="slat", type="int", help="Latitude for Skew-t")
parser.add_option("--slon", dest="slon", type="int", help="Longitude for Skew-t")
parser.add_option("--getij", dest="getij", action="store_true", help="Get i,j and nearest Lat/Lon for entered Lat/Lon")
parser.add_option("--skewt2", dest="skewt2", action="store_true", help="Plot Skew-t for a location using SHARPpy")
parser.add_option("--uh25", dest="uh25", action="store_true", help="Plot 2-5km Updraft Helicity")
(opt, arg) = parser.parse_args()

indir = opt.indir # dir of input file
filein = opt.infile
if opt.auto: # for auto file input for daily runs
    run = opt.run # model init time
    filein = 'wrfout_d01_'+datetime.datetime.utcnow().strftime('%Y-%m-%d')+'_'+run+':00:00' # auto filename for current days run
while os.path.isfile(indir+filein) is False and not opt.auto: #if file doesnt exist get filename
    print "File", filein, "not found! in directory:", indir
    indir = raw_input("Please enter a directory (blank for current dir): ")
    filein = raw_input("Please enter a filename: ")
try: #check if file exists and read in
    print "Reading in file: ", indir+filein
    #nc = Dataset(indir+filein) # for netcdf 4
    nc = netcdf.netcdf_file(indir+filein,'r') # for scipy
except: # quit if cant read file
    print "Something went wrong reading in the file"
    print "QUITTING"
    sys.exit()

outdir = opt.outdir # output image dir

## BASEMAP STUFF

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
# dimensions of domain
x_dim = len(xlat[0,0,:])
y_dim = len(xlong[0,:,0])

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
    projname = 'Lambert Conformal'
elif map_proj == 2: # polar stereographic
    proj = 'npstere'
    projname = 'Polar Stereographic'
elif map_proj == 3: # mercator
    proj = 'merc'
    projname = 'Mercator'
else: # not supported and quit
    print "Projection ", map_proj, "unknown"
    print "QUITTING"
    sys.exit()

# make map
m = Basemap(resolution='i',projection=proj,width=width_meters,height=height_meters,lat_0=cen_lat,lon_0=cen_lon,lat_1=truelat1,lat_2=truelat2) 
#m = Basemap(resolution='i',projection=proj,llcrnrlon=xlong[0,0,0],llcrnrlat=xlat[0,0,0],urcrnrlon=xlong[0,-1,-1],urcrnrlat=xlat[0,-1,-1],lat_0=cen_lat,lon_0=cen_lon)    
#x, y = m(xlong[0,:,:],xlat[0,:,:])            
# get lat/lons of ny by nx evenly space grid 
# make lons, lats and x, y co ordinates
lons, lats = m.makegrid(x_dim, y_dim)
x, y = m(lons, lats) # compute map proj coordinates.
print "Using map projection: ", projname

## GET THIS DATA FOR NOW
times = nc.variables['Times'] #each time output in wrf nc file
t2 = nc.variables['T2'] #temp at 2m / Kelvin
u10 = nc.variables['U10'] #u10 wind / ms/s
v10 = nc.variables['V10'] #v10 wind / ms/s
psfc = nc.variables['PSFC'] #surface pressure / Pascals
rainc = nc.variables['RAINC'] # accumulated total cumulus precip
rainnc = nc.variables['RAINNC'] # accumulated total grid scale precip
thgt = nc.variables['HGT'] #terrain height

# general info
init = str(''.join(times[0])).replace('_',' ') # model init time
alltimes = [] #list to hold all times

### BEGIN PLOT FUNCTIONS ###
# savefile and makeplot and the functions for putting data on maps may stay here for now #

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
    #ax.set_title(title+currtime)
    ax.text(0,1.01*height_meters,title+'\nValid:'+currtime,fontsize=14) 
    ax.text(0.65*width_meters,1.01*height_meters,'Init: '+init, fontsize=12)
    #fig.suptitle('Init: '+init+'', fontsize=12) #init title
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
        
def t2wind(): # plot t2 and wind barbs
    # create figure
    plt.figure(figsize=(8,8))
    temps = t2[time] # temps in K
    if opt.tunit == 'F':        
        t2f = conv.k_to_f(temps) # convert to F
        clevs = np.arange(-30,115,5) # levels / degF
        cs = m.contourf(x,y,t2f,clevs,cmap=cm.get_cmap('gist_ncar'))
    elif opt.tunit == 'C':
        t2c = conv.k_to_c(temps) # convert to C
        clevs = np.arange(-40,55,5) # levels / degC
        cs = m.contourf(x,y,t2c,clevs,cmap=cm.get_cmap('gist_ncar'))
        
    #make x and y grid points for barbs
    #yy = np.arange(0, len(y), 8)
    #xx = np.arange(0, len(x), 8)
    #gp = np.meshgrid(yy, xx)

    #print x[::thin,::thin].shape #check x co-ord thinning
    #print u[time,::thin,::thin].shape #check u10 thinning
    #x_th,y_th = m(xlong[0,::thin,::thin],xlat[0,::thin,::thin]) #another method to thin barbs
    
    #convert wind to kts        
    u10kts = conv.ms_to_kts(u10[time])
    v10kts = conv.ms_to_kts(v10[time])    
    m.barbs(x[::thin,::thin], y[::thin,::thin], u10kts[::thin,::thin], v10kts[::thin,::thin],length=opt.barbsize) #plot barbs
    title = "2m Temperature and Wind Barbs (kts)" 
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
    x, y = m(lons, lats)
    psfchpa = conv.pa_to_hpa(psfc[time]) #convert Pa to hPa     
    mslp = calc.calc_mslp(psfchpa, thgt[0], t2[time]) # get mslp
    mslp = gaussian_filter(mslp, sigma=3) #smooth wiggles
    #find local min and local max
    local_min, local_max = funcs.extrema(mslp, mode='wrap', window=50)       
    clevs = np.arange(900,1055,2.)        
    cs = m.contour(x,y,mslp,clevs,colors='k',linewidths=2.)
    plt.clabel(cs, inline=True, fmt='%1.0f', fontsize=12, colors='k')
    xlows = x[local_min];   xhighs = x[local_max]
    ylows = y[local_min];   yhighs = y[local_max]
    lowvals = mslp[local_min]; highvals = mslp[local_max]
    # plot lows as blue L's, with min pressure value underneath.
    xyplotted = []
    # don't plot if there is already a L or H within dmin meters.
    yoffset = 0.022*(m.ymax-m.ymin)
    dmin = yoffset
    for x,y,p in zip(xlows, ylows, lowvals):
        if x < m.xmax and x > m.xmin and y < m.ymax and y > m.ymin:
            dist = [np.sqrt((x-x0)**2+(y-y0)**2) for x0,y0 in xyplotted]
            if not dist or min(dist) > dmin:
                plt.text(x,y,'L',fontsize=14,fontweight='bold', ha='center',va='center',color='b')
                plt.text(x,y-yoffset,repr(int(p)),fontsize=12, ha='center',va='top',color='b', bbox = dict(boxstyle="square",ec='None',fc=(1,1,1,0.5)))
                xyplotted.append((x,y))
    # plot highs as red H's, with max pressure value underneath.
    xyplotted = []
    for x,y,p in zip(xhighs, yhighs, highvals):
        if x < m.xmax and x > m.xmin and y < m.ymax and y > m.ymin:
            dist = [np.sqrt((x-x0)**2+(y-y0)**2) for x0,y0 in xyplotted]
            if not dist or min(dist) > dmin:
                plt.text(x,y,'H',fontsize=14,fontweight='bold', ha='center',va='center',color='r')
                plt.text(x,y-yoffset,repr(int(p)),fontsize=12, ha='center',va='top',color='r', bbox = dict(boxstyle="square",ec='None',fc=(1,1,1,0.5)))
                xyplotted.append((x,y))
    title = "MSLP (hPa)"
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
        ppn = conv.mm_to_in(ppn) # convert ppn to inches
    norm = matplotlib.colors.BoundaryNorm(clevs, 15) # set boundary of data by normalizing (0,1)
    cs = m.contourf(x,y,ppn,clevs,norm=norm,cmap=cmap.precip_colormap) #plot total
    title = "Precipitation Accumulation"
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
        currppn = conv.mm_to_in(currppn) # convert ppn to inches
    norm = matplotlib.colors.BoundaryNorm(clevs, 15) # set boundary of data by normalizing (0,1)
    cs = m.contourf(x,y,currppn,clevs,norm=norm,cmap=cmap.precip_colormap) #plot total
    title = "Precipitation"
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
        currppn = conv.mm_to_in(currppn) # convert ppn to inches
    norm = matplotlib.colors.BoundaryNorm(clevs, 15) # set boundary of data by normalizing (0,1)
    cs = m.contourf(x,y,currppn,clevs,norm=norm,cmap=cmap.precip_colormap) #plot total
    title = "Convective Precipitation"
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
    q2 = nc.variables['Q2'][time] # water vapour mixing ratio at 2m
    t2c = conv.k_to_c(t2[time]) #convert temp to celcius
    psfchpa = conv.pa_to_hpa(psfc[time]) # pres to hPa
    es = calc.calc_es(t2c[time]) # calc es      
    ws = calc.calc_ws(es, psfchpa) # calc ws
    u10kts = conv.ms_to_kts(u10[time])
    v10kts = conv.ms_to_kts(v10[time])

    if opt.rh:
        rh = calc.calc_rh(q2, ws) #calc rh      
        clevs = np.arange(0,105,5)     
        cs = m.contourf(x,y,rh,clevs,cmap=cm.get_cmap('jet')) #plot RH
        cblabel='RH \ %'
        title = "Relative Humidity \n Valid: "
        ftitle = 'rh-'
        cbticks = True
    elif opt.td:
        rh = calc.calc_rh(q2, ws) # calc rh
        td = calc.calc_dewpoint(es, rh) # calc td (deg C)
        title = "2m Dew Point"
        ftitle = 'td-'
        if opt.tunit == 'C':
            clevs = np.arange(-30,65,5) # levels / degC
            cblabel = r'$\degree$C'
        elif opt.tunit == 'F':
            clevs = np.arange(-20,125,5) # levels / degF
            td = conv.c_to_f(td) #convert celcius to fahrenheit
            cblabel = r'$\degree$F'
        cs = m.contourf(x,y,td,clevs,cmap=cm.get_cmap('gist_ncar')) #plot Td
        m.barbs(x[::thin,::thin], y[::thin,::thin], u10kts[::thin,::thin], v10kts[::thin,::thin],length=opt.barbsize) #plot barbs
        cbticks=True
    
    makeplot(cs,title,cblabel,clevs,cbticks,ftitle)
    
def upperair(): # plot upper air chart for given level. geopotential height, wind bards and temp
    pb = nc.variables['PB'][time] #base state pressure, Pa
    p = nc.variables['P'][time] # perturbation pressure, Pa
    totalp = pb + p # total pressure in Pa
    U = nc.variables['U'][time] # U wind component
    V = nc.variables['V'][time] # V wind component
    Unew = funcs.unstagger(U,'U') # unstagger u
    Vnew = funcs.unstagger(V,'V') # unstagger v
    ph = nc.variables['PH'][time] #perturbation geopotential
    phb = nc.variables['PHB'][time] #base state geopotential
    totalgp = phb + ph # total geopotential
    totalgp = funcs.unstagger(totalgp,'Z') #total geopotential unstaggered
    theta = nc.variables['T'][time] #perturbation potential temperature (theta-t0)
    theta0 = nc.variables['T00'][0] #base state theta
    totalTheta = theta + theta0 # total potential temp
    totalT = conv.k_to_c(calc.theta_to_temp(totalTheta, totalp)) # calc temps in C
    
    levels = opt.lvl.split(',') # get list of levels
    for level in levels: 
        plt.figure(figsize=(8,8)) #create fig for each plot
        level = int(level) # make it int
        #interp data for level
        gphgt = funcs.linear_interp(totalgp,totalp,level)
        totalTfinal = funcs.linear_interp(totalT,totalp,level)     
        uinterp = funcs.linear_interp(Unew,totalp,level)
        vinterp = funcs.linear_interp(Vnew,totalp,level) 
        Ufinal = conv.ms_to_kts(uinterp) #convert to kts
        Vfinal = conv.ms_to_kts(vinterp)
        #speed = calc.calc_wspeed(Ufinal, Vfinal)
        gphgt = conv.gphgt_to_hgt(gphgt) # convert to height (m)
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
        title = level+'mb Height (m), Temp (C), Wind Barbs (kts)'
        ftitle = level+'mb-' 
        cblabel = 'kts'
        clevs = False
        cbticks = False
        makeplot(cs,title,cblabel,clevs,cbticks,ftitle)
        
def surface(): # plot surface chart. t2, wind barbs and mslp
    # create figure
    plt.figure(figsize=(8,8))
    x, y = m(lons, lats)
    t2c = conv.k_to_c(t2[time]) #convert temp to celcius
    if opt.tunit == 'F':
        t2f = conv.c_to_f(t2c) #convert celcius to fahrenheit
        clevs = np.arange(-30,115,5) # levels / degF
        cs = m.contourf(x,y,t2f,clevs,cmap=cm.get_cmap('gist_ncar'))
        cblabel = r'$\degree$F'
    elif opt.tunit == 'C':
        clevs = np.arange(-40,55,5) # levels / degC
        cs = m.contourf(x,y,t2c,clevs,cmap=cm.get_cmap('gist_ncar'))
        cblabel = r'$\degree$C'
    cbticks = True    
    
    psfchpa = conv.pa_to_hpa(psfc[time]) #convert Pa to hPa     
    mslp = calc.calc_mslp(psfchpa, thgt[0], t2[time]) # get mslp
    mslp = gaussian_filter(mslp, sigma=3) # smooth wiggles
    local_min, local_max = funcs.extrema(mslp, mode='wrap', window=50)
    #make x and y grid points for barbs
    #yy = np.arange(0, len(y), 8)
    #xx = np.arange(0, len(x), 8)
    #gp = np.meshgrid(yy, xx)

    #print x[::thin,::thin].shape #check x co-ord thinning
    #print u[time,::thin,::thin].shape #check u10 thinning
    #x_th,y_th = m(xlong[0,::thin,::thin],xlat[0,::thin,::thin]) #another method to thin barbs
        
    #convert wind to kts        
    u10kts = conv.ms_to_kts(u10[time])
    v10kts = conv.ms_to_kts(v10[time]) 
    m.barbs(x[::thin,::thin], y[::thin,::thin], u10kts[::thin,::thin], v10kts[::thin,::thin],length=opt.barbsize) #plot barbs
    title = "2m Temp, Wind Barbs (kts), MSLP (hPa)"
    ftitle = 'sfc-'       
    pclevs = np.arange(900,1055,2.)        
    pcs = m.contour(x,y,mslp,pclevs,colors='k',linewidths=2.)
    plt.clabel(pcs, inline=True, fmt='%1.0f', fontsize=12, colors='k')
    xlows = x[local_min];   xhighs = x[local_max]
    ylows = y[local_min];   yhighs = y[local_max]
    lowvals = mslp[local_min]; highvals = mslp[local_max]
    # plot lows as blue L's, with min pressure value underneath.
    xyplotted = []
    # don't plot if there is already a L or H within dmin meters.
    yoffset = 0.022*(m.ymax-m.ymin)
    dmin = yoffset
    for x,y,p in zip(xlows, ylows, lowvals):
        if x < m.xmax and x > m.xmin and y < m.ymax and y > m.ymin:
            dist = [np.sqrt((x-x0)**2+(y-y0)**2) for x0,y0 in xyplotted]
            if not dist or min(dist) > dmin:
                plt.text(x,y,'L',fontsize=14,fontweight='bold', ha='center',va='center',color='b')
                plt.text(x,y-yoffset,repr(int(p)),fontsize=12, ha='center',va='top',color='b', bbox = dict(boxstyle="square",ec='None',fc=(1,1,1,0.5)))
                xyplotted.append((x,y))
    # plot highs as red H's, with max pressure value underneath.
    xyplotted = []
    for x,y,p in zip(xhighs, yhighs, highvals):
        if x < m.xmax and x > m.xmin and y < m.ymax and y > m.ymin:
            dist = [np.sqrt((x-x0)**2+(y-y0)**2) for x0,y0 in xyplotted]
            if not dist or min(dist) > dmin:
                plt.text(x,y,'H',fontsize=14,fontweight='bold',
                        ha='center',va='center',color='r')
                plt.text(x,y-yoffset,repr(int(p)),fontsize=12, ha='center',va='top',color='r', bbox = dict(boxstyle="square",ec='None',fc=(1,1,1,0.5)))
                xyplotted.append((x,y))
                
    makeplot(cs,title,cblabel,clevs,cbticks,ftitle)
    
def snowaccum(): # plot snow accumulation
    # create figure
    plt.figure(figsize=(8,8))
    snow = nc.variables['SNOWNC'][time] # total accumulated grid scale snow and ice / mm at each time
    if opt.punit == 'mm':
        clevs = [0,0.5,1,2.5,3,4,5,8,10,15,20,30,40,50,80,100,150,200,250,500]
        cblabel = 'mm'
    elif opt.punit == 'in':
        snow = conv.mm_to_in(snow) # convert to inches
        clevs = [0.25,0.5,0.75,1,1.5,2,2.5,3,4,5,6,8,10,12,14,16,18,20,22,24]
        cblabel = 'inches'
    cbticks = True
    norm = matplotlib.colors.BoundaryNorm(clevs, 19) # set boundary of data by normalizing (0,1)    
    cs = m.contourf(x,y,snow,clevs,norm=norm,cmap=cmap.snow_colormap)
    title = "Snow Accumulation"
    ftitle = 'snow-'

    makeplot(cs,title,cblabel,clevs,cbticks,ftitle)
    
def hailaccum(): # plot hail accumulation
    # create figure
    plt.figure(figsize=(8,8))
    hail = nc.variables['HAILNC'][time] # accimulated total grid scale hail / mm at each time
    if opt.punit == 'mm':
        clevs = [0.5,1.,1.5,2.,2.5,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.]
        cblabel = 'mm'
    elif opt.punit == 'in':
        hail = conv.mm_to_in(hail) # convert to inches
        clevs = [0.01,0.02,0.04,0.06,0.08,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55]
        cblabel = 'inches'
    cbticks = True
    norm = matplotlib.colors.BoundaryNorm(clevs, 14) # set boundary of data by normalizing (0,1)    
    cs = m.contourf(x,y,hail,clevs,norm=norm,cmap=cmap.hail_colormap)
    title = "Hail Accumulation"
    ftitle = 'hail-'

    makeplot(cs,title,cblabel,clevs,cbticks,ftitle)
    
def simudbz(): # plot simulated reflectivity, mp_physics dependent
    # create figure
    plt.figure(figsize=(8,8))
    qrain = nc.variables['QRAIN'] # rain water mixing ratio
    t2c = conv.k_to_c(t2[time]) #convert temp to celcius
    rhoa = calc.calc_rhoa(psfc[time], t2[time])
    Qrain = qrain[time,1] # rain mixing ratio
    Qrain = np.nan_to_num(Qrain) # change NaN to zeroes, changge infs to nums
    try: #depends on MP scheme
        Qsn = nc.variables['QSNOW'] # try to get snow mixing ratio
    except:
        Qsn = np.zeros(np.shape(qrain)) # else create zeros array same shape as qrain
    Qsnow = Qsn[time,1] # snow mixing ratio
    Qsnow = np.nan_to_num(Qsnow) # change NaN to zeros
    dBZ = calc.calc_dbz(t2c, rhoa, Qrain, Qsnow)
    clevs = np.arange(0,85,5)
    norm = matplotlib.colors.BoundaryNorm(clevs, 17) # normalize levels
    cs = m.contourf(x,y,dBZ,clevs,norm=norm,cmap=cmap.dbz_colormap)
    title = "Simulated Reflectivity"
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
        Qrainall = nc.variables['QRAIN'][time] # rain water mixing ratio at all levels
        t2c = conv.k_to_c(t2[time]) #convert temp to celcius
        rhoa = calc.calc_rhoa(psfc[time], t2[time])
        try: # depends on MP scheme
            Qsn = nc.variables['QSNOW'] # try to get snow mixing ratio
        except:
            Qsn = np.zeros(np.shape(Qrainall)) # else create zeros array same shape as qrain
        Qsnowall = Qsn[time] # get all Qsnow values at all levels for each time
        Qrainmax = np.max(Qrainall, axis=0) #max rain QV
        Qsnowmax = np.max(Qsnowall, axis=0) #max snow QV
        dBZ = calc.calc_dbz(t2c, rhoa, Qrainmax, Qsnowmax)
    clevs = np.arange(0,85,5)
    norm = matplotlib.colors.BoundaryNorm(clevs, 17) # normalize levels
    cs = m.contourf(x,y,dBZ,clevs,norm=norm,cmap=cmap.dbz_colormap)
    title = "Composite Reflectivity"
    ftitle = 'compdbz-'
    cblabel = 'dBZ'
    cbticks = True
    makeplot(cs,title,cblabel,clevs,cbticks,ftitle)
    
def lclhgt(): # plot lcl height
    # create figure
    plt.figure(figsize=(8,8))
    q2 = nc.variables['Q2'][time] # water vapour mixing ratio at 2m
    t2c = conv.k_to_c(t2[time]) #convert temp to celcius
    psfchpa = conv.pa_to_hpa(psfc[time])
    es = calc.calc_es(t2c)
    ws = calc.calc_ws(es, psfchpa)
    rh = calc.calc_rh(q2, ws)
    td = calc.calc_dewpoint(es, rh)
    lcl = calc.calc_lcl(t2c, td)
    clevs = np.arange(0,6000,500)
    cs = m.contourf(x,y,lcl,clevs,cmap=cmap.lcl_colormap)
    title = "LCL Height"
    ftitle = 'lcl-'
    cblabel = 'm'
    cbticks = True
    makeplot(cs,title,cblabel,clevs,cbticks,ftitle)
    
def thetaE(): # plot theta-e
    # create figure
    plt.figure(figsize=(8,8))
    theta = nc.variables['T'][time] #perturbation potential temperature (theta-t0)
    theta0 = nc.variables['T00'][0] #base state theta
    theta = theta[0] + theta0 # total theta
    psfchpa = conv.pa_to_hpa(psfc[time])
    t2c = conv.k_to_c(t2[time]) #convert temp to celcius
    es = calc.calc_es(t2c)
    ws = calc.calc_ws(es, psfchpa)     
    thetae = calc.calc_thetae(theta, t2[time], ws)
    clevs = np.arange(260,372,4) # set by max and min of data    
    cs = m.contourf(x,y,thetae,clevs,cmap=cm.get_cmap('gist_ncar'))
    title = "Theta-e"
    ftitle = 'thetae-'
    cblabel = 'K'
    cbticks = True
    makeplot(cs,title,cblabel,clevs,cbticks,ftitle)
    
def h75lr(): # 700-500mb lapse rates
    # create figure
    plt.figure(figsize=(8,8))
    pb = nc.variables['PB'][time] #base state pressure, Pa
    p = nc.variables['P'][time] # perturbation pressure, Pa
    totalp = pb + p # total pressure in Pa
    theta = nc.variables['T'][time] #perturbation potential temperature (theta-t0)
    theta0 = nc.variables['T00'][0] #base state theta
    totalTheta = theta + theta0 # total potential temp
    totalT= conv.k_to_c(calc.theta_to_temp(totalTheta, totalp)) # calc temp in deg C
    # interp temps to levels
    totalT700 = funcs.linear_interp(totalT,totalp,700)
    totalT500 = funcs.linear_interp(totalT,totalp,500)
    # calc h7-h5 lapse rates
    lr = totalT700 - totalT500
    clevs = np.arange(5,10.5,.5) # conditionally unstable levels
    cs = m.contourf(x,y,lr,clevs,cmap=cm.get_cmap('gist_ncar'))
    title = "H7-H5 Lapse Rates"
    ftitle = 'h75lr-'
    cblabel = r'$\degree$C'
    cbticks = True
    makeplot(cs,title,cblabel,clevs,cbticks,ftitle)
    
def absvort500(): # plot 500mb absolute vorticity
    # create figure
    plt.figure(figsize=(8,8))
    pb = nc.variables['PB'][time] #base state pressure, Pa
    p = nc.variables['P'][time] # perturbation pressure, Pa
    totalp = pb + p # total pressure in Pa
    U = funcs.unstagger(nc.variables['U'][time],'U') # U wind component UNSTAGGERED
    V = funcs.unstagger(nc.variables['V'][time],'V') # V wind component
    fcoriolis = calc.calc_fcoriolis(xlat[0])
    uinterp = funcs.linear_interp(U,totalp,500) #interp to 500mb
    vinterp = funcs.linear_interp(V,totalp,500) 
    vertvort = calc.calc_vertvort(uinterp, vinterp, dx)
    avort = vertvort + fcoriolis # absolute vorticity 
    avort = np.multiply(avort, 1e5) # scale up for levels
    clevs = np.arange(-6, 52, 2)    
    cs = m.contourf(x,y,avort,clevs,cmap=cm.get_cmap('gist_ncar'))
    title = '500mb Absolute Vorticity'
    ftitle = '500absvort-' 
    cblabel = r'$10^{-5} s^{-1}$'
    cbticks = True
    makeplot(cs,title,cblabel,clevs,cbticks,ftitle)
    
def shr06(): # plot the 0-6km shear vector
    # create figure
    plt.figure(figsize=(8,8))
    ph = nc.variables['PH'][time] #perturbation geopotential
    phb = nc.variables['PHB'][time] #base state geopotential
    totalgp = phb + ph # total geopotential
    totalgp = funcs.unstagger(totalgp,'Z') #total geopotential unstaggered
    U = funcs.unstagger(nc.variables['U'][time],'U') # U wind component # UNSTAGGERED
    V = funcs.unstagger(nc.variables['V'][time],'V') # V wind component
    u10kts = conv.ms_to_kts(u10[time]) # sfc wind in kts
    v10kts = conv.ms_to_kts(v10[time]) 
    u6 = funcs.interp_generic(6000, (totalgp/9.81), U) # interp to 6km
    v6 = funcs.interp_generic(6000, (totalgp/9.81), V)
    u6kts = conv.ms_to_kts(u6) # convert 6km wind to kts
    v6kts = conv.ms_to_kts(v6)
    #using 10m wind as sfc wind
    ushr = u6kts - u10kts # calc 0-6 shr in kts
    vshr = v6kts - v10kts
    speed = calc.calc_wspeed(ushr, vshr)
    # plot data
    clevs = np.arange(20,145,5)
    cs = m.contourf(x, y, speed, clevs, cmap=cm.get_cmap('gist_ncar'))
    m.barbs(x[::thin,::thin], y[::thin,::thin], ushr[::thin,::thin], vshr[::thin,::thin],length=opt.barbsize) #plot barbs
    title = '0-6km Shear'
    ftitle = 'shr06-' 
    cblabel = 'kts'
    cbticks = True
    makeplot(cs,title,cblabel,clevs,cbticks,ftitle)
    
def vertvol(): # plot the vertical velocity at levels. NEEDS CORRECTING TO VERTICAL MOTION OMEGA EQUATION
    W = funcs.unstagger(nc.variables['W'][time],'W') # unstaggered vertical velocity
    pb = nc.variables['PB'][time] #base state pressure, Pa
    p = nc.variables['P'][time] # perturbation pressure, Pa
    totalp = pb + p # total pressure in Pa
    
    levels = opt.lvl.split(',') # get list of levels
    for level in levels: 
        plt.figure(figsize=(8,8)) #create fig for each plot
        level = int(level) # make it int
        Wfinal = funcs.linear_interp(W,totalp,level) # interpolate W to levels
        clevs = np.arange(-2.0,2.2,0.2)        
        cs = m.contourf(x,y,Wfinal,clevs,cmap=cm.get_cmap('gist_ncar'))
        level = str(level)
        title = level+'mb Vertical Velocity'
        ftitle = level+'mbvv-' 
        cblabel = r'$ms^{-1}$'
        cbticks = True
        makeplot(cs,title,cblabel,clevs,cbticks,ftitle)
        
def olr_to_temp(): # convert OLR to IR temp
    plt.figure(figsize=(8,8))
    olr = nc.variables['OLR'][time]
    olrtemp = np.power(olr / 5.67e-8, 0.25) - 273.15 # calc temp using Stefan-Boltzman law and convert to deg C
    clevs = np.arange(-80, 36 ,4)
    cs = m.contourf(x,y,olrtemp,clevs,cmap=cmap.irsat_colormap)
    title = 'IR Brightness Temp'
    ftitle = 'irtemp-' 
    cblabel = r'$\degree$C'
    cbticks = True
    makeplot(cs,title,cblabel,clevs,cbticks,ftitle)
    
def pymeteo_skewt(): # uses pyMeteo package (https://github.com/cwebster2/pyMeteo) to plot skew-t for lat/lon. Credit Casey Webster
    import pymeteo.skewt as skewt
    try:
        skewt.plot_wrf(filein,opt.slat,opt.slon,time,'skewt'+str(time)+'.png')
    except:
        print "LAT/LON NOT IN DOMAIN. QUITTING"
        sys.exit()
        
def plot_skewt(): # plot skew-t by writing data to file and use SHARPpy available at: https://github.com/sharppy/SHARPpy
    i, j = funcs.latlon_ij(opt.slat, opt.slon, xlat, xlong)
    inlat = xlat[0,i,j]
    inlon = xlong[0,i,j]
    pb = nc.variables['PB'][time,:,i,j] #base state pressure, Pa
    p = nc.variables['P'][time,:,i,j] # perturbation pressure, Pa
    totalp = p + pb # total pressure
    ph = nc.variables['PH'][time,:,i,j] #perturbation geopotential
    phb = nc.variables['PHB'][time,:,i,j] #base state geopotential
    totalgp = phb + ph # total geopotential
    totalgp = funcs.unstagger(totalgp,'Z') #total geopotential unstaggered
    U = nc.variables['U'][time,:,i,j] # U wind component
    V = nc.variables['V'][time,:,i,j] # V wind component
    theta = nc.variables['T'][time,:,i,j] #perturbation potential temperature (theta-t0)
    theta0 = nc.variables['T00'][0] #base state theta
    totaltheta = theta+theta0 # total potential temp
    qvapor = nc.variables['QVAPOR'][time,:,i,j] #water vapor mixing ratio kg/kg
    #need to calc these variables for skewt
    level = conv.pa_to_hpa(totalp) # levels in hPa
    height = conv.gphgt_to_hgt(totalgp) # heights in m
    temps = calc.theta_to_temp(totaltheta, totalp) # temps in degK
    tempc = conv.k_to_c(temps) # temps in degC
    es = calc.calc_es(tempc) # calc es
    ws = calc.calc_ws(es, level) # calc ws
    rh = calc.calc_rh(qvapor, ws) # calc rh
    dwpt = calc.calc_dewpoint(es, rh) # calc dewpoint in degC
    winddir = calc.calc_wdir(U, V) # calc wind dir
    wspd = conv.ms_to_kts(calc.calc_wspeed(U, V)) # calc wind spd
    skewt_data = funcs.skewt_data(timestamp, level, height, tempc, dwpt, winddir, wspd, inlat, inlon) # write the data to SPC file format
    pltfuncs.do_sharppy(skewt_data) # use SHARPpy to plot skew-t
    
def updraft_hel(): # plot the 2-5km updraft helicity
    plt.figure(figsize=(8,8))
    U = funcs.unstagger(nc.variables['U'][time],'U') # U wind component # UNSTAGGERED
    V = funcs.unstagger(nc.variables['V'][time],'V') # V wind component
    W = funcs.unstagger(nc.variables['W'][time],'W') # unstaggered vertical velocity
    ph = nc.variables['PH'][time] #perturbation geopotential
    phb = nc.variables['PHB'][time] #base state geopotential
    totalgp = phb + ph # total geopotential
    totalgp = funcs.unstagger(totalgp,'Z') #total geopotential unstaggered
    heights = totalgp / 9.81
    levels = 6 # no of levels in between bottom and top of a layer (add extra one to get to very top of layer)
    depth = 1000 # depth of layer
    dz = depth / (levels-1) # increment / m
    #create arrays to hold all the values at each level
    u2km = np.zeros((levels, np.shape(U)[1], np.shape(U)[2]))
    v2km = np.zeros((levels, np.shape(V)[1], np.shape(V)[2]))
    u3km = np.zeros((levels, np.shape(U)[1], np.shape(U)[2]))
    v3km = np.zeros((levels, np.shape(V)[1], np.shape(V)[2]))
    u4km = np.zeros((levels, np.shape(U)[1], np.shape(U)[2]))
    v4km = np.zeros((levels, np.shape(V)[1], np.shape(V)[2]))
    #u5km = np.zeros((levels, np.shape(U)[1], np.shape(U)[2]))
    #v5km = np.zeros((levels, np.shape(V)[1], np.shape(V)[2]))
    w2km = np.zeros((levels, np.shape(W)[1], np.shape(W)[2]))
    w3km = np.zeros((levels, np.shape(W)[1], np.shape(W)[2]))
    w4km = np.zeros((levels, np.shape(W)[1], np.shape(W)[2]))
    #w5km = np.zeros((levels, np.shape(W)[1], np.shape(W)[2]))
    zeta2km = np.zeros((levels, np.shape(U)[1], np.shape(U)[2]))
    zeta3km = np.zeros((levels, np.shape(U)[1], np.shape(U)[2]))
    zeta4km = np.zeros((levels, np.shape(U)[1], np.shape(U)[2]))
    #zeta5km = np.zeros((levels, np.shape(U)[1], np.shape(U)[2]))
    for i in range(0,levels): # loop through to interpolate to levels and store in array
        print "Interpolating...doing loop ", i, "of ", (levels-1)
        increment = i*dz
        u2km[i] = funcs.interp_generic(2000+increment, heights, U)
        v2km[i] = funcs.interp_generic(2000+increment, heights, V)
        u3km[i] = funcs.interp_generic(3000+increment, heights, U)
        v3km[i] = funcs.interp_generic(3000+increment, heights, V)
        u4km[i] = funcs.interp_generic(4000+increment, heights, U)
        v4km[i] = funcs.interp_generic(4000+increment, heights, V)
        #u5km[i] = funcs.interp_generic(5000+increment, heights, U)
        #v5km[i] = funcs.interp_generic(5000+increment, heights, V)
        w2km[i] = funcs.interp_generic(2000+increment, heights, W)
        w3km[i] = funcs.interp_generic(2000+increment, heights, W)
        w4km[i] = funcs.interp_generic(2000+increment, heights, W)
        #w5km[i] = funcs.interp_generic(2000+increment, heights, W)
        zeta2km[i] = calc.calc_vertvort(u2km[i], v2km[i], dx)
        zeta3km[i] = calc.calc_vertvort(u3km[i], v3km[i], dx)
        zeta4km[i] = calc.calc_vertvort(u4km[i], v4km[i], dx)
        #zeta5km[i] = calc.calc_vertvort(u5km[i], v5km[i], dx)
    # calc the layer mean
    w2to3 = np.mean(w2km, axis=0)
    w3to4 = np.mean(w3km, axis=0)
    w4to5 = np.mean(w4km, axis=0)
    zeta2to3 = np.mean(zeta2km, axis=0)
    zeta3to4 = np.mean(zeta3km, axis=0)
    zeta4to5 = np.mean(zeta4km, axis=0)
    # calc the 2-5km UH
    UH = ( w2to3*zeta2to3 + w3to4*zeta3to4 + w4to5*zeta4to5 ) * 1000
    UH = funcs.nine_point_smooth(UH)
    #u2km = funcs.interp_generic(2000, heights, U)
    #v2km = funcs.interp_generic(2000, heights, V)
    #u3km = funcs.interp_generic(3000, heights, U)
    #v3km = funcs.interp_generic(3000, heights, V)
    #u4km = funcs.interp_generic(4000, heights, U)
    #v4km = funcs.interp_generic(4000, heights, V)
    #u5km = funcs.interp_generic(5000, heights, U)
    #v5km = funcs.interp_generic(5000, heights, V)
    #w2km = funcs.interp_generic(2000, heights, W)
    #w3km = funcs.interp_generic(2000, heights, W)
    #w4km = funcs.interp_generic(2000, heights, W)
    #w5km = funcs.interp_generic(2000, heights, W)
    #w2to3 = 0.5 * ( w2km + w3km )
    #w3to4 = 0.5 * ( w3km + w4km )
    #w4to5 = 0.5 * ( w4km + w5km )
    #zeta2km = calc.calc_vertvort(u2km, v2km, dx)
    #zeta3km = calc.calc_vertvort(u3km, v3km, dx)
    #zeta4km = calc.calc_vertvort(u4km, v4km, dx)
    #zeta5km = calc.calc_vertvort(u5km, v5km, dx)
    #zeta2to3 = 0.5 * ( zeta2km + zeta3km )
    #zeta3to4 = 0.5 * ( zeta3km + zeta4km )
    #zeta4to5 = 0.5 * ( zeta4km + zeta5km )
    #UH = ( w2to3*zeta2to3 + w3to4*zeta3to4 + w4to5*zeta4to5 ) * 1000
    clevs = np.arange(0,210,10) 
    cs = m.contourf(x,y,UH,cmap=cmap.uh_colormap)
    title = '2-5km Updraft Helicity'
    ftitle = 'uh-' 
    cblabel = r'$m^{2}s^{-2}$'
    cbticks = True
    makeplot(cs,title,cblabel,clevs,cbticks,ftitle)
    
### END PLOT FUNCTIONS ###
flag = False # to check for plotting options
#### BEGIN TIME LOOP ####    
for time in range(times.shape[0]):
    currtime = str(''.join(times[time])).replace('_', ' ') #get current model time    
    filetime = currtime.translate(None, ':').replace(' ', '_') # time for filename
    alltimes.append(currtime) # all times in output 
    timestamp = currtime[8:10]+currtime[5:7]+currtime[2:4]+'/'+currtime[11:13]+currtime[14:16]
    
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
        
    if opt.shear06:
        print "Plotting 0-6km Shear for time: ", currtime
        shr06()
        flag = True
        
    if opt.vv:
        print "Plotting vertical velocity for time: ", currtime
        vertvol()
        flag = True
        
    if opt.irtemp:
        print "Plotting IR Brightness Temp for time: ", currtime
        olr_to_temp()
        flag = True
        
    if opt.skewt:
        print "Plotting Skew-t for time: ", currtime
        pymeteo_skewt()
        flag = True
        
    if opt.getij:
        print "Getting i, j for lat=",opt.slat, ', lon=',opt.slon
        funcs.latlon_ij(opt.slat, opt.slon, xlat, xlong)
        #print "A less accurate method:"
        #funcs.latlon_ij2(opt.slat, opt.slon, xlat, xlong)
        flag = True
        sys.exit()
        
    if opt.skewt2:
        print "Plotting Skew-t for time: ", currtime
        plot_skewt()
        flag = True
        
    if opt.uh25:
        print "Plotting 2-5km Updraft Helicity for time: ", currtime
        updraft_hel()
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
    print "Map projection: ", proj, '-' , projname

nc.close() # close netcdf file