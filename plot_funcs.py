# -*- coding: utf-8 -*-
"""
Liam Till
University of Oklahoma / Unversity of Reading
plot_funcs.py

Functions for plotting for wrfplot package
Last modified: 23/04/16
"""

import numpy as np

# do the SHARPpy sounding plot    
def do_sharppy(spc_file):
    """
    Based on the tutorial which can be found here: http://nbviewer.ipython.org/github/sharppy/SHARPpy/blob/master/tutorials/SHARPpy_basics.ipynb
    SHARPpy can be found here: https://github.com/sharppy/SHARPpy
    Credit goes to:
    Patrick Marsh (SPC)
    Kelton Halbert (OU School of Meteorology)
    Greg Blumberg (OU/CIMMS)
    Tim Supinie (OU School of Meteorology)
    
    """
    import sharppy
    import sharppy.sharptab.profile as profile
    import sharppy.sharptab.interp as interp
    import sharppy.sharptab.winds as winds
    import sharppy.sharptab.utils as utils
    import sharppy.sharptab.params as params
    import sharppy.sharptab.thermo as thermo
    import matplotlib.pyplot as plt
    from StringIO import StringIO
    from matplotlib.axes import Axes
    import matplotlib.transforms as transforms
    import matplotlib.axis as maxis
    import matplotlib.spines as mspines
    import matplotlib.path as mpath
    from matplotlib.projections import register_projection
    
    spc_file = open('skewt_data', 'r').read()


    def parseSPC(spc_file):
        ## read in the file
        data = np.array([l.strip() for l in spc_file.split('\n')])

        ## necessary index points
        title_idx = np.where( data == '%TITLE%')[0][0]
        start_idx = np.where( data == '%RAW%' )[0] + 1
        finish_idx = np.where( data == '%END%')[0]
    
        ## create the plot title
        data_header = data[title_idx + 1].split()
        location = data_header[0]+' '+data_header[1]
        time = data_header[2]
        title = location+' '+time
        ## put it all together for StringIO
        full_data = '\n'.join(data[start_idx : finish_idx][:])
        sound_data = StringIO( full_data )
    
        ## read the data into arrays
        p, h, T, Td, wdir, wspd = np.genfromtxt( sound_data, delimiter=',', comments="%", unpack=True )
    
        return p, h, T, Td, wdir, wspd, title
        
    pres, hght, tmpc, dwpc, wdir, wspd, title = parseSPC(spc_file)
    prof = profile.create_profile(profile='default', pres=pres, hght=hght, tmpc=tmpc, \
    dwpc=dwpc, wspd=wspd, wdir=wdir, missing=-9999, strictQC=True)
    
    sfcpcl = params.parcelx( prof, flag=1 ) # Surface Parcel
    fcstpcl = params.parcelx( prof, flag=2 ) # Forecast Parcel
    mupcl = params.parcelx( prof, flag=3 ) # Most-Unstable Parcel
    mlpcl = params.parcelx( prof, flag=4 ) # 100 mb Mean Layer Parcel
         
    msl_hght = prof.hght[prof.sfc] # Grab the surface height value
    print "SURFACE HEIGHT (m MSL):",msl_hght
    agl_hght = interp.to_agl(prof, msl_hght) # Converts to AGL
    print "SURFACE HEIGHT (m AGL):", agl_hght
    msl_hght = interp.to_msl(prof, agl_hght) # Converts to MSL
    print "SURFACE HEIGHT (m MSL):",msl_hght
    print "Most-Unstable CAPE:", mupcl.bplus # J/kg
    print "Most-Unstable CIN:", mupcl.bminus # J/kg
    print "Most-Unstable LCL:", mupcl.lclhght # meters AGL
    print "Most-Unstable LFC:", mupcl.lfchght # meters AGL
    print "Most-Unstable EL:", mupcl.elhght # meters AGL
    print "Most-Unstable LI:", mupcl.li5 # C
    
    class SkewXTick(maxis.XTick):
        def draw(self, renderer):
            if not self.get_visible(): return
            renderer.open_group(self.__name__)
    
            lower_interval = self.axes.xaxis.lower_interval
            upper_interval = self.axes.xaxis.upper_interval
    
            if self.gridOn and transforms.interval_contains(
                    self.axes.xaxis.get_view_interval(), self.get_loc()):
                self.gridline.draw(renderer)
    
            if transforms.interval_contains(lower_interval, self.get_loc()):
                if self.tick1On:
                    self.tick1line.draw(renderer)
                if self.label1On:
                    self.label1.draw(renderer)
    
            if transforms.interval_contains(upper_interval, self.get_loc()):
                if self.tick2On:
                    self.tick2line.draw(renderer)
                if self.label2On:
                    self.label2.draw(renderer)
    
            renderer.close_group(self.__name__)
    
    
    # This class exists to provide two separate sets of intervals to the tick,
    # as well as create instances of the custom tick
    class SkewXAxis(maxis.XAxis):
        def __init__(self, *args, **kwargs):
            maxis.XAxis.__init__(self, *args, **kwargs)
            self.upper_interval = 0.0, 1.0
    
        def _get_tick(self, major):
            return SkewXTick(self.axes, 0, '', major=major)
    
        @property
        def lower_interval(self):
            return self.axes.viewLim.intervalx
    
        def get_view_interval(self):
            return self.upper_interval[0], self.axes.viewLim.intervalx[1]
    
    
    # This class exists to calculate the separate data range of the
    # upper X-axis and draw the spine there. It also provides this range
    # to the X-axis artist for ticking and gridlines
    class SkewSpine(mspines.Spine):
        def _adjust_location(self):
            trans = self.axes.transDataToAxes.inverted()
            if self.spine_type == 'top':
                yloc = 1.0
            else:
                yloc = 0.0
            left = trans.transform_point((0.0, yloc))[0]
            right = trans.transform_point((1.0, yloc))[0]
    
            pts  = self._path.vertices
            pts[0, 0] = left
            pts[1, 0] = right
            self.axis.upper_interval = (left, right)
    
    
    # This class handles registration of the skew-xaxes as a projection as well
    # as setting up the appropriate transformations. It also overrides standard
    # spines and axes instances as appropriate.
    class SkewXAxes(Axes):
        # The projection must specify a name.  This will be used be the
        # user to select the projection, i.e. ``subplot(111,
        # projection='skewx')``.
        name = 'skewx'
    
        def _init_axis(self):
            #Taken from Axes and modified to use our modified X-axis
            self.xaxis = SkewXAxis(self)
            self.spines['top'].register_axis(self.xaxis)
            self.spines['bottom'].register_axis(self.xaxis)
            self.yaxis = maxis.YAxis(self)
            self.spines['left'].register_axis(self.yaxis)
            self.spines['right'].register_axis(self.yaxis)
    
        def _gen_axes_spines(self):
            spines = {'top':SkewSpine.linear_spine(self, 'top'),
                      'bottom':mspines.Spine.linear_spine(self, 'bottom'),
                      'left':mspines.Spine.linear_spine(self, 'left'),
                      'right':mspines.Spine.linear_spine(self, 'right')}
            return spines
    
        def _set_lim_and_transforms(self):
            """
            This is called once when the plot is created to set up all the
            transforms for the data, text and grids.
            """
            rot = 30
    
            #Get the standard transform setup from the Axes base class
            Axes._set_lim_and_transforms(self)
    
            # Need to put the skew in the middle, after the scale and limits,
            # but before the transAxes. This way, the skew is done in Axes
            # coordinates thus performing the transform around the proper origin
            # We keep the pre-transAxes transform around for other users, like the
            # spines for finding bounds
            self.transDataToAxes = self.transScale + (self.transLimits +
                    transforms.Affine2D().skew_deg(rot, 0))
    
            # Create the full transform from Data to Pixels
            self.transData = self.transDataToAxes + self.transAxes
    
            # Blended transforms like this need to have the skewing applied using
            # both axes, in axes coords like before.
            self._xaxis_transform = (transforms.blended_transform_factory(
                        self.transScale + self.transLimits,
                        transforms.IdentityTransform()) +
                    transforms.Affine2D().skew_deg(rot, 0)) + self.transAxes
    
    # Now register the projection with matplotlib so the user can select
    # it.
    register_projection(SkewXAxes)
    
    pcl = mupcl
    # Create a new figure. The dimensions here give a good aspect ratio
    fig = plt.figure(figsize=(6.5875, 6.2125))
    ax = fig.add_subplot(111, projection='skewx')
    ax.grid(True)
    
    pmax = 1000
    pmin = 10
    dp = -10
    presvals = np.arange(int(pmax), int(pmin)+dp, dp)
    
    # plot the moist-adiabats
    for t in np.arange(-10,45,5):
        tw = []
        for p in presvals:
            tw.append(thermo.wetlift(1000., t, p))
        ax.semilogy(tw, presvals, 'k-', alpha=.2)
    
    def thetas(theta, presvals):
        return ((theta + thermo.ZEROCNK) / (np.power((1000. / presvals),thermo.ROCP))) - thermo.ZEROCNK
    
    # plot the dry adiabats
    for t in np.arange(-50,110,10):
        ax.semilogy(thetas(t, presvals), presvals, 'r-', alpha=.2)
    
    plt.title(title, fontsize=14, loc='left')
    # Plot the data using normal plotting functions, in this case using
    # log scaling in Y, as dicatated by the typical meteorological plot
    ax.semilogy(prof.tmpc, prof.pres, 'r', lw=2)
    ax.semilogy(prof.dwpc, prof.pres, 'g', lw=2)
    ax.semilogy(pcl.ttrace, pcl.ptrace, 'k-.', lw=2)
    
    # An example of a slanted line at constant X
    l = ax.axvline(0, color='b', linestyle='--')
    l = ax.axvline(-20, color='b', linestyle='--')
    
    # Disables the log-formatting that comes with semilogy
    ax.yaxis.set_major_formatter(plt.ScalarFormatter())
    ax.set_yticks(np.linspace(100,1000,10))
    ax.set_ylim(1050,100)
    
    ax.xaxis.set_major_locator(plt.MultipleLocator(10))
    ax.set_xlim(-50,50)
    plt.show()
    
    ##PLOTS SKEWT OK ABOVE HERE ##
    """
    sfc = prof.pres[prof.sfc]
    p3km = interp.pres(prof, interp.to_msl(prof, 3000.))
    p6km = interp.pres(prof, interp.to_msl(prof, 6000.))
    p1km = interp.pres(prof, interp.to_msl(prof, 1000.))
    mean_3km = winds.mean_wind(prof, pbot=sfc, ptop=p3km)
    sfc_6km_shear = winds.wind_shear(prof, pbot=sfc, ptop=p6km)
    sfc_3km_shear = winds.wind_shear(prof, pbot=sfc, ptop=p3km)
    sfc_1km_shear = winds.wind_shear(prof, pbot=sfc, ptop=p1km)
    print "0-3 km Pressure-Weighted Mean Wind (kt):", utils.comp2vec(mean_3km[0], mean_3km[1])[1]
    print "0-6 km Shear (kt):", utils.comp2vec(sfc_6km_shear[0], sfc_6km_shear[1])[1]
    srwind = params.bunkers_storm_motion(prof)
    print "Bunker's Storm Motion (right-mover) [deg,kts]:", utils.comp2vec(srwind[0], srwind[1])
    print "Bunker's Storm Motion (left-mover) [deg,kts]:", utils.comp2vec(srwind[2], srwind[3])
    srh3km = winds.helicity(prof, 0, 3000., stu = srwind[0], stv = srwind[1])
    srh1km = winds.helicity(prof, 0, 1000., stu = srwind[0], stv = srwind[1])
    print "0-3 km Storm Relative Helicity [m2/s2]:",srh3km[0]
    
    stp_fixed = params.stp_fixed(sfcpcl.bplus, sfcpcl.lclhght, srh1km[0], utils.comp2vec(sfc_6km_shear[0], sfc_6km_shear[1])[1])
    ship = params.ship(prof)
    eff_inflow = params.effective_inflow_layer(prof)
    ebot_hght = interp.to_agl(prof, interp.hght(prof, eff_inflow[0]))
    etop_hght = interp.to_agl(prof, interp.hght(prof, eff_inflow[1]))
    print "Effective Inflow Layer Bottom Height (m AGL):", ebot_hght
    print "Effective Inflow Layer Top Height (m AGL):", etop_hght
    effective_srh = winds.helicity(prof, ebot_hght, etop_hght, stu = srwind[0], stv = srwind[1])
    print "Effective Inflow Layer SRH (m2/s2):", effective_srh[0]
    ebwd = winds.wind_shear(prof, pbot=eff_inflow[0], ptop=eff_inflow[1])
    ebwspd = utils.mag( ebwd[0], ebwd[1] )
    print "Effective Bulk Wind Difference:", ebwspd
    scp = params.scp(mupcl.bplus, effective_srh[0], ebwspd)
    stp_cin = params.stp_cin(mlpcl.bplus, effective_srh[0], ebwspd, mlpcl.lclhght, mlpcl.bminus)
    print "Supercell Composite Parameter:", scp
    print "Significant Tornado Parameter (w/CIN):", stp_cin
    print "Significant Tornado Parameter (fixed):", stp_fixed
    
    indices = {'SBCAPE': [int(sfcpcl.bplus), 'J/kg'],\
           'SBCIN': [int(sfcpcl.bminus), 'J/kg'],\
           'SBLCL': [int(sfcpcl.lclhght), 'm AGL'],\
           'SBLFC': [int(sfcpcl.lfchght), 'm AGL'],\
           'SBEL': [int(sfcpcl.elhght), 'm AGL'],\
           'SBLI': [int(sfcpcl.li5), 'C'],\
           'MLCAPE': [int(mlpcl.bplus), 'J/kg'],\
           'MLCIN': [int(mlpcl.bminus), 'J/kg'],\
           'MLLCL': [int(mlpcl.lclhght), 'm AGL'],\
           'MLLFC': [int(mlpcl.lfchght), 'm AGL'],\
           'MLEL': [int(mlpcl.elhght), 'm AGL'],\
           'MLLI': [int(mlpcl.li5), 'C'],\
           'MUCAPE': [int(mupcl.bplus), 'J/kg'],\
           'MUCIN': [int(mupcl.bminus), 'J/kg'],\
           'MULCL': [int(mupcl.lclhght), 'm AGL'],\
           'MULFC': [int(mupcl.lfchght), 'm AGL'],\
           'MUEL': [int(mupcl.elhght), 'm AGL'],\
           'MULI': [int(mupcl.li5), 'C'],\
           '0-1 km SRH': [int(srh1km[0]), 'm2/s2'],\
           '0-1 km Shear': [int(utils.comp2vec(sfc_1km_shear[0], sfc_1km_shear[1])[1]), 'kts'],\
           '0-3 km SRH': [int(srh3km[0]), 'm2/s2'],\
           'Eff. SRH': [int(effective_srh[0]), 'm2/s2'],\
           'EBWD': [int(ebwspd), 'kts'],\
           'PWV': [round(params.precip_water(prof), 2), 'inch'],\
           'K-index': [int(params.k_index(prof)), ''],\
           'STP(fix)': [round(stp_fixed, 1), ''],\
           'SHIP': [round(ship, 1), ''],\
           'SCP': [round(scp, 1), ''],\
           'STP(cin)': [round(stp_cin, 1), '']}
    
    # Set the parcel trace to be plotted as the Most-Unstable parcel.
    pcl = mupcl
    
    # Create a new figure. The dimensions here give a good aspect ratio
    fig = plt.figure(figsize=(6.5875, 6.2125))
    ax = fig.add_subplot(111, projection='skewx')
    ax.grid(True)
    
    pmax = 1000
    pmin = 10
    dp = -10
    presvals = np.arange(int(pmax), int(pmin)+dp, dp)
    
    # plot the moist-adiabats
    for t in np.arange(-10,45,5):
        tw = []
        for p in presvals:
            tw.append(thermo.wetlift(1000., t, p))
        ax.semilogy(tw, presvals, 'k-', alpha=.2)
    
    def thetas(theta, presvals):
        return ((theta + thermo.ZEROCNK) / (np.power((1000. / presvals),thermo.ROCP))) - thermo.ZEROCNK
    
    # plot the dry adiabats
    for t in np.arange(-50,110,10):
        ax.semilogy(thetas(t, presvals), presvals, 'r-', alpha=.2)
    
    plt.title(' OAX 140616/1900 (Observed)', fontsize=12, loc='left')
    # Plot the data using normal plotting functions, in this case using
    # log scaling in Y, as dicatated by the typical meteorological plot
    ax.semilogy(prof.tmpc, prof.pres, 'r', lw=2) # Plot the temperature profile
    ax.semilogy(prof.wetbulb, prof.pres, 'c-') # Plot the wetbulb profile
    ax.semilogy(prof.dwpc, prof.pres, 'g', lw=2) # plot the dewpoint profile
    ax.semilogy(pcl.ttrace, pcl.ptrace, 'k-.', lw=2) # plot the parcel trace 
    # An example of a slanted line at constant X
    l = ax.axvline(0, color='b', linestyle='--')
    l = ax.axvline(-20, color='b', linestyle='--')
    
    # Plot the effective inflow layer using blue horizontal lines
    ax.axhline(eff_inflow[0], color='b')
    ax.axhline(eff_inflow[1], color='b')
    
    #plt.barbs(10*np.ones(len(prof.pres)), prof.pres, prof.u, prof.v)
    # Disables the log-formatting that comes with semilogy
    ax.yaxis.set_major_formatter(plt.ScalarFormatter())
    ax.set_yticks(np.linspace(100,1000,10))
    ax.set_ylim(1050,100)
    ax.xaxis.set_major_locator(plt.MultipleLocator(10))
    ax.set_xlim(-50,50)
    
    # List the indices within the indices dictionary on the side of the plot.
    string = ''
    for key in np.sort(indices.keys()):
        string = string + key + ': ' + str(indices[key][0]) + ' ' + indices[key][1] + '\n'
    plt.text(1.02, 1, string, verticalalignment='top', transform=plt.gca().transAxes)
    
    # Draw the hodograph on the Skew-T.
    # TAS 2015-4-16: hodograph doesn't plot for some reason ...
    ax2 = plt.axes([.625,.625,.25,.25])
    below_12km = np.where(interp.to_agl(prof, prof.hght) < 12000)[0]
    u_prof = prof.u[below_12km]
    v_prof = prof.v[below_12km]
    ax2.plot(u_prof[~u_prof.mask], v_prof[~u_prof.mask], 'k-', lw=2)
    ax2.get_xaxis().set_visible(False)
    ax2.get_yaxis().set_visible(False)
    for i in range(10,90,10):
        # Draw the range rings around the hodograph.
        circle = plt.Circle((0,0),i,color='k',alpha=.3, fill=False)
        ax2.add_artist(circle)
    ax2.plot(srwind[0], srwind[1], 'ro') # Plot Bunker's Storm motion right mover as a red dot
    ax2.plot(srwind[2], srwind[3], 'bo') # Plot Bunker's Storm motion left mover as a blue dot
    
    ax2.set_xlim(-60,60)
    ax2.set_ylim(-60,60)
    ax2.axhline(y=0, color='k')
    ax2.axvline(x=0, color='k')
    plt.show()
    """