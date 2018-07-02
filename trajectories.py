import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpl_toolkits.basemap
from mpl_toolkits.basemap import Basemap, maskoceans
import pygrib, os, sys
from netCDF4 import Dataset
from numpy import *
import numpy as np
from pylab import *
import time
from datetime import date, timedelta
from matplotlib import animation
import matplotlib.animation as animation
import types
import matplotlib.lines as mlines
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap
import nclcmaps
import pandas as pd
import xarray as xr
from scipy import interpolate
import operator


#Read in with xarray
ds = xr.open_dataset('/uufs/chpc.utah.edu/common/home/steenburgh-group7/tom/cm1/cm1r19B/run/cm1run_20ms_1500m_tug.nc')


#######################   Calculate Trajectories ####################


num_x = ds.uinterp[0,0,0,:].size
num_y = ds.uinterp[0,0,:,0].size
num_z = ds.uinterp[0,:,0,0].size

x = np.arange(0,num_x,1)
y = np.arange(0,num_y,1)
z = np.arange(0,num_z,1)


###############################################################################
#############  These are variables the user changes  ##########################
###############################################################################
#Backward trajectories
num_seeds = 800
time_steps = 60
sets_seeds = 1
start_time_step = 260
hor_resolution = 200 #meters
vert_resolution = 100 #meters
time_step_length = 120 #seconds
###############################################################################
###############################################################################


#x and y location
xpos = np.zeros((sets_seeds, time_steps, num_seeds))
ypos = np.zeros((sets_seeds, time_steps, num_seeds))
zpos = np.zeros((sets_seeds, time_steps, num_seeds))

#Color variable (theta)
color = np.zeros((sets_seeds, time_steps, num_seeds))



###############################################################################
###############  Initial locations (important to user)  #######################
###############################################################################
#Set 1
xpos[0,0,:] = np.arange(800,1600,1) #Start on right side of domain
ypos[0,0,:] = 252 #Evenly spaced in y
zpos[0,0,:] = 20 #Height to start trajectories at

##Set 2
#xpos[1,0,:] = 1350 #Start on right side of domain
#ypos[1,0,:] = np.linspace(162, 342, num_seeds) #Evenly spaced in y
#zpos[1,0,:] = 8 #Height to start trajectories at
#
##Set 3
#xpos[2,0,:] = 1050 #Start on right side of domain
#ypos[2,0,:] = np.linspace(162, 342, num_seeds) #Evenly spaced in y
#zpos[2,0,:] = 1 #Height to start trajectories at
###############################################################################
###############################################################################






#Loop over all time spteps and compute trajectory
for s in range(sets_seeds):
    for t in range(1, time_steps):
        
        print t
    
        #Create list of current coordinates
        yloc = np.array(ypos[s,t-1,:])
        xloc = np.array(xpos[s,t-1,:])
        zloc = np.array(zpos[s,t-1,:])
        
        coor = []
        for i in range(num_seeds):
            coor.append((zloc[i], yloc[i], xloc[i]))
        coor_array = np.array(coor)
    
        start = time.time()
        
        
          
        #Interpolation function only accepts numpy arrays
        #(This takes a long time.  Just use for testing)
    #    u = ds.uinterp[start_time_step-t,:,:,:].values
    #    v = ds.vinterp[start_time_step-t,:,:,:].values
    #    w = ds.winterp[start_time_step-t,:,:,:].values
    #    th = ds.th[start_time_step-t,:,:,:].values
    
    
    
    ###########  Interpolation function only accepts rrays with all values ########
    # All code code in this chuck was developed to speed up placing the 
    # huge data arrays into numpy arrays. The main idea is to extract 
    # the chunkc of data that is necessary to interpolate at all points
    # along the trajectory
    
        #Create arrays of zeros to place data chunks in to maintain locations
        u = np.zeros((ds.uinterp[0,:,0,0].size, ds.uinterp[0,0,:,0].size, ds.uinterp[0,0,0,:].size))
        v = np.zeros((ds.uinterp[0,:,0,0].size, ds.uinterp[0,0,:,0].size, ds.uinterp[0,0,0,:].size))
        w = np.zeros((ds.uinterp[0,:,0,0].size, ds.uinterp[0,0,:,0].size, ds.uinterp[0,0,0,:].size))
        th = np.zeros((ds.uinterp[0,:,0,0].size, ds.uinterp[0,0,:,0].size, ds.uinterp[0,0,0,:].size))
        
        #Number of surrounding points to use.  2 is min necessary for linear interp
        size = 2
        
        #Check if they've all left the domain
        if np.all(np.isnan(coor_array)):
            xmin, xmax, ymin, ymax = (0,0,0,0)
        else:
            #Find max and min in each direction for entire set of coordinates
            xmin = np.int(np.nanmin(coor_array[:,2]))-size
            xmax = np.int(np.nanmax(coor_array[:,2]))+size
            ymin = np.int(np.nanmin(coor_array[:,1]))-size
            ymax = np.int(np.nanmax(coor_array[:,1]))+size
            zmin = np.int(np.nanmin(coor_array[:,0]))-size
            zmax = np.int(np.nanmax(coor_array[:,0]))+size
    
        #Make sure points don't leave domain
        if xmin < 0:
            xmin = 0
        if ymin < 0:
            ymin = 0
        if zmin < 0:
            zmin = 0
    
        if xmax > ds.uinterp[0,0,0,:].size:
            xmax = ds.uinterp[0,0,0,:].size
        if ymax > ds.uinterp[0,0,:,0].size:
            ymax = ds.uinterp[0,0,:,0].size
        if zmax > ds.uinterp[0,:,0,0].size:
            zmax = ds.uinterp[0,:,0,0].size
    
        #Specify and add necessary chunk of data to arrays with zeros
        u[zmin:zmax,ymin:ymax,xmin:xmax] = ds.uinterp[start_time_step-t,zmin:zmax,ymin:ymax,xmin:xmax].values
        v[zmin:zmax,ymin:ymax,xmin:xmax] = ds.vinterp[start_time_step-t,zmin:zmax,ymin:ymax,xmin:xmax].values
        w[zmin:zmax,ymin:ymax,xmin:xmax] = ds.winterp[start_time_step-t,zmin:zmax,ymin:ymax,xmin:xmax].values
        
        #Color variable (theta)
        th[zmin:zmax,ymin:ymax,xmin:xmax] = ds.th[start_time_step-t,zmin:zmax,ymin:ymax,xmin:xmax].values
    
        #Timer
        stop = time.time()
        print(stop-start)
    
    
                
        #Calc new xpos
        xpos[s,t,:] = xpos[s,t-1,:] - interpolate.interpn((z,y,x), u, coor, method='linear', bounds_error=False, fill_value=np.nan)*time_step_length/hor_resolution
    
        #Calc new ypos
        ypos[s,t,:] = ypos[s,t-1,:] - interpolate.interpn((z,y,x), v, coor, method='linear', bounds_error=False, fill_value=np.nan)*time_step_length/hor_resolution
    
        #Calc new zpos
        zpos[s,t,:] = zpos[s,t-1,:] - interpolate.interpn((z,y,x), w, coor, method='linear', bounds_error=False, fill_value= np.nan)*time_step_length/vert_resolution
    
        #Color variable (theta)
        color[s,t-1,:] = interpolate.interpn((z,y,x), th, coor, bounds_error=False, fill_value= np.nan)
    
    
        
    ########  To get the last variable for colors ###############
        
    #Create list of current coordinates
    t = time_steps-1
    yloc = np.array(ypos[s,t,:])
    xloc = np.array(xpos[s,t,:])
    zloc = np.array(zpos[s,t,:])
    
    coor = []
    for i in range(num_seeds):
        coor.append((zloc[i], yloc[i], xloc[i]))
      
    #Color variable (theta)
    th = ds.th[start_time_step-t,:,:,:].values
    
    #Color variable (theta)
    color[s,t,:] = interpolate.interpn((z,y,x), th, coor, bounds_error=False, fill_value= np.nan)





    ###############################################################################
    ##############Add zs to terrain hright in model has terrain####################
    ###############################################################################
    
    #Create new array with terrain height
    zpos_terrain = np.array(zpos)*vert_resolution
    
    #check if run has terrain (if not, zs is zero and zpos should not be affected)
    try:
        zs = np.array(ds.zs[0,:,:])
    except:
        zs = np.zeros((ds.zs[0,:,0].size, ds.zs[0,0,:].size))  
      
    #loop over all xy parcel location, interpolate to grid to get terrain height,
    #then add terrain height to zpos
    
    #Create list of current coordinates
    yloc = np.array(ypos[s,:,:]).flatten()
    xloc = np.array(xpos[s,:,:]).flatten()
    
    coor = []
    for i in range(num_seeds*time_steps):
        coor.append((yloc[i], xloc[i]))
    coor_array = np.array(coor)
    
    zpos_terrain[s,:,:] = zpos_terrain[s,:,:] + np.reshape(interpolate.interpn((y,x), zs, coor, method='linear', bounds_error=False, fill_value= np.nan), (time_steps, num_seeds))






#%%




###############################################################################
############## Set ncl_cmap as the colormap you want ##########################
###############################################################################

### 1) In order for this to work the files "nclcmaps.py" and "__init__.py"
### must be present in the dirctory.
### 2) You must "import nclcmaps"
### 3) The path to nclcmaps.py must be added to tools -> PYTHONPATH manager in SPyder
### 4) Then click "Update module names list" in tolls in Spyder and restart Spyder
                
## The steps above describe the general steps for adding "non-built in" modules to Spyder

###############################################################################
###############################################################################


#Read in colormap and put in proper format
colors1 = np.array(nclcmaps.colors['MPL_YlGnBu'])#perc2_9lev'])
colors_int = colors1.astype(int)
colors = list(colors_int)
cmap_var = nclcmaps.make_cmap(colors, bit=True)

colors1_t = np.array(nclcmaps.colors['MPL_Greys'])
colors_int_t = colors1_t.astype(int)
colors_t = list(colors_int_t)
cmap_terrain = nclcmaps.make_cmap(colors_t, bit=True)


colors1 = np.array(nclcmaps.colors['prcp_1'])#perc2_9lev'])
colors_int = colors1.astype(int)
colors = list(colors_int)
cmap_precip = nclcmaps.make_cmap(colors, bit=True)

colors1_t = np.array(nclcmaps.colors['OceanLakeLandSnow'])
colors_int_t = colors1_t.astype(int)
colors_t = list(colors_int_t)
cmap_terrain = nclcmaps.make_cmap(colors_t, bit=True)








##############################   Plots ########################################
    
    
for i in range(150,151):
#for i in range(0,ds.dbz[:,0,0,0].size-10):
    
    secs = (i*120)+1200
 
    fig = plt.figure(num=None, figsize=(12,8), facecolor='w', edgecolor='k')
    
    
    ###################  Plot Trajectories in xy  ###############################
    
    #Plot characteristics
    ax = plt.subplot(211,aspect = 'equal')
    plt.title("{0} Set[s] of {1} Backward Trajectories from {2}00m".format(sets_seeds, num_seeds, np.int(zpos[0,0,0])), fontsize = 18, y = 1.05) 
    plt.subplots_adjust(left=0.07, bottom=0.11, right=0.87, top=0.93, wspace=0, hspace=0.15)
    #plt.axis('equal')
    #plt.axis('off')
    
    
    #Levels for theta
    levels = np.linspace(np.nanpercentile(color, 5), np.nanpercentile(color,95),100)
    levels_ticks = np.linspace(np.nanpercentile(color, 5), np.nanpercentile(color,95), 8)
    levels_ticks_labels = np.around(np.linspace(np.nanpercentile(color, 5), np.nanpercentile(color,95), 8),2)
    
    #Plot
    cmap = cm.viridis
    trajectories_xy = plt.scatter(xpos, ypos, c = color, cmap = cmap, norm=matplotlib.colors.BoundaryNorm(levels,cmap.N), s = 2.5, zorder = 5)
    plt.xticks(np.arange(0,ds.uinterp[0,0,0,:].size,200))
    plt.yticks(np.arange(0,ds.uinterp[0,0,:,0].size,100))
    ax.set_xticklabels(np.arange(0,ds.uinterp[0,0,0,:].size*200/100,40), fontsize = 10)
    ax.set_yticklabels(np.arange(0,ds.uinterp[0,0,0,:].size*200/100,40), fontsize = 10)
    plt.ylabel('Distance (km)', fontsize = 15)
    plt.xlim([0,ds.uinterp[0,0,0,:].size])
    plt.ylim([0,ds.uinterp[0,0,:,0].size])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    #Water and land details
    levels_water = [1.5, 2.5]
    levels_terrain = [0,1.5]
    xland_plt = plt.contourf(ds.xland[0,:,:], levels_water, alpha = 1, colors = ('lightsteelblue'), zorder = 2)
    xland_plt = plt.contourf(ds.xland[0,:,:], levels_terrain, alpha = 1, colors = ('gainsboro'), zorder = 1)
    
    #Terrain
    terrain_levels = np.arange(-1, 3000.1, 200)
    terrain = plt.contourf(ds.zs[0,:,:], terrain_levels, alpha = 1, cmap = cm.Greys, vmin = -800, vmax = 3000, zorder = 1)
    



    ###################  Plot Trajectories in xz  ###############################


    #Plot characteristics
    ax = plt.subplot(212)
    
    #Plot
    trajectories_xz = plt.scatter(xpos.flatten(), zpos.flatten()*100, c = color.flatten(), cmap = cmap, norm=matplotlib.colors.BoundaryNorm(levels,cmap.N), s = 2.5, zorder = 3)


    
    plt.grid(True, color = 'white', )
    plt.ylabel('Height (m)', fontsize = 15)
    plt.xlabel('Distance (km)', fontsize = 15)
    plt.xticks(np.arange(0,ds.uinterp[0,0,0,:].size,200))
    ytick = np.arange(0,4001,250)
    plt.yticks(ytick)
    ax.set_xticklabels(np.arange(0,ds.uinterp[0,0,0,:].size*200/100,40), fontsize = 10)
    plt.xlim([0,ds.uinterp[0,0,0,:].size])
    plt.ylim([0,ytick[-1]])
    ax.tick_params(axis='y', labelsize=11)
    plt.axvspan(0,ds.uinterp[0,0,0,:].size,color='gainsboro',lw=0)
    
    #Labels
    sub_title = ['[Run: 20 $\mathregular{ms^{-1}}$ and No Terrain]']
    props = dict(boxstyle='square', facecolor='white', alpha=1)
    ax.text(1350, -450, sub_title[0], fontsize = 15, zorder = 5)


    #Plot Terrain
    x1d = np.arange(0,ds.zs[0,0,:].size,1)
    ymid = np.int(ds.zs[0,:,0].size/2)
    terrain = plt.plot(x1d, ds.zs[0,ymid,:], c = 'slategrey', linewidth = 3, zorder = 4)
    
    #Plot Lake
    lake = np.array(ds.xland[0,ymid,:])
    lake[lake == 1] = np.nan
    lake_plt = plt.plot(x1d, lake-3, c = 'blue', linewidth = 4, zorder = 5)
    
    #Colorbar
    cbaxes = fig.add_axes([0.9, 0.27, 0.035, 0.52])             
    cbar = plt.colorbar(trajectories_xy, cax = cbaxes, ticks = levels_ticks)
    cbar.ax.set_yticklabels(levels_ticks_labels)
    cbar.ax.tick_params(labelsize=10)
    plt.text(-0.30, -0.08, 'Theta (K)', fontsize = 15)
    
    
    plt.savefig("/uufs/chpc.utah.edu/common/home/u1013082/public_html/phd_plots/cm1/plots/backward_trajectory_20ms_1500m_height{0}00m_ydist{1}grdpts_{2}sets.png".format(np.int(zpos[0,0,0]), np.int(xpos[0,0,0]), np.int(sets_seeds)), dpi = 600)
    plt.close(fig)
    

##Build GIF
#os.system('module load imagemagick')
#os.system('convert -delay 12 -quality 100 /uufs/chpc.utah.edu/common/home/u1013082/public_html/phd_plots/cm1/png_for_gifs/single_run_cref_th_div_*.png  /uufs/chpc.utah.edu/common/home/u1013082/public_html/phd_plots/cm1/gifs/single_run_cref_th_div.gif')

###Delete PNGs





