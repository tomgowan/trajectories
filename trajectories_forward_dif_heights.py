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


###############################################################################
###############################################################################   
#######################   Calculate Trajectories ##############################
###############################################################################
###############################################################################

#Dimension size variables
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
num_seeds = 48
time_steps = 220
resolution = 200 #meters
time_step_length = 120 #seconds
###############################################################################
###############################################################################


#Create arrays for location and variable color of parcels
xpos = np.zeros((time_steps, time_steps, num_seeds))
ypos = np.zeros((time_steps, time_steps, num_seeds))
zpos = np.zeros((time_steps, time_steps, num_seeds))

#Color variable (height)
color = np.zeros((time_steps, time_steps, num_seeds))


    
#Initialize arrays of parcels for each time step
for j in range(time_steps-1):
    
    ###############################################################################
    #############  These may also be important to the user  #######################
    ###############################################################################
    #Choose starting point (left side of domain, evenly spaced in y, in gridpoints)
    xpos[j,j,:] = 0 #Start on left side of domain
    ypos[j,j,:] = np.linspace(0, num_y-1, num_seeds) #Evenly spaced in y
    for i in np.arange(0,num_seeds,4):
        zpos[j,j,i:i+4] = np.arange(5,20.1,5) #Height to start trajectories at
    #Colors for heights
    for i in np.arange(0,num_seeds,4):
        color[j,j:,i:i+4] = np.arange(5,20.1,5) #Height to start trajectories at
    
    ###############################################################################
    ###############################################################################
     
  
    
#Outside loop goes over all time steps.  This was done to reduce time and memory
#issues by saving only one time step at a time to an attay
for t in range(time_steps-1):
    
    #Timer
    print t
    start = time.time()
    

    #Interpolation function only accepts numpy arrays. This method is quickest for animations
    start = time.time()
    u = ds.uinterp[t,:,:,:].values
    v = ds.vinterp[t,:,:,:].values
    w = ds.winterp[t,:,:,:].values


    
    #Loop over each set of parcels (one for each time step)
    for j in range(t+1): 

        #Create list of current coordinates
        yloc = np.array(ypos[j,t,:])
        xloc = np.array(xpos[j,t,:])
        zloc = np.array(zpos[j,t,:])
        
        coor = []
        for i in range(num_seeds):
            coor.append((zloc[i], yloc[i], xloc[i]))
        coor_array = np.array(coor)
        
        
        #Calc new xpos
        xpos[j,t+1,:] = xpos[j,t,:] + interpolate.interpn((z,y,x), u, coor, method='linear', bounds_error=False, fill_value=np.nan)*time_step_length/resolution
    
        #Calc new ypos
        ypos[j,t+1,:] = ypos[j,t,:] + interpolate.interpn((z,y,x), v, coor, method='linear', bounds_error=False, fill_value=np.nan)*time_step_length/resolution
    
        #Calc new zpos
        zpos[j,t+1,:] = zpos[j,t,:] + interpolate.interpn((z,y,x), w, coor, method='linear', bounds_error=False, fill_value= np.nan)*time_step_length/resolution
    


    #Timer
    stop = time.time()
    print(stop-start)






#Replace all zeros with nan
xpos[xpos == 0] = 'nan'
ypos[ypos == 0] = 'nan'
zpos[zpos == 0] = 'nan'
color[color == 0] = 'nan'

#Save arrays
#np.save('xpos_1000m', xpos)
#np.save('ypos_1000m', ypos)
#np.save('zpos_1000m', zpos)
#np.save('color_1000m', color)



###############################################################################
##############Add zs to terrain hright in model has terrain####################
###############################################################################

#Create new array with terrain height
zpos_terrain = np.array(zpos)*100

#check if run has terrain (if not, zs is zero and zpos should not be affected)
try:
    zs = np.array(ds.zs[0,:,:])
except:
    zs = np.zeros((ds.zs[0,:,0].size, ds.zs[0,0,:].size))  
  
#loop over all xy parcel location, interpolate to grid to get terrain height,
#then add terrain height to zpos

for t in range(time_steps-1):
    for j in range(t+1):    
        yloc = np.array(ypos[j,t,:])
        xloc = np.array(xpos[j,t,:])
        
        coor = []  
        for i in range(num_seeds):
            coor.append((yloc[i], xloc[i]))

        zpos_terrain[j,t,:] = zpos_terrain[j,t,:] + interpolate.interpn((y,x), zs, coor, method='linear', bounds_error=False, fill_value= np.nan)

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
    
    
#for i in range(0,time_steps):
for i in range(165,168):

    
    secs = (i*120)
 
    fig = plt.figure(num=None, figsize=(12,8), facecolor='w', edgecolor='k')
    
    
    ###################  Plot Trajectories in xy  ###############################
    
    #Plot characteristics
    ax = plt.subplot(211,aspect = 'equal')
    plt.title("Forward Trajectories from %d00m" % zpos[0,0,0], fontsize = 18, y = 1.05) 
    plt.subplots_adjust(left=0.07, bottom=0.11, right=0.87, top=0.93, wspace=0, hspace=0.15)
    #plt.axis('equal')
    #plt.axis('off')
    
    
    #Levels for theta
#    levels = np.linspace(np.int(np.nanmin(color[0,:,:]))+0.25,np.int(np.nanmax(color[0,:,:]))-0.24999,100)
#    levels_ticks = np.arange(np.int(np.nanmin(color[0,:,:]))+0.25,np.int(np.nanmax(color[0,:,:]))-0.24999,0.5)

    #levels = np.linspace(265,272.01,100)
    #levels_ticks = np.arange(265,272.01,0.5)
   
    #levels = np.linspace(268,279.01,100)
    #levels_ticks = np.arange(268,279.01,0.5)
    
    levels = np.linspace(5,20.1,5)
    levels_ticks = np.arange(5,20.1,5) 


    #Plot
    cmap = cm.tab10
    trajectories_xy = plt.scatter(xpos[:,i,:], ypos[:,i,:], c = color[:,i,:], cmap = cmap, norm=matplotlib.colors.BoundaryNorm(levels,cmap.N), s = 2.5, zorder = 5)
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
    trajectories_xz = plt.scatter(xpos[:,i,:], zpos_terrain[:,i,:], c = color[:,i,:], cmap = cmap, norm=matplotlib.colors.BoundaryNorm(levels,cmap.N), s = 2.5, zorder = 5)
    plt.grid(True, color = 'white', )
    plt.ylabel('Height (m)', fontsize = 15)
    plt.xlabel('Distance (km)', fontsize = 15)
    plt.xticks(np.arange(0,ds.uinterp[0,0,0,:].size,200))
    plt.yticks(np.arange(0,3501,500))
    ax.set_xticklabels(np.arange(0,ds.uinterp[0,0,0,:].size*200/100,40), fontsize = 10)
    plt.xlim([0,ds.uinterp[0,0,0,:].size])
    plt.ylim([0,3500])
    ax.tick_params(axis='y', labelsize=11)
    plt.axvspan(0,ds.uinterp[0,0,0,:].size,color='gainsboro',lw=0)
    
    #Plot Terrain
    x1d = np.arange(0,ds.zs[0,0,:].size,1)
    ymid = np.int(ds.zs[0,:,0].size/2)
    terrain = plt.plot(x1d, ds.zs[0,ymid,:], c = 'slategrey', linewidth = 3, zorder = 4)
    
    #Plot Lake
    lake = np.array(ds.xland[0,ymid,:])
    lake[lake == 1] = np.nan
    lake_plt = plt.plot(x1d, lake-3, c = 'blue', linewidth = 4, zorder = 5)
    
    #Labels
    sub_title = ['[Run: 20 $\mathregular{ms^{-1}}$ and 1500m]']
    props = dict(boxstyle='square', facecolor='white', alpha=1)
    ax.text(1350, -650, sub_title[0], fontsize = 15, zorder = 5)



    
    #Colorbar
    cbaxes = fig.add_axes([0.9, 0.27, 0.035, 0.52])             
    cbar = plt.colorbar(trajectories_xy, cax = cbaxes, ticks = levels_ticks)
    cbar.ax.tick_params(labelsize=10)
    plt.text(-0.30, -0.08, 'Theta [K]', fontsize = 15)
    

    
    plt.savefig("/uufs/chpc.utah.edu/common/home/u1013082/public_html/phd_plots/cm1/png_for_gifs/forward_trajectory_1500m_diff_heights%d" % zpos[0,0,0] + "00m%03d.png" % i, dpi = 200)
    plt.close(fig)
    print i

##Build GIF
#os.system('module load imagemagick')
#os.system('convert -delay 10 -quality 100 /uufs/chpc.utah.edu/common/home/u1013082/public_html/phd_plots/cm1/png_for_gifs/forward_trajectory_1500m_500m_release*.png  /uufs/chpc.utah.edu/common/home/u1013082/public_html/phd_plots/cm1/gifs/forward_trajectory_1500m_20ms_500m_release.gif')

###Delete PNGs





