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
import multiprocessing


#Read in with xarray
ds = xr.open_dataset('/uufs/chpc.utah.edu/common/home/steenburgh-group7/tom/cm1/cm1r19B/run/cm1run_20ms_1500m_tug_wide_45sec.nc')

#xpos = np.load('xpos_traj_500m.npy')
#ypos = np.load('ypos_traj_500m.npy')
#zpos_terrain = np.load('zpos_traj_500m.npy')
#color = np.load('color_traj_500m.npy')


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
num_seeds_y = 80
num_seeds_x = 475
time_steps = 50 #ds.th[:,0,0,0].size-2
resolution = 200 #meters
time_step_length = 45 #seconds
###############################################################################
###############################################################################


#Create arrays for location and variable color of parcels
xpos = np.zeros((time_steps, num_seeds_x, num_seeds_y))
ypos = np.zeros((time_steps, num_seeds_x, num_seeds_y))
zpos = np.zeros((time_steps, num_seeds_x, num_seeds_y))

#Color variable (theta)
color = np.zeros((time_steps, num_seeds_x, num_seeds_y))

    

###############################################################################
#############  These may also be important to the user  #######################
###############################################################################
#Choose starting point (left side of domain, evenly spaced in y, in gridpoints)
#xpos
for i in range(num_seeds_x):
    xpos[0,i,:] = i*10

#ypos
for i in range(num_seeds_y):
    ypos[0,:,i] = i*10

#zpos   
zpos[0,:,:] = 5 #Height to start trajectories at
###############################################################################
###############################################################################

##Add chunck of nan's to array to parcels always fill domain
#n = np.empty((len(xpos[:,0,0]), len(xpos[0,:,0]), len(xpos[0,0,:])))
#n[:] = np.nan
#
#xpos = np.append(xpos, n, axis = 1)
#ypos = np.append(ypos, n, axis = 1)
#zpos = np.append(zpos, n, axis = 1)
#color = np.append(color, n, axis = 1)



     

    
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
    th = ds.th[t,:,:,:].values

    

    #Create list of current coordinates
    yloc = np.array(ypos[t,:,:]).flatten()
    xloc = np.array(xpos[t,:,:]).flatten()
    zloc = np.array(zpos[t,:,:]).flatten()
    
    coor = []
    for i in range(len(xloc)):
        coor.append((zloc[i], yloc[i], xloc[i]))
    coor_array = np.array(coor)

    #Calc new xpos
    xpos[t+1,1:,:] = xpos[t,:-1,:] + np.reshape(interpolate.interpn((z,y,x), u, coor[:-num_seeds_y], method='linear', bounds_error=False, fill_value=np.nan)*time_step_length/resolution, (num_seeds_x-1,num_seeds_y))

    #Calc new ypos
    ypos[t+1,1:,:] = ypos[t,:-1,:] + np.reshape(interpolate.interpn((z,y,x), v, coor[:-num_seeds_y], method='linear', bounds_error=False, fill_value=np.nan)*time_step_length/resolution, (num_seeds_x-1,num_seeds_y))

    #Calc new zpos
    zpos[t+1,1:,:] = zpos[t,:-1,:] + np.reshape(interpolate.interpn((z,y,x), w, coor[:-num_seeds_y], method='linear', bounds_error=False, fill_value=np.nan)*time_step_length/resolution, (num_seeds_x-1,num_seeds_y))

    #Color variable (theta)
    color[t,:,:] = np.reshape(interpolate.interpn((z,y,x), th, coor, bounds_error=False, fill_value= np.nan), (num_seeds_x,num_seeds_y))
    
#    #Initialize new set of parcels for each timestep
#    xpos[t+1,0,:] = np.array(xpos[0,0,:])
#    ypos[t+1,0,:] = np.array(ypos[0,0,:])
#    zpos[t+1,0,:] = np.array(zpos[0,0,:])

    
#    #Initialize new set of parcels for each timestep if on an odd time step
#    if t % 2 == 0:
#        xpos[t+1,0,:] = np.nan
#        ypos[t+1,0,:] = np.nan
#        zpos[t+1,0,:] = np.nan
#        
#    else:
#        xpos[t+1,0,:] = np.array(xpos[0,0,:])
#        ypos[t+1,0,:] = np.array(ypos[0,0,:])
#        zpos[t+1,0,:] = np.array(zpos[0,0,:])
    

    #Timer
    stop = time.time()
    print(stop-start)



########  To get the last variable for colors ###############
t = time_steps-1
#Create list of current coordinates
yloc = np.array(ypos[t,:,:]).flatten()
xloc = np.array(xpos[t,:,:]).flatten()
zloc = np.array(zpos[t,:,:]).flatten()

coor = []
for i in range(len(xloc)):
    coor.append((zloc[i], yloc[i], xloc[i]))
  
#Color variable (theta)
th = ds.th[t,:,:,:].values

#Color variable (theta)
color[t,:,:] = np.reshape(interpolate.interpn((z,y,x), th, coor, bounds_error=False, fill_value= np.nan), (num_seeds_x,num_seeds_y))






#Replace all zeros with nan
xpos[xpos == 0] = 'nan'
ypos[ypos == 0] = 'nan'
zpos[zpos == 0] = 'nan'
color[color == 0] = 'nan'

#Save arrays
np.save('xpos_traj_%d00m' % zpos[0,0,0], xpos)
np.save('ypos_traj_%d00m' % zpos[0,0,0], ypos)
np.save('zpos_traj_%d00m' % zpos[0,0,0], zpos)
np.save('color_traj_%d00m' % zpos[0,0,0], color)



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
    #Create list of current coordinates
    yloc = np.array(ypos[t,:,:]).flatten()
    xloc = np.array(xpos[t,:,:]).flatten()

    coor = []
    for i in range(len(xloc)):
        coor.append((yloc[i], xloc[i]))
    coor_array = np.array(coor)

    zpos_terrain[t,:,:] = zpos_terrain[t,:,:] + np.reshape(interpolate.interpn((y,x), zs, coor, method='linear', bounds_error=False, fill_value= np.nan), (num_seeds_x, num_seeds_y))
    
#Save arrays
np.save('xpos_traj_%d00m' % zpos[0,0,0], xpos)
np.save('ypos_traj_%d00m' % zpos[0,0,0], ypos)
np.save('zpos_traj_%d00m' % zpos[0,0,0], zpos_terrain)
np.save('color_traj_%d00m' % zpos[0,0,0], color)

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




##############################   Plots ########################################
    
#Use multiple processors to create images
def Plotting_Loop(i): 


#for i in range(0,time_steps):    
    secs = (i*120)
 
    fig = plt.figure(num=None, figsize=(12,8), facecolor='w', edgecolor='k')
    
    
    ###################  Plot Trajectories in xy  ###############################
    
    #Plot characteristics
    ax = plt.subplot(211,aspect = 'equal')
    plt.title("Forward Trajectories from 500m", fontsize = 18, y = 1.05) 
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
    
    #levels = np.linspace(262.5,271.01,100)
    #levels_ticks = np.arange(262.5,271.01,0.5) 
    
    #Levels for theta
    levels = np.linspace(np.nanpercentile(color, 5), np.nanpercentile(color,95),100)
    levels_ticks = np.linspace(np.nanpercentile(color, 5), np.nanpercentile(color,95), 8)
    levels_ticks_labels = np.around(np.linspace(np.nanpercentile(color, 5), np.nanpercentile(color,95), 8),2)


    #Plot
    cmap = cm.viridis
    trajectories_xy = plt.scatter(xpos[i,:,:], ypos[i,:,:], c = color[i,:,:], cmap = cmap, norm=matplotlib.colors.BoundaryNorm(levels,cmap.N), s = 2.5, zorder = 5)
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
    trajectories_xz = plt.scatter(xpos[i,:,:], zpos_terrain[i,:,:], c = color[i,:,:], cmap = cmap, norm=matplotlib.colors.BoundaryNorm(levels,cmap.N), s = 2.5, zorder = 5)
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
    

    
    plt.savefig("/uufs/chpc.utah.edu/common/home/u1013082/public_html/phd_plots/cm1/png_for_gifs/forward_trajectory_1500_everywhere45sec5500m%03d.png" % i, dpi = 100)
    plt.close(fig)
    print i

#run function to create images
pool = multiprocessing.Pool(30) #number of processors
pool.map(Plotting_Loop, range(0,time_steps))

##Build GIF
#os.system('module load imagemagick')
#os.system('convert -delay 10 -quality 100 /uufs/chpc.utah.edu/common/home/u1013082/public_html/phd_plots/cm1/png_for_gifs/forward_trajectory_1500m_everywhere500m*.png  /uufs/chpc.utah.edu/common/home/u1013082/public_html/phd_plots/cm1/gifs/forward_trajectory_1500m_20ms_everywhere_500m.gif')

###Delete PNGs





