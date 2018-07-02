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
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.gridspec as gridspec
import multiprocessing

#Read in with xarray
ds = xr.open_dataset('/uufs/chpc.utah.edu/common/home/steenburgh-group7/tom/cm1/output/netcdf_output/tug/cm1run_150m_15ms_2000m_90sec.nc')

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
#Forward trajectories
num_seeds_y = 130
num_seeds_x = 400
time_steps = ds.th[:,0,0,0].size-2
hor_resolution = 200 #meters
vert_resolution = 100 #meters
time_step_length = 45 #seconds
###############################################################################
###############################################################################
#%%
#
##Create arrays for location and variable color of parcels
#xpos = np.zeros((time_steps, num_seeds_x, num_seeds_y))
#ypos = np.zeros((time_steps, num_seeds_x, num_seeds_y))
#zpos = np.zeros((time_steps, num_seeds_x, num_seeds_y))
#
##Color variable (theta)
#color = np.zeros((time_steps, num_seeds_x, num_seeds_y))
#
#    
#
################################################################################
##############  These may also be important to the user  #######################
################################################################################
##Choose starting point (left side of domain, evenly spaced in y, in gridpoints)
##xpos
#for i in range(num_seeds_x):
#    xpos[0,i,:] = i*10
#
##ypos
#for i in range(num_seeds_y):
#    ypos[0,:,i] = i*10
#
##zpos   
#zpos[0,:,:] = 5 #Height to start trajectories at
################################################################################
################################################################################
#     
##Add chunck of nan's to array to parcels always fill domain
#n = np.empty((len(xpos[:,0,0]), len(xpos[0,:,0]), len(xpos[0,0,:])))
#n[:] = np.nan
#
#xpos = np.append(xpos, n, axis = 1)
#ypos = np.append(ypos, n, axis = 1)
#zpos = np.append(zpos, n, axis = 1)
#color = np.append(color, n, axis = 1)
#
#
#
#     
#
#    
##Outside loop goes over all time steps.  This was done to reduce time and memory
##issues by saving only one time step at a time to an attay
#for t in range(time_steps-1):
#    
#    #Timer
#    print t
#    start = time.time()
#    
#
#    #Interpolation function only accepts numpy arrays. This method is quickest for animations
#    start = time.time()
#    u = ds.uinterp[t,:,:,:].values
#    v = ds.vinterp[t,:,:,:].values
#    w = ds.winterp[t,:,:,:].values
#    th = ds.th[t,:,:,:].values
#
#    
#
#    #Create list of current coordinates
#    yloc = np.array(ypos[t,:,:]).flatten()
#    xloc = np.array(xpos[t,:,:]).flatten()
#    zloc = np.array(zpos[t,:,:]).flatten()
#    
#    coor = []
#    for i in range(len(xloc)):
#        coor.append((zloc[i], yloc[i], xloc[i]))
#    coor_array = np.array(coor)
#
#    #Calc new xpos
#    xpos[t+1,1:,:] = xpos[t,:-1,:] + np.reshape(interpolate.interpn((z,y,x), u, coor[:-num_seeds_y], method='linear', bounds_error=False, fill_value=np.nan)*time_step_length/hor_resolution, (2*num_seeds_x-1,num_seeds_y))
#
#    #Calc new ypos
#    ypos[t+1,1:,:] = ypos[t,:-1,:] + np.reshape(interpolate.interpn((z,y,x), v, coor[:-num_seeds_y], method='linear', bounds_error=False, fill_value=np.nan)*time_step_length/hor_resolution, (2*num_seeds_x-1,num_seeds_y))
#
#    #Calc new zpos
#    zpos[t+1,1:,:] = zpos[t,:-1,:] + np.reshape(interpolate.interpn((z,y,x), w, coor[:-num_seeds_y], method='linear', bounds_error=False, fill_value=np.nan)*time_step_length/vert_resolution, (2*num_seeds_x-1,num_seeds_y))
#
#    #Color variable (theta)
#    color[t,:,:] = np.reshape(interpolate.interpn((z,y,x), th, coor, bounds_error=False, fill_value= np.nan), (2*num_seeds_x,num_seeds_y))
#    
##    #Initialize new set of parcels for each timestep
##    xpos[t+1,0,:] = np.array(xpos[0,0,:])
##    ypos[t+1,0,:] = np.array(ypos[0,0,:])
##    zpos[t+1,0,:] = np.array(zpos[0,0,:])
#
#    
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
#    
#
#    #Timer
#    stop = time.time()
#    print(stop-start)
#
#
#########  To get the last variable for colors ###############
#t = time_steps-1
##Create list of current coordinates
#yloc = np.array(ypos[t,:,:]).flatten()
#xloc = np.array(xpos[t,:,:]).flatten()
#zloc = np.array(zpos[t,:,:]).flatten()
#
#coor = []
#for i in range(len(xloc)):
#    coor.append((zloc[i], yloc[i], xloc[i]))
#  
##Color variable (theta)
#th = ds.th[t,:,:,:].values
#
##Color variable (theta)
#color[t,:,:] = np.reshape(interpolate.interpn((z,y,x), th, coor, bounds_error=False, fill_value= np.nan), (2*num_seeds_x,num_seeds_y))
#
#
#
#
#
#
##Replace all zeros with nan
#xpos[xpos == 0] = 'nan'
#ypos[ypos == 0] = 'nan'
#zpos[zpos == 0] = 'nan'
#color[color == 0] = 'nan'
#
##Save arrays
#np.save('xpos_traj_%d00m' % zpos[0,0,0], xpos)
#np.save('ypos_traj_%d00m' % zpos[0,0,0], ypos)
#np.save('zpos_traj_%d00m' % zpos[0,0,0], zpos)
#np.save('color_traj_%d00m' % zpos[0,0,0], color)
#
#
#
################################################################################
###############Add zs to terrain hright in model has terrain####################
################################################################################
#
##Create new array with terrain height
#zpos_terrain = np.array(zpos)*100
#
##check if run has terrain (if not, zs is zero and zpos should not be affected)
#try:
#    zs = np.array(ds.zs[0,:,:])
#except:
#    zs = np.zeros((ds.zs[0,:,0].size, ds.zs[0,0,:].size))  
#  
##loop over all xy parcel location, interpolate to grid to get terrain height,
##then add terrain height to zpos
#for t in range(time_steps-1):
#    #Create list of current coordinates
#    yloc = np.array(ypos[t,:,:]).flatten()
#    xloc = np.array(xpos[t,:,:]).flatten()
#
#    coor = []
#    for i in range(len(xloc)):
#        coor.append((yloc[i], xloc[i]))
#    coor_array = np.array(coor)
#
#    zpos_terrain[t,:,:] = zpos_terrain[t,:,:] + np.reshape(interpolate.interpn((y,x), zs, coor, method='linear', bounds_error=False, fill_value= np.nan)/vert_resolution, (2*num_seeds_x, num_seeds_y))
#    
##Save arrays
#np.save('xpos_traj_%d00m' % zpos[0,0,0], xpos)
#np.save('ypos_traj_%d00m' % zpos[0,0,0], ypos)
#np.save('zpos_traj_%d00m' % zpos[0,0,0], zpos_terrain)
#np.save('color_traj_%d00m' % zpos[0,0,0], color)

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
colors1 = np.array(nclcmaps.colors['temp1'])#perc2_9lev'])
colors_int = colors1.astype(int)
colors = list(colors_int)
cmap_var = nclcmaps.make_cmap(colors, bit=True)




##############################   Plots ########################################
    
    
#Use multiple processors to create images
#def Plotting_Loop(i):
for i in range(200,201):    
    secs = (i*45)
 
    fig = plt.figure(num=None, figsize=(14,7), facecolor='w', edgecolor='k')

    ax = plt.axes([-0.2, 0, 1.2, 1], projection='3d')
    ax.set_aspect(1)
    plt.title("Trajectories from 500m [elapsed time: {0} seconds]".format(np.int(int(secs))), fontsize = 20, y = 1) 

    
#    #color levels
#    levels = np.linspace(262.5,271.01,100)
#    levels_ticks = np.arange(262.5,271.01,1) 
    
    #Levels for theta
    levels = np.linspace(np.nanpercentile(color, 5), np.nanpercentile(color,95),100)
    levels_ticks = np.linspace(np.nanpercentile(color, 5), np.nanpercentile(color,95), 8)
    levels_ticks_labels = np.around(np.linspace(np.nanpercentile(color, 5), np.nanpercentile(color,95), 8),2)
    
    #trajectory plot
    cmap = cm.viridis
    traj = ax.scatter(xpos[i,:,:], ypos[i,:,:], zpos_terrain[i,:,:], c = color[i,:,:].flatten(), cmap = cmap, norm=matplotlib.colors.BoundaryNorm(levels,cmap.N), s = 2.5, zorder = 4)
    ax.view_init(25, -90)
    
    #terrain plot
    X, Y = np.meshgrid(x, y)
    Z = ds.zs[0,:,:]
#    ax.plot_surface(X, Y, Z, rstride=10, cstride=10,  alpha = 1, cmap = cm.Greys, vmin = -800, vmax = 2000, zorder = 1)
    
    #2m temp plot
    
    # fourth dimention - colormap
    # create colormap according to x-value (can use any 50x50 array)
    t2 = ds.t2[i,:,:]
    color_dimension =  t2# change to desired fourth dimension
    minn, maxx = color_dimension.min(), color_dimension.max()
    norm = matplotlib.colors.Normalize(minn+6, maxx)
    m = plt.cm.ScalarMappable(norm=norm, cmap='jet')
    m.set_array([])
    fcolors = m.to_rgba(color_dimension)

    ax.plot_surface(X, Y, Z, rstride=5, cstride=5,  alpha = 0.75, facecolors=fcolors, vmin=minn, vmax=maxx, shade=False, zorder = 3)
    
    #water and land details
    lake = np.array(ds.xland[0,:,:])
    lake[lake == 1] = 'nan'
    ax.plot_surface(X, Y, lake, rstride=5, cstride=5, alpha = 0.5,  color = ('lightskyblue'), zorder = 2)
    
    #xtick
    xtick = np.arange(0,1700,200)
    ax.set_xticks(xtick)
    ax.set_xticklabels(xtick*200/1000)
    
    #yticks
    ytick = np.arange(0,600,100)
    ax.set_yticks(ytick)
    ax.set_yticklabels(ytick*200/1000, horizontalalignment='left', verticalalignment='baseline')
    
    #zticks
    ztick = np.arange(500,4000,500)
    ax.set_zticks(ztick)
    ax.set_zticklabels(ztick, horizontalalignment='left')
    
    #tick labels
    ax.set_xlabel('Distance (km)', labelpad = 30, fontsize = 14)
    ax.set_ylabel('Distance (km)            ', labelpad = 15, fontsize = 12)
    ax.set_zlabel('Height (m)', labelpad = 22, fontsize = 14, rotation = -45)

    #Limits
    ax.set_xlim3d(30,1664)
    ax.set_zlim3d(0,3500)
    ax.set_ylim3d(0,504)
    
    #Colorbar
    cbaxes = fig.add_axes([0.86, 0.27, 0.035, 0.52])             
    cbar = plt.colorbar(traj, cax = cbaxes, ticks = levels_ticks)
    cbar.ax.set_yticklabels(levels_ticks_labels)
    cbar.ax.tick_params(labelsize=10)
    plt.text(-0.30, -0.1, 'Theta [K]', fontsize = 15)
    
    plt.savefig("/uufs/chpc.utah.edu/common/home/u1013082/public_html/phd_plots/cm1/png_for_gifs/forward_trajectory_1500_3deverywhere500m00m%03d_45sec_25deg.png" % i, dpi = 100)
    plt.close(fig)
    print i

#run function to create images
#if __name__ == '__main__':
#    pool = multiprocessing.Pool(15) #number of processors
#    pool.map(Plotting_Loop, range(0,time_steps))
    


##Build GIF
#os.system('module load imagemagick')
#os.system('convert -delay 10 -quality 100 /uufs/chpc.utah.edu/common/home/u1013082/public_html/phd_plots/cm1/png_for_gifs/forward_trajectory_1500_3deverywhere500m*.png  /uufs/chpc.utah.edu/common/home/u1013082/public_html/phd_plots/cm1/gifs/forward_trajectory_1500m_20ms_everywhere_500m2.gif')

###Delete PNGs





