import matplotlib
matplotlib.use('Agg')
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
import multiprocessing as mp
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.patches as patches

#Read in with xarray
ds = xr.open_dataset('/uufs/chpc.utah.edu/common/home/steenburgh-group7/tom/cm1/output/netcdf_output/tug/cm1run_150m_15ms_0000m_60sec.nc')



###############################################################################
###############################################################################   
##########################   Set up Trajectories ##############################
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
##################### INFO TO CALCULATE SEEDS #################################
#############  These are variables the user changes  ##########################
###############################################################################
#Backward trajectories
parcel_spread = 13
num_seeds_z = 151 #Up to 5000m (3 seeds every vertical grid point)
num_seeds_y = np.int(ds.ny/parcel_spread) 
num_seeds_x = np.int(ds.nx/parcel_spread) 
time_steps = 450#ds.uinterp[:,0,0,0].size-2 #Number of time steps to run trajectories forward
start_time_step = 10 #Starting time step
hor_resolution = 150 #[meters]
vert_resolution = 100 #[meters] (resolution away from terrain)
time_step_length = 60.0 #[seconds]
###############################################################################
###############################################################################


#Create arrays for location a(dim =d variable color of parcels
xpos = np.zeros((time_steps, num_seeds_y, num_seeds_x))
ypos = np.zeros((time_steps, num_seeds_y, num_seeds_x))
zpos = np.zeros((time_steps, num_seeds_y, num_seeds_x))
zpos_grid = np.zeros((time_steps, num_seeds_y, num_seeds_x))
change_zs = np.zeros((time_steps, num_seeds_y, num_seeds_x))


#Variable
var_name1 = 'th'
variable1 = np.zeros((time_steps, num_seeds_y, num_seeds_x))






###############################################################################
##################### INITIAL LOCATION OF SEEDS ###############################
#########  (These may also be important to the user)  #########################
###############################################################################
#Choose starting point 
#xpos
for i in range(num_seeds_x):
    xpos[0,:,i] = i*parcel_spread #Seeds at every 10th x-gridpoint
xpos_init = np.array(xpos[0,:,0])

#ypos   
for i in range(num_seeds_y):
    ypos[0,i,:] = i*parcel_spread #Seeds at every 10th y-gridpoint
ypos_init = np.array(ypos[0,:,0])

#zpos
zpos[0,:,:] =  5 #Height to start trajectories at
zpos_init = np.array(zpos[0,:,0])


###############################################################################
###############################################################################



#%%
###############################################################################
##################### Initilize surface terrain height ########################
###############################################################################
#####  Use model terrain to correctly calculate trajectpries.  Essentially ####
#####  apply parcel movement due to u, v, and w to surface terrain height  ####
####   modified array                                                      ####
###############################################################################
 
#    
##check if run has terrain (if not, zs is zero and zpos should not be affected)
#try:
#    zh = np.array(ds.zh[0,:,:,:])
#    
#    #Add terrain height
#    #Create list of current coordinates for terrain addition
#    xloc = np.array(xpos[0,:,:]).flatten()
#    yloc = np.array(ypos[0,:,:]).flatten()
#    zloc = np.array(zpos[0,:,:]).flatten()
#    coor_terrain = []
#    for i in range(len(xloc)):
#        coor_terrain.append((zloc[i], yloc[i], xloc[i]))
#    
#    zpos_terrain = np.array(zpos)
#    zpos_terrain[0,:,:] = np.reshape(interpolate.interpn((z,y,x), zh, coor_terrain, method='linear', bounds_error=False, fill_value= 0), (num_seeds_y, num_seeds_x))
#    zpos_terrain_init = np.array(zpos_terrain[0,:,0])
#    
#    #This gets us the grid spacing for the vertical grid
#    z_grid = zh[1:,:,:]-zh[:-1,:,:]
#
#except:
#    zh = np.zeros((ds.nz,ds.ny, ds.nx))
#    zpos_terrain = np.copy(zpos)*vert_resolution
#    zpos_terrain_init = np.copy(zpos_terrain[0,:,0])
#    z_grid = np.ones((ds.nz-1,ds.ny, ds.nx))*vert_resolution
#    
#
#
################################################################################
################################################################################   
###########################  Calculate Trajectories ############################
################################################################################
################################################################################
#
##Loop over all time spteps and compute trajectory
#for t in range(time_steps-1):
#    
#    print t
#    start = time.time() #Timer
#    
#    #Get model data
#    u = ds.uinterp[start_time_step + t,:,:,:].values
#    v = ds.vinterp[start_time_step + t,:,:,:].values
#    w = ds.winterp[start_time_step + t,:,:,:].values
#    var1 = getattr(ds,var_name1)[start_time_step + t,:,:,:].values
#
#                                
#
#    #####################   Get x,y,z for variables ###########################
#    ##  Because the model output is terrain following, the change in surface ##
#    ##  height at each timestep must be added to the height (zpos) of the    ##
#    ##  parcels.                                                             ##
#    ###########################################################################
#        
#    #Get surface height grid
#    try:
#        zs = np.array(ds.zs[0,:,:])
#    except:
#        zs = np.zeros((ds.ny, ds.nx))  
#    
#    #Need because length of arrays change
#    num_seeds_x = len(xpos[0,0,:])
#    num_seeds_y = len(xpos[0,:,0])
#
#    #Get previous and current x and y positions to calc change in surface height
#    if t > 0:
#        
#        #x and y for one time step back
#        xloc = np.copy(xpos[t-1,:,:]).flatten()
#        yloc = np.copy(ypos[t-1,:,:]).flatten()
#        coor_xminus1 = []
#        for i in range(len(xloc)):
#            coor_xminus1.append((yloc[i], xloc[i]))
#            
#            
#        #xand y for current time step
#        xloc = np.copy(xpos[t,:,:]).flatten()
#        yloc = np.copy(ypos[t,:,:]).flatten()
#        coor_x = []
#        for i in range(len(xloc)):
#            coor_x.append((yloc[i], xloc[i]))
#            
#        #Calc surface height at each time
#        xminus1_h = interpolate.interpn((y,x), zs, coor_xminus1, method='linear', bounds_error=False, fill_value= 0)
#        x_h = interpolate.interpn((y,x), zs, coor_x, method='linear', bounds_error=False, fill_value= 0)
#        
#        #Calc change in surface height during last timestep, which will be added to zpos
#        change_zs[t,:,:] =  np.reshape(x_h - xminus1_h, (num_seeds_y, num_seeds_x))
#    
#    
#    #Get get x, y, and z positions to calc parcel movement
#    xloc = np.copy(xpos[t,:,:]).flatten()
#    yloc = np.copy(ypos[t,:,:]).flatten()
#    zloc = np.copy(zpos[t,:,:]).flatten()
#    coor_var = []
#    for i in range(len(xloc)):
#        coor_var.append((zloc[i], yloc[i], xloc[i])) 
#    
#    
#    #####################   Calc new xpos #####################################
#    xpos[t+1,:,:] = xpos[t,:,:] + np.reshape(interpolate.interpn((z,y,x), u, coor_var, method='linear', bounds_error=False, fill_value=np.nan)*time_step_length/hor_resolution, (num_seeds_y, num_seeds_x))
#
#    #####################   Calc new ypos #####################################
#    ypos[t+1,:,:]  = ypos[t,:,:] + np.reshape(interpolate.interpn((z,y,x), v, coor_var, method='linear', bounds_error=False, fill_value=np.nan)*time_step_length/hor_resolution, (num_seeds_y, num_seeds_x))
#
#    #####################   Calc new zpos #####################################
#    #zpos grid spacing
#    zpos_grid[t,:,:] = np.reshape(interpolate.interpn((z[:-1],y,x), z_grid, coor_var, method='linear', bounds_error=False, fill_value= 0), (num_seeds_y, num_seeds_x))
#    #terrain-following (includes change in surface height)
#    zpos[t+1,:,:]  = zpos[t,:,:] - change_zs[t,:,:]/zpos_grid[t,:,:] + np.reshape(interpolate.interpn((z,y,x), w, coor_var, method='linear', bounds_error=False, fill_value= 0), (num_seeds_y, num_seeds_x))*time_step_length/zpos_grid[t,:,:]
#    #terrain-height coordinates
#    zpos_terrain[t+1,:,:]  = zpos_terrain[t,:,:] + np.reshape(interpolate.interpn((z,y,x), w, coor_var, method='linear', bounds_error=False, fill_value= 0)*time_step_length, (num_seeds_y, num_seeds_x))
#
#    zpos = zpos.clip(min=0) #Prevents parcel from going into the ground
#    
#        
#    #Prevents parcel from going into ground
#    xloc = np.copy(xpos[t,:,:]).flatten()
#    yloc = np.copy(ypos[t,:,:]).flatten()
#    coor_terrain = []
#    for i in range(len(xloc)):
#        coor_terrain.append((yloc[i], xloc[i])) 
#    surface_height = np.reshape(interpolate.interpn((y,x), zs, coor_terrain, method='linear', bounds_error=False, fill_value= 0)/vert_resolution, (num_seeds_y, num_seeds_x))
#    zpos_terrain[t,:,:] = zpos_terrain[t,:,:].clip(min=surface_height) 
#    
#    
#    #Variables
#    variable1[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var1, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan), (num_seeds_y, num_seeds_x))
#    
#    
#    
#    ##################   Initialize new set of parcels #######################
#    #Initalize new set of parcels to enter the domain every other time step. #
#    #Idea is to create new 2d array and append it to the originals          ##
#    ##########################################################################
#
#    if t % 2 == 0: #Only send in new parcels every other timestep
#        #Arrays for new set of trajectories as that enter the domain
#        xpos_new = np.zeros((time_steps, num_seeds_y, 1))
#        ypos_new = np.zeros((time_steps, num_seeds_y, 1))
#        zpos_new = np.zeros((time_steps, num_seeds_y, 1))
#        zpos_terrain_new = np.zeros((time_steps, num_seeds_y, 1))
#        zpos_grid_new = np.zeros((time_steps, num_seeds_y, 1))
#        variable1_new = np.zeros((time_steps, num_seeds_y, 1))
#        change_zs_new = np.zeros((time_steps, num_seeds_y, 1))
#        
#        
#        #Location of new parcels that enter domain
#        xpos_new[t+1,:,0] = np.array(xpos_init)
#        ypos_new[t+1,:,0] = np.array(ypos_init)
#        zpos_new[t+1,:,0] = np.array(zpos_init)
#        zpos_terrain_new[t+1,:,0] = np.array(zpos_terrain_init)
#        
#        #Set past time steps to nan
#        xpos_new[:t,:,:] = np.nan
#        ypos_new[:t,:,:] = np.nan
#        zpos_new[:t,:,:] = np.nan
#        zpos_terrain_new[:t,:,:] = np.nan
#    
#        
#        #Append new arrays to old arrays
#        xpos = np.append(xpos_new, xpos, axis = 2)
#        ypos = np.append(ypos_new, ypos, axis = 2)
#        zpos = np.append(zpos_new, zpos, axis = 2)
#        zpos_terrain = np.append(zpos_terrain_new, zpos_terrain, axis = 2)
#        zpos_grid = np.append(zpos_grid_new, zpos_grid, axis = 2)
#        variable1 = np.append(variable1_new, variable1, axis = 2)
#        change_zs = np.append(change_zs_new, change_zs, axis = 2)
#        
#
#
#
#    #Timer
#    stop = time.time()
#    print(stop-start)
# 
#
#
#
#
###############  Get variable data for final time step #########################
#t = time_steps-1
#
##Need because length of arrays change
#num_seeds_x = len(xpos[0,0,:])
#num_seeds_y = len(xpos[0,:,0])
#
##Get get x, y, and z positions to calc parcel movement
#xloc = np.copy(xpos[t,:,:]).flatten()
#yloc = np.copy(ypos[t,:,:]).flatten()
#zloc = np.copy(zpos[t,:,:]).flatten()
#coor_var = []
#for i in range(len(xloc)):
#    coor_var.append((zloc[i], yloc[i], xloc[i])) 
#
#
##Color variable (theta)
#variable1[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var1, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan), (num_seeds_y, num_seeds_x))
#
##%%
##Save arrays
#np.save('xpos_traj_disp_15ms_0000m_everywhere', xpos)
#np.save('ypos_traj_disp_15ms_0000m_everywhere', ypos)
#np.save('zpos_traj_disp_15ms_0000m_everywhere', zpos_terrain)
#np.save('%s_traj_disp_15ms_0000m_everywhere' %var_name1, variable1)




#%%
xpos = np.load('xpos_traj_disp_15ms_0000m_everywhere.npy')
ypos = np.load('ypos_traj_disp_15ms_0000m_everywhere.npy')
zpos_terrain = np.load('zpos_traj_disp_15ms_0000m_everywhere.npy')
variable1 = np.load('%s_traj_disp_15ms_0000m_everywhere.npy' % var_name1)



#%%





###############################################################################
#############################   PLOTS   #######################################
###############################################################################



#########  Create Grid ########
## The code below makes the data terrain following 
ymid = np.int(ds.ny/2)
x1d = np.arange(0,ds.nx,1)
y1d = np.arange(0,ds.nz,1)
try:
    z = np.array(ds.zs[0,ymid,:ds.nx])/1000*10 #Div by 1000 to go to m and mult by 30 to match y dim
except:
    z = np.zeros((ds.nx))


#Create 2D arrays for plotting data
x2d = np.zeros((ds.nz, ds.nx))
y2d = np.zeros((ds.nz, ds.nx))

for i in range(ds.nz):
    x2d[i,:] = x1d
for j in range(ds.nx):
    y2d[:,j] = y1d+z[j]
        
#Variables for plotting
xmin = 0
xmax = ds.nx
xlen = xmax-xmin

zmin = 0
zmax = ds.nz

#Variables from output to build plot
land = ds.xland[0,:,:].values
try:
    zs = ds.zs[0,:,:].values
except:
    zs = np.zeros((ds.ny,ds.nx))
    

dbz_ymid = ds.dbz[:,:,ymid,:].values
print('hi')
cref = ds.cref.values



#Read in colormap and put in proper format
colors1 = np.array(nclcmaps.colors['amwg256'])#perc2_9lev'])
colors_int = colors1.astype(int)
colors = list(colors_int)
cmap_dth = nclcmaps.make_cmap(colors, bit=True)


#Read in colormap and put in proper format
colors1 = np.array(nclcmaps.colors['WhiteBlueGreenYellowRed'])#perc2_9lev'])
colors_int = colors1.astype(int)
colors = list(colors_int)
cmap_th = nclcmaps.make_cmap(colors, bit=True)

colors1 = np.array(nclcmaps.colors['prcp_1'])#perc2_9lev'])
colors_int = colors1.astype(int)
colors = list(colors_int)
cmap_precip = nclcmaps.make_cmap(colors, bit=True)


###############################################################################
###############################################################################

    
def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap
###############################################################################
###############################################################################


    shifted_cmapz = shiftedColorMap(cmap_dth, midpoint=1 - zlmax/(zlmax + abs(zlmin)), name='shifted')
    shifted_cmapx = shiftedColorMap(cmap_dth, midpoint=1 - xlmax/(xlmax + abs(xlmin)), name='shifted')
    

###############################################################################
###################  Backward Trajectory Plot  ################################
###############################################################################
#Use multiple processors to create images
def plotting(i): 


#for i in range(time_steps):   
#for i in range(0,1): 
    secs = (i*time_step_length)
    fig = plt.figure(num=None, figsize=(21,9), facecolor='w', edgecolor='k')
    
    
#####################  Subplot 1 (Trajectories in xy)  ########################
    
    #Plot characteristics
    ax = plt.subplot(221,aspect = 'equal')
    plt.title("Parcel Trajectories [elapsed time = %d seconds]"  % secs, fontsize = 18, y = 1.05) 
    plt.subplots_adjust(left=0.07, bottom=0.22, right=0.94, top=0.9, wspace=0.12, hspace=0.2)
    
    xlmin = 263
    xlmax = 273.01
    xlevels = np.arange(xlmin,xlmax, 0.05)
    xlevels_ticks = np.arange(xlmin,xlmax,1)
    xlevels_ticks_labels = np.arange(xlmin,xlmax, 1).astype(int)
    


    
    cmap = cmap_th
    trajectories_xy = plt.scatter(xpos[i,:,:], ypos[i,:,:], c = variable1[i,:,:], cmap = cmap_th, norm=matplotlib.colors.BoundaryNorm(xlevels,cmap.N), s = 2.5, zorder = 7)
    
    plt.xticks(np.arange(0,num_x,40000/hor_resolution))
    plt.yticks(np.arange(0,num_y,20000/hor_resolution))
    ax.set_xticklabels(np.arange(0,num_x*hor_resolution/1000,40), fontsize = 12)
    ax.set_yticklabels(np.arange(0,num_x*hor_resolution/1000,20), fontsize = 12)
    plt.ylabel('Distance (km)', fontsize = 16)
    plt.xlim([0,num_x])
    plt.ylim([0,num_y])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    #Water and land details
    levels_water = [1.5, 2.5]
    levels_terrain = [0,1.5]
    xland_plt = plt.contourf(land, levels_water, alpha = 1, colors = ('lightsteelblue'), zorder = 2)
    xland_plt = plt.contourf(land, levels_terrain, alpha = 1, colors = ('gainsboro'), zorder = 1)
    
    #Terrain
    terrain_levels = np.arange(-1, 3000.1, 200)
    terrain = plt.contourf(zs, terrain_levels, alpha = 1, cmap = cm.Greys, vmin = -800, vmax = 3000, zorder = 1)
    
    
    
#################### Subplot 2 (Trajectories in xz)  ##########################
    
    
    #Plot characteristics
    ax = plt.subplot(223)
    
    #Plot
    trajectories_xz = plt.scatter(xpos[i,:,:].flatten(), zpos_terrain[i,:,:].flatten(), c = variable1[i,:,:].flatten(), cmap = cmap_th, norm=matplotlib.colors.BoundaryNorm(xlevels,cmap.N), s = 2.5, zorder = 3)
    
    
    plt.grid(True, color = 'white', )
    plt.xticks(np.arange(0,num_x,40000/hor_resolution))
    ytick = np.arange(0,4001,500)
    plt.yticks(ytick)
    ax.set_xticklabels(np.arange(0,num_x*hor_resolution/1000,40), fontsize = 12)
    plt.xlim([0,num_x])
    plt.ylim([-50,ytick[-1]])
    ax.tick_params(axis='y', labelsize=14)
    plt.axvspan(0,num_x,color='gainsboro',lw=0)
    plt.xlabel('Distance (km)', fontsize = 16)
    plt.ylabel('Height (m)', fontsize = 18)

    
    #Plot Terrain
    x1d = np.arange(0,num_x,1)
    terrain = plt.plot(x1d, zs[ymid,:]-40, c = 'slategrey', linewidth = 4, zorder = 4)

    
    #Plot Lake
    lake = np.array(land[ymid,:])
    lake[lake == 1] = np.nan
    lake_plt = plt.plot(x1d, lake-40, c = 'blue', linewidth = 4, zorder = 5)
    
    #Colorbar  
    cbaxes = fig.add_axes([0.1, 0.09, 0.35, 0.035])         
    cbar = plt.colorbar(trajectories_xz, orientation='horizontal', cax = cbaxes, ticks = xlevels_ticks)
    cbar.ax.set_yticklabels(xlevels_ticks_labels)
    cbar.ax.tick_params(labelsize=12)
    plt.text(0.42, -1.7, 'Theta (K)', fontsize = 16)
    
    
    
    
################################# Subplot 3 (cref)  ###########################
    
    creflmin = 10
    creflmax = 55.01
    creflevels = np.arange(creflmin,creflmax, 0.05)
    creflevels_ticks = np.arange(creflmin,creflmax,5)
    creflevels_ticks_labels = np.arange(creflmin,creflmax, 5).astype(int)
    
    
    ax = plt.subplot(222,aspect = 'equal')
    plt.title("Reflectivity [elapsed time = %d seconds]"  % secs, fontsize = 18, y = 1.05) 

    
    #Plot reflectivity
    cref_plot = plt.contourf(cref[i+start_time_step,:,:], creflevels, cmap = cmap_precip, extend = 'max', alpha = 1,  zorder = 3)
    
    plt.xticks(np.arange(0,num_x,40000/hor_resolution))
    plt.yticks(np.arange(0,num_y,20000/hor_resolution))
    ax.set_xticklabels(np.arange(0,num_x*hor_resolution/1000,40), fontsize = 12)
    ax.set_yticklabels(np.arange(0,num_x*hor_resolution/1000,20), fontsize = 12)
    plt.xlim([0,num_x])
    plt.ylim([0,num_y])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

        
    #Water and land details
    levels_water = [1.5, 2.5]
    levels_terrain = [0,1.5]
    xland_plt = plt.contourf(land, levels_water, alpha = 1, colors = ('lightsteelblue'), zorder = 2)
    xland_plt = plt.contourf(land, levels_terrain, alpha = 1, colors = ('gainsboro'), zorder = 1)
    
    #Terrain
    terrain_levels = np.arange(-1, 3000.1, 200)
    terrain = plt.contourf(zs, terrain_levels, alpha = 1, cmap = cm.Greys, vmin = -800, vmax = 3000, zorder = 1)

    
    #Colorbar
    cbaxes = fig.add_axes([0.565, 0.09, 0.35, 0.035])
    cbar = plt.colorbar(cref_plot, orientation='horizontal', cax = cbaxes, ticks = creflevels_ticks)
    cbar.ax.set_yticklabels(creflevels_ticks_labels)
    cbar.ax.tick_params(labelsize=12)
    plt.text(0.49, -1.7, 'dBZ', fontsize = 16)
    
    
    #Labels
    sub_title = ['[Run: 15 $\mathregular{ms^{-1}}$]']
    props = dict(boxstyle='square', facecolor='white', alpha=1)
    ax.text(1.9, -2.08, sub_title[0], fontsize = 18, zorder = 5, transform=ax.transAxes)
    
    
    
    #################### Subplot 4 (ref in xz)  ##########################
    
    
    #Plot characteristics
    ax = plt.subplot(224)
    
    #Plot
    cref_plot = plt.contourf(x2d, y2d*vert_resolution, dbz_ymid[i+start_time_step,:,:], creflevels, cmap = cmap_precip, extend = 'max', alpha = 1,  zorder = 3)
    
    
    plt.grid(True, color = 'white', )
    plt.xticks(np.arange(0,num_x,40000/hor_resolution))
    ytick = np.arange(0,4001,500)
    plt.yticks(ytick)
    ax.set_xticklabels(np.arange(0,num_x*hor_resolution/1000,40), fontsize = 12)
    plt.xlim([0,num_x])
    plt.ylim([-50,ytick[-1]])
    ax.tick_params(axis='y', labelsize=14)
    plt.axvspan(0,num_x,color='gainsboro',lw=0)
    plt.xlabel('Distance (km)', fontsize = 16)

    
    
    #Plot Terrain
    x1d = np.arange(0,num_x,1)
    terrain = plt.plot(x1d, zs[0,:]-40, c = 'slategrey', linewidth = 4, zorder = 4)

    
    #Plot Lake
    lake = np.array(land[ymid,:])
    lake[lake == 1] = np.nan
    lake_plt = plt.plot(x1d, lake-40, c = 'blue', linewidth = 4, zorder = 5)
    

    
    plt.savefig("/uufs/chpc.utah.edu/common/home/u1013082/public_html/phd_plots/cm1/png_for_gifs/backward_trajectory_0000m_15ms_60s_%d" % np.int(zpos_init[0]) + "00m_lowres_%03d.png" % i, dpi = 65)
    plt.close(fig)
    print(i)

plt.switch_backend('Agg')
  
  
#run function to create images
pool = mp.Pool(processes = 20)
pool.map(plotting, range(time_steps))#number of processors
pool.close()
pool.join()
        



