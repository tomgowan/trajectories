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
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.patches as patches


#Read in with xarray
#ds = xr.open_dataset('/uufs/chpc.utah.edu/common/home/steenburgh-group8/tom/cm1/output/tug/cm1run_150m_15ms_0000m_90sec.nc')
#ds = xr.open_dataset('/uufs/chpc.utah.edu/common/home/steenburgh-group7/tom/cm1/output/netcdf_output/tug/cm1run_150m_15ms_0500m_90sec.nc')
#ds = xr.open_dataset('/uufs/chpc.utah.edu/common/home/steenburgh-group7/tom/cm1/output/netcdf_output/tug/cm1run_150m_15ms_2000m_90sec.nc')

ds = xr.open_dataset('/uufs/chpc.utah.edu/common/home/steenburgh-group7/tom/cm1/output/netcdf_output/tug/cm1run_150m_25ms_0000m_90sec.nc')
#ds = xr.open_dataset('//uufs/chpc.utah.edu/common/home/steenburgh-group8/tom/cm1/output/tug/cm1run_150m_25ms_0500m_90sec.nc')
#ds = xr.open_dataset('/uufs/chpc.utah.edu/common/home/steenburgh-group8/tom/cm1/output/tug/cm1run_150m_25ms_2000m_90sec.nc')



###############################################################################
###############################################################################   
##########################   Set up Trajectories ##############################
###############################################################################
###############################################################################



###############################################################################
##################### INFO TO CALCULATE SEEDS #################################
#############  These are variables the user changes  ##########################
###############################################################################
#Backward trajectories
num_seeds_z = 1 
time_steps = 100 #Number of time steps to run trajectories back
start_time_step = 210 #Starting time step
hor_resolution = 150 #[meters]
vert_resolution = 100 #[meters] (resolution away from terrain)
time_step_length = 90.0 #[seconds]
left = 1495 #Starting xpos of traj and left side of data chunks
near = 160 #close side of ypos for traj
far = ds.ny-160 #far side of ypos for traj
num_seeds_y = far-near
bottom = 1 #starting zpos
top = 40 #top of data chucks


#Dimension size variables
num_x = ds.nx - left
num_y = ds.ny
num_z = top

x = np.arange(0,num_x,1)
y = np.arange(0,num_y,1)
z = np.arange(0,num_z,1)
###############################################################################
###############################################################################


#Create arrays for location 
xpos = np.zeros((time_steps, num_seeds_z, num_seeds_y))
ypos = np.zeros((time_steps, num_seeds_z, num_seeds_y))
zpos = np.zeros((time_steps, num_seeds_z, num_seeds_y))
zpos_grid = np.zeros((time_steps, num_seeds_z, num_seeds_y))
change_zs = np.zeros((time_steps, num_seeds_z, num_seeds_y))


#Variable
var_name1 = 'th'
variable1 = np.zeros((time_steps, num_seeds_z, num_seeds_y))


#Budgets
var_name2 = 'ptb_mp'#Coded to be a budget variable
variable2 = np.zeros((time_steps, num_seeds_z, num_seeds_y))

var_name6 = 'ptb_vturb'#Coded to be a budget variable
variable6 = np.zeros((time_steps, num_seeds_z, num_seeds_y))

var_name8 = 'ptb_rad'#Coded to be a budget variable
variable8 = np.zeros((time_steps, num_seeds_z, num_seeds_y))


#Tendancies
var_name11 = 'tt_cond'#Coded to be a budget variable
variable11 = np.zeros((time_steps, num_seeds_z, num_seeds_y))

var_name15 = 'tt_subl'#Coded to be a budget variable
variable15 = np.zeros((time_steps, num_seeds_z, num_seeds_y))




###############################################################################
##################### INITIAL LOCATION OF SEEDS ###############################
#########  (These may also be important to the user)  #########################
###############################################################################
#Choose starting point 
#xpos
xpos[0,:,:] = 5#Seeds at every x-gridpoint
#For each variable I only pull out the 1495 x-gridpoint on so this 
#represent starting at the 1500 gridpoint

#ypos   
for i in range(0,far-near):
    ypos[0,:,i] = i + near #Seeds all at middle y-gridpoint

#zpos
zpos[0,:,:] = bottom

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
 
    
#check if run has terrain (if not, zs is zero and zpos should not be affected)
try:
    zh = np.array(ds.zh[0,:top,:,left:])
    #Add terrain height
    #Create list of current coordinates for terrain addition
    
    xloc = np.array(xpos[0,:,:]).flatten()
    yloc = np.array(ypos[0,:,:]).flatten()
    zloc = np.array(zpos[0,:,:]).flatten()
    coor_terrain = []
    for i in range(len(xloc)):
        coor_terrain.append((zloc[i], yloc[i], xloc[i]))
    
    zpos_terrain = np.array(zpos)
    zpos_terrain[0,:,:] = np.reshape(interpolate.interpn((z,y,x), zh, coor_terrain, method='linear', bounds_error=False, fill_value= 0), (num_seeds_z, num_seeds_y))
    
    #This gets us the grid spacing for the vertical grid
    z_grid = zh[1:,:,:]-zh[:-1,:,:]


except:
    zh = np.zeros((top,ds.ny, ds.nx-left))
    zpos_terrain = np.copy(zpos)*vert_resolution
    z_grid = np.ones((top-1,ds.ny, ds.nx-left))*vert_resolution
    


###############################################################################
###############################################################################   
##########################  Calculate Trajectories ############################
###############################################################################
###############################################################################

#Loop over all time spteps and compute trajectory
for t in range(time_steps-1):
    
    print t
    start = time.time() #Timer
    
    #Get model data
    u = ds.uinterp[start_time_step+t,:top,:,left:].values
    v = ds.vinterp[start_time_step+t,:top,:,left:].values
    w = ds.winterp[start_time_step+t,:top,:,left:].values
    var1 = getattr(ds,var_name1)[start_time_step+t,:top,:,left:].values
    var2 = getattr(ds,var_name2)[start_time_step+t,:top,:,left:].values
    var6 = getattr(ds,var_name6)[start_time_step+t,:top,:,left:].values
    var8 = getattr(ds,var_name8)[start_time_step+t,:top,:,left:].values
    var11 = getattr(ds,var_name11)[start_time_step+t,:top,:,left:].values
    var15 = getattr(ds,var_name15)[start_time_step+t,:top,:,left:].values

                                




    #####################   Get x,y,z for variables ###########################
    ##  Because the model output is terrain following, the change in surface ##
    ##  height at each timestep must be added to the height (zpos) of the    ##
    ##  parcels.                                                             ##
    ###########################################################################
        
    #Get surface height grid
    try:
        zs = np.array(ds.zs[0,:,left:])
    except:
        zs = np.zeros((ds.ny, ds.nx-left))  


    #Get previous and current x and y positions to calc change in surface height
    if t > 0:
        
        #x and y for one time step back
        xloc = np.copy(xpos[t-1,:,:]).flatten()
        yloc = np.copy(ypos[t-1,:,:]).flatten()
        coor_xminus1 = []
        for i in range(len(xloc)):
            coor_xminus1.append((yloc[i], xloc[i]))
            
            
        #xand y for current time step
        xloc = np.copy(xpos[t,:,:]).flatten()
        yloc = np.copy(ypos[t,:,:]).flatten()
        coor_x = []
        for i in range(len(xloc)):
            coor_x.append((yloc[i], xloc[i]))
            
        #Calc surface height at each time
        xminus1_h = interpolate.interpn((y,x), zs, coor_xminus1, method='linear', bounds_error=False, fill_value= 0)#/vert_resolution
        x_h = interpolate.interpn((y,x), zs, coor_x, method='linear', bounds_error=False, fill_value= 0)#/vert_resolution
        
        #Calc change in surface height during last timestep, which will be added to zpos
        change_zs[t,:,:] =  np.reshape(x_h - xminus1_h, (num_seeds_z, num_seeds_y))
    
    
    #Get get x, y, and z positions to calc parcel movement
    xloc = np.copy(xpos[t,:,:]).flatten()
    yloc = np.copy(ypos[t,:,:]).flatten()
    zloc = np.copy(zpos[t,:,:]).flatten()
    coor_var = []
    for i in range(len(xloc)):
        coor_var.append((zloc[i], yloc[i], xloc[i])) 
    
    
    #####################   Calc new xpos #####################################
    xpos[t+1,:,:] = xpos[t,:,:] + np.reshape(interpolate.interpn((z,y,x), u, coor_var, method='linear', bounds_error=False, fill_value=np.nan)*time_step_length/hor_resolution, (num_seeds_z, num_seeds_y))

    #####################   Calc new ypos #####################################
    ypos[t+1,:,:]  = ypos[t,:,:] + np.reshape(interpolate.interpn((z,y,x), v, coor_var, method='linear', bounds_error=False, fill_value=np.nan)*time_step_length/hor_resolution, (num_seeds_z, num_seeds_y))

    #####################   Calc new zpos #####################################
    #zpos grid spacing
    zpos_grid[t,:,:] = np.reshape(interpolate.interpn((z[:-1],y,x), z_grid, coor_var, method='linear', bounds_error=False, fill_value= 0), (num_seeds_z, num_seeds_y))
    #terrain-following (includes change in surface height)
    zpos[t+1,:,:]  = zpos[t,:,:] + change_zs[t,:,:]/zpos_grid[t,:,:] + np.reshape(interpolate.interpn((z,y,x), w, coor_var, method='linear', bounds_error=False, fill_value= 0), (num_seeds_z, num_seeds_y))*time_step_length/zpos_grid[t,:,:]
    #terrain-height coordinates
    zpos_terrain[t+1,:,:]  = zpos_terrain[t,:,:] + np.reshape(interpolate.interpn((z,y,x), w, coor_var, method='linear', bounds_error=False, fill_value= 0)*time_step_length, (num_seeds_z, num_seeds_y))

    zpos = zpos.clip(min=0) #Prevents parcel from going into the ground
    
        
    #Prevents parcel from going into ground
    xloc = np.copy(xpos[t,:,:]).flatten()
    yloc = np.copy(ypos[t,:,:]).flatten()
    coor_terrain = []
    for i in range(len(xloc)):
        coor_terrain.append((yloc[i], xloc[i])) 
    surface_height = np.reshape(interpolate.interpn((y,x), zs, coor_terrain, method='linear', bounds_error=False, fill_value= 0)/vert_resolution, (num_seeds_z, num_seeds_y))
    zpos_terrain[t,:,:] = zpos_terrain[t,:,:].clip(min=surface_height) 
    
    
    
    #Variables
    variable1[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var1, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan), (num_seeds_z, num_seeds_y))
    
    #Budget variable
    variable2[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var2, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_y))
    variable6[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var6, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_y))
    variable8[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var8, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_y))
    
    
    #Tendancy variable
    variable11[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var11, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_y))
    variable15[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var15, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_y))




    #Timer
    stop = time.time()
    print(stop-start)


##############  Get variable data for final time step #########################
t = time_steps-1

#Get get x, y, and z positions to calc parcel movement
xloc = np.copy(xpos[t,:,:]).flatten()
yloc = np.copy(ypos[t,:,:]).flatten()
zloc = np.copy(zpos[t,:,:]).flatten()
coor_var = []
for i in range(len(xloc)):
    coor_var.append((zloc[i], yloc[i], xloc[i])) 


#Color variable (theta)
variable1[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var1, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan), (num_seeds_z, num_seeds_y))
variable2[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var2, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_y))
variable6[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var6, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_y))
variable8[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var8, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_y))

#Tendancy variable
variable11[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var11, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_y))
variable15[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var15, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_y))
#%%
    
#Save arrays
np.save('xpos_terrain_zoom_traj_disp_25ms_0000m_%dtime_steps' %time_steps, xpos)
np.save('ypos_terrain_zoom_traj_disp_25ms_0000m_%dtime_steps' %time_steps, ypos)
np.save('zpos_terrain_zoom_traj_disp_25ms_0000m_%dtime_steps' %time_steps, zpos_terrain)
np.save('%s_terrain_zoom_traj_disp_25ms_0000m' %var_name1 + '_%dtime_steps' %time_steps, variable1)
np.save('%s_terrain_zoom_traj_disp_25ms_0000m' %var_name2 + '_%dtime_steps' %time_steps, variable2)
np.save('%s_terrain_zoom_traj_disp_25ms_0000m' %var_name6 + '_%dtime_steps' %time_steps, variable6)
np.save('%s_terrain_zoom_traj_disp_25ms_0000m' %var_name8 + '_%dtime_steps' %time_steps, variable8)
np.save('%s_terrain_zoom_traj_disp_25ms_0000m' %var_name11 + '_%dtime_steps' %time_steps, variable11)
np.save('%s_terrain_zoom_traj_disp_25ms_0000m' %var_name15 + '_%dtime_steps' %time_steps, variable15)


#%%


    
    
    
    


