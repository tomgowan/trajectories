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
ds = xr.open_dataset('/uufs/chpc.utah.edu/common/home/steenburgh-group8/tom/cm1/output/tug/cm1run_150m_15ms_0000m_90sec.nc')
#ds = xr.open_dataset('//uufs/chpc.utah.edu/common/home/steenburgh-group8/tom/cm1/output/tug/cm1run_150m_25ms_0500m_90sec.nc')
#ds = xr.open_dataset('/uufs/chpc.utah.edu/common/home/steenburgh-group8/tom/cm1/output/tug/cm1run_150m_25ms_2000m_90sec.nc')

#ds = xr.open_dataset('/uufs/chpc.utah.edu/common/home/steenburgh-group7/tom/cm1/output/netcdf_output/tug/cm1run_150m_15ms_0000m_60sec.nc')
#ds = xr.open_dataset('/uufs/chpc.utah.edu/common/home/steenburgh-group7/tom/cm1/output/netcdf_output/tug/cm1run_150m_15ms_0500m_90sec.nc')
#ds = xr.open_dataset('/uufs/chpc.utah.edu/common/home/steenburgh-group7/tom/cm1/output/netcdf_output/tug/cm1run_150m_15ms_2000m_90sec.nc')

#%%

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
num_seeds_z = 151 #Up to 5000m (3 seeds every vertical grid point)
num_seeds_x = ds.nx #One for each x gridpoint
time_steps = 250 #Number of time steps to run trajectories back
start_time_step = 310 #Starting time step
hor_resolution = 150 #[meters]
vert_resolution = 100 #[meters] (resolution away from terrain)
time_step_length = 90.0 #[seconds]
###############################################################################
###############################################################################


#Create arrays for location a(dim =d variable color of parcels
xpos = np.zeros((time_steps, num_seeds_z, num_seeds_x))
ypos = np.zeros((time_steps, num_seeds_z, num_seeds_x))
zpos = np.zeros((time_steps, num_seeds_z, num_seeds_x))
zpos_grid = np.zeros((time_steps, num_seeds_z, num_seeds_x))
change_zs = np.zeros((time_steps, num_seeds_z, num_seeds_x))


#Variable
var_name1 = 'th'
variable1 = np.zeros((time_steps, num_seeds_z, num_seeds_x))


#Budgets
var_name2 = 'ptb_mp'#Coded to be a budget variable
variable2 = np.zeros((time_steps, num_seeds_z, num_seeds_x))

#var_name3 = 'ptb_hadv'#Coded to be a budget variable
#variable3 = np.zeros((time_steps, num_seeds_z, num_seeds_x))
#
#var_name4 = 'ptb_vadv'#Coded to be a budget variable
#variable4 = np.zeros((time_steps, num_seeds_z, num_seeds_x))
#
#var_name5 = 'ptb_hturb'#Coded to be a budget variable
#variable5 = np.zeros((time_steps, num_seeds_z, num_seeds_x))

var_name6 = 'ptb_vturb'#Coded to be a budget variable
variable6 = np.zeros((time_steps, num_seeds_z, num_seeds_x))

#var_name7 = 'ptb_rdamp'#Coded to be a budget variable
#variable7 = np.zeros((time_steps, num_seeds_z, num_seeds_x))

var_name8 = 'ptb_rad'#Coded to be a budget variable
variable8 = np.zeros((time_steps, num_seeds_z, num_seeds_x))

#var_name9 = 'ptb_div'#Coded to be a budget variable
#variable9 = np.zeros((time_steps, num_seeds_z, num_seeds_x))
#
#var_name10 = 'ptb_diss'#Coded to be a budget variable
#variable10 = np.zeros((time_steps, num_seeds_z, num_seeds_x))


#Tendancies
var_name11 = 'tt_cond'#Coded to be a budget variable
variable11 = np.zeros((time_steps, num_seeds_z, num_seeds_x))

#var_name12 = 'tt_evac'#Coded to be a budget variable
#variable12 = np.zeros((time_steps, num_seeds_z, num_seeds_x))
#
#var_name13 = 'tt_evar'#Coded to be a budget variable
#variable13 = np.zeros((time_steps, num_seeds_z, num_seeds_x))
#
#var_name14 = 'tt_dep'#Coded to be a budget variable
#variable14 = np.zeros((time_steps, num_seeds_z, num_seeds_x))

var_name15 = 'tt_subl'#Coded to be a budget variable
variable15 = np.zeros((time_steps, num_seeds_z, num_seeds_x))

#var_name16 = 'tt_melt'#Coded to be a budget variable
#variable16 = np.zeros((time_steps, num_seeds_z, num_seeds_x))
#
#var_name17 = 'tt_frz'#Coded to be a budget variable
#variable17 = np.zeros((time_steps, num_seeds_z, num_seeds_x))


###############################################################################
##################### INITIAL LOCATION OF SEEDS ###############################
#########  (These may also be important to the user)  #########################
###############################################################################
#Choose starting point 
#xpos
for i in range(0, num_seeds_x):
    xpos[0,:,i] = i#Seeds at every x-gridpoint

#ypos   
ymid = np.int(ds.ny/2)
ypos[0,:,:] = ymid #Seeds all at middle y-gridpoint

#zpos
for i in range(num_seeds_z):
    zpos[0,i,:] = i/3.0 #Surface to top of desired regions (3 per z-gridpoint)

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
    zh = np.array(ds.zh[0,:,:,:])
    #Add terrain height
    #Create list of current coordinates for terrain addition
    
    xloc = np.array(xpos[0,:,:]).flatten()
    yloc = np.array(ypos[0,:,:]).flatten()
    zloc = np.array(zpos[0,:,:]).flatten()
    coor_terrain = []
    for i in range(len(xloc)):
        coor_terrain.append((zloc[i], yloc[i], xloc[i]))
    
    zpos_terrain = np.array(zpos)
    zpos_terrain[0,:,:] = np.reshape(interpolate.interpn((z,y,x), zh, coor_terrain, method='linear', bounds_error=False, fill_value= 0), (num_seeds_z, num_seeds_x))
    
    #This gets us the grid spacing for the vertical grid
    z_grid = zh[1:,:,:]-zh[:-1,:,:]


except:
    zh = np.zeros((ds.nz,ds.ny, ds.nx))
    zpos_terrain = np.copy(zpos)*vert_resolution
    z_grid = np.ones((ds.nz-1,ds.ny, ds.nx))*vert_resolution
    


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
    u = ds.uinterp[start_time_step-t,:,:,:].values
    v = ds.vinterp[start_time_step-t,:,:,:].values
    w = ds.winterp[start_time_step-t,:,:,:].values
    var1 = getattr(ds,var_name1)[start_time_step-t,:,:,:].values
    var2 = getattr(ds,var_name2)[start_time_step-t,:,:,:].values
#    var3 = getattr(ds,var_name3)[start_time_step-t,:,:,:].values
#    var4 = getattr(ds,var_name4)[start_time_step-t,:,:,:].values
#    var5 = getattr(ds,var_name5)[start_time_step-t,:,:,:].values
    var6 = getattr(ds,var_name6)[start_time_step-t,:,:,:].values
#    var7 = getattr(ds,var_name7)[start_time_step-t,:,:,:].values
    var8 = getattr(ds,var_name8)[start_time_step-t,:,:,:].values
#    var9 = getattr(ds,var_name9)[start_time_step-t,:,:,:].values
#    var10 = getattr(ds,var_name10)[start_time_step-t,:,:,:].values
    var11 = getattr(ds,var_name11)[start_time_step-t,:,:,:].values
#    var12 = getattr(ds,var_name12)[start_time_step-t,:,:,:].values
#    var13 = getattr(ds,var_name13)[start_time_step-t,:,:,:].values
#    var14 = getattr(ds,var_name14)[start_time_step-t,:,:,:].values
    var15 = getattr(ds,var_name15)[start_time_step-t,:,:,:].values
#    var16 = getattr(ds,var_name16)[start_time_step-t,:,:,:].values
#    var17 = getattr(ds,var_name17)[start_time_step-t,:,:,:].values
                                




    #####################   Get x,y,z for variables ###########################
    ##  Because the model output is terrain following, the change in surface ##
    ##  height at each timestep must be added to the height (zpos) of the    ##
    ##  parcels.                                                             ##
    ###########################################################################
        
    #Get surface height grid
    try:
        zs = np.array(ds.zs[0,:,:])
    except:
        zs = np.zeros((ds.ny, ds.nx))  


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
        change_zs[t,:,:] =  np.reshape(x_h - xminus1_h, (num_seeds_z, num_seeds_x))
    
    
    #Get get x, y, and z positions to calc parcel movement
    xloc = np.copy(xpos[t,:,:]).flatten()
    yloc = np.copy(ypos[t,:,:]).flatten()
    zloc = np.copy(zpos[t,:,:]).flatten()
    coor_var = []
    for i in range(len(xloc)):
        coor_var.append((zloc[i], yloc[i], xloc[i])) 
    
    
    #####################   Calc new xpos #####################################
    xpos[t+1,:,:] = xpos[t,:,:] - np.reshape(interpolate.interpn((z,y,x), u, coor_var, method='linear', bounds_error=False, fill_value=np.nan)*time_step_length/hor_resolution, (num_seeds_z, num_seeds_x))

    #####################   Calc new ypos #####################################
    ypos[t+1,:,:]  = ypos[t,:,:] - np.reshape(interpolate.interpn((z,y,x), v, coor_var, method='linear', bounds_error=False, fill_value=np.nan)*time_step_length/hor_resolution, (num_seeds_z, num_seeds_x))

    #####################   Calc new zpos #####################################
    #zpos grid spacing
    zpos_grid[t,:,:] = np.reshape(interpolate.interpn((z[:-1],y,x), z_grid, coor_var, method='linear', bounds_error=False, fill_value= 0), (num_seeds_z, num_seeds_x))
    #terrain-following (includes change in surface height)
    zpos[t+1,:,:]  = zpos[t,:,:] - change_zs[t,:,:]/zpos_grid[t,:,:] - np.reshape(interpolate.interpn((z,y,x), w, coor_var, method='linear', bounds_error=False, fill_value= 0), (num_seeds_z, num_seeds_x))*time_step_length/zpos_grid[t,:,:]
    #terrain-height coordinates
    zpos_terrain[t+1,:,:]  = zpos_terrain[t,:,:] - np.reshape(interpolate.interpn((z,y,x), w, coor_var, method='linear', bounds_error=False, fill_value= 0)*time_step_length, (num_seeds_z, num_seeds_x))

    zpos = zpos.clip(min=0) #Prevents parcel from going into the ground
    
        
    #Prevents parcel from going into ground
    xloc = np.copy(xpos[t,:,:]).flatten()
    yloc = np.copy(ypos[t,:,:]).flatten()
    coor_terrain = []
    for i in range(len(xloc)):
        coor_terrain.append((yloc[i], xloc[i])) 
    surface_height = np.reshape(interpolate.interpn((y,x), zs, coor_terrain, method='linear', bounds_error=False, fill_value= 0)/vert_resolution, (num_seeds_z, num_seeds_x))
    zpos_terrain[t,:,:] = zpos_terrain[t,:,:].clip(min=surface_height) 
    
    
    
    #Variables
    variable1[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var1, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan), (num_seeds_z, num_seeds_x))
    
    #Budget variable
    variable2[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var2, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))
#    variable3[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var3, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))
#    variable4[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var4, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))
#    variable5[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var5, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))
    variable6[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var6, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))
#    variable7[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var7, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))
    variable8[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var8, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))
#    variable9[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var9, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))
#    variable10[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var10, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))
    
    
    #Tendancy variable
    variable11[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var11, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))
#    variable12[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var12, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))
#    variable13[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var13, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))
#    variable14[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var14, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))
    variable15[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var15, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))
#    variable16[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var16, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))
#    variable17[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var17, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))



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
variable1[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var1, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan), (num_seeds_z, num_seeds_x))
variable2[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var2, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))
#variable3[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var3, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))
#variable4[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var4, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))
#variable5[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var5, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))
variable6[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var6, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))
#variable7[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var7, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))
variable8[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var8, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))
#variable9[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var9, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))
#variable10[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var10, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))

#Tendancy variable
variable11[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var11, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))
#variable12[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var12, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))
#variable13[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var13, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))
#variable14[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var14, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))
variable15[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var15, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))
#variable16[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var16, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))
#variable17[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var17, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))
#%%
    
#Save arrays
np.save('xpos_traj_disp_15ms_0000m_%dtime_steps' %time_steps, xpos)
np.save('ypos_traj_disp_15ms_0000m_%dtime_steps' %time_steps, ypos)
np.save('zpos_traj_disp_15ms_0000m_%dtime_steps' %time_steps, zpos_terrain)
np.save('%s_traj_disp_15ms_0000m' %var_name1 + '_%dtime_steps' %time_steps, variable1)
np.save('%s_traj_disp_15ms_0000m' %var_name2 + '_%dtime_steps' %time_steps, variable2)
#np.save('%s_traj_disp_15ms_0500m' %var_name3 + '_%dtime_steps' %time_steps, variable3)
#np.save('%s_traj_disp_15ms_0500m' %var_name4 + '_%dtime_steps' %time_steps, variable4)
#np.save('%s_traj_disp_15ms_0500m' %var_name5 + '_%dtime_steps' %time_steps, variable5)
np.save('%s_traj_disp_15ms_0000m' %var_name6 + '_%dtime_steps' %time_steps, variable6)
#np.save('%s_traj_disp_15ms_0500m' %var_name7 + '_%dtime_steps' %time_steps, variable7)
np.save('%s_traj_disp_15ms_0000m' %var_name8 + '_%dtime_steps' %time_steps, variable8)
#np.save('%s_traj_disp_15ms_0500m' %var_name9 + '_%dtime_steps' %time_steps, variable9)
#np.save('%s_traj_disp_15ms_0500m' %var_name10 + '_%dtime_steps' %time_steps, variable10)
np.save('%s_traj_disp_15ms_0000m' %var_name11 + '_%dtime_steps' %time_steps, variable11)
#np.save('%s_traj_disp_15ms_0500m' %var_name12 + '_%dtime_steps' %time_steps, variable12)
#np.save('%s_traj_disp_15ms_0500m' %var_name13 + '_%dtime_steps' %time_steps, variable13)
#np.save('%s_traj_disp_15ms_0500m' %var_name14 + '_%dtime_steps' %time_steps, variable14)
np.save('%s_traj_disp_15ms_0000m' %var_name15 + '_%dtime_steps' %time_steps, variable15)
#np.save('%s_traj_disp_15ms_0500m' %var_name16 + '_%dtime_steps' %time_steps, variable16)
#np.save('%s_traj_disp_15ms_0500m' %var_name17 + '_%dtime_steps' %time_steps, variable17)



#%%
xpos = np.load('xpos_traj_disp_25ms_0000m_250time_steps.npy')
ypos = np.load('ypos_traj_disp_25ms_0000m_250time_steps.npy')
zpos_terrain = np.load('zpos_traj_disp_25ms_0000m_250time_steps.npy')
variable1 = np.load('%s_traj_disp_25ms_0000m_250time_steps.npy' % var_name1)
variable2 = np.load('%s_traj_disp_25ms_0000m_250time_steps.npy' %var_name2)
#variable3 = np.load('%s_traj_disp_15ms_0500m_300time_steps.npy' %var_name3)
#variable4 = np.load('%s_traj_disp_15ms_0500m_300time_steps.npy' %var_name4)
#variable5 = np.load('%s_traj_disp_15ms_0500m_300time_steps.npy' %var_name5)
variable6 = np.load('%s_traj_disp_25ms_0000m_250time_steps.npy' %var_name6)
#variable7 = np.load('%s_traj_disp_15ms_0500m_300time_steps.npy' %var_name7)
variable8 = np.load('%s_traj_disp_25ms_0000m_250time_steps.npy' %var_name8)
#variable9 = np.load('%s_traj_disp_15ms_0500m_300time_steps.npy' %var_name9)
#variable10 = np.load('%s_traj_disp_15ms_0500m_300time_steps.npy' %var_name10)
variable11 = np.load('%s_traj_disp_25ms_0000m_250time_steps.npy' %var_name11)
#variable12 = np.load('%s_traj_disp_15ms_0500m_300time_steps.npy' %var_name12)
#variable13 = np.load('%s_traj_disp_15ms_0500m_300time_steps.npy' %var_name13)
#variable14 = np.load('%s_traj_disp_15ms_0500m_300time_steps.npy' %var_name14)
variable15 = np.load('%s_traj_disp_25ms_0000m_250time_steps.npy' %var_name15)
#variable16 = np.load('%s_traj_disp_15ms_0500m_300time_steps.npy' %var_name16)
#variable17 = np.load('%s_traj_disp_15ms_0500m_300time_steps.npy' %var_name17)

#%%
###############################################################################
######################  Calculate displacement  ###############################
###############################################################################

##### ENTER VALUE FOR NUMBER OF TIME STEPS BACK TO CALCULATE TIME-MEAN AND PLOT
ts_plot = 90

z_disp = np.zeros((time_steps-2, num_seeds_z, num_seeds_x))
x_disp = np.zeros((time_steps-2, num_seeds_z, num_seeds_x))
y_disp = np.zeros((time_steps-2, num_seeds_z, num_seeds_x))
var1_disp = np.zeros((time_steps-2, num_seeds_z, num_seeds_x))
var2_disp = np.zeros((time_steps-2, num_seeds_z, num_seeds_x))
#var3_disp = np.zeros((time_steps-2, num_seeds_z, num_seeds_x))
#var4_disp = np.zeros((time_steps-2, num_seeds_z, num_seeds_x))
#var5_disp = np.zeros((time_steps-2, num_seeds_z, num_seeds_x))
var6_disp = np.zeros((time_steps-2, num_seeds_z, num_seeds_x))
#var7_disp = np.zeros((time_steps-2, num_seeds_z, num_seeds_x))
var8_disp = np.zeros((time_steps-2, num_seeds_z, num_seeds_x))
#var9_disp = np.zeros((time_steps-2, num_seeds_z, num_seeds_x))
#var10_disp = np.zeros((time_steps-2, num_seeds_z, num_seeds_x))
var11_disp = np.zeros((time_steps-2, num_seeds_z, num_seeds_x))
#var12_disp = np.zeros((time_steps-2, num_seeds_z, num_seeds_x))
#var13_disp = np.zeros((time_steps-2, num_seeds_z, num_seeds_x))
#var14_disp = np.zeros((time_steps-2, num_seeds_z, num_seeds_x))
var15_disp = np.zeros((time_steps-2, num_seeds_z, num_seeds_x))
#var16_disp = np.zeros((time_steps-2, num_seeds_z, num_seeds_x))
#var17_disp = np.zeros((time_steps-2, num_seeds_z, num_seeds_x))


for t in range(time_steps-2):
    z_disp[t,:,:] = zpos_terrain[0,:,:] - zpos_terrain[t+1,:,:]
    x_disp[t,:,:] = xpos[0,:,:] - xpos[t+1,:,:]
    y_disp[t,:,:] = ypos[0,:,:] - ypos[t+1,:,:]
    var1_disp[t,:,:] = variable1[0,:,:] - variable1[t+1,:,:]

    #Budget variables
    var2_disp[t,:,:] =  np.sum(variable2[:t+1,:,:], axis = 0)
#    var3_disp[t,:,:] =  np.sum(variable3[:t+1,:,:], axis = 0)
#    var4_disp[t,:,:] =  np.sum(variable4[:t+1,:,:], axis = 0)
#    var5_disp[t,:,:] =  np.sum(variable5[:t+1,:,:], axis = 0)
    var6_disp[t,:,:] =  np.sum(variable6[:t+1,:,:], axis = 0)
#    var7_disp[t,:,:] =  np.sum(variable7[:t+1,:,:], axis = 0)
    var8_disp[t,:,:] =  np.sum(variable8[:t+1,:,:], axis = 0)
#    var9_disp[t,:,:] =  np.sum(variable9[:t+1,:,:], axis = 0)
#    var10_disp[t,:,:] =  np.sum(variable10[:t+1,:,:], axis = 0)
    
    #Tendancy variables
    var11_disp[t,:,:] =  np.sum(variable11[:t+1,:,:], axis = 0)
#    var12_disp[t,:,:] =  np.sum(variable12[:t+1,:,:], axis = 0)
#    var13_disp[t,:,:] =  np.sum(variable13[:t+1,:,:], axis = 0)
#    var14_disp[t,:,:] =  np.sum(variable14[:t+1,:,:], axis = 0)
    var15_disp[t,:,:] =  np.sum(variable15[:t+1,:,:], axis = 0)
#    var16_disp[t,:,:] =  np.sum(variable16[:t+1,:,:], axis = 0)
#    var17_disp[t,:,:] =  np.sum(variable17[:t+1,:,:], axis = 0)




mean_z_disp = np.zeros((num_seeds_z, num_seeds_x))
mean_z_disp = np.mean(z_disp[:ts_plot,:,:], axis = 0)

mean_x_disp = np.zeros((num_seeds_z, num_seeds_x))
mean_x_disp = np.mean(x_disp[:ts_plot,:,:], axis = 0)

mean_y_disp = np.zeros((num_seeds_z, num_seeds_x))
mean_y_disp = np.mean(y_disp[:ts_plot,:,:], axis = 0)

#Variable
mean_var1_disp = np.zeros((num_seeds_z, num_seeds_x))
mean_var1_disp = np.mean(var1_disp[:ts_plot,:,:], axis = 0)


##Budget Variables
mean_var2_disp = np.zeros((num_seeds_z, num_seeds_x))
mean_var2_disp = np.mean(var2_disp[:ts_plot,:,:], axis = 0)

#mean_var3_disp = np.zeros((num_seeds_z, num_seeds_x))
#mean_var3_disp = np.mean(var3_disp[:ts_plot,:,:], axis = 0)
#
#mean_var4_disp = np.zeros((num_seeds_z, num_seeds_x))
#mean_var4_disp = np.mean(var4_disp[:ts_plot,:,:], axis = 0)
#
#mean_var5_disp = np.zeros((num_seeds_z, num_seeds_x))
#mean_var5_disp = np.mean(var5_disp[:ts_plot,:,:], axis = 0)

mean_var6_disp = np.zeros((num_seeds_z, num_seeds_x))
mean_var6_disp = np.mean(var6_disp[:ts_plot,:,:], axis = 0)

#mean_var7_disp = np.zeros((num_seeds_z, num_seeds_x))
#mean_var7_disp = np.mean(var7_disp[:ts_plot,:,:], axis = 0)

mean_var8_disp = np.zeros((num_seeds_z, num_seeds_x))
mean_var8_disp = np.mean(var8_disp[:ts_plot,:,:], axis = 0)
#
#mean_var9_disp = np.zeros((num_seeds_z, num_seeds_x))
#mean_var9_disp = np.mean(var9_disp[:ts_plot,:,:], axis = 0)
#
#mean_var10_disp = np.zeros((num_seeds_z, num_seeds_x))
#mean_var10_disp = np.mean(var10_disp[:ts_plot,:,:], axis = 0)


#Tendencay Variable
mean_var11_disp = np.zeros((num_seeds_z, num_seeds_x))
mean_var11_disp = np.mean(var11_disp[:ts_plot,:,:], axis = 0)

#mean_var12_disp = np.zeros((num_seeds_z, num_seeds_x))
#mean_var12_disp = np.mean(var12_disp[:ts_plot,:,:], axis = 0)
#
#mean_var13_disp = np.zeros((num_seeds_z, num_seeds_x))
#mean_var13_disp = np.mean(var13_disp[:ts_plot,:,:], axis = 0)
#
#mean_var14_disp = np.zeros((num_seeds_z, num_seeds_x))
#mean_var14_disp = np.mean(var14_disp[:ts_plot,:,:], axis = 0)

mean_var15_disp = np.zeros((num_seeds_z, num_seeds_x))
mean_var15_disp = np.mean(var15_disp[:ts_plot,:,:], axis = 0)
#
#mean_var16_disp = np.zeros((num_seeds_z, num_seeds_x))
#mean_var16_disp = np.mean(var16_disp[:ts_plot,:,:], axis = 0)
#
#mean_var17_disp = np.zeros((num_seeds_z, num_seeds_x))
#mean_var17_disp = np.mean(var17_disp[:ts_plot,:,:], axis = 0)







#%%

###############################################################################
#############################   PLOTS   #######################################
###############################################################################



############## Set ncl_cmap as the colormap you want ##########################

### 1) In order for this to work the files "nclcmaps.py" and "__init__.py"
### must be present in the dirctory.
### 2) You must "import nclcmaps"
### 3) The path to nclcmaps.py must be added to tools -> PYTHONPATH manager in SPyder
### 4) Then click "Update module names list" in tolls in Spyder and restart Spyder
                
## The steps above describe the general steps for adding "non-built in" modules to Spyder

###############################################################################
###############################################################################


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



##########  Create Grid ########
### The code below makes the data terrain following 
x1d = np.arange(0,num_seeds_x,1)
y1d = np.arange(0,num_seeds_z,1)
try:
    z = np.array(ds.zs[0,ymid,:num_seeds_x])/1000*30 #Div by 1000 to go to m and mult by 30 to match y dim
except:
    z = np.zeros((num_seeds_x))

#Create 2D arrays for plotting data
x2d = np.zeros((num_seeds_z, num_seeds_x))
y2d = np.zeros((num_seeds_z, num_seeds_x))

for i in range(num_seeds_z):
    x2d[i,:] = x1d
for j in range(num_seeds_x):
    y2d[:,j] = y1d+z[j]
        
#Variables for plotting
xmin = 0
xmax = num_seeds_x
xlen = xmax-xmin

zmin = 0
zmax = num_seeds_z






###############################################################################
############ FILL OUT TO CHOOSE WHICH AREA OF GRID TO PLOT   ##################
###############################################################################

xleft = 165000/hor_resolution #160km in on left side
xright = 320000/hor_resolution #320km on right side

z_scaling = 4

###############################################################################


##############################   Plot ########################################

    
fig = plt.figure(num=None, figsize=(20,8),  facecolor='w', edgecolor='k')
for j in range(1,3):
    subplot = 210 + j
    ax = plt.subplot(subplot,aspect = 'equal')
    plt.subplots_adjust(left=0.04, bottom=0.1, right=0.9, top=0.95, wspace=0, hspace=0)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    
    ##Levels variable 1
    zlmin = -1000
    zlmax = 2000.01
    zlevels = np.arange(zlmin,zlmax, 10)
    zlevels_ticks = np.arange(zlmin,zlmax,500)
    zlevels_ticks_labels = np.arange(zlmin,zlmax, 500).astype(int)
    
    #Levels variable 2
    xlmin = -20
    xlmax = 80.01
    xlevels = np.arange(xlmin,xlmax, 1)
    xlevels_ticks = np.arange(xlmin,xlmax,20)
    xlevels_ticks_labels = np.arange(xlmin,xlmax, 20).astype(int)
    
#    ##Levels variable 1
#    zlmin = -4
#    zlmax = 6.01
#    zlevels = np.arange(zlmin,zlmax, 0.1)
#    zlevels_ticks = np.arange(zlmin,zlmax,2)
#    zlevels_ticks_labels = np.arange(zlmin,zlmax, 2).astype(int)
#    
#    #Levels variable 2
#    xlmin = -4
#    xlmax = 6.01
#    xlevels = np.arange(xlmin,xlmax, 0.1)
#    xlevels_ticks = np.arange(xlmin,xlmax,2)
#    xlevels_ticks_labels = np.arange(xlmin,xlmax, 2).astype(int)
#    
    
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
    
    z = zpos_terrain[0,:,xleft:xright]/100*z_scaling
    #Plot variables
    if j == 1:
        z_disp_plot = plt.contourf(xpos[0,:,xleft:xright], z, mean_z_disp[:,xleft:xright], zlevels, cmap = shifted_cmapz, alpha = 1, zorder = 3, extend = 'both')
        #z_disp_plot = plt.contourf(x1d, y1d, mean_var1_disp, zlevels,  cmap = cmap_dth, vmin = -zlmax, alpha = 1, zorder = 3)

    if j == 2:
        #total = mean_var2_disp[:,xleft:xright]+mean_var5_disp[:,xleft:xright]+mean_var6_disp[:,xleft:xright]+mean_var7_disp[:,xleft:xright]+mean_var8_disp[:,xleft:xright]+mean_var9_disp[:,xleft:xright]+mean_var10_disp[:,xleft:xright]
        #total_disp = var5_disp + var6_disp + var7_disp + var8_disp + var9_disp + var10_disp        #x_disp_plot = plt.contourf(x2d, y2d, mean_var2_disp, xlevels,  cmap = cmap_var,vmin = -xlmax, alpha = 1, zorder = 3)
        x_disp_plot = plt.contourf(xpos[0,:,xleft:xright], z, mean_x_disp[:,xleft:xright]*hor_resolution/1000, xlevels,  cmap = shifted_cmapx, alpha = 1, zorder = 3, extend = 'both') #Percentage
        #x_disp_plot = plt.contourf(x1d, y1d, variable1[0,:,:], zlevels,  cmap = cmap_dth, alpha = 0.5, zorder = 4) #Percentage
    
    #Plot Terrain
    z_terrain  = zpos_terrain[0,0,xleft:xright]/100*z_scaling-2
    terrain = plt.plot(x1d[xleft:xright], z_terrain, c = 'slategrey', linewidth = 5)
    
    #Plot Lake
    lake = ds.xland[0,ymid,:].values
    lake[lake == 1] = np.nan
    lake_plt = plt.plot(x1d[xleft:xright], lake[xleft:xright]-2, c = 'blue', linewidth = 4, zorder = 6)
    
    
    #Title
    if j == 2:
        sub_title = '[Run: Low Wind]'
        ax.text(0.85, -0.28, sub_title, fontsize = 20, transform=ax.transAxes)
    
    
    ### Label and define grid area
    #y-axis
    plt.yticks(np.arange(0, np.max(z[:,0]+1),z_scaling*10))
    ax.set_yticklabels(np.arange(0,np.max(z[:,0]+1),1).astype(int), fontsize = 15, zorder = 6)
    plt.ylim([-2,np.max(z[:,0])-z_scaling*10+1])
    #x-axis
    plt.xticks(np.arange(0,ds.nx, 100))
    ax.set_xticklabels(np.arange(0, ds.nx*hor_resolution/1000, 100*hor_resolution/1000).astype(int), fontsize = 15)
    plt.xlim([xleft,xright])
    #axis labels
    plt.ylabel('Height (km)', fontsize = 20, labelpad = 8)
    if j == 2:
        plt.xlabel('Distance within Domain (km)', fontsize = 20, labelpad = 9)
                
    #Colorbar
    if j == 1:
        zcbaxes = fig.add_axes([0.92, 0.59, 0.03, 0.3])             
        zcbar = plt.colorbar(z_disp_plot, cax = zcbaxes, ticks = zlevels_ticks)
        zcbar.ax.set_yticklabels(zlevels_ticks_labels)
        zcbar.ax.tick_params(labelsize=15)
        plt.text(0.28, -0.17, 'm', fontsize = 21)
        
    if j == 2:
        xcbaxes = fig.add_axes([0.92, 0.16, 0.03, 0.3])             
        xcbar = plt.colorbar(x_disp_plot, cax = xcbaxes, ticks = xlevels_ticks)
        xcbar.ax.set_yticklabels(xlevels_ticks_labels)
        xcbar.ax.tick_params(labelsize=15)
        plt.text(0.2, -0.17, 'km', fontsize = 21)
    
    #Labels
    if j == 1:
        #sub_title = "Time-mean change in ${\Theta}$"
        sub_title = "Time-mean vertical displacement"
    if j == 2:
        #sub_title = 'Time-mean change in ${\Theta}$ due to MP'
        sub_title = 'Time-mean streamwise displacement'
    props = dict(boxstyle='square', facecolor='white', alpha=1)
    ax.text(0.01, 0.845, sub_title, fontsize = 20, bbox = props, zorder = 5, transform=ax.transAxes)
    

plt.savefig("/uufs/chpc.utah.edu/common/home/u1013082/public_html/phd_plots/cm1/plots/time_xy_0000m_15ms_60sec.png", dpi=200)
plt.close(fig)


#%%

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

###############################################################################
###################  Backward Trajectory Plot  ################################
###############################################################################

fig = plt.figure(num=None, figsize=(12,8), facecolor='w', edgecolor='k')


###################  Plot Trajectories in xy  ###############################

#Plot characteristics
ax = plt.subplot(211,aspect = 'equal')
#plt.title("Backward Trajectories", fontsize = 18, y = 1.05) 
plt.subplots_adjust(left=0.1, bottom=0.11, right=0.87, top=0.93, wspace=0, hspace=0.15)

xlmin = 262
xlmax = 278.01
xlevels = np.arange(xlmin,xlmax, 0.05)
xlevels_ticks = np.arange(xlmin,xlmax,2)
xlevels_ticks_labels = np.arange(xlmin,xlmax, 2).astype(int)

#xlmin = -0.3
#xlmax = 0.3
#xlevels = np.arange(xlmin,xlmax, 0.005)
#xlevels_ticks = np.arange(xlmin,xlmax,0.5)
#xlevels_ticks_labels = np.arange(xlmin,xlmax, 0.5).astype(int)

###################   Pick wich trajecotires to plot ##########################

##Convective region ###
#colmin = 1100
#colmax = 1190
#rowmin = 45
#rowmax = 65
#tmin = 0
#tmax = 175
#step = 4
#tstep = 2

#Windward cold pool ###
colmin = 1650
colmax = 2060
rowmin = 0
rowmax = 30
tmin = 0
tmax = 250
step = 4
tstep = 3

### end ###
#colmin = 1120#1080
#colmax = 1160#1140
#rowmin = 10
#rowmax = 25
#tmin = 0
#tmax = 150
#step = 2
#tstep = 2

### mid ###
#colmin = 880
#colmax = 940
#rowmin = 10
#rowmax = 25
#tmin = 0
#tmax = 200
#step = 2
#tstep = 3

###############################################################################


#Plot
xpos_reg = xpos[tmin:tmax:tstep,rowmin:rowmax:step, colmin:colmax:step]
ypos_reg = ypos[tmin:tmax:tstep,rowmin:rowmax:step, colmin:colmax:step]
zpos_terrain_reg = zpos_terrain[tmin:tmax:tstep,rowmin:rowmax:step, colmin:colmax:step]
color_reg = variable1[tmin:tmax:tstep,rowmin:rowmax:step, colmin:colmax:step]
color = variable1[tmin:tmax,rowmin:rowmax:step, colmin:colmax:step]

#Left half trajectories
rowsl,colsl = np.where(ypos_reg[-1,:,:] > ymid)
#Right half trajectories
rowsr,colsr = np.where(ypos_reg[-1,:,:] < ymid)

#Left
meanxl = np.nanmean(xpos_reg[:,rowsl,colsl], axis = 1)
meanyl = np.nanmean(ypos_reg[:,rowsl,colsl], axis = 1)
meanzl = np.nanmean(zpos_terrain_reg[:,rowsl,colsl],axis = 1)
meancl = np.nanmean(color_reg[:,rowsl,colsl], axis = 1)

#Right
meanxr = np.nanmean(xpos_reg[:,rowsr,colsr], axis = 1)
meanyr = np.nanmean(ypos_reg[:,rowsr,colsr], axis = 1)
meanzr = np.nanmean(zpos_terrain_reg[:,rowsr,colsr],axis = 1)
meancr = np.nanmean(color_reg[:,rowsr,colsr], axis = 1)

#All
meanx = np.nanmean(np.nanmean(xpos_reg, axis = 1), axis = 1)
meany = np.nanmean(np.nanmean(ypos_reg, axis = 1), axis = 1)
meanz = np.nanmean(np.nanmean(zpos_terrain_reg, axis = 1), axis = 1)
meanc = np.nanmean(np.nanmean(color_reg, axis = 1), axis = 1)


#Flip variables so highest parcels are on top
flipx = np.flip(xpos,0)
flipy = np.flip(ypos,0)
flipcolor = np.flip(color,0)
fliptmin = len(xpos[:,0,0]) - tmax
fliptmax = len(xpos[:,0,0]) - tmin

cmap = cmap_th
trajectories_xy = plt.scatter(flipx[fliptmin:fliptmax,rowmin:rowmax:step, colmin:colmax:step], flipy[fliptmin:fliptmax,rowmin:rowmax:step, colmin:colmax:step], c = flipcolor, cmap = cmap_th, norm=matplotlib.colors.BoundaryNorm(xlevels,cmap.N), s = 2.5, zorder = 7)
#trajectories_xy = plt.scatter(xpos[tmin:tmax,rowmin:rowmax:step, colmin:colmax:step], ypos[tmin:tmax,rowmin:rowmax:step, colmin:colmax:step], c = variable1[tmin:tmax,rowmin:rowmax:step, colmin:colmax:step], cmap = cmap_th, norm=matplotlib.colors.BoundaryNorm(xlevels,cmap.N), s = 2.5, zorder = 7)

trajectories_xy = plt.scatter(meanxl, meanyl, c = meancl, cmap = cmap_th, norm=matplotlib.colors.BoundaryNorm(xlevels,cmap.N), edgecolors='k',s = 40, zorder = 10)
trajectories_xy = plt.scatter(meanxr, meanyr, c = meancr, cmap = cmap_th, norm=matplotlib.colors.BoundaryNorm(xlevels,cmap.N), edgecolors='k',s = 40, zorder = 10)


plt.xticks(np.arange(0,ds.uinterp[0,0,0,:].size,200))
plt.yticks(np.arange(0,ds.uinterp[0,0,:,0].size,100))
ax.set_xticklabels(np.arange(0,ds.uinterp[0,0,0,:].size*200/100,40), fontsize = 14)
ax.set_yticklabels(np.arange(0,ds.uinterp[0,0,0,:].size*200/100,40), fontsize = 14)
plt.ylabel('Distance (km)', fontsize = 20)
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
xland_plt = plt.contour(ds.xland[0,:,:], levels_water, alpha = 1, colors = ('k'), zorder = 10)
xland_plt = plt.contourf(ds.xland[0,:,:], levels_terrain, alpha = 1, colors = ('gainsboro'), zorder = 1)

#Terrain
terrain_levels = np.arange(-1, 3000.1, 200)
try:
    zs = ds.zs[0,:,:].values
except:
    zs = np.zeros((ds.ny,ds.nx))

terrain = plt.contourf(zs[:,:], terrain_levels, alpha = 1, cmap = cm.Greys, vmin = -800, vmax = 3000, zorder = 1)




###################  Plot Trajectories in xz  ###############################


#Plot characteristics
ax = plt.subplot(212)

#Plot
#total_disp = var2_disp + var5_disp + var6_disp + var7_disp + var8_disp + var9_disp + var10_disp
#total_variable = variable2 + variable5 + variable6 + variable7 + variable8 + variable9 + variable10


trajectories_xz = plt.scatter(xpos[tmin:tmax,rowmin:rowmax:step, colmin:colmax:step], zpos_terrain[tmin:tmax,rowmin:rowmax:step, colmin:colmax:step], c = color, cmap = cmap_th, norm=matplotlib.colors.BoundaryNorm(xlevels,cmap.N), s = 2.5, zorder = 3)

trajectories_xz = plt.scatter(meanx, meanz, c = meanc, cmap = cmap_th, norm=matplotlib.colors.BoundaryNorm(xlevels,cmap.N), edgecolors='k', s = 40, zorder = 10)



plt.grid(True, color = 'white', )
plt.ylabel('Height (m)', fontsize = 20)
plt.xlabel('Distance (km)', fontsize = 20)
plt.xticks(np.arange(0,ds.nx,40*1000/200))
ytick = np.arange(0,4001,500)
plt.yticks(ytick)
ax.set_xticklabels(np.arange(0,ds.nx*200/1000,40), fontsize = 14)
plt.xlim([0,ds.uinterp[0,0,0,:].size])
plt.ylim([-50,ytick[-1]])
ax.tick_params(axis='y', labelsize=14)
plt.axvspan(0,ds.uinterp[0,0,0,:].size,color='gainsboro',lw=0)

#Labels
#sub_title = ['[Run: Low Wind]']
#props = dict(boxstyle='square', facecolor='white', alpha=1)
#ax.text(1800, -850, sub_title[0], fontsize = 15, zorder = 5)


#Plot Terrain
x1d = np.arange(0,ds.nx,1)
ymid = np.int(ds.ny/2)
terrain = plt.plot(x1d, zs[ymid,:]-40, c = 'slategrey', linewidth = 4, zorder = 4)

#Plot Lake
lake = np.array(ds.xland[0,ymid,:])
lake[lake == 1] = np.nan
lake_plt = plt.plot(x1d, lake-40, c = 'blue', linewidth = 4, zorder = 5)

#Colorbar
cbaxes = fig.add_axes([0.9, 0.27, 0.035, 0.52])             
cbar = plt.colorbar(trajectories_xz, cax = cbaxes, ticks = xlevels_ticks)
cbar.ax.set_yticklabels(xlevels_ticks_labels)
cbar.ax.tick_params(labelsize=14)
plt.text(-0.6, -0.1, 'Theta (K)', fontsize = 20)

zs_patch = np.int((colmin+colmax)/2)
try:
    ax.add_patch(patches.Rectangle((colmin,rowmin*33.333+ds.zs[0,ymid,zs_patch]),colmax-colmin,(rowmax-rowmin)*33.333,fill=False, zorder = 10, linewidth = 1.5))
except:
    ax.add_patch(patches.Rectangle((colmin,rowmin*33.333),colmax-colmin,(rowmax-rowmin)*33.333,fill=False, zorder = 10, linewidth = 1.5))


plt.savefig("/uufs/chpc.utah.edu/common/home/u1013082/public_html/phd_plots/cm1/plots/backward_trajectory_0000m_25ms_90sec.png", dpi = 250)
plt.close(fig)
    



