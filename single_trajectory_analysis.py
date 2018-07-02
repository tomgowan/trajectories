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
#ds = xr.open_dataset('/uufs/chpc.utah.edu/common/home/steenburgh-group7/tom/cm1/cm1r19B/run/cm1run_20ms_1500m_tug.nc')
ds = xr.open_dataset('/uufs/chpc.utah.edu/common/home/steenburgh-group7/tom/cm1/output/netcdf_output/tug/cm1run_15ms_1500m_90sec.nc')



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
time_steps = 300 #Number of time steps to run trajectories back
start_time_step = 330 #Starting time step
hor_resolution = 200 #[meters]
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

var_name3 = 'ptb_hadv'#Coded to be a budget variable
variable3 = np.zeros((time_steps, num_seeds_z, num_seeds_x))

var_name4 = 'ptb_vadv'#Coded to be a budget variable
variable4 = np.zeros((time_steps, num_seeds_z, num_seeds_x))

var_name5 = 'ptb_hturb'#Coded to be a budget variable
variable5 = np.zeros((time_steps, num_seeds_z, num_seeds_x))

var_name6 = 'ptb_vturb'#Coded to be a budget variable
variable6 = np.zeros((time_steps, num_seeds_z, num_seeds_x))

var_name7 = 'ptb_rdamp'#Coded to be a budget variable
variable7 = np.zeros((time_steps, num_seeds_z, num_seeds_x))

var_name8 = 'ptb_rad'#Coded to be a budget variable
variable8 = np.zeros((time_steps, num_seeds_z, num_seeds_x))

var_name9 = 'ptb_div'#Coded to be a budget variable
variable9 = np.zeros((time_steps, num_seeds_z, num_seeds_x))

var_name10 = 'ptb_diss'#Coded to be a budget variable
variable10 = np.zeros((time_steps, num_seeds_z, num_seeds_x))


#Tendancies
var_name11 = 'tt_cond'#Coded to be a budget variable
variable11 = np.zeros((time_steps, num_seeds_z, num_seeds_x))

var_name12 = 'tt_evac'#Coded to be a budget variable
variable12 = np.zeros((time_steps, num_seeds_z, num_seeds_x))

var_name13 = 'tt_evar'#Coded to be a budget variable
variable13 = np.zeros((time_steps, num_seeds_z, num_seeds_x))

var_name14 = 'tt_dep'#Coded to be a budget variable
variable14 = np.zeros((time_steps, num_seeds_z, num_seeds_x))

var_name15 = 'tt_subl'#Coded to be a budget variable
variable15 = np.zeros((time_steps, num_seeds_z, num_seeds_x))

var_name16 = 'tt_melt'#Coded to be a budget variable
variable16 = np.zeros((time_steps, num_seeds_z, num_seeds_x))

var_name17 = 'tt_frz'#Coded to be a budget variable
variable17 = np.zeros((time_steps, num_seeds_z, num_seeds_x))


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
xpos = np.load('xpos_traj_disp_15ms_1500m_300time_steps.npy')
ypos = np.load('ypos_traj_disp_15ms_1500m_300time_steps.npy')
zpos_terrain = np.load('zpos_traj_disp_15ms_1500m_300time_steps.npy')
variable1 = np.load('%s_traj_disp_15ms_1500m_300time_steps.npy' % var_name1)
variable2 = np.load('%s_traj_disp_15ms_1500m_300time_steps.npy' %var_name2)
variable3 = np.load('%s_traj_disp_15ms_1500m_300time_steps.npy' %var_name3)
variable4 = np.load('%s_traj_disp_15ms_1500m_300time_steps.npy' %var_name4)
variable5 = np.load('%s_traj_disp_15ms_1500m_300time_steps.npy' %var_name5)
variable6 = np.load('%s_traj_disp_15ms_1500m_300time_steps.npy' %var_name6)
variable7 = np.load('%s_traj_disp_15ms_1500m_300time_steps.npy' %var_name7)
variable8 = np.load('%s_traj_disp_15ms_1500m_300time_steps.npy' %var_name8)
variable9 = np.load('%s_traj_disp_15ms_1500m_300time_steps.npy' %var_name9)
variable10 = np.load('%s_traj_disp_15ms_1500m_300time_steps.npy' %var_name10)
variable11 = np.load('%s_traj_disp_15ms_1500m_300time_steps.npy' %var_name11)
variable12 = np.load('%s_traj_disp_15ms_1500m_300time_steps.npy' %var_name12)
variable13 = np.load('%s_traj_disp_15ms_1500m_300time_steps.npy' %var_name13)
variable14 = np.load('%s_traj_disp_15ms_1500m_300time_steps.npy' %var_name14)
variable15 = np.load('%s_traj_disp_15ms_1500m_300time_steps.npy' %var_name15)
variable16 = np.load('%s_traj_disp_15ms_1500m_300time_steps.npy' %var_name16)
variable17 = np.load('%s_traj_disp_15ms_1500m_300time_steps.npy' %var_name17)


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
#colmin = 1200
#colmax = 1250
#rowmin = 5
#rowmax = 20
#tmin = 0
#tmax = 200
#step = 2
#tstep = 3

### end ###
colmin = 1120#1080
colmax = 1160#1140
rowmin = 10
rowmax = 25
tmin = 0
tmax = 150
step = 2
tstep = 2

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

var2_reg = variable2[tmin:tmax:tstep,rowmin:rowmax:step, colmin:colmax:step]
var6_reg = variable6[tmin:tmax:tstep,rowmin:rowmax:step, colmin:colmax:step]

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
meanvar2 = np.nanmean(np.nanmean(var2_reg, axis = 1), axis = 1)
meanvar6 = np.nanmean(np.nanmean(var6_reg, axis = 1), axis = 1)







###################  Plot Trajectories in xz  ###############################
###############################################################################
###################  Backward Trajectory Plot  ################################
###############################################################################

fig = plt.figure(num=None, figsize=(10,4), facecolor='w', edgecolor='k')


###################  Plot Trajectories in xy  ###############################

#Plot characteristics
ax = plt.subplot(111,aspect = 'equal')
plt.title("Backward Trajectories", fontsize = 18, y = 1.05) 
plt.subplots_adjust(left=0.1, bottom=0.11, right=0.9, top=0.93, wspace=0, hspace=0.15)

xlmin = 262
xlmax = 278.01
xlevels = np.arange(xlmin,xlmax, 0.05)
xlevels_ticks = np.arange(xlmin,xlmax,2)
xlevels_ticks_labels = np.arange(xlmin,xlmax, 2).astype(int)



#Left Axis
ax1 = plt.subplot(111)
plt.grid(True, color = 'white', )
trajectories_xz = ax1.plot(meanx, meanc, 'k', zorder = 10)

ax1.set_ylabel('Potential Temperature (K)', fontsize = 11, labelpad = 10)
ax1.set_xlabel('Distance within Domain (km)', fontsize = 11)
ax1.set_xticks(np.arange(0,ds.nx,20*1000/200))
ytick = np.arange(264,269.01,0.5)
ax1.set_yticks(ytick)
ax1.set_xticklabels(np.arange(0,ds.nx*200/1000,20), fontsize = 8)
ax1.set_xlim([200,ds.uinterp[0,0,0,:].size-500])
ax1.set_ylim([ytick[0],ytick[-1]])
ax1.tick_params(axis='y', labelsize=8)
ax1.axvspan(0,ds.uinterp[0,0,0,:].size,color='gainsboro',lw=0)

#Right Axis
ax2 = ax1.twinx()
ax2.set_ylabel('Diabatic Heating Rate (K/s)', fontsize = 11, rotation = 270, labelpad = 20)
ytick = np.arange(-0.5,0.5001,0.1)
ax2.set_yticks(ytick)
ax2.set_ylim([ytick[0],ytick[-1]])
ax2.tick_params(axis='y', labelsize=8)

width = 6
plt.bar(meanx, meanvar2, width = width, color='r', zorder = 0)
plt.bar(meanx, meanvar6, width = width, bottom = meanvar2, color='b', zorder = 0)

#Labels
sub_title = ['[Run: 15 $\mathregular{ms^{-1}}$ and 1500m]']
props = dict(boxstyle='square', facecolor='white', alpha=1)
ax1.text(1250, -850, sub_title[0], fontsize = 11, zorder = 0)


##Plot Terrain
#x1d = np.arange(0,ds.nx,1)
#ymid = np.int(ds.zs[0,:,0].size/2)
#terrain = ax1.plot(x1d, ds.zs[0,ymid,:]-40, c = 'slategrey', linewidth = 4, zorder = 4)

#Plot Lake
lake = np.array(ds.xland[0,ymid,:])
lake[lake == 1] = np.nan
lake[lake == 2] = ytick[0]
x1d = np.arange(0,num_seeds_x,1)
lake_plt = ax2.plot(x1d, lake, c = 'blue', linewidth = 4, zorder = 5)





plt.savefig("/uufs/chpc.utah.edu/common/home/u1013082/public_html/phd_plots/cm1/plots/individual_trajectory_analysis.png", dpi = 200)
plt.close(fig)
    




