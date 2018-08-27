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
import shiftedcolormap
import matplotlib.patches as mpatches

#Read in with xarray
dsl_no = xr.open_dataset('/uufs/chpc.utah.edu/common/home/steenburgh-group8/tom/cm1/output/tug/cm1run_150m_15ms_0000m_90sec.nc')
dsl_small = xr.open_dataset('/uufs/chpc.utah.edu/common/home/steenburgh-group7/tom/cm1/output/netcdf_output/tug/cm1run_150m_15ms_0500m_90sec.nc')
dsl_tall = xr.open_dataset('/uufs/chpc.utah.edu/common/home/steenburgh-group7/tom/cm1/output/netcdf_output/tug/cm1run_150m_15ms_2000m_90sec.nc')

dsh_no = xr.open_dataset('/uufs/chpc.utah.edu/common/home/steenburgh-group7/tom/cm1/output/netcdf_output/tug/cm1run_150m_25ms_0000m_90sec.nc')
dsh_small = xr.open_dataset('//uufs/chpc.utah.edu/common/home/steenburgh-group8/tom/cm1/output/tug/cm1run_150m_25ms_0500m_90sec.nc')
dsh_tall = xr.open_dataset('/uufs/chpc.utah.edu/common/home/steenburgh-group8/tom/cm1/output/tug/cm1run_150m_25ms_2000m_90sec.nc')



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
far = dsl_no.ny-160 #far side of ypos for traj
num_seeds_y = far-near
bottom = 1 #starting zpos
top = 40 #top of data chucks
ymid = np.int(dsl_no.ny/2)


#Dimension size variables
num_x = dsl_no.nx-left
num_y = dsl_no.ny
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






############## Read in trajectory Data for high wind ###########################

## Low wind
xposl_no = np.load('xpos_terrain_zoom_traj_disp_15ms_0000m_100time_steps.npy')
yposl_no = np.load('ypos_terrain_zoom_traj_disp_15ms_0000m_100time_steps.npy')
zpos_terrainl_no = np.load('zpos_terrain_zoom_traj_disp_15ms_0000m_100time_steps.npy')

xposl_small = np.load('xpos_terrain_zoom_traj_disp_15ms_0500m_100time_steps.npy')
yposl_small = np.load('ypos_terrain_zoom_traj_disp_15ms_0500m_100time_steps.npy')
zpos_terrainl_small = np.load('zpos_terrain_zoom_traj_disp_15ms_0500m_100time_steps.npy')

xposl_tall = np.load('xpos_terrain_zoom_traj_disp_15ms_2000m_100time_steps.npy')
yposl_tall = np.load('ypos_terrain_zoom_traj_disp_15ms_2000m_100time_steps.npy')
zpos_terrainl_tall = np.load('zpos_terrain_zoom_traj_disp_15ms_2000m_100time_steps.npy')


#High wind
xposh_no = np.load('xpos_terrain_zoom_traj_disp_25ms_0000m_100time_steps.npy')
yposh_no = np.load('ypos_terrain_zoom_traj_disp_25ms_0000m_100time_steps.npy')
zpos_terrainh_no = np.load('zpos_terrain_zoom_traj_disp_25ms_0000m_100time_steps.npy')

xposh_small = np.load('xpos_terrain_zoom_traj_disp_25ms_0500m_100time_steps.npy')
yposh_small = np.load('ypos_terrain_zoom_traj_disp_25ms_0500m_100time_steps.npy')
zpos_terrainh_small = np.load('zpos_terrain_zoom_traj_disp_25ms_0500m_100time_steps.npy')

xposh_tall = np.load('xpos_terrain_zoom_traj_disp_25ms_2000m_100time_steps.npy')
yposh_tall = np.load('ypos_terrain_zoom_traj_disp_25ms_2000m_100time_steps.npy')
zpos_terrainh_tall = np.load('zpos_terrain_zoom_traj_disp_25ms_2000m_100time_steps.npy')





### Calulate array to label by starting over water or not
water_mark = np.zeros((6,time_steps, num_seeds_z, num_seeds_y))
land_mark = np.zeros((6,time_steps, num_seeds_z, num_seeds_y))

xland = dsl_no.xland[0,:,:].values #1500 is starting point of trajectories




    #Water and Land Labels
for run in range(6):
    
    #xpos
    run_name_x = ['xposl_no', 'xposl_small', 'xposl_tall','xposh_no', 'xposh_small', 'xposh_tall']
    xpos = eval(run_name_x[run])
    
    #ypos
    run_name_y = ['yposl_no', 'yposl_small', 'yposl_tall','yposh_no', 'yposh_small', 'yposh_tall']
    ypos = eval(run_name_y[run])
    
    for i in range(num_seeds_y):
        for j in range(time_steps):
            
            try:
                x_ind = np.int(xpos[j,0,i])+1500
                xwater = np.where(xland[:,x_ind] == 2)
                water_max = np.max(xwater)
                water_min = np.min(xwater)
        
        
                if ypos[j,0,i] <= water_max and ypos[j,0,i] >= water_min:
                    water_mark[run,:,0,i] = 1
            except:
                pass
    
land_mark[water_mark == 0] = 1

water_mark[water_mark == 0] = np.nan
land_mark[land_mark == 0] = np.nan






#x = xpos[:,0,100]
#y = ypos[:,0,100]
#z = zpos_terrain[:,0,100]
#x_new = np.linspace(x.min(), x.max(),10000)
#f = interpolate.interp1d(x, z, kind='linear')
#z_smooth=f(x_new)
#f = interpolate.interp1d(x, y, kind='linear')
#y_smooth=f(x_new)

#%%
###############################################################################
############################## Plot XY ########################################
###############################################################################


left = 1500
right = 2050
bottom = 0
top = 50
t_xy = 1
    
#Colormap
colors1_t = np.array(nclcmaps.colors['BlueRed'])#amwg256'])
colors_int_t = colors1_t.astype(int)
colors_t = list(colors_int_t)
cmap_th = nclcmaps.make_cmap(colors_t, bit = True)

#cmap_wind = nclcmaps.cmap('ncl_default')

#Levels
lmin = 0.4
lmax = 0.6
levels = np.arange(lmin,lmax, 0.25)
levels_ticks = np.arange(lmin,lmax,1)
levels_ticks_labels = np.arange(lmin,lmax, 1).astype(int)

#Colors
cland = 'red'
cwater = 'blue'
skip = 3

#shifted_cmap = shiftedcolormap.shiftedColorMap(cmap_wind, midpoint=1 - lmax/(lmax + abs(lmin)), name='shifted')

    
#########################  Create Fig  ########################################

    


#Create Fig
fig = plt.figure(num=None, figsize=(18.2,9), facecolor='w', edgecolor='k')
#Loop over subplots
for run in range(1,7):
    

    #loop over runs
    #For run info
    run_name = ['dsl_no', 'dsl_small', 'dsl_tall','dsh_no', 'dsh_small', 'dsh_tall']
    model_run = eval(run_name[run-1])
    
    #xpos
    run_name_x = ['xposl_no', 'xposl_small', 'xposl_tall','xposh_no', 'xposh_small', 'xposh_tall']
    xpos = eval(run_name_x[run-1])
    
    #ypos
    run_name_y = ['yposl_no', 'yposl_small', 'yposl_tall','yposh_no', 'yposh_small', 'yposh_tall']
    ypos = eval(run_name_y[run-1])
    
        
    #############################  Plot xy  ###################################
    ax = plt.subplot(230 + run,aspect = 'equal')
    plt.subplots_adjust(left=0.06, bottom=0.15, right=0.96, top=0.95, wspace=0.06, hspace=0.01)
    #ax.axis('equal')

    
    #Plot Trajectories
    traj_xy_land = plt.scatter(xpos[:,0,::skip]*land_mark[run-1,:,0,::skip], ypos[:,0,::skip]*land_mark[run-1,:,0,::skip]-near, color = cland,  s = 3, zorder = 7, alpha  = 0.6)
    traj_xy_water = plt.scatter(xpos[:,0,::skip]*water_mark[run-1,:,0,::skip], ypos[:,0,::skip]*water_mark[run-1,:,0,::skip]-near, color = cwater,  s = 3, zorder = 7, alpha  = 0.6)


    #Limits
    ax.set_xlim([0,right-left])
    ax.set_ylim([0,far-near])

    #Land, water, terrain
    levels_water = [1.5, 2.5]
    terrain_levels_tall = np.arange(0, 2000.1, 200)
    terrain_levels_small = np.arange(0, 800.1, 200)
    terrain_levels = np.arange(0, 300.1, 200)

    
    #Plot Lake
    water = plt.contourf(model_run.xland[0,near:far,left:right], levels_water, alpha = 1, colors = ('skyblue'), zorder = 5, linewidths = 4)
    water = plt.contour(model_run.xland[0,near:far,left:right], levels_water, alpha = 1, colors = ('skyblue'), zorder = 20, linewidths = 4)
    
    #Plot Terrain
    try:
        if run == 3 or run == 6:
            terrain = plt.contourf(model_run.zs[0,near:far,left:right], terrain_levels_tall, alpha = 1, cmap = cm.Greys, vmin = -800, vmax = 2500, zorder = 4, linewidths = 1.75)
            terrain = plt.contour(model_run.zs[0,near:far,left:right], terrain_levels, alpha = 1, colors = 'w', zorder = 4, linewidths = 1.75)
        if run == 2 or run == 5:
            terrain = plt.contourf(model_run.zs[0,near:far,left:right], terrain_levels_small, alpha = 1, cmap = cm.Greys, vmin = -800, vmax = 2500, zorder = 4, linewidths = 1.75)
            terrain = plt.contour(model_run.zs[0,near:far,left:right], terrain_levels, alpha = 1, colors = 'w', zorder = 4, linewidths = 1.75)
        
        #Plot Point at peak
        peak = np.where(model_run.zs[0,near:far,left:right].values == np.max(model_run.zs[0,near:far,left:right].values))
        ax.add_patch(patches.Circle((peak[1][0], peak[0][0]), 7, color = 'w', zorder = 20))

    except:
        zs = np.zeros((far-near, right-left))
        terrain = plt.contourf(zs, terrain_levels_tall, alpha = 1, cmap = cm.Greys, vmin = -800, vmax = 2500, zorder = 4, linewidths = 1.75)
    
    
    #Plot Characteristics

    ax.set_xlim([0,right-left])
    ax.set_ylim([0,far-near])

    ytick = np.arange(0,far-near, 10000/hor_resolution+1)
    ax.set_yticks(ytick)
    xtick =np.arange(0,right-left,10000/hor_resolution+1)
    ax.set_xticks(xtick)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    

    
    if run == 1 or run == 4:
        ax.set_yticks(ytick)
        ax.set_yticklabels(ytick*hor_resolution/1000, fontsize = 16)
        ax.set_ylabel("Distance (km)", fontsize = 22)
    else:
        for tic in ax.yaxis.get_major_ticks():
            tic.tick1On = tic.tick2On = False
    
    if run == 4 or run == 5 or run == 6:
        ax.set_xticks(xtick)
        ax.set_xticklabels(xtick*hor_resolution/1000, fontsize = 16)
        ax.set_xlabel("Distance (km)", fontsize = 22)
    else:
        for tic in ax.xaxis.get_major_ticks():
            tic.tick1On = tic.tick2On = False


    #Draw border
    ax.add_patch(patches.Rectangle((1.5,1),right-left-3,far-near-2,fill=False, zorder = 10, linewidth = 2))
    
    #Titles
    if run == 1:
        plt.title("No Mountain", fontsize = 26, y = 1)
        ax.text(3.15, 0.67, 'Low Wind', transform=ax.transAxes, fontsize = 26, rotation = -90)
    if run == 2: 
        plt.title("500m Mountain", fontsize = 26, y = 1)
    if run == 3: 
        plt.title("2000m Mountain", fontsize = 26, y = 1)
    if run == 4:
        ax.text(3.15, 0.67, 'High Wind', transform=ax.transAxes, fontsize = 26, rotation = -90)
    if run == 6:
        props = dict(boxstyle='square', facecolor='white', alpha=1, linewidth = 2)

        
#Legend
land_lab = mpatches.Patch(color=cland, label='Over Land Trajectories',alpha = 1)
water_lab = mpatches.Patch(color=cwater, label='Over Water Trajectories', alpha = 1)
plt.legend(bbox_to_anchor=(0.36, -0.23), handles=[water_lab, land_lab], fontsize = 21, ncol = 2)
    
##Colorbars
#cbaxes = fig.add_axes([0.27, 0.09, 0.5, 0.045])           
#cbar = plt.colorbar(traj_xy, orientation='horizontal', cax = cbaxes, ticks = levels_ticks)
#cbar.ax.tick_params(labelsize=19)
#plt.text(0.35, -1.85, 'Potential Temperature (K)', fontsize = 24)



#Save and Close
plt.savefig("/uufs/chpc.utah.edu/common/home/u1013082/public_html/phd_plots/cm1/plots/trajectories_terrain.png", dpi = 150)
plt.close(fig)  
plt.switch_backend('Agg')


#%%

###############################################################################
############################## Plot XZ ########################################
###############################################################################

    
################################################################################
########################### Set up coordinates ################################
################################################################################

vert_resolution = 100
hor_resolution = 150
z_scale = 5
ymid = np.int(dsl_no.ny/2)
t = 270
left = 1500
right = 2050
bottom = 0
top = 35
near = 160

## The code below makes the data terrain following 
x1d = np.arange(0,right-left,1)
z1d = np.arange(0,top,1)

#Create 2D arrays for plotting data (first demnsion for each run)
x2d = np.zeros((6,top, right-left))
z2d = np.zeros((6,top, right-left))
run = ['dsl_no', 'dsh_no', 'dsl_small','dsh_small', 'dsl_tall', 'dsh_tall']


for j in range(6):
    model_run = eval(run[j])
    try:
        z = np.array(model_run.zs[0,ymid,left:right])/vert_resolution #Convert to gridpoints
    except:
        z = np.zeros((right-left))+0.4

    for i in range(top):
        x2d[j,i,:] = x1d
    for k in range(right-left):
        z2d[j,:,k] = z1d+z[k]
        
        
        
fig = plt.figure(num=None, figsize=(16,8.5), facecolor='w', edgecolor='k')
#Loop over subplots
for run in range(1,7):
    

    #loop over runs
    #loop over runs
    #For run info
    run_name = ['dsl_no', 'dsh_no', 'dsl_small','dsh_small', 'dsl_tall', 'dsh_tall']
    model_run = eval(run_name[run-1])
    
    #xpos
    run_name_x = ['xposl_no', 'xposh_no', 'xposl_small','xposh_small', 'xposl_tall', 'xposh_tall']
    xpos = eval(run_name_x[run-1])
    
    #ypos
    run_name_z = ['zpos_terrainl_no', 'zpos_terrainh_no', 'zpos_terrainl_small','zpos_terrainh_small', 'zpos_terrainl_tall', 'zpos_terrainh_tall']
    zpos_terrain = eval(run_name_z[run-1])
    
        
    #############################  Plot xz ###################################
    ax = plt.subplot(320 + run,aspect = 'equal')
    plt.subplots_adjust(left=0.05, bottom=0.17, right=0.95, top=0.94, wspace=0.02, hspace=0.02)
    #plt.axis('equal')
    
    #Plot Trajectories
    ind = np.array([0,3,1,4,2,5]) #Conversion becuase of different ordering of subplots
    traj_xy_land = plt.scatter(xpos[:,0,::skip]*land_mark[ind[run-1],:,0,::skip], zpos_terrain[:,0,::skip]/vert_resolution*z_scale*land_mark[ind[run-1],:,0,::skip], color = cland,  s = 2.5, zorder = 7, alpha  = 0.6)
    traj_xy_water = plt.scatter(xpos[:,0,::skip]*water_mark[ind[run-1],:,0,::skip], zpos_terrain[:,0,::skip]/vert_resolution*z_scale*water_mark[ind[run-1],:,0,::skip], color = cwater,  s = 2.5, zorder = 7, alpha  = 0.6)




    
    print(run)
    
    #Plot Terrain
    terrain = plt.plot(x1d, z2d[run-1,0,:]*z_scale+0.5, c = 'grey', linewidth = 4, zorder = 14)
    
    #Plot Lake
    lake = model_run.xland[0,ymid,left:right].values
    lake[lake == 1] = np.nan
    lake_plt = plt.plot(x1d, lake+1, c = 'skyblue', linewidth = 4, zorder = 15)
    
    #Plot Characteristics
    plt.grid(True, color = 'white', )
    ax.set_xlim([0,right-left])
    ax.set_ylim([0,top*z_scale-5*z_scale+1])
    plt.axvspan(0,dsl_no.nx,color='gainsboro',lw=0)

    ytick = np.arange(0,top*z_scale,5*z_scale)
    ax.set_yticks(ytick)
    xtick =np.arange(0,right-left,10000/hor_resolution+1)
    ax.set_xticks(xtick)
    ax.set_xticklabels([])
    ax.set_yticklabels([])


    
    if run == 1 or run == 3 or run == 5:
        ax.set_yticks(ytick)
        yticklabs = ytick*vert_resolution/z_scale/1000.
        ax.set_yticklabels(yticklabs, fontsize = 13)
        ax.set_ylabel("Height (km)", fontsize = 16)
    else:
        for tic in ax.yaxis.get_major_ticks():
            tic.tick1On = tic.tick2On = False
    
    
    if run == 5 or run == 6:
        ax.set_xticks(xtick)
        ax.set_xticklabels(xtick*hor_resolution/1000, fontsize = 13)
        ax.set_xlabel("Distance (km)", fontsize = 16)
    else:
        for tic in ax.xaxis.get_major_ticks():
            tic.tick1On = tic.tick2On = False


    
    #Titles
    props = dict(boxstyle='square', facecolor='white', alpha=1, linewidth = 2)
    if run == 1:
        plt.title("Low Wind", fontsize = 23, y = 1.03)
    if run == 2: 
        plt.title("High Wind", fontsize = 23, y = 1.03)
        ax.text(1.01, 0.95, 'No Mountain'.center(15), transform=ax.transAxes,fontsize = 20, rotation = -90)
    if run == 4: 
        ax.text(1.01, 0.58, '500m\nMountain'.center(20), transform=ax.transAxes, fontsize = 20, rotation = -90)
    if run == 6:
        #ax.text(0.5, -0.5, 'Elapsed time: {:,} seconds'.format(t*90), transform=ax.transAxes, bbox = props, fontsize = 16)
        ax.text(1.01, 0.58, '2000m\nMountain'.center(18), transform=ax.transAxes,fontsize = 20, rotation = -90)

        
#Legend
land_lab = mpatches.Patch(color=cland, label='Over Land Trajectories',alpha = 1)
water_lab = mpatches.Patch(color=cwater, label='Over Water Trajectories', alpha = 1)
plt.legend(bbox_to_anchor=(0.66, -0.35), handles=[water_lab, land_lab], fontsize = 21, ncol = 2)
    
##Colorbars
#cbaxes = fig.add_axes([0.3, 0.075, 0.4, 0.03])           
#cbar = plt.colorbar(traj_xz, orientation='horizontal', cax = cbaxes, ticks = levels_ticks)
#cbar.ax.tick_params(labelsize=15)
#plt.text(0.35, -2.2, 'U-wind  $\mathregular{(ms^{-1})}$', fontsize = 18)

#Quiver
#    plt.quiverkey(quiv, X=-1.65, Y=-0.2, U=20, linewidth = 0.75, color = 'k', label='20 $\mathregular{ms^{-1}}$', labelpos='W', fontproperties={'size': '28'})

#Save and Close
plt.savefig("/uufs/chpc.utah.edu/common/home/u1013082/public_html/phd_plots/cm1/plots/trajectories_terrain_xz.png", dpi = 150)
plt.close(fig)  
plt.switch_backend('Agg')