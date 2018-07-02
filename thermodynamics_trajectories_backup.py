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
ds = xr.open_dataset('/uufs/chpc.utah.edu/common/home/steenburgh-group7/tom/cm1/cm1r19B/run/cm1run_20ms_1500m_tug.nc')
#ds = xr.open_dataset('/uufs/chpc.utah.edu/common/home/steenburgh-group7/tom/cm1/cm1r19B/run/cm1run_20ms_1500m_tug_wide_45sec.nc')


#xpos = np.load('xpos_traj_500m.npy')
#ypos = np.load('ypos_traj_500m.npy')
#zpos_terrain = np.load('zpos_traj_500m.npy')
#color = np.load('color_traj_500m.npy')
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
num_seeds_z = 151 #Up to 5000m (one seed every grid point so they can be contoured)
num_seeds_x = ds.nx #About right half of domain (one seed every fgridpoint)
time_steps = 60 #Run trajectories back 100 time steps (all during steady-state)
start_time_step = 260 #Start near end of simulation
hor_resolution = 200 #meters
vert_resolution = 100.0 #meters
time_step_length = 120.0 #seconds
###############################################################################
###############################################################################

#%%
#Create arrays for location and variable color of parcels
xpos = np.zeros((time_steps, num_seeds_z, num_seeds_x))
ypos = np.zeros((time_steps, num_seeds_z, num_seeds_x))
zpos = np.zeros((time_steps, num_seeds_z, num_seeds_x))

#Variable
var_name1 = 'th'
variable1 = np.zeros((time_steps, num_seeds_z, num_seeds_x))

var_name2 = 'ptb_mp'#Coded to be a budget variable
variable2 = np.zeros((time_steps, num_seeds_z, num_seeds_x))



###############################################################################
##################### INITIAL LOCATION OF SEEDS ###############################
#########  (These may also be important to the user)  #########################
###############################################################################
#Choose starting point 
#xpos
rdist = 50 #distance from back
for i in range(0, num_seeds_x):
    xpos[0,:,i] = ds.th[0,0,0,:].size-num_seeds_x-rdist + i #800 gridpoints from right side of domain

#ypos   
ymid = np.int(ds.th[0,0,:,0].size/2)
ypos[0,:,:] = ymid #Middle y-gridpoint

#zpos
for i in range(num_seeds_z):
    zpos[0,i,:] = i/3.0 #Surface to top of desired regions

###############################################################################
###############################################################################


#####  Use model terrain to correct errors from computing trajectories with
#####  with terrain following coordinates
    
#check if run has terrain (if not, zs is zero and zpos should not be affected)
try:
    zs = np.array(ds.zs[0,:,:])
except:
    zs = np.zeros((ds.yh.size, ds.xh.size))  
    
#Add terrain height
#Create list of current coordinates for terrain addition

xloc = np.array(xpos[0,:,:]).flatten()
yloc = np.array(ypos[0,:,:]).flatten()
coor_terrain = []
for i in range(len(xloc)):
    coor_terrain.append((yloc[i], xloc[i]))

zpos_terrain = np.array(zpos)
zpos_terrain[0,:,:] = zpos_terrain[0,:,:] + np.reshape(interpolate.interpn((y,x), zs, coor_terrain, method='linear', bounds_error=False, fill_value= 0)/vert_resolution, (num_seeds_z, num_seeds_x))



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
    
    #Create list of current coordinates
    xloc = np.array(xpos[t,:,:]).flatten()
    yloc = np.array(ypos[t,:,:]).flatten()
    zloc = np.array(zpos[t,:,:]).flatten()
    
    coor = []
    for i in range(len(xloc)):
        coor.append((zloc[i], yloc[i], xloc[i]))
    coor_array = np.array(coor)




    #Timer
    stop = time.time()
    print(stop-start)


            
    #####################   Calc new xpos #####################################
    xpos[t+1,:,:] = xpos[t,:,:] - np.reshape(interpolate.interpn((z,y,x), u, coor, method='linear', bounds_error=False, fill_value=np.nan)*time_step_length/hor_resolution, (num_seeds_z, num_seeds_x))

    #####################   Calc new ypos #####################################
    ypos[t+1,:,:]  = ypos[t,:,:] - np.reshape(interpolate.interpn((z,y,x), v, coor, method='linear', bounds_error=False, fill_value=np.nan)*time_step_length/hor_resolution, (num_seeds_z, num_seeds_x))

    #####################   Calc new zpos #####################################
    #terrain-following
    zpos[t+1,:,:]  = zpos[t,:,:] - np.reshape(interpolate.interpn((z,y,x), w, coor, method='linear', bounds_error=False, fill_value= 0)*time_step_length/vert_resolution, (num_seeds_z, num_seeds_x))
    zpos = zpos.clip(min=0) #Prevents z from being negative
    #terrain-height coordinates
    zpos_terrain[t+1,:,:]  = zpos_terrain[t,:,:] - np.reshape(interpolate.interpn((z,y,x), w, coor, method='linear', bounds_error=False, fill_value= 0)*time_step_length/vert_resolution, (num_seeds_z, num_seeds_x))
    zpos_terrain = zpos_terrain.clip(min=0) #Prevents z from being negative


    #####################   Get x,y,z for variables ###########################
    ##  This convoluted code gets the right location to match variable data  ##
    ##  with parcels.  I'm still not totally sure what's going on but it     ##
    ## appears to work.  The key issue is that the variable data is all      ##
    ## terrain-following.  u, v, and w is also terrain following but the     ##
    ## only modification need to be made are adding the surface elvation     ##
    ## before starting calculations.  Variable data must have surface height ##
    ## accounted for at every time step, hence this method.                  ##
    
    #Get surface height
    try:
        zs = np.array(ds.zs[0,:,:])
    except:
        zs = np.zeros((ds.yh.size, ds.xh.size))  
    
    #Get get x and y positions
    xloc = np.copy(xpos[t,:,:]).flatten()
    yloc = np.copy(ypos[t,:,:]).flatten()

    coor_terrain = []
    for i in range(len(xloc)):
        coor_terrain.append((yloc[i], xloc[i]))
    
    #Add the difference between the max surface height and the current surface 
    #height of all parcels to the terrain-following zpos
    zpos_var = np.copy(zpos)
    zpos_var[t,:,:] = zpos_var[t,:,:] + (np.max(zpos_terrain[0,0,:]-np.reshape(interpolate.interpn((y,x), zs, coor_terrain, method='linear', bounds_error=False, fill_value= 0)/vert_resolution, (num_seeds_z, num_seeds_x)))#-(15-height))
    
    #Add difference between max surface height and inital surface elevation added
    #to z
    xloc = np.copy(xpos[t,:,:]).flatten()
    yloc = np.copy(ypos[t,:,:]).flatten()
    zloc = np.copy(zpos_var[t,:,:]).flatten()-(np.max(zpos_terrain[0,0,:])-(zpos_terrain[0,:,:].flatten()-zpos[0,:,:].flatten()))

    coor_var = []
    for i in range(len(xloc)):
        coor_var.append((zloc[i], yloc[i], xloc[i]))  
    
    #Variables
    variable1[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var1, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan), (num_seeds_z, num_seeds_x))
    variable2[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var2, coor_var, method = 'linear', bounds_error=False, fill_value= np.nan), (num_seeds_z, num_seeds_x))
    
    #Budget variable
    #variable2[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var2, coor, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))

    
    
    

########  To get the last variable for colors ###############
t = time_steps-1

#Get surface height
try:
    zs = np.array(ds.zs[0,:,:])
except:
    zs = np.zeros((ds.yh.size, ds.xh.size))  

#Get get x and y positions
xloc = np.copy(xpos[t,:,:]).flatten()
yloc = np.copy(ypos[t,:,:]).flatten()

coor_terrain = []
for i in range(len(xloc)):
    coor_terrain.append((yloc[i], xloc[i]))

#Add the difference between the max surface height and the current surface 
#height of all parcels to the terrain-following zpos
zpos_var = np.copy(zpos)
zpos_var[t,:,:] = zpos_var[t,:,:] + (np.max(zpos_terrain[0,0,:]-np.reshape(interpolate.interpn((y,x), zs, coor_terrain, method='linear', bounds_error=False, fill_value= 0)/vert_resolution, (num_seeds_z, num_seeds_x)))#-(15-height))

#Add difference between max surface height and inital surface elevation added
#to z
xloc = np.copy(xpos[t,:,:]).flatten()
yloc = np.copy(ypos[t,:,:]).flatten()
zloc = np.copy(zpos_var[t,:,:]).flatten()-(np.max(zpos_terrain[0,0,:])-(zpos_terrain[0,:,:].flatten()-zpos[0,:,:].flatten()))

coor_var = []
for i in range(len(xloc)):
    coor_var.append((zloc[i], yloc[i], xloc[i]))  

#Color variable (theta)
variable1[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var1, coor, method = 'linear', bounds_error=False, fill_value= np.nan), (num_seeds_z, num_seeds_x))
variable2[t,:,:] = np.reshape(interpolate.interpn((z,y,x), var2, coor, method = 'linear', bounds_error=False, fill_value= np.nan)*time_step_length, (num_seeds_z, num_seeds_x))



 
    
#Save arrays
np.save('xpos_traj_disp', xpos)
np.save('ypos_traj_disp', ypos)
np.save('zpos_traj_disp', zpos_terrain)
np.save('%s_traj_disp' % var_name1, variable1)
np.save('%s_traj_disp' %var_name2, variable2)


###############################################################################
######################  Calculate displacement  ###############################
###############################################################################

z_disp = np.zeros((time_steps-2, num_seeds_z, num_seeds_x))
x_disp = np.zeros((time_steps-2, num_seeds_z, num_seeds_x))
y_disp = np.zeros((time_steps-2, num_seeds_z, num_seeds_x))
var1_disp = np.zeros((time_steps-2, num_seeds_z, num_seeds_x))
var2_disp = np.zeros((time_steps-2, num_seeds_z, num_seeds_x))


for t in range(time_steps-2):
    #z_disp[t,:,:] = zpos_terrain[0,:,:] - zpos_terrain[t+1,:,:]
    z_disp[t,:,:] = zpos[0,:,:] - zpos[t+1,:,:]
    x_disp[t,:,:] = xpos[0,:,:] - xpos[t+1,:,:]
    y_disp[t,:,:] = ypos[0,:,:] - ypos[t+1,:,:]
    var1_disp[t,:,:] = variable1[0,:,:] - variable1[t+1,:,:]
    var2_disp[t,:,:] = variable2[0,:,:] - variable2[t+1,:,:]

    #var2 is different bc it's a budget variable. It must be summed
    #at each time step to acquire total change
    #var2_disp[t,:,:] =  np.sum(variable2[:t+1,:,:], axis = 0)


mean_z_disp = np.zeros((num_seeds_z, num_seeds_x))
mean_z_disp = np.mean(z_disp, axis = 0)

mean_x_disp = np.zeros((num_seeds_z, num_seeds_x))
mean_x_disp = np.mean(x_disp, axis = 0)

mean_y_disp = np.zeros((num_seeds_z, num_seeds_x))
mean_y_disp = np.mean(y_disp, axis = 0)

mean_var1_disp = np.zeros((num_seeds_z, num_seeds_x))
mean_var1_disp = np.mean(var1_disp, axis = 0)

mean_var2_disp = np.zeros((num_seeds_z, num_seeds_x))
mean_var2_disp = np.mean(var2_disp, axis = 0)



################################################################################
###############Add zs to terrain hright in model has terrain####################
################################################################################

#Create new array with terrain height
zpos_terrain2 = np.copy(zpos)

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

    zpos_terrain2[t,:,:] = zpos_terrain2[t,:,:] + np.reshape(interpolate.interpn((y,x), zs, coor, method='linear', bounds_error=False, fill_value= 0)/vert_resolution, (num_seeds_z, num_seeds_x))

 



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
z = np.array(ds.zs[0,ymid,-num_seeds_x-rdist:-rdist])/1000*30 #Div by 1000 to go to m and mult by 30 to match y dim

x2d = np.zeros((num_seeds_z, num_seeds_x))
y2d = np.zeros((num_seeds_z, num_seeds_x))

for i in range(num_seeds_z):
    x2d[i,:] = x1d
for j in range(num_seeds_x):
    y2d[:,j] = y1d+z[j]
        

xmin = ds.th[0,0,0,:].size-num_seeds_x-rdist
xmax = ds.th[0,0,0,:].size-rdist
xlen = xmax-xmin

zmin = 0
zmax = num_seeds_z



##############################   Plot ########################################
    
fig = plt.figure(num=None, figsize=(18,9),  facecolor='w', edgecolor='k')
for j in range(1,3):
    subplot = 210 + j
    ax = plt.subplot(subplot,aspect = 'equal')
    plt.subplots_adjust(left=0.04, bottom=0.1, right=0.9, top=0.95, wspace=0, hspace=0)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    

    
    ##Levels
    zlmin = -10
    zlmax = 16.01
    zlevels = np.arange(zlmin,zlmax, 0.5)
    zlevels_ticks = np.arange(zlmin,zlmax,2)
    zlevels_ticks_labels = np.arange(zlmin,zlmax, 2).astype(int)
    
    xlmin = 260
    xlmax = 290
    xlevels = np.arange(xlmin,xlmax, 0.05)
    xlevels_ticks = np.arange(xlmin,xlmax,2)
    xlevels_ticks_labels = np.arange(xlmin,xlmax, 2).astype(int)
        
    #Plot mean vert disp
    if j == 1:
        z_disp_plot = plt.contourf(xpos[0,:,:]-xpos[0,0,0], zpos_terrain[0,:,:]*3, mean_var1_disp, zlevels,  cmap = cmap_dth, vmin = -zlmax, alpha = 1, zorder = 3)
        #z_disp_plot = plt.contourf(x1d, y1d, mean_var1_disp, zlevels,  cmap = cmap_dth, vmin = -zlmax, alpha = 1, zorder = 3)

    if j == 2:
        #x_disp_plot = plt.contourf(x2d, y2d, mean_var2_disp, xlevels,  cmap = cmap_var,vmin = -xlmax, alpha = 1, zorder = 3)
        x_disp_plot = plt.contourf(xpos[0,:,:]-xpos[0,0,0], zpos_terrain[0,:,:]*3, variable1[0,:,:], xlevels,  cmap = cmap_dth, alpha = 1, zorder = 3) #Percentage
        #x_disp_plot = plt.contourf(x1d, y1d, variable1[0,:,:], zlevels,  cmap = cmap_dth, alpha = 0.5, zorder = 4) #Percentage
    
    #Plot Terrain
    terrain = plt.plot(x1d, z, c = 'slategrey', linewidth = 10)
    
    #Plot Lake
    lake = np.array(ds.xland[0,ymid,-num_seeds_x-rdist:-rdist])
    lake[lake == 1] = np.nan
    lake_plt = plt.plot(x1d, lake-3, c = 'blue', linewidth = 6, zorder = 5)
    
    
    #Title
    if j == 2:
        sub_title = '[Run: 20 $\mathregular{ms^{-1}}$ and 1500m]'
        ax.text(590, -50, sub_title, fontsize = 20)
    
    
    ### Label and define grid area
    #y-axis
    plt.ylim([-3,np.max(y2d[:,0])+1])
    plt.yticks(np.arange(0, np.max(y2d[:,0]+1),30))
    ax.set_yticklabels(np.arange(0,np.max(y2d[:,0]+1),1).astype(int), fontsize = 15, zorder = 6)
    #x-axis
    plt.xticks(np.arange(12*5,xlen,100))
    ax.set_xticklabels(np.arange(xmin*0.2+12,xmax*0.2+12.1,20).astype(int), fontsize = 15)
    plt.xlim([-1,xlen-1])
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
        plt.text(0.25, -0.12, 'K', fontsize = 21)
        
    if j == 2:
        xcbaxes = fig.add_axes([0.92, 0.16, 0.03, 0.3])             
        xcbar = plt.colorbar(x_disp_plot, cax = xcbaxes, ticks = xlevels_ticks)
        xcbar.ax.set_yticklabels(xlevels_ticks_labels)
        xcbar.ax.tick_params(labelsize=15)
        plt.text(0.25, -0.12, 'm/s', fontsize = 21)
    
    #Labels
    if j == 1:
        sub_title = 'Time-mean change in Potential Temperature'
    if j == 2:
        sub_title = 'Potential Temperature'
    props = dict(boxstyle='square', facecolor='white', alpha=1)
    ax.text(10, 130, sub_title, fontsize = 20, bbox = props, zorder = 5)


plt.savefig("/uufs/chpc.utah.edu/common/home/u1013082/public_html/phd_plots/cm1/plots/time_mean_theta_change_theta_1500m_20ms_120sec_testt.png", dpi=350)
plt.close(fig)


#%%



###############################################################################
###################  Backward Trajectory Plot  ################################
###############################################################################

fig = plt.figure(num=None, figsize=(12,8), facecolor='w', edgecolor='k')


###################  Plot Trajectories in xy  ###############################

#Plot characteristics
ax = plt.subplot(211,aspect = 'equal')
plt.title("Backward Trajectories", fontsize = 18, y = 1.05) 
plt.subplots_adjust(left=0.07, bottom=0.11, right=0.87, top=0.93, wspace=0, hspace=0.15)

xlmin = 260
xlmax =288
xlevels = np.arange(xlmin,xlmax, 0.05)
xlevels_ticks = np.arange(xlmin,xlmax,2)
xlevels_ticks_labels = np.arange(xlmin,xlmax, 2).astype(int)

#Region of arrya to plot
colmin = 695#570
colmax = 696#620
rowmin = 0#30
rowmax = 40
tmin = 0
tmax =58
step = 1

#Plot
cmap = cm.viridis
trajectories_xy = plt.scatter(xpos[tmin:tmax,rowmin:rowmax:step, colmin:colmax:step], ypos[tmin:tmax,rowmin:rowmax:step, colmin:colmax:step], c = variable1[tmin:tmax,rowmin:rowmax:step, colmin:colmax:step], cmap = cmap_th, norm=matplotlib.colors.BoundaryNorm(xlevels,cmap.N), s = 2.5, zorder = 5)
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
trajectories_xz = plt.scatter(xpos[tmin:tmax,rowmin:rowmax:step, colmin:colmax:step].flatten(), zpos_terrain[tmin:tmax,rowmin:rowmax:step, colmin:colmax:step].flatten()*100, c = variable1[tmin:tmax,rowmin:rowmax:step, colmin:colmax:step].flatten(), cmap = cmap_th, norm=matplotlib.colors.BoundaryNorm(xlevels,cmap.N), s = 2.5, zorder = 3)



plt.grid(True, color = 'white', )
plt.ylabel('Height (m)', fontsize = 15)
plt.xlabel('Distance (km)', fontsize = 15)
plt.xticks(np.arange(0,ds.uinterp[0,0,0,:].size,200))
ytick = np.arange(0,5001,500)
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
cbar = plt.colorbar(trajectories_xy, cax = cbaxes, ticks = xlevels_ticks)
cbar.ax.set_yticklabels(xlevels_ticks_labels)
cbar.ax.tick_params(labelsize=10)
plt.text(-0.30, -0.08, 'Theta (K)', fontsize = 15)


plt.savefig("/uufs/chpc.utah.edu/common/home/u1013082/public_html/phd_plots/cm1/plots/backward_trajectory.png", dpi = 600)
plt.close(fig)
    



