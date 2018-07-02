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
time_step_length = 90 #seconds
###############################################################################
###############################################################################
#%%
 

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
    
    

for i in range(200,201):    
    secs = (i*45)
 
    fig = plt.figure(num=None, figsize=(14,7), facecolor='w', edgecolor='k')

    ax = plt.axes([-0.2, 0, 1.2, 1], projection='3d')
    ax.set_aspect(1)
    plt.title("CM1 Model Domain", fontsize = 20, y = 1) 

    
    ax.view_init(25, -90)
    
    #terrain plot
    X, Y = np.meshgrid(x, y)
    Z = ds.zs[0,:,:]-110
    ax.plot_surface(X, Y, Z, rstride=10, cstride=10,  alpha = 1, cmap = cm.Greys, vmin = -800, vmax = 2000, zorder = 1)
    
    #water and land details
    lake = np.array(ds.xland[0,:,:])
    lake[lake == 1] = 'nan'
    ax.plot_surface(X, Y, lake-110, rstride=5, cstride=5, alpha = 0.5,  color = ('lightskyblue'), zorder = 2)
    
    #xtick
    xtick = np.arange(0,2208,200).astype(int)
    ax.set_xticks(xtick)
    ax.set_xticklabels(xtick*150/1000)
    
    #yticks
    ytick = np.arange(0,661,100).astype(int)
    ax.set_yticks(ytick)
    ax.set_yticklabels(ytick*150/1000, horizontalalignment='left', verticalalignment='baseline')
    
    #zticks
    ztick = np.arange(500,5000.01,500).astype(int)
    ax.set_zticks(ztick)
    ax.set_zticklabels(ztick, horizontalalignment='left')
    
    #tick labels
    ax.set_xlabel('Distance (km)', labelpad = 30, fontsize = 14)
    ax.set_ylabel('Distance (km)            ', labelpad = 15, fontsize = 12)
    ax.set_zlabel('Height (m)', labelpad = 22, fontsize = 14, rotation = -2)

    #Limits
    ax.set_xlim3d(30,2160)
    ax.set_zlim3d(50,5000.01)
    ax.set_ylim3d(0,661)
    
    
    plt.savefig("/uufs/chpc.utah.edu/common/home/u1013082/public_html/phd_plots/cm1/plots/cm1_model_domain", dpi = 200)
    plt.close(fig)
    print i






