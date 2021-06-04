#!/usr/bin/env python
# coding: utf-8

# # ECCOv4 Loading llc binary files in the 'compact' format
# 
# This section demonstrates loading ECCO binary files written in the so-called llc 'compact' format. Compact binary files have a non-intuitive layout of the 13 lat-lon-cap tiles required by the MITgcm.  This tutorial demonstrates using the 'read_llc_to_tiles' routine to read and re-organize the llc compact binary files into a more familiar 13-tile layout.

# ## Objectives
# 
# By the end of the tutorial you will be able to read llc compact binary files of any dimension, plot them, and convert them into DataArrays.

# In[1]:


## Import the ecco_v4_py library into Python
## =========================================

## -- If ecco_v4_py is not installed in your local Python library, 
##    tell Python where to find it.  For example, if your ecco_v4_py
##    files are in /Users/ifenty/ECCOv4-py/ecco_v4_py, then use:
import sys

sys.path.append('/Users/ifenty/git_repos/my_forks/ECCOv4-py')
import ecco_v4_py as ecco

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr


# ## The *read_llc_to_tiles* subroutine
# 
# *read_llc_to_tiles* reads a llc compact format binary file and converts to a numpy ndarray of dimension:
# [N_recs, N_z, N_tiles, llc, llc]
# 
# For ECCOv4 our convention is:
# ```
# 'N_recs' = number of time levels
# 'N_z' = number of depth levels
# 'N_tiles' = 13
# 'llc' = 90
# ```
# 
# By default the routine will try to load a single 2D slice of a llc90 compact binary file: (N_recs = 1, N_z =1, N_tiles = 13, and llc=90).  
# 
# There are several other options which you can learn about using the 'help' command:

# In[2]:


help(ecco.read_llc_to_tiles)


# ## Related routines
# 
# Two related routines which you might find useful: 
# 
# 1. *read_llc_to_compact*: Loads an MITgcm binary file in the 'compact' format of the 
#     lat-lon-cap (LLC) grids and preserves its original dimension
#     
# 2.  *read_llc_to_faces* : Loads an MITgcm binary file in the 'compact' format of the 
#     lat-lon-cap (LLC) grids and converts it to the '5 faces' dictionary. 
#     
# For the remainder of the tutorial we will only use *read_llc_to_tiles*.

# ## Example 1: Load a 2D llc 'compact' binary file
# 
# The file 'bathy_eccollc_90x50_min2pts.bin' contains the 2D array of bathymetry for the model. 

# In[3]:


input_dir = '/Users/ifenty/tmp/input_init/'
input_file = 'bathy_eccollc_90x50_min2pts.bin'


# *read_llc_to_tiles* actually runs several other subroutines: *load_binary_array* which does the lower level reading of the binary file, *llc_compact_to_faces* which converts the array to 5 'faces', and finally *llc_faces_to_tiles* which extracts the 13 tiles from the 5 faces:

# In[4]:


bathy = ecco.read_llc_to_tiles(input_dir, input_file)


# *bathy* is a numpy float32 array with dimension [13, 90, 90]

# ### Plot the 13 tiles bathymetry data 

# In[5]:


# Use plot_tiles to make a quick plot of the 13 tiles.  See the tutorial on plotting for more examples.
ecco.plot_tiles(bathy, layout='latlon',rotate_to_latlon=True,show_tile_labels=False, show_colorbar=True);


# ## Load ecco-grid information to make a fancier lat-lon plot
# 
# With the longitudes (XC) and latitudes (YC) and the 13 tile ndarray we can plot the field in different geographic projections.  See the tutorial on plotting for more examples.

# In[6]:


ecco_grid_dir = '/Users/ifenty/tmp/nctiles_grid/'
ecco_grid = ecco.load_ecco_grid_nc(input_dir, 'ECCO-GRID.nc')


# In[7]:


plt.figure(figsize=(15,6));
ecco.plot_proj_to_latlon_grid(ecco_grid.XC, ecco_grid.YC, bathy, show_colorbar=True, user_lon_0=-66);


# In[8]:


plt.figure(figsize=(12,8));
ecco.plot_proj_to_latlon_grid(ecco_grid.XC, ecco_grid.YC, bathy, show_colorbar=True,                              projection_type='stereo', lat_lim =-45, less_output=True,dx=.25,dy=.25);


# ## Convert the ndarray into a DataArray
# 
# Converting the ndarray to a DataArray can be useful for broadcasting operations.  One approach is to create the DataArray manually by specifying the names and values of the new dimension coordinates: 

# In[9]:


tile = range(13)
i = range(90)
j = range(90)


# In[10]:


# Convert numpy array to an xarray DataArray with matching dimensions as the monthly mean fields
bathy_DA = xr.DataArray(bathy,coords={'tile': tile,
                                         'j': j,
                                         'i': i}, dims=['tile','j','i'])


# In[11]:


print(bathy_DA.dims)


# Another approach is to use the routine *llc_tiles_to_xda*.  *llc_tiles_to_xda* uses 

# In[12]:


bathy.shape


# In[13]:


bathy_DA2 = ecco.llc_tiles_to_xda(bathy, var_type='c',grid_da=ecco_grid.XC)
print(bathy_DA2.dims)
print(bathy_DA2.coords)


# ## Example 2: Load a 3D 'compact' llc binary file with 3rd dimension = Time
# 
# 
# The file 'runoff-2d-Fekete-1deg-mon-V4-SMOOTH.bin' contains the 12 month climatology of river runoff, dimensions of [time, j, i].  

# In[14]:


input_file = 'runoff-2d-Fekete-1deg-mon-V4-SMOOTH.bin'


# specify the length of the n_recs dimension, nl, as 12. By default, nk = 1

# In[15]:


runoff = ecco.read_llc_to_tiles(input_dir, input_file, nl=12)


# Shape is n_recs (12), n_z (1), n_tiles (13), llc (90), llc (90)

# In[16]:


print(runoff.shape)


# ### Plot the November runoff climatology

# In[17]:


plt.figure(figsize=(15,6));
ecco.plot_proj_to_latlon_grid(ecco_grid.XC, ecco_grid.YC, runoff[10,:], cmin=0,cmax=1e-7,                              show_colorbar=True, cmap='bwr',user_lon_0=-66);


# ## Convert the ndarray into a DataArray
# 
# Converting the ndarray to a DataArray can be useful for broadcasting operations.  Two methods: 'manual' and using *llc_tiles_to_xda*.
# 
# ### Method 1: Manual

# In[18]:


tile = range(1,14)
i = range(90)
j = range(90)
time = range(12) # months
k = [0];


# In[19]:


# Convert numpy array to an xarray DataArray with matching dimensions as the monthly mean fields
runoff_DA = xr.DataArray(runoff,coords={'time': time,
                                        'k': k,
                                        'tile': tile,
                                        'j': j,
                                        'i': i}, dims=['time','k','tile','j','i'])


# In[20]:


print(runoff_DA.dims)
print(runoff_DA.shape)


# ### Method 2: *llc_tiles_to_xda*
# 
# runoff is a 3D array with the 3rd dimension being time and therefore we need to pass *llc_tiles_to_xda* a similarly-dimensioned DataArray or tell the subroutine that the 3rd dimension is depth. The 'ecco_grid' DataSet doesn't have any similary-dimensioned DataArrays (no DataArrays with a time dimension). Therefore we will tell the routine that the new 4th dimension should be time:

# In[21]:


##### specify that the 5th dimension should be time
runoff_DA2 = ecco.llc_tiles_to_xda(runoff, var_type='c',dim4='depth', dim5='time')
print(runoff_DA2.dims)
print(runoff_DA2.coords)


# ## Example 3: Load a 3D 'compact' llc binary file with 3rd dimension = Depth
# 
# 
# The file 'total_kapredi_r009bit11.bin' is a 50 depth level array of the adjusted Redi parameter from Release 1 (first guess + adjustments), dimensions of [depth, j, i]
# Release 4's adjusted Redi parameter can be found at https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-protected/ECCO_L4_OCEAN_3D_MIX_COEFFS_LLC0090GRID_V4R4/OCEAN_3D_MIXING_COEFFS_ECCO_V4r4_native_llc0090.nc
#
# In[22]:


input_file = 'total_kapredi_r009bit11.bin'


# specify the number of depth levels as 50.  n_recs defaults to 1 and is dropped by default.

# In[23]:


kapredi = ecco.read_llc_to_tiles(input_dir, input_file, nk=50)


# In[24]:


print(kapredi.shape)


# ### Plot log10 of the parameter at the 10th depth level (105m)

# In[25]:


plt.figure(figsize=(15,6));
ecco.plot_proj_to_latlon_grid(ecco_grid.XC, ecco_grid.YC, np.log10(kapredi[10,:]),                              cmin=2,cmax=4,show_colorbar=True);


# ## Convert the ndarray into a DataArray
# 
# Converting the ndarray to a DataArray can be useful for broadcasting operations.  Two methods: 'manual' and using *llc_tiles_to_xda*.
# 
# ### Method 1: Manual

# In[26]:


tile = range(1,14)
i = range(90)
j = range(90)
k = range(50)


# In[27]:


# Convert numpy array to an xarray DataArray with matching dimensions as the monthly mean fields
kapredi_DA = xr.DataArray(kapredi,coords={'k': k,
                                          'tile': tile,
                                          'j': j,
                                          'i': i},dims=['k','tile','j','i'])


# In[28]:


print(kapredi_DA.dims)
print(kapredi_DA.shape)


# ### Method 2: *llc_tiles_to_xda*
# 
# kapredi is a 4D array (depth, tile, j, i) and therefore we need to pass *llc_tiles_to_xda* a similarly-dimensioned DataArray or tell the subroutine that the 4th dimension is depth

# In[29]:


# use similarly-dimensioned DataArray: hFacC
kapredi_DA2 = ecco.llc_tiles_to_xda(kapredi, var_type='c',grid_da=ecco_grid.hFacC)
print(kapredi_DA2.dims)
print(kapredi_DA2.coords)


# In[30]:


# specify that the 4th dimension should be depth
kapredi_DA3 = ecco.llc_tiles_to_xda(kapredi, var_type='c', dim4='depth')
print(kapredi_DA3.dims)
print(kapredi_DA3.coords)


# ## Parting thoughts
# 
# 1. *read_llc_to_tiles* can also be used to read ECCO '*.data'* generated when re-running the ECCO model.
# 
# 1. Converting from numpy ndarrays or xarray DataArrays in the 'tile' format to a 5-faces format or compact format can be made with routines like: *llc_tiles_to_compact* and *llc_tiles_to_faces*
# 
# ### *llc_tiles_to_faces*

# In[31]:


help(ecco.llc_tiles_to_faces)


# ### *llc_tiles_to_compact*

# In[32]:


help(ecco.llc_tiles_to_compact)

