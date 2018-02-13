
# coding: utf-8

# # Combining multiple `Datasets` 
# 
# ## Objectives:
# 
# Show how to combine the multiple `Datasets` that are created each time an ECCO v4 state estimate variable is loaded.
# 
# ## Loading multiple `Datasets`  
# 
# In previous tutorials we've loaded lat-lon-cap NetCDF tile files for different state estimate variables and model grid parameters.  Here we will show you how to merge the resulting `Datasets` together.  Some benefits of merging `Datasets` include having a tidier workspace and simplying subsetting operations, the subject of a future tutorial.  
# 
# First, load three `Datasets` for state estimate variables corresponding to 'c','u', and 'v' c-grid points.

# In[1]:


import matplotlib.pylab as plt
import numpy as np
import sys
import xarray as xr
from copy import deepcopy 
import ecco_v4_py as ecco


# ### 'c' point: ``SSH``

# In[2]:


data_dir='/Volumes/ECCO_BASE/ECCO_v4r3/nctiles_monthly/SSH/'    
var = 'SSH'
var_type = 'c'
ssh_all_tiles = ecco.load_all_tiles_from_netcdf(data_dir, 
                                                var, var_type)
ecco.minimal_metadata(ssh_all_tiles)


# ### 'u' point: ``ADVxSNOW``
# 
# ``ADVxSNOW`` is the horizontal advective flux of snow in each tile's $x$ direction.

# In[3]:


data_dir='/Volumes/ECCO_BASE/ECCO_v4r3/nctiles_monthly/ADVxSNOW/'    
var = 'ADVxSNOW'
var_type = 'u'
advxsnow_all_tiles= ecco.load_all_tiles_from_netcdf(data_dir, 
                                                    var, var_type)
ecco.minimal_metadata(advxsnow_all_tiles)


# ### 'v' point: ``ADVySNOW``
# 
# ``ADVySNOW`` is the horizontal advective flux of snow in each tile's $y$ direction.

# In[4]:


data_dir='/Volumes/ECCO_BASE/ECCO_v4r3/nctiles_monthly/ADVySNOW/'    
var = 'ADVySNOW'
var_type = 'v'
advysnow_all_tiles= ecco.load_all_tiles_from_netcdf(data_dir, 
                                                    var, var_type)
ecco.minimal_metadata(advysnow_all_tiles)


# ### Examining the dimensions and coordinates of these `Datasets`

# Recall that these three `Datasets` are each storing just one `DataArray` and that each  `DataArray` uses a different coordinate system.

# In[5]:


ssh_all_tiles.SSH.dims


# In[6]:


advysnow_all_tiles.ADVySNOW.dims


# In[7]:


advxsnow_all_tiles.ADVxSNOW.dims


# ## Merging multiple `Datasets`  from state estimate variables
# 
# Using the `xarray` method `merge` we can create a single `Dataset` that stores multiple `DataArrays`.  

# In[8]:


output_merged = xr.merge([ssh_all_tiles, advxsnow_all_tiles, advysnow_all_tiles])


# In[9]:


output_merged


# ### Examining the merged `Dataset`
# 
# As before, let's look at the contents of the new merged `Dataset`
# 
# #### 1. Dimensions
# `Dimensions:   (i: 90, i_g: 90, j: 90, j_g: 90, tile: 13, time: 288)`
# 
# The new *Dimensions* list is the union of the dimensions of the three merged variables. 
# 
# #### 2. Coordinates
# ```
# Coordinates:
#   * time      (time) float64 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 ...
#   * j         (j) float64 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 ...
#   * i         (i) float64 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 ...
#     tim       (time) datetime64[ns] 1992-01-16 1992-02-16 1992-03-16 ...
#     timestep  (time) float64 732.0 1.428e+03 2.172e+03 2.892e+03 3.636e+03 ...
#     lon_c     (tile, j, i) float64 -111.6 -111.3 -110.9 -110.5 -110.0 -109.3 ...
#     lat_c     (tile, j, i) float64 -88.24 -88.38 -88.52 -88.66 -88.8 -88.94 ...
#   * tile      (tile) int64 1 2 3 4 5 6 7 8 9 10 11 12 13
#   * i_g       (i_g) float64 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 ...
#     lon_u     (tile, j, i_g) float64 -115.0 -115.0 -115.0 -115.0 -115.0 ...
#     lat_u     (tile, j, i_g) float64 -88.24 -88.38 -88.52 -88.66 -88.8 ...
#   * j_g       (j_g) float64 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 ...
#     lon_v     (tile, j_g, i) float64 -111.6 -111.3 -110.9 -110.5 -110.0 ...
#     lat_v     (tile, j_g, i) float64 -88.18 -88.32 -88.46 -88.6 -88.74 ...
# ``` 
# 
# Notice the three ``lon_*`` and ``lat_*`` *non-dimensional* coordinates are also merged.
# 
# 
# #### 3. Data Variables
# ```
# Data variables:
#     SSH       (time, tile, j, i) float64 nan nan nan nan nan nan nan nan nan ...
#     ADVxSNOW  (time, tile, j, i_g) float64 nan nan nan nan nan nan nan nan ...
#     ADVySNOW  (time, tile, j_g, i) float64 nan nan nan nan nan nan nan nan ...
# ```
# 
# The *Data variables* are now the three `Data Arrays`.  
# 
# #### 4. Attributes
# 
# Notice that all of the high level `Dataset` attributes are gone!  Each of the three `Datasets` had different attributes and the `merge` routine simply drops them.  The attributes of the *Data variables* remain intact:

# In[10]:


print output_merged.SSH.attrs['long_name']


# In[11]:


print output_merged.SSH.attrs['units']


# ## Merging ``GRID`` parameter  `Dataset`
# 
# Let's use the `merge` routine to combine a 13-tile `Dataset` holding the model grid parameters with out merged `output_merged`.
# 
# ### Load the model grid `Dataset`

# In[12]:


# specify the location of your nctiles_grid directory
grid_dir='/Volumes/ECCO_BASE/ECCO_v4r3/nctiles_grid/'
var = 'GRID'
var_type = 'grid'
grid_all_tiles = ecco.load_all_tiles_from_netcdf(grid_dir, var, var_type)

# minimize the metadata (optional)
ecco.minimal_metadata(grid_all_tiles)


# ### Merge ``grid_all_tiles`` with ``output_merged``

# In[13]:


output_merged = xr.merge([output_merged, grid_all_tiles])


# In[14]:


output_merged


# ### Examining the merged `Dataset`
# 
# The result of this final merging is a single `Dataset` with 20+ variables, a combination of state estimate output variables and grid parameters.
# 
# ## Merging and memory
# 
# Merging `Datasets` together does not make copies of the data in memory.  Instead, merged `Datasets` are in fact just a reorganized collection of pointers.  You may want to delete the original variables to clear your name, but it is not necessary.

# In[15]:


del ssh_all_tiles
del advxsnow_all_tiles
del advysnow_all_tiles
del grid_all_tiles


# In[16]:


whos


# ## Conclusion
# 
# Now you know how to merge multiple `Datasets` using the `merge` command.  We demonstrated merging of `Datasets` constructed from all 13 tile files from three different variables and the model grid parameters.  `merge` can just as easily combine `Dataset` variables loaded from individual tiles or even the ECCO v4 R3 interpolated lat-lon output fields.
