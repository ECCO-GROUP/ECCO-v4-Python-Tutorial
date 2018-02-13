
# coding: utf-8

# # Loading all 13 lat-lon-cap NetCDF tile files at once
# 
# ## Objectives
# 
# Introduce an easy method for loading all 13 llc NetCDF tile files at once.  
# 
# ## Loading all 13 lat-lon-cap NetCDF grid tile files using  `load_all_tiles_from_netcdf`
# 
# So far we've loaded NetCDF tile files one tile at a time.  The `load_all_tiles_from_netcdf` method automates the loading of all 13 tiles.
# 
# Let's jump right in and use `load_all_tiles_from_netcdf` to load all 13 GRID tile files.  Call the new `Dataset` object `grid_all_tiles`.   Because we are loading GRID tile files, we specify 'grid' as the *var_type*.

# In[1]:


import matplotlib.pylab as plt
import numpy as np
import sys
import xarray as xr
from copy import deepcopy 
import ecco_v4_py as ecco

# specify the locaiotn of your nctiles_grid directory
grid_dir='/Volumes/ECCO_BASE/ECCO_v4r3/nctiles_grid/'
var = 'GRID'
var_type = 'grid'
grid_all_tiles = ecco.load_all_tiles_from_netcdf(grid_dir, var, var_type)

# minimize the metadata (optional)
ecco.minimal_metadata(grid_all_tiles)


# Let's look at `grid_all_tiles`:

# In[2]:


grid_all_tiles


# ### Examining the Dataset object contents
# 
# #### 1. Dimensions
# `Dimensions:  (i: 90, i_g: 90, j: 90, j_g: 90, k: 50, k_g: 50, tile: 13)`
# 
# The *Dimensions* list now inculdes a new **tile** dimension.
# 
# #### 2. Coordinates
# ```
# Coordinates:
#   * k        (k) float64 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 ...
#   * i        (i) float64 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 ...
#   * j        (j) float64 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 ...
#   * i_g      (i_g) float64 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 ...
#   * j_g      (j_g) float64 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 ...
#   * k_g      (k_g) float64 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 ...
#   * tile     (tile) int64 1 2 3 4 5 6 7 8 9 10 11 12 13
# ``` 
# 
# The **tile** coordinate now appears as an array: 1 .. 13
# 
# 
# #### 3. Data Variables
# ```
# Data variables:
#     XC       (tile, j, i) float64 -111.6 -111.3 -110.9 -110.5 -110.0 -109.3 ...
#     hFacC    (tile, k, j, i) float64 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ...
#     XG       (tile, j_g, i_g) float64 -115.0 -115.0 -115.0 -115.0 -115.0 ...
#     DXC      (tile, j, i_g) float64 1.558e+04 1.559e+04 1.559e+04 1.56e+04 ...
#     hFacW    (tile, k, j, i_g) float64 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ...
#     DYC      (tile, j_g, i) float64 1.156e+04 1.159e+04 1.162e+04 1.165e+04 ...
#     hFacS    (tile, k, j_g, i) float64 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ...
#     RC       (k) float64 -5.0 -15.0 -25.0 -35.0 -45.0 -55.0 -65.0 -75.0 ...
#     RF       (k_g) float64 0.0 -10.0 -20.0 -30.0 -40.0 -50.0 -60.0 -70.0 ...
# ```
# 
# Each *Data variable* that varies in the horizontal now has a new dimension, **tile**. 
# 
# > **Note:** The ordering of the arrays is (tile, k, j, i)

# ### A closer look at the DataArray after loading with `load_all_tiles_from_netcdf`
# 
# Let's look at one of the grid variables loaded with the `load_all_tiles_from_netcdf` routine:

# In[3]:


grid_all_tiles.XC


# The ``XC`` `DataArray` object now has a **tile** dimension of length 13.

# ## Loading non-grid variables using `load_all_tiles_from_netcdf`
# 
# As with `load_tile_from_netcdf` we can use `load_all_tiles_from_netcdf` to give variables better coordinate labels by specifying their c-grid point category.
# 
# We'll demonstrate by loading the ``SSH`` variable that is on c-grid 'c' points.  Loading variables on 'g', 'u',or 'v' points only requires changing the argument of **var_type**.
# 
# ### A 'c' point variable: SSH
# 
# Let's load all tiles of the 'c' point sea surface height (SSH) variable.  

# In[4]:


data_dir='/Volumes/ECCO_BASE/ECCO_v4r3/nctiles_monthly/SSH/'    
var = 'SSH'
var_type = 'c'
ssh_all_tiles = ecco.load_all_tiles_from_netcdf(data_dir, var, var_type)
ecco.minimal_metadata(ssh_all_tiles)


# ### Examining the Dataset object contents
# 
# Let's examine `ssh_all_tiles` since this is first time we're loading all 13 tiles of an output variable.

# In[5]:


ssh_all_tiles


# #### 1. Dimensions
# `Dimensions:   (i: 90, j: 90, tile: 13, time: 288)`
# 
# The *Dimensions* list now inculdes a new **tile** dimension of 13.
# 
# #### 2. Coordinates
# ```
# Coordinates:
#     lon_c     (tile, j, i) float64 -111.6 -111.3 -110.9 -110.5 -110.0 -109.3 ...
#     lat_c     (tile, j, i) float64 -88.24 -88.38 -88.52 -88.66 -88.8 -88.94 ...
#   * tile      (tile) int64 1 2 3 4 5 6 7 8 9 10 11 12 13
# ``` 
# 
# The **tile** coordinate now appears as an array: 1 .. 13 and **lon_c** and **lat_c** have a new **tile** dimension.
# 
# 
# #### 3. Data Variables
# ```
# Data variables:
#     SSH       (time, tile, j, i) float64 nan nan nan nan nan nan nan nan nan ...
# ```
# 
# The SSH *Data variable* now has a new **tile** dimension.

# ## Conclusion
# 
# Now you know how to load all 13 NetCDF tile files for grid parameters and output variables using `load_all_tiles_from_netcdf`.
