#!/usr/bin/env python
# coding: utf-8

# # Loading the ECCOv4 native model grid parameters 
# 
# ## Objectives
# 
# Introduce two methods of loading the ECCOv4 native model grid parameters.
# 
# ## Introduction
# 
# The ECCOv4 model grid parameters are provided as a single NetCDF file.  The file you have may look a little different than the one shown here because we have been working hard at improving how what exactly goes into our NetCDF files.  
# 
# This tutorial document is current as of Sep 2019 with the ECCOv4 NetCDF grid files provided in the following directories:
# 
# https://ecco.jpl.nasa.gov/drive/files/Version4/Release3_alt (po.daac drive, recommended)
# 
# https://web.corral.tacc.utexas.edu/OceanProjects/ECCO/ECCOv4/Release3_alt/  (mirror at U. Texas, Austin)

# ## Two methods to load the ECCOv4 model grid parameter NetCDF file 
# 
# Because the ECCOv4 model grid parameter data is provided in a single file you can use the ``open_dataset`` routine from ``xarray`` to open it. 
# 
# Alternatively, our subroutine ``load_ecco_grid_nc`` allows you to (optionally) specify a subset of vertical levels or a subset of tiles to load.
# 
# We'll show both methods.  Let's start with ``open_dataset``.
# 
# First set up your environment.

# In[1]:


import numpy as np
import xarray as xr
import sys
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')


# In[2]:


## Import the ecco_v4_py library into Python
## =========================================

## -- If ecco_v4_py is not installed in your local Python library, 
##    tell Python where to find it.  For example, if your ecco_v4_py
##    files are in /Users/ifenty/ECCOv4-py/ecco_v4_py, then use:

sys.path.append('/Users/ifenty/ECCOv4-py')
import ecco_v4_py as ecco


# In[3]:


## Set top-level file directory for the ECCO NetCDF files
## =================================================================

## -- If files are on a local machine, use something like 
# base_dir = '/Users/ifenty/'
base_dir = '/Users/ifenty/ECCOv4-release/'

## define a high-level directory for ECCO fields
ECCO_dir = base_dir + '/Release3_alt/'


# ### Method 1: Loading the model grid parameters using ``load_ecco_grid_nc``
# 
# Method 2 is super simple, just use ``open_dataset``:

# In[4]:


grid_dir = ECCO_dir + 'nctiles_grid/'

## load the grid
grid = xr.open_dataset(grid_dir + 'ECCOv4r3_grid.nc')
grid


# Let's plot two of the model grid parameter fields ``hFacC`` (tracer cell thickness factor) and ``rA`` (grid cell surface area)
# 
# First we plot ``hFac``:

# In[5]:


ecco.plot_tiles(grid.hFacC.sel(k=0), show_colorbar=True, cmap='gray');


# In[6]:


ecco.plot_tiles(grid.rA, show_colorbar=True);
'Model grid cell surface area [m^2]'


# ### Method 2: Loading the model grid parameters using ``load_ecco_grid_nc``
# 
# A more advanced routine, ``load_ecco_grid_nc``, allows you to load only a subset of tiles and vertical levels.  If no optional parameters are given, the entire grid object is loaded, just like ``open_dataset``

# In[7]:


grid_dir = ECCO_dir + 'nctiles_grid'

grid = ecco.load_ecco_grid_nc(grid_dir, 'ECCOv4r3_grid.nc')
grid


# Alternatively we can load just a subset of tiles and vertical levels.

# In[8]:


grid_subset = ecco.load_ecco_grid_nc(grid_dir, 'ECCOv4r3_grid.nc', tiles_to_load = [1, 10, 12], k_subset=[0,1,2,3])
grid_subset


# notice that ``grid_subset`` only has 3 tiles (9,10, 11) and 4 depth levels (0, 1, 2, 3), as expected.
# 
# Let's plot ``hFacC`` and ``rA`` again

# In[9]:


ecco.plot_tiles(grid_subset.hFacC.sel(k=0), show_colorbar=True, cmap='gray');
'Model grid cell surface area [m^2] in tiles 1, 10, and 12 '


# Notice that 10 of the 13 tiles are blank because they were not loaded.

# In[10]:


ecco.plot_tiles(grid_subset.rA, show_colorbar=True);
'Model grid cell surface area [m^2]'


# ## Summary
# 
# Now you know two ways to load the ECCOv4 grid parameter file!  
