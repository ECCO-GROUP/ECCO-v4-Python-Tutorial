#!/usr/bin/env python
# coding: utf-8

# # Combining multiple `Datasets` 
# 
# ## Objectives:
# 
# Show how to combine multiple ECCO v4 state estimate `Datasets` after loading.
# 
# ## Loading multiple `Datasets`  
# 
# In previous tutorials we've loaded single lat-lon-cap NetCDF tile files for ECCO state estimate variables and model grid parameters.  Here we will demonstrate merging the resulting `Datasets` together.  Some benefits of merging `Datasets` include having a tidier workspace and simplying subsetting operations, the subject of a future tutorial.  
# 
# First, we'll load three ECCOv4 NetCDF state estimate variables and the model grid.  For this exercise use the ECCOv4 NetCDF files in Format 2 (one NetCDF file corresponds to one variable and one year).  
# 
# First let's define our environment

# In[1]:


import numpy as np
import xarray as xr
import sys
import matplotlib.pyplot as plt
import json


# In[2]:


## Import the ecco_v4_py library into Python
## =========================================
#import ecco_v4_py as ecco

## -- If ecco_v4_py is not installed in your local Python library, 
##    tell Python where to find it.  For example, if your ecco_v4_py
##    files are in /Users/ifenty/ECCOv4-py/ecco_v4_py, then use:

sys.path.append('/Users/ifenty/ECCOv4-py')
import ecco_v4_py as ecco


# In[3]:


## Set top-level file directory for the ECCO NetCDF files
## =================================================================
# base_dir = '/home/username/'
base_dir = '/Users/ifenty/ECCOv4-release'

## define a high-level directory for ECCO fields
ECCO_dir = base_dir + '/Release3_alt'


# ### *c* point: ``SSH``

# In[4]:


SSH_dir= ECCO_dir + '/nctiles_monthly/SSH/'  
var = 'SSH'

# make use of the ``load_ecco_var_from_years_nc`` routine to 
# load a single variable 
ecco_dataset_A = ecco.load_ecco_var_from_years_nc(SSH_dir,                                                   var,                                                   tiles_to_load=[0,1,2,3],                                                  years_to_load=[2010])


# to see the data variables in a dataset, use ``.data_vars``:

# In[5]:


ecco_dataset_A.data_vars


# As expected, *ecco_dataset_A* has one ``data variable``, *SSH*, which has dimensions **i**, **j**, **tile**, and **time**.
# Its coordinates:

# In[6]:


ecco_dataset_A.SSH.coords


# ### *u* point: ``ADVx_TH``
# 
# ``ADVx_TH`` is the horizontal advective flux of potential temperature in each tile's $x$ direction.

# In[7]:


ADVx_dir= ECCO_dir + '/nctiles_monthly//ADVx_TH/'  
var = 'ADVx_TH'

# make use of the ``load_ecco_var_from_years_nc`` routine to 
# load a single variable 
ecco_dataset_B= ecco.load_ecco_var_from_years_nc(ADVx_dir,                                                  var,                                                  years_to_load=2010).load();
ecco_dataset_B.attrs=[]

ecco_dataset_B.data_vars


# *ecco_dataset_A* has one ``data variable``, *ADVx_TH*, which has dimensions **i_g**, **j**, **k** **tile**, and **time**.  *ADVx_TH* is a three-dimensional field.
# 
# Notice that the coordinates of *ecco_dataset_B* is distinct from *ecco_dataset_A*.  Specifically, *ecco_dataset_B* has a **k** dimension, and several new non-dimension coords such as **drF** and **dyG**

# In[8]:


ecco_dataset_B.coords


# ### *v* point: ``ADVy_TH``
# 
# ``ADVy_TH`` is the horizontal advective flux of potential temperature in each tile's $y$ direction.

# In[9]:


ADVy_dir= ECCO_dir + '/nctiles_monthly/ADVy_TH/'  
var = 'ADVy_TH'

ecco_dataset_C= ecco.load_ecco_var_from_years_nc(ADVy_dir,                                                  var,                                                  years_to_load=2010).load();
ecco_dataset_C.attrs=[]
ecco_dataset_C.data_vars


# *ecco_dataset_C* has one ``data variable``, *ADVy_TH*, which has dimensions **i**, **j_g**, **k** **tile**, and **time**.  *ADVy_TH* is a three-dimensional field.

# ### Examining the dimensions and coordinates of these `Datasets`

# Each of our three `Datasets` contain a single  `DataArray`.  Each of these `DataArray` objects has different horizontal dimension labels.  
# 
# * **i** and **j** for *SSH*
# * **i_g** and **j** for *ADVx_TH*
# * **i** and **j_g** for *ADVy_TH*
# 

# In[10]:


print(ecco_dataset_A.data_vars)
print(ecco_dataset_B.data_vars)
print(ecco_dataset_C.data_vars)


# ## Merging multiple `Datasets`  from state estimate variables
# 
# Using the `xarray` method ``merge`` we can create a single `Dataset` with multiple `DataArrays`.  

# In[11]:


# merge together
ecco_dataset_ABC = xr.merge([ecco_dataset_A, ecco_dataset_B, ecco_dataset_C]).load()


# ### Examining the merged `Dataset`
# 
# As before, let's look at the contents of the new merged `Dataset`

# In[12]:


ecco_dataset_ABC


# 
# #### 1. Dimensions
# `Dimensions:    (i: 90, i_g: 90, j: 90, j_g: 90, k: 50, nv: 2, tile: 13, time: 12)`
# 
# *ecco_dataset_merged* is a container of ``Data Arrays`` and as such it lists all of the unique dimensions its ``Data Arrays``. In other words, *Dimensions* shows all of the dimensions used by its variables. 
# 
# #### 2. Dimension Coordinates
# ```
# Coordinates:
# Coordinates:
#   * j          (j) int32 0 1 2 3 4 5 6 7 8 9 ... 80 81 82 83 84 85 86 87 88 89
#   * i          (i) int32 0 1 2 3 4 5 6 7 8 9 ... 80 81 82 83 84 85 86 87 88 89
#   * tile       (tile) int32 0 1 2 3 4 5 6 7 8 9 10 11 12
#   * time       (time) datetime64[ns] 2010-01-16T12:00:00 ... 2010-12-16T12:00:00
#   * i_g        (i_g) int32 0 1 2 3 4 5 6 7 8 9 ... 80 81 82 83 84 85 86 87 88 89
#   * k          (k) int32 0 1 2 3 4 5 6 7 8 9 ... 40 41 42 43 44 45 46 47 48 49
#   * j_g        (j_g) int32 0 1 2 3 4 5 6 7 8 9 ... 80 81 82 83 84 85 86 87 88 89
# ``` 
# 
# Notice that the **tile** and **time** coordinates are unchanged.  ``merge`` recognizes identical coordiantes and keeps them.
# 
# #### 3. Non-Dimension Coordinates
# ```
# Coordinates:
#     XC         (tile, j, i) float32 -111.60647 -111.303 -110.94285 ... nan nan
#     YC         (tile, j, i) float32 -88.24259 -88.382515 -88.52242 ... nan nan
#     rA         (tile, j, i) float32 362256450.0 363300960.0 ... nan nan
#     time_bnds  (time, nv) datetime64[ns] 2010-01-01 2010-02-01 ... 2011-01-01
#     iter       (time) int32 158532 159204 159948 160668 ... 165084 165804 166548
#     Z          (k) float32 -5.0 -15.0 -25.0 -35.0 ... -5039.25 -5461.25 -5906.25
#     PHrefC     (k) float32 49.05 147.15 245.25 ... 49435.043 53574.863 57940.312
#     drF        (k) float32 10.0 10.0 10.0 10.0 10.0 ... 387.5 410.5 433.5 456.5
#     dxC        (tile, j, i_g) float32 15583.418 15588.104 ... 23406.256
#     rAw        (tile, j, i_g) float32 361699460.0 362790240.0 ... 364760350.0
#     dyG        (tile, j, i_g) float32 23210.262 23273.26 ... 15595.26 15583.685
#     hFacW      (tile, k, j, i_g) float32 0.0 0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0
#     rAs        (tile, j_g, i) float32 179944260.0 180486990.0 ... 364150620.0
#     dxG        (tile, j_g, i) float32 15584.907 15589.316 ... 23142.107
#     dyC        (tile, j_g, i) float32 11563.718 11593.785 ... 15578.138
#     hFacS      (tile, k, j_g, i) float32 0.0 0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0
# ```
# 
# The list of non-dimension coordinates is now much longer.  Like *Dimensions*, the non-dimension coordinates of the merged ``Dataset`` contain all of the non-dimension coordinates of the ``Data Arrays``.
# 
# 
# #### 4. Attributes
# 
# Notice that all of the high level `Dataset` attributes are gone!  Each of the three `Datasets` had different attributes and the `merge` routine simply drops them.  The attributes of the *Data variables* remain intact:

# In[13]:


# (this json command makes Python dictionaries easier to read)
print(json.dumps(ecco_dataset_ABC.SSH.attrs, indent=2,sort_keys=True))


# ## Adding the model grid  `Dataset`
# 
# Let's use the ``merge`` routine to combine a `Dataset` of the model grid parameters with `output_merged`.
# 
# ### Load the model grid parameters

# In[14]:


# Load the llc90 grid parameters
grid_dir= ECCO_dir + '/nctiles_grid/'
grid_dataset = ecco.load_ecco_grid_nc(grid_dir, 'ECCOv4r3_grid.nc')
grid_dataset.coords


# ### Merge ``grid_all_tiles`` with ``output_merged``

# In[15]:


ecco_dataset_ABCG= xr.merge([ecco_dataset_ABC, grid_dataset])
ecco_dataset_ABCG


# ### Examining the merged `Dataset`
# 
# The result of this last merge is a single `Dataset` with 3 *Data variables*, and a complete set of model grid parameters (distances and areas).
# 
# ## Merging and memory
# 
# Merging `Datasets` together does not make copies of the data in memory.  Instead, merged `Datasets` are in fact just a reorganized collection of pointers.  You may want to delete the original variables to clear your namespace, but it is not necessary.

# ## Summary
# 
# Now you know how to merge multiple `Datasets` using the `merge` command.  We demonstrated merging of `Datasets` constructed from three different variables types and the model grid parameters.
