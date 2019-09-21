#!/usr/bin/env python
# coding: utf-8

# # Saving `Datasets` and `DataArrays` to NetCDF
# 
# ## Objectives
# 
# Introduce an easy method for saving `Datasets` and `DataArrays` objects to NetCDF
# 
# ## Introduction
# 
# Saving your `Datasets` and `DataArrays` objects to NetCDF files couldn't be simpler.  The `xarray` module that we've been using to load NetCDF files provides methods for saving your `Datasets` and `DataArrays` as NetCDF files.
# 
# Here is the manual page on the subjet: http://xarray.pydata.org/en/stable/generated/xarray.Dataset.to_netcdf.html
# 
# The method `._to_netcdf( )` is available to **both** `Datasets` and `DataArrays` objects.  So useful!
# 
# ### Syntax
# ``
# your_dataset.to_netcdf('/your_filepath/your_netcdf_filename.nc')
# ``
# 
# ## Saving an existing `Dataset` to NetCDF
# 
# First, let's set up the environment and load a `Dataset`

# In[1]:


import numpy as np
import xarray as xr
import sys
import matplotlib.pyplot as plt
import json
import sys


# In[2]:


## Import the ecco_v4_py library into Python
## =========================================

#import ecco_v4_py as ecco

## -- If ecco_v4_py is not installed in your local Python library, 
##    tell Python where to find it.  For example, if your ecco_v4_py
##    files are in /Users/ifenty/ECCOv4-py/ecco_v4_py, then use:
sys.path.append('/Users/ifenty/ECCOv4-py')
import ecco_v4_py as ecco


# load a single tile, monthly averaged *THETA* for March 2010 for model tile 2. 

# In[3]:


## Set top-level file directory for the ECCO NetCDF files
## =================================================================
# base_dir = '/home/username/'
base_dir = '/Users/ifenty/ECCOv4-release'

## define a high-level directory for ECCO fields
ECCO_dir = base_dir + '/Release3_alt'


# In[4]:


## LOAD NETCDF FILE
## ================

# directory of the file
data_dir= ECCO_dir + '/nctiles_monthly/THETA/'

# filename
fname = 'THETA_2010.nc'

# load the dataset file using xarray
theta_dataset = xr.open_dataset(data_dir + fname).load()


# Now that we've loaded *theta_dataset*, let's save it in the **/tmp** file directory with a new name.

# In[5]:


new_filename_1 = './test_output.nc'
print ('saving to ', new_filename_1)
theta_dataset.to_netcdf(path=new_filename_1)
print ('finished saving')


# *It's really that simple!*
# 
# ## Saving a new custom ``Dataset`` to NetCDF
# 
# 
# Now let's create a new custom `Dataset` that with *THETA*, *SSH* and model grid parameter variables for a few tiles and depth level 10.

# In[6]:


data_dir= ECCO_dir + '/nctiles_monthly/'

SSH_THETA_201003 =     ecco.recursive_load_ecco_var_from_years_nc(data_dir,                                               ['SSH', 'THETA'],                                               tiles_to_load = [0,1,2],
                                              years_to_load = 2010)
grid_dir = ECCO_dir + '/nctiles_grid/'    
grid = ecco.load_ecco_grid_nc(grid_dir)
grid.close()

custom_dataset = xr.merge([SSH_THETA_201003, grid])


# and now we can easily save it:

# In[7]:


new_filename_2 = './test_output_2.nc'
print ('saving to ', new_filename_2)
custom_dataset.to_netcdf(path=new_filename_2)
custom_dataset.close()
print ('finished saving')


# In[8]:


custom_dataset


# ## Verifying our new NetCDF files
# 
# To verify that ``to_netcdf()`` worked, load them and compare with the originals.
# 
# ### Compare *theta_dataset* with *dataset_1*

# In[9]:


# the first test dataset
dataset_1 = xr.open_dataset(new_filename_1)

# release the file handle (not necessary but generally a good idea)
dataset_1.close()


# The `np.allclose` method does element-by-element comparison of variables

# In[10]:


# loop through the data variables in dataset_1
for key in dataset_1.keys():
    print ('checking %s ' % key)
    print ('-- identical in dataset_1 and theta_dataset : %s' %            np.allclose(dataset_1[key], theta_dataset[key], equal_nan=True))
    
# note: ``equal_nan`` means nan==nan (default nan != nan)


# *THETA* is the same in both datasets.
# 
# ### Compare *custom_dataset* with *dataset_2*

# In[11]:


# our custom dataset
dataset_2 = xr.open_dataset(new_filename_2)
dataset_2.close()
print ('finished loading')


# In[12]:


for key in dataset_2.keys():
    print ('checking %s ' % key)
    print ('-- identical in dataset_2 and custom_dataset : %s'           % np.allclose(dataset_2[key], custom_dataset[key], equal_nan=True))


# *THETA* and *SSH* are the same in both datasets.

# So nice to hear!

# ## Saving the results of calculations
# 
# ### Calculations in the form of `DataArrays`
# Often we would like to store the results of our calculations to disk.  If your operations are made at the level of `DataArray` objects (and not the lower `ndarray` level) then you can use these same methods to save your results.  All of the coordinates will be preserved (although attributes be lost).  Let's demonstrate by making a dummy calculation on SSH
# 
# $$SSH_{sq}(i) = SSH(i)^2$$

# In[13]:


SSH_sq = custom_dataset.SSH * custom_dataset.SSH

SSH_sq


# *SSH_sq* is itself a `DataArray`.

# Before saving, let's give our new *SSH_sq* variable a better name and descriptive attributes. 

# In[14]:


SSH_sq.name = 'SSH^2'
SSH_sq.attrs['long_name'] = 'Square of Surface Height Anomaly'
SSH_sq.attrs['units'] = 'm^2'

# Let's see the result
SSH_sq


# much better!  Now we'll save.

# In[15]:


new_filename_3 = './ssh_sq_DataArray.nc'
print ('saving to ', new_filename_3)

SSH_sq.to_netcdf(path=new_filename_3)
print ('finished saving')


# 
# ### Calculations in the form of `numpy ndarrays`
# 
# If calculations are made at the `ndarray` level then the results will also be `ndarrays`.

# In[16]:


SSH_dummy_ndarray = custom_dataset.SSH.values *  custom_dataset.SSH.values

type(SSH_dummy_ndarray)


# You'll need to use different methods to save these results to NetCDF files, one of which is described here: http://pyhogs.github.io/intro_netcdf4.html

# ## Summary
# 
# Saving `Datasets` and `DataArrays` to disk as NetCDF files is fun and easy with ``xarray``!
