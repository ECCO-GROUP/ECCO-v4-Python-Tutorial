
# coding: utf-8

# # Saving `Datasets` and `DataArrays` to NetCDF
# 
# ## Objectives
# 
# Introduce an easy method for saving `Datasets` and `DataArrays` objects to NetCDF
# 
# ## Introduction
# 
# Saving your `Datasets` and `DataArrays` objects to NetCDF couldn't be simpler.  The `xarray` module that we've been using to load NetCDF files from disk provides methods for direct writing to NetCDF.
# 
# Here is the manual page on the subjet: http://xarray.pydata.org/en/stable/generated/xarray.Dataset.to_netcdf.html
# 
# The method is called `._to_netcdf( )` and it available to both `Datasets` and `DataArrays` objects.
# 
# ### Syntax
# ``
# your_dataset.to_netcdf('/your_filepath/your_netcdf_filename.nc')
# ``
# 
# ## Saving a `Dataset`
# 
# First, let's load a `Dataset`

# In[1]:


import matplotlib.pylab as plt
import numpy as np
import sys
import xarray as xr
from copy import deepcopy 

get_ipython().magic(u'matplotlib inline')
sys.path.append('/Users/ifenty/Documents/Work/git_repo/ECCO-GROUP/ECCOv4-py')
sys.path.append('/Users/ifenty/Documents/Work/git_repo/ECCO-GROUP/ECCOv4-py/ecco_v4_py/')
import ecco_v4_py as ecco

ECCO_dir = '/Users/ifenty/ECCOv4/R3'

# Load all tiles of the LLC90 Grid    
data_dir= ECCO_dir + '/nctiles_grid/'    
var = 'GRID'
var_type = 'grid'

grid_all_tiles = ecco.load_all_tiles_from_netcdf(data_dir, 
                                                 var, var_type,
                                                 less_output=True)

# Load all tiles of SSH
data_dir= ECCO_dir + '/nctiles_monthly/SSH/'    
var = 'SSH'
var_type = 'c'
ssh_all_tiles = ecco.load_all_tiles_from_netcdf(data_dir, 
                                                var, var_type,
                                                less_output=True)

# minimize the metadata (optional)
data = xr.merge([ssh_all_tiles, grid_all_tiles])


# ### Saving a `Dataset`
# 
# Now that we've loaded *ssh_all_tiles*, let's save it in the *SSH* file directory.

# In[2]:


new_filename = data_dir + 'data_all_tiles.nc'
print 'saving to ', new_filename

data.to_netcdf(path=new_filename)
print 'finished saving'


# Now let's create a new `Dataset` that only including *SSH* and some grid parameter variables that are on the same 'c' grid points as *SSH*.
# 
# First, make three new `Datasets`, one for each variable.

# In[3]:


# Extract these three DataArrays from the data object
# and simultaneously convert them to one-variable Dataset objects
SSH = data.SSH.to_dataset(name = 'SSH')
XC  = data.XC.to_dataset(name = 'XC')
YC  = data.YC.to_dataset(name = 'YC')

# Merge these Datasets
data_subset = xr.merge([SSH, XC , YC])

# Give data_subset the same metadata attributes as data
data_subset.attrs = data.attrs

# Examine the results
data_subset


# and now we can easily save it:

# In[4]:


new_filename = data_dir + 'data_subset_all_tiles.nc'
print 'saving to ', new_filename

data_subset.to_netcdf(path=new_filename)
print 'finished saving'


# ### Loading our saved  `Dataset`
# 
# To verify that our worked let's load it up and compare with *data_subset*

# In[5]:


new_filename = data_dir + 'data_subset_all_tiles.nc'

data_subset_loaded = xr.open_dataset(new_filename)

print 'finished loading'


# In[6]:


data_subset_loaded


# The `xarray` equals method does element by element comparison between two `Dataset` and `DataArray` objects

# In[7]:


data_subset.equals(data_subset_loaded)


# So nice to hear!

# ### Saving a `DataArray`
# 
# Let's save the *SSH* `DataArray` object inside *ssh_all_tiles*

# In[8]:


new_filename = data_dir + 'ssh_dataArray.nc'
print 'saving to ', new_filename

ssh_all_tiles.to_netcdf(path=new_filename)
print 'finished saving'


# ## Saving the results of calculations
# 
# ### Calculations in the form of `DataArrays`
# Often we would like to store the results of our calculations to disk.  If your operations are made at the level of `DataArray` objects (and not the lower `ndarray` level) then you can use these same methods to save your results.  All of the coordinates will be preserved (although attributes be lost).  Let's demonstrate by making a dummy calculation on SSH
# 
# $$SSH_{sq}(i) = SSH(i)^2$$

# In[9]:


SSH_sq = data.SSH * data.SSH

SSH_sq


# *SSH_sq* is itself a `DataArray`.  If instead we had operated on the `numpy` arrays stored in data.SSH directly we would have been returned a `ndarray`

# Before saving, let's give our new *SSH_sq* variable a better name and descriptive attributes. 

# In[10]:


SSH_sq.name = 'SSH^2'
SSH_sq.attrs['long_name'] = 'Square of Surface Height Anomaly'
SSH_sq.attrs['units'] = 'm^2'
SSH_sq.attrs['grid_layout'] = 'original llc'

# Let's see the result
SSH_sq


# much better!  Now we'll save.

# In[11]:


new_filename = data_dir + 'ssh_sq_DataArray.nc'
print 'saving to ', new_filename

SSH_sq.to_netcdf(path=new_filename)
print 'finished saving'


# 
# ### Calculations in the form of `numpy ndarrays`
# 
# If calculations are made at the `ndarray` level then you'll need a different method to save the results to NetCDF.  If you want to learn how to save `ndarrays` to NetCDF, here is some documentation: http://pyhogs.github.io/intro_netcdf4.html

# In[12]:


SSH_dummy_ndarray = data.SSH.values *  data.SSH.values

type(SSH_dummy_ndarray)


# ## Conclusion
# 
# Now you know how to save `Datasets` and `DataArrays` to disk as NetCDF. This is a very useful operation for saving the results and intermediate results of your work.
