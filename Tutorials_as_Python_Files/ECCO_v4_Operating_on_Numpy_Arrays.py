#!/usr/bin/env python
# coding: utf-8

# # Operating on Numpy arrays
# 
# ## Objectives
# 
# Introduce numpy's **pass-by-reference** approach handling numpy arrays and methods for avoiding pitfalls when operating on numpy arrays.
# 
# ## Introduction
# 
# From: http://span.ece.utah.edu/common-mistakes-in-moving-from-matlab-to-python:
# 
# "Whenever you reference an array in Python, the computer will provide the memory address for the thing you are accessing, not the actual value. This is called **pass-by-reference**. This saves memory and makes your programs faster, but it is also harder to keep straight."  
# 
# From: https://docs.python.org/2/library/copy.html
# 
# "Assignment statements in Python do not copy objects, they create bindings [pointers] between a target and an object." "... a copy is sometimes needed so one can change one copy without changing the other. The 'copy' module provides generic ... copy operations."
# 
# If you are not familiar with the **pass-by-reference** aspect of Python then I strongly suggest you read this short, informative essay on "Python Names and Values": https://nedbatchelder.com/text/names.html
# 
# We've briefly touched on this important subject in earlier tutorials.  Now we'll go into a bit more detail.

# ## Variable assignments
# 
# Unlike some other languages, creating a new variable with an assignment statement in Python such as
# `
# x = some_numpy_array
# `
# 
# does not make a copy of ``some_numpy_array``.  Instead, the assignment statement makes ``x`` and ``some_numpy_array`` both point to the same `numpy` array in memory.  Because ``x`` and ``some_numpy_array`` are both refer (or pointer) to the same `numpy` array in memory, the `numpy` array can be changed by operations on either ``x`` or ``some_numpy_array``.  If you aren't aware of this behavior then you may run into very difficult to identify bugs in your calculations!
# 
# ### A simple demonstration
# 
# Let's demonstrate this issue with a very simple `numpy` array

# In[1]:


import numpy as np
import xarray as xr
import sys
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import json
from copy import deepcopy 
import warnings
warnings.filterwarnings('ignore')


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
# base_dir = '/home/username/'
base_dir = '/Users/ifenty/ECCOv4-release'

## define a high-level directory for ECCO fields
ECCO_dir = base_dir + '/Release3_alt'


# Create a simple numpy array

# In[4]:


a=np.array([1, 2, 3, 4, 5])

# Assign 'b' to point to the same numpy array
b=a

# Test to see if b and a point to the same thing
b is a


# Now change the fourth element of ``b`` and print both ``a`` and ``b``

# In[5]:


b[3] = 10
print (a)
print (b)


# ### A fancier demonstration
# 
# Let's now demonstrate with a `numpy` array that stores ``SSH`` output.

# In[6]:


## LOAD NETCDF SSH FILE

# directory of the file
data_dir= ECCO_dir + '/nctiles_monthly/SSH/'
# filename
fname = 'SSH_2010.nc'
# load the dataset file
ssh_dataset = xr.open_dataset(data_dir + fname).load()

## LOAD NETCDF GRID FILE
grid_dir= ECCO_dir + '/nctiles_grid/'
# filename
fname = 'ECCOv4r3_grid.nc'
# load the dataset file
grid_dataset = xr.open_dataset(grid_dir + fname).load()

## Merge SSH and GRID
output_all = xr.merge((ssh_dataset, grid_dataset))


# Recall the dimensions of our ``SSH`` `DataArray`:

# In[7]:


output_all.SSH.dims


# Show the first four SSH values in **j** and **i** for the fifth month (May 1992) and second tile:

# In[8]:


output_all.SSH[4,1,0:4,0:4].values


# Assign the variable `ssh_tmp` to this *subset* of the `numpy` array that ``SSH`` points to:

# In[9]:


ssh_tmp = output_all.SSH[4,1,0:2,0:2].values
ssh_tmp


# Now change the values of all elements of ``ssh_tmp`` to 10

# In[10]:


ssh_tmp[:] = 10
ssh_tmp


# And see that yes, in fact, this change is reflected in our ``SSH`` `DataArray`:

# In[11]:


output_all.SSH[4,1,0:4,0:4].values


# ## Dealing with *pass-by-reference*: right hand side operations
# 
# One way to have a new variable assignment not point to the original variable is to *perform an operation on the right hand side of the assignment statement*.  
# 
# "Python evaluates expressions from left to right. Notice that while evaluating an assignment, the right-hand side is evaluated before the left-hand side."
# https://docs.python.org/2/reference/expressions.html#evaluation-order
# 
# Performing an operation on the right hand side creates new values in memory.  The new variable assignment will then point to these new values, leaving the original untouched.
# 
# ### Simple demonstration 1
# Operate on ``a`` by adding 1 before the assigment statement

# In[12]:


# Create a simple numpy array
a=np.array([1, 2, 3, 4, 5])

b = a + 1

print (a)
print (b)


# Now change the fourth element of ``b`` and print both ``a`` and ``b``

# In[13]:


b[3] = 10
print (a)
print (b)


# ``a`` and ``b`` do indeed point to different values in memory.

# ### Simple demonstration 2
# 
# Operate on ``a`` by adding 0 before the assigment statement.  This is a kind of dummy operation.

# In[14]:


# Create a simple numpy array
a=np.array([1, 2, 3, 4, 5])

# Add 0 to `a`:
b = a + 0

print (a)
print (b)


# In[15]:


# Test to see if b and a point to the same thing
b is a


# Now change the fourth element of ``b`` and print both ``a`` and ``b``

# In[16]:


b[3] = 10
print (a)
print (b)


# Once again we see that ``a`` and ``b`` do indeed point to different values in memory.

# ### A fancier demonstration
# 
# Let's now demonstrate with a `numpy` array that stores ``SSH`` output.

# In[17]:


output_all.SSH[4,1,5:9,5:9].values


# In[18]:


ssh_tmp = output_all.SSH[4,1,5:9,5:9].values * output_all.rA[1,5:9,5:9].values
ssh_tmp[:] = 10
ssh_tmp


# In[19]:


output_all.SSH[4,1,5:9,5:9].values


# Operating on the right hand side of the assignment does indeed new arrays in memory leaving the original SSH `numpy` array untouched.

# ## Dealing with *pass-by-reference*: copy and deepcopy
# 
# A second way to have a new variable assignment not point to the original variable is to *use the copy or deepcopy command*.
# 
# ### Simple demonstration
# Use the `numpy` command.

# In[20]:


# Create a simple numpy array
a=np.array([1, 2, 3, 4, 5])
b=np.copy(a)

print (a)
print (b)


# Now change the fourth element of ``b`` and print both ``a`` and ``b``

# In[21]:


b[3] = 10
print (a)
print (b)


# In[22]:


output_all.SSH


# ### Fancier demonstration
# 
# `Dataset` and `DataArray` objects are too complicated for `numpy`'s `copy` command.  For complex objects such as these use the `deepcopy` command.

# In[23]:


ssh_tmp = deepcopy(output_all.SSH)
ssh_tmp[:] = 10
ssh_tmp[4,1,5:9,5:9].values


# In[24]:


output_all.SSH[4,1,5:9,5:9].values


# Using `deepcopy` gives us an entirely new array in memory.  Operations on ``ssh_tmp`` do not affect the original fields that we found in the `output_all_SSH` `DataArray`.
# 
# #### alternative to `deepcopy`
# `xarray` give us another way to deepcopy `DataArrays` and `Datasets`:
# 
# ``
# ssh_tmp = output_all.copy(deep=True)
# ``

# ## Conclusion
# 
# You now know about the possible pitfalls for dealing with Python's **pass-by-reference** way of handling assignment statements and different methods for making copies of `numpy` arrays and `Datasets` and `DataArrays`.  
