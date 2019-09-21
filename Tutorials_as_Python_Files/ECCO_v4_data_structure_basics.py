#!/usr/bin/env python
# coding: utf-8

# # The Dataset and DataArray objects used in the ECCOv4 Python package.
# 
# ## Objectives
# 
# To introduce the two high-level data structures, `Dataset` and `DataArray`, that are used in by the `ecco_v4_py` Python package to load and store the ECCO v4 model grid parameters and state estimate variables.
# 
# ## Introduction
# 
# The ECCOv4 files are provided as NetCDF files.  The file you have may look a little different than the ones shown here because we have been working hard at improving how what exactly goes into our NetCDF files.  
# 
# This tutorial document is current as of Sep 2019 with the ECCOv4 NetCDF grid files provided in the following directories:
# 
# https://ecco.jpl.nasa.gov/drive/files/Version4/Release3_alt (po.daac drive, recommended)
# 
# https://web.corral.tacc.utexas.edu/OceanProjects/ECCO/ECCOv4/Release3_alt/  (mirror at U. Texas, Austin)
# 
# In this first tutorial we will start slowly, providing detail at every step.  Later tutorials will be assume knowledge of some basic operations introduced here.
# 
# Let's get started.
# 
# ## Import external packages and modules
# 
# Before using Python libraries we must import them.  Usually this is done at the beginning of every Python program or interactive Juypter notebook instance but one can import a library at any point in the code.  Python libraries, called **packages**, contain subroutines and/or define data structures that provide useful functionality.
# 
# Before we go further, let's import some packages needed for this tutorial:

# In[1]:


# NumPy is the fundamental package for scientific computing with Python. 
# It contains among other things:
#    a powerful N-dimensional array object
#    sophisticated (broadcasting) functions
#    tools for integrating C/C++ and Fortran code
#    useful linear algebra, Fourier transform, and random number capabilities
# http://www.numpy.org/
#
# make all functions from the 'numpy' module available with the prefix 'np'
import numpy as np

# xarray is an open source project and Python package that aims to bring the 
# labeled data power of pandas to the physical sciences, by providing
# N-dimensional variants of the core pandas data structures.
# Our approach adopts the Common Data Model for self- describing scientific 
# data in widespread use in the Earth sciences: xarray.Dataset is an in-memory
# representation of a netCDF file.
# http://xarray.pydata.org/en/stable/
#
# import all function from the 'xarray' module available with the prefix 'xr'
import xarray as xr


# ### Load the ECCO Version 4 Python package
# 
# The *ecco_v4_py* is a Python package written specifically for working with the NetCDF output provided in the [nctiles_monthly](ftp://ecco.jpl.nasa.gov/Version4/Release3/nctiles_monthly/) directory of the [ECCO v4 release](ftp://ecco.jpl.nasa.gov/Version4/Release3/)
# 
# See the "Getting Started" page in the tutorial for instructions about installing the *ecco_v4_py* module on your machine.

# In[2]:


## Import the ecco_v4_py library into Python
## =========================================
## -- If ecco_v4_py is not installed in your local Python library, 
##    tell Python where to find it.  For example, if your ecco_v4_py
##    files are in /home/username/ECCOv4-py/ecco_v4_py, then use:
import sys
sys.path.append('/Users/ifenty/ECCOv4-py')
import ecco_v4_py as ecco


# The syntax 
# 
# ```Python
#   import XYZ package as ABC
# ```
# 
# allows you to access all of the subroutines and/or objects in a package with perhaps a long complicated name with a shorter, easier name.
# 
# Here, we import `ecco_v4_py` as `ecco` because typing `ecco` is easier than `ecco_v4_py` every time.  Also, `ecco_v4_py` is actually comprised of multiple python modules and by importing just `ecco_v4_py` we can actually access all of the subroutines in those modules as well.  Fancy.

# ## Load a single state estimate variable NetCDF tile file
# 
# To load ECCO v4's NetCDF files we will use the *open_dataset* command from the Python package [xarray](http://xarray.pydata.org/en/stable/index.html). The *open_dataset* routine creates a `Dataset` object and loads the contents of the NetCDF file, including its metadata, into a data structure.    
# 
# Let's open one monthly mean THETA file associated with *tile 2* (the North East Atlantic Ocean).

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

# directory containing the file
data_dir= ECCO_dir + '/nctiles_monthly/THETA/'

# filename
fname = 'THETA_2010.nc'

# load the file
ds = xr.open_dataset(data_dir + fname).load()


# What is *ds*?  It is a `Dataset` object which is defined somewhere deep in the `xarray` package:

# In[5]:


type(ds)


# ## The Dataset object 
# 
# According to the xarray documentation, a [Dataset](http://xarray.pydata.org/en/stable/generated/xarray.Dataset.html) is a Python object designed as an "in-memory representation of the data model from the NetCDF file format."
# 
# What does that mean?  NetCDF files are *self-describing* in the sense that they [include information about the data they contain](https://www.unidata.ucar.edu/software/netcdf/docs/faq.html).  When `Datasets` are created by loading a NetCDF file they load all of the same data and metadata.
# 
# Just as a NetCDF file can contain many variables, a `Dataset` can contain many variables.  These variables are referred to as `Data Variables` in the `xarray` nomenclature.
# 
# `Datasets` contain three main classes of fields:
# 
# 1. **Coordinates**   : arrays identifying the coordinates of the data variables
# 2. **Data Variables**: the data variable arrays and their associated coordinates
# 3. **Attributes**    : metadata describing the dataset
# 
# Now that we've loaded `GRID.0003.nc` as the *ds* `Dataset` object let's examine its contents.  
# 
# > **Note:** *You can get information about objects and their contents by typing the name of the variable and hitting **enter** in an interactive session of an IDE such as Spyder or by executing the cell of a Jupyter notebook.*
# 
# 

# In[6]:


ds


# ### Examining the Dataset object contents
# 
# Let's go through *ds* piece by piece, starting from the top.
# 
# #### 1. Object type
# `<xarray.Dataset>`
# 
# The top line tells us what type of object the variable is.  *ds* is an instance of a`Dataset` defined in `xarray`.
# 
# #### 2. Dimensions
# ```Dimensions:    (i: 90, j: 90, k: 50, nv: 2, tile: 13, time: 12)```
# 
# The *Dimensions* list shows all of the different dimensions used by all of the different arrays stored in the NetCDF file (and now loaded in the `Dataset` object).
#   
# Arrays may use any combination of these dimensions.  In the case of this *grid* datasets, we find 1D (e.g., depth), 2D (e.g., lat/lon), and 3D (e.g., mask) arrays.
#   
# The lengths of these dimensions are next to their name: `(i: 90, j: 90, k: 50, nv: 2, tile: 13, time: 12)`.  There are 50 vertical levels in the ECCO v4 model grid so the `k` corresponds to vertical dimension.  `i` and `j` correspond to the horizontal dimensions.  The lat-lon-cap grid has 13 tiles.  This THETA file has 12 monthly-mean records for 2010.  The dimension `nv` is a time dimension that corresponds to the start and end times of the monthly-mean averaging periods.  In other words, for every 1 month, there are 2 (nv = 2) time records, one describing when the month started and the other when the month ended.
# 
# > **Note:** Each tile in the llc90 grid used by ECCO v4 has 90x90 horizontal grid points.  That's where the 90 in llc**90** comes from!  
# 
# #### 3. Coordinates
# 
# Some coordinates have an asterix **"\*"** in front of their names.  They are known as *dimension coordinates* and are always one-dimensional arrays of length $n$ which specify the length of arrays in the dataset in different dimensions.
# 
# ```
# Coordinates:
#   * j          (j) int32 0 1 2 3 4 5 6 7 8 9 ... 80 81 82 83 84 85 86 87 88 89
#   * i          (i) int32 0 1 2 3 4 5 6 7 8 9 ... 80 81 82 83 84 85 86 87 88 89
#   * k          (k) int32 0 1 2 3 4 5 6 7 8 9 ... 40 41 42 43 44 45 46 47 48 49
#   * tile       (tile) int32 0 1 2 3 4 5 6 7 8 9 10 11 12
#   * time       (time) datetime64[ns] 2010-01-16T12:00:00 ... 2010-12-16T12:00:00
# ``` 
#   
# These [coordinates](http://xarray.pydata.org/en/stable/data-structures.html#coordinates) are arrays whose values *label* each grid cell in the arrays.  They are used for *label-based indexing* and *alignment*.
# 
# Let's look at the three primary spatial coordiates, i, j, k. 

# In[7]:


print(ds.i.long_name)
print(ds.j.long_name)
print(ds.k.long_name)


# `i` indexes (or labels) the tracer grid cells in the `x` direction,  `j` indexes the tracer grid cells in the `y` direction, and similarly `k` indexes the tracer grid cells in the `z` direction.

# #### 4. Data Variables
# ```
# Data variables:
#     THETA      (time, tile, k, j, i) float32 ...
# ``` 
# The *Data Variables* are one or more `xarray.DataArray` objects.  `DataArray` objects are labeled, multi-dimensional arrays that may also contain metadata (attributes).  `DataArray` objects are very important to understand because they are container objects which store the  numerical arrays of the state estimate fields.  We'll investigate these objects in more detail after completing our survey of this `Dataset`.
# 
# In this NetCDF file there is one *Data variables*, `THETA`, which is stored as a five dimensional array (**time, tile, k,j,i**) field of average potential temperature. The llc grid has 13 tiles.  Each tile has two horizontal dimensions (i,j) and one vertical dimension (k).
#   
# `THETA` is stored here as a 32 bit floating point precision.
#   
# > **Note:** The meaning of all MITgcm grid parameters can be found [here](https://mitgcm.readthedocs.io/en/latest/algorithm/horiz-grid.html).
# 

# #### 5. Attributes
# ```
# Attributes:
#     product_time_coverage_start:  1992-01-01T12:00:00
#     author:                       Ian Fenty and Ou Wang
#     Insitution:                   JPL
#     product_version:              ECCO Version 4 Release 3 (ECCOv4r3) 1992-2015
#     time_units:                   days since 1992-01-01 00:00:00
#     Conventions:                  CF-1.6
#     Project:                      Estimating the Circulation and Climate of t...
#     cdm_data_type:                Grid
#     geospatial_lon_units:         degrees_east
#     Metadata_Conventions:         CF-1.6, Unidata Dataset Discovery v1.0, GDS...
#     no_data:                      NaNf
#     geospatial_lat_units:         degrees_north
#     product_time_coverage_end:    2015-12-31T12:00:00
#     geospatial_vertical_min:      0
#     nz:                           50
#     geospatial_vertical_units:    meter
#     geospatial_vertical_max:      6134.5
#     date_created:                 Mon May 13 02:40:10 2019
#     geospatial_lat_max:           89.739395
#     geospatial_lat_min:           -89.873055
#     nx:                           90
#     ny:                           90
#     geospatial_lon_max:           179.98691
#     geospatial_lon_min:           -179.98895
#     time_coverage_start:          2010-01-01T00:00:00
#     time_coverage_end:            2011-01-01T00:00:00
# ```
#   
# The `attrs` variable is a Python [dictionary object](https://www.python-course.eu/dictionaries.php) containing metadata or any auxilliary information.
#   
# Metadata is presented as a set of dictionary `key-value` pairs.  Here the `keys` are *description, A, B,  ... missing_value.* while the `values` are the corresponding text and non-text values.  
#   
# To see the metadata `value` associated with the metadata `key` called "Conventions" we can print the value as follows:

# In[8]:


print (ds.attrs['Conventions'])


# "CF-1.6" tells us that ECCO NetCDF output conforms to the [**Climate and Forecast Conventions version 1.6**](http://cfconventions.org/).  How convenient.  

# ### Map of the `Dataset` object
# 
# Now that we've completed our survey, we see that a `Dataset` is a really a kind of *container* comprised of (actually pointing to) many other objects.  
# 
# + dims: A `dict` that maps dimension names (keys) with dimension lengths (values)
# + coords: A `dict` that maps dimension names (keys such as **k, j, i**) with arrays that label each point in the dimension (values) 
# + One or more *Data Variables* that are pointers to `DataArray` objects 
# + attrs A `dict` that maps different attribute names (keys) with the attributes themselves (values).
# 
# ![Dataset-diagram](../figures/Dataset-diagram.png)

# ## The `DataArray` Object
# 
# It is worth looking at the `DataArray` object in more detail because `DataArrays` store the arrays that store the ECCO output.  Please see the [xarray documentation on the DataArray object](http://xarray.pydata.org/en/stable/data-structures.html#dataarray) for more information.
# 
# `DataArrays` are actually very similar to `Datasets`.  They also contain dimensions, coordinates, and attributes.  The two main differences between `Datasets` and `DataArrays` is that `DataArrays` have a **name** (a string) and an array of **values**.  The **values** array is a [numpy n-dimensional array](https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.array.html), an `ndarray`.
# 
# ### Examining the contents of a `DataArray` 
# 
# Let's examine the contents of one of the coordinates `DataArrays` found in *ds*, *XC*.  

# In[9]:


ds.XC


# ### Examining the `DataArray`
# 
# The layout of `DataArrays` is very similar to those of `Datasets`.  Let's examine each part of *ds.XC*, starting from the top.
# 
# #### 1. Object type
# `<xarray.DataArray>`
# 
# This is indeed a `DataArray` object from the `xarray` package.
# 
# > Note: You can also find the type of an object with the `type` command: `print type(ds.XC)`

# In[10]:


print (type(ds.XC))


# #### 2. Object Name
# `XC`
# 
# The top line shows `DataArray` name, `XC`.

# #### 3. Dimensions
# `(tile: 13, j: 90, i: 90)`  
# 
# Unlike *THETA*, *XC* does not have time or depth dimensions which makes sense since the longitude of the grid cell centers do not vary with time or depth.

# #### 4. The `numpy` Array
# ````
# array([[[-111.60647 , -111.303   , -110.94285 , ...,   64.791115,
#            64.80521 ,   64.81917 ],
#         [-104.8196  , -103.928444, -102.87706 , ...,   64.36745 ,
#            64.41012 ,   64.4524  ],
#         [ -98.198784,  -96.788055,  -95.14185 , ...,   63.936497,
#            64.008224,   64.0793  ],
#         ...,
# ````
# 
# In `Dataset` objects there are *Data variables*.  In `DataArray` objects we find `numpy` **arrays**.  Python prints out a subset of the entire array.  
# 
# > **Note**: `DataArrays` store **only one** array while `DataSets` can store **one or more** `DataArrays`.
# 
# We access the `numpy` array by invoking the `.values` command on the `DataArray`.

# In[11]:


ds.XC.values


# The array that is returned is a numpy n-dimensional array:

# In[12]:


type(ds.XC.values)


# Being a numpy array, one can use all of the numerical operations provided by the numpy module on it.
# 
# 
# > ** Note: ** You may find it useful to learn about the operations that can be made on numpy arrays. Here is a quickstart guide: 
# https://docs.scipy.org/doc/numpy-dev/user/quickstart.html
# 
# We'll learn more about how to access the values of this array in a later tutorial.  For now it is sufficient to know how to access the arrays!

# #### 4. Coordinates
# 
# The dimensional coordinates (with the asterixes) are
# ```
# Coordinates:
#   * j        (j) int32 0 1 2 3 4 5 6 7 8 9 10 ... 80 81 82 83 84 85 86 87 88 89
#   * i        (i) int32 0 1 2 3 4 5 6 7 8 9 10 ... 80 81 82 83 84 85 86 87 88 89
#   * tile     (tile) int32 0 1 2 3 4 5 6 7 8 9 10 11 12
# ```
# 
# We find three 1D arrays with coordinate labels for **j**, **i**, and **tile**.

# In[13]:


ds.XC.coords


# two other important coordinates here are `tile` and `time`

# In[14]:


print('tile: ')
print(ds.tile.values)
print('time: ')
print(ds.time.values)


# The file we loaded was `/nctiles_mean/THETA/THETA_2010.nc`, the 2010 monthly-mean potential temperature field.  Here the time coordinates are the center of the averaging periods.

# #### 5. Attributes
# ```
# Attributes:
#     units:          degrees_east
#     long_name:      longitude at center of tracer cell
#     standard_name:  longitude_at_c_location
#     valid_range:    -180., 180.
# ```
# 
# The `XC` variable has a `long_name` (longitude at center of tracer cell) and units (degrees_east) and other information.  This metadata was loaded from the NetCDF file.  The entire attribute dictionary is accessed using `.attrs`.

# In[15]:


ds.XC.attrs


# In[16]:


ds.XC.attrs['units']


# In[17]:


ds.XC.attrs['valid_range']


# ### Map of the `DataArray` Object
# 
# The `DataArray` can be mapped out with the following diagram:
# 
# ![DataArray-diagram](../figures/DataArray-diagram.png)

# ## Summary
# 
# Now you know the basics of the `Dataset` and `DataArray` objects that will store the ECCO v4 model grid parameters and state estimate output variables.  Go back and take a look athe grid $ds$ object that we originally loaded.  It should make a lot more sense now!
