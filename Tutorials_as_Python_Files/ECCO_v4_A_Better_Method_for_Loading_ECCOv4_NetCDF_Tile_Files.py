
# coding: utf-8

# # A better method for loading ECCOv4 NetCDF tile files
# 
# ## Objectives
# 
# Introduce an alternative method for loading ECCO v4 NetCDF tile files that returns `Dataset` and `DataArray` objects with better labelling of variable coordinates with respect to *where* they are situated on the Arakawa-C grid.
# 
# ## Introduction
# 
# As we showed in the first tutorial, we can use the `open_dataset` method from `xarray` to load a NetCDF tile file into Python as a `Dataset` object.  `open_dataset` is very convienent because it automatically parses the NetCDF file and constructs a `Dataset` object using all of the  dimensions, coordinates, variables, and metadata information.  However, by default the names of the coordinates are pretty generic: **i1, i2, i3**, etc. We can do a lot better. 
# 
# In the last tutorial we loaded a single ECCOv4 grid tile file and examined its contents.  Let's load it up again and take another look at its coordinates.  This time we'll name the new `Dataset` object  `grid_3_od` since we are loading the file using `open_dataset`.

# In[1]:


import matplotlib.pylab as plt
import numpy as np
import sys
import xarray as xr
from copy import deepcopy 
import ecco_v4_py as ecco


# In[2]:


# point to your local directory holding the nctiles_grid files
grid_dir='/Volumes/ECCO_BASE/ECCO_v4r3/nctiles_grid/'
fname = 'GRID.0003.nc'
grid_3_od = xr.open_dataset(grid_dir + fname)


# In[3]:


grid_3_od


# ### Examining the Dataset object contents
# 
# We see that all of the Data variables in `grid_3_od` use one of three dimensions, **i1**, **i2**, and **i3**.  As we saw before, some variables are 3D, ``hFacC(i1,i2,i3)``; others are 2D, ``XC(i2,i3)``; and others are 1D, ``RF(i1)``.
#     
# While this `Dataset` object is useful, its coordinate names are somewhat ambiguous.  The MITgcm uses the staggered Arakawa-C grid (hereafter c-grid) to discretize the model equations.  Model variables are not co-located in c-grid models; they can fall on one of four different categories of point.  From the above, one cannot distinguish *where* on the *c-grid* these different variables are situated.  
# 
# To understand this issue better, let's review the c-grid coordinate system.
# 
# ## The four horizontal points of the c-grid
# 
# ![C-grid-points.png](../figures/C-grid-points.png)
# **The four different categories of point used in the staggered Arakawa-C grid (C-grid)**
# 
# ### *c* points
# 
# Scalar variables (e.g., ``THETA, SALT, SSH, OBP, SIarea``) and variables associated with vertical velocity (e.g., ``WVEL``) are situated at the center of the tracer grid cell in the horizontal plane.  These are $c$ points
# 
# Define the coordinates $(i,j)$ for the discrete indices of $c$ points in the $x$ and $y$ directions, respectively.
# 
# In the ECCO v4 NetCDF files, $c(0,0)$ is the $-x$ most and $-y$ most tracer grid cell.
# 
# * In the +$y$ direction, the next $c$ point is $c(0,1)$.
# * In the +$x$ direction, the next $c$ point is $c(1,0)$ 
# 
# ### *u* points
# 
# Vector variables related to horizontal velocity in the $x$ direction are  staggered along the edges of tracer cells between $c$ points in the horizontal plane. Examples include horizontal velocity in the $x$ direction (``UVEL``) and horizontal advective flux of snow in the $x$ direction (``ADVxSNOW``).  They are situated along the edges (if 2D) or faces (if 3D) of the tracer grid cells in the $x$ direction.     
# 
# Define the coordinates $(i_g,j)$ for the discrete indices of $u$ points in the $x$ and $y$ directions, respectively.
# 
# We use $i_g$ as the coordinate in the $x$ direction because $u$ points are situated along the tracer grid cell ed***G***es.  We use $j$ for its $y$ coordinate because $u$ points and $c$ points fall along the same lines in $y$.
# 
# In the ECCO v4 netCDF files, $u(0,0)$ is the $-x$ most and $-y$ most $u$ point.
# 
# ### *v* points
# 
# Vector variables related to horizontal velocity in the $y$ direction are  staggered along the edges of tracer cells between $c$ points in the horizontal plane. Examples include horizontal velocity in the $y$ direction (``VVEL``) and horizontal advective flux of snow in the $y$ direction (``ADVySNOW``).  They are situated along the edges (if 2D) or faces (if 3D) of the tracer grid cells in the $y$ direction.     
# 
# Define the coordinates $(i,j_g)$ for the discrete indices of $v$ points in the $x$ and $y$ directions, respectively.
# 
# We use $j_g$ as the coordinate in the $y$ direction because $v$ points are situated along the tracer grid cell ed**G**es.  We use $i$ for its $x$ coordinate because $v$ points and $c$ points fall along the same lines in $x$.  
# 
# In the ECCO v4 NetCDF files, $v(0,0)$ is the $-x$ most and $-y$ most $v$ point.
# 
# ### *g* points
# 
# Variables that are explictly related to horizontal velocities in the model in both the $x$ and $y$ direction are situated at $g$ points in the horizontal plane.  $g$ points are situated at the corners of tracer grid cells.  
# 
# Define the coordinates $(i_g,j_g)$ for the discrete indices of $g$ points in the $x$ and $y$ directions, respectively.  
# 
# We use $i_g$ and $j_g$ because $g$ points are on the ed**G**es (corners) of tracer grid cells.
# 
# In the ECCO v4 NetCDF files, $g(0,0)$ is the $-x$ most and $-y$ most $g$ point.
# 
# ## The two vertical points of the c-grid
# 
# There are two different coordinates in the vertical $z$ dimension:
# 
# ### *w* points
# 
# Variables related to vertical velocity or vertical fluxes are situated at $w$ points in the vertical direction.  These variables are situated on the upper and lower faces of the tracer grid cell.   
# 
# Define the coordinate $k_g$ for $w$ points using the same reasoning as above: $w$ points fall along the the ed**G**es of tracer grid cells in the $z$ direction (top and bottom of the grid cells).
# 
# In the ECCO v4 NetCDF files, $k_g(0)$ is the sea surface.
# 
# ### *k* points
# 
# Variables that are not related to vertical velocity or vertical fluxes are situated at $k$ points in the vertical direction.  These variables are situated on the upper and lower faces of the tracer grid cell.   
# 
# Define the coordinate $k$ for points situated in the center of the grid cells.
# 
# In the ECCO v4 NetCDF files, $k(0)$ is the middle of the uppermost tracer grid cell.
# 
# 
# ## Applying the C-grid coordinates to the variables
# 
# The default coordinate names in the ECCO v4 NetCDF tile files do not distinguish between the four horizontal coordinates, $i, i_g, j, j_g$ and the two vertical coordinates, $k_g$ and $k$, used by our c-grid model.
# 
# To apply these more descriptive coordinates to the `Dataset` objects that are created when we load netCDF files, we provide a special routine, `load_tile_from_netcdf`. 
# 
# `load_tile_form_netcdf` takes four arguments,
# 
# 1. *data_dir*: the directory of the NetCDF file
# 
# 2. *var*: the name of the variable stored in the NetCDF file (without the tile number).  Examples include *THETA*, *VVEL*, *DFxEHFF*
# 
# 3. *var_type*: one of 'c', 'g', 'u', 'v', or 'grid' corresponding with the variable's c-grid point type.  'grid' is a special case because unlike every other variable NetCDF file, grid files actually contain several variables that are on different c-grid points.
# 
# 4. *tile_index*: the tile number [1 .. 13]

# ### Loading an ECCO v4 netCDF grid tile file using  `load_tile_from_netcdf`
# 
# Let's use `load_tile_from_netcdf` to load grid tile 3 again. This time we'll call the new `Dataset` object `grid_3_new`

# In[4]:


var = 'GRID'
var_type = 'grid'
tile_index = 3
grid_3_new = ecco.load_tile_from_netcdf(grid_dir, var, var_type, tile_index)


# In[5]:


grid_3_new


# ### Examining the Dataset object contents
# 
# 
# #### 1. Dimensions
# `Dimensions:  (i: 90, i_g: 90, j: 90, j_g: 90, k: 50, k_g: 50)`
# 
# The *Dimensions* list now lists the six different coordinates and their dimension used by variables stored in this new grid tile `Dataset` object.  
# 
# #### 2. Coordinates
# ```
# Coordinates:
#     tile     int64 3
#   * k        (k) float64 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 ...
#   * i        (i) float64 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 ...
#   * j        (j) float64 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 ...
#   * i_g      (i_g) float64 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 ...
#   * j_g      (j_g) float64 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 ...
#   * k_g      (k_g) float64 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 ...
# ``` 
# 
# Each of the 6 dimensions has a corresponding set of coordinate labels.
# 
# 
# #### 3. Data Variables
# ```
# Data variables:
#     XC       (j, i) float64 ...
#     hFacC    (k, j, i) float64 ...
#     XG       (j_g, i_g) float64 ...
#     DXC      (j, i_g) float64 ...
#     hFacW    (k, j, i_g) float64 ...
#     DYC      (j_g, i) float64 ...
#     hFacS    (k, j_g, i) float64 ...
#     RF       (k_g) float64 ...
#     RC       (k) float64 ...
# ```
# 
# We now see that each *data variable* is now described with its proper coordinates.
# For example, ``XC``, the longitude of the tracer grid cell center uses $i,j$ coordiantes while ``XG``, the longitude grid cell corners, now has $i_g, j_g$ coordinates.  ``DXC``, the distance between adjacent cell centers, uses $i_g, j$ coordiantes: consistent with the notion that distances between cell centers in $x$ are situated between cell centers.
# 
# > **Note:** The ordering of arrays is not (x,y,z) but (z,y,x).  
# 
# In addition to the variables now having more descriptive coordinates, the `load_tile_from_netcdf` routine also adds 3 new three-dimensional land masks, one each for grid cell at the 'u' (land_u), 'c' (land_c), and 'v' (land_v) points.

# ## A closer look at the DataArray after loading with `load_tile_from_netcdf`
# 
# Let's take a quick look at one of the variables loaded using the `load_tile_from_netcdf` routine:

# In[6]:


grid_3_new.XC


# The `load_tile_from_netcdf` routine adds a new coordinate to the variables, *tile*.   This will be helpful when loading and combining multiple tiles.
# 
# Finally, we see a new attribute: *rotated_to_latlon*.  This attribute is a flag telling us that the arrays are in the original lat-lon-cap tile layout.  We'll discuss this later.

# ## Loading non-grid ECCO variables using `load_tile_from_netcdf`
# 
# So far we've only looked at ECCO v4 grid tile files.  With the `load_tile_from_netcdf` routine you can assign the proper coordinates to any variable by specifying its c-grid point category.  
# 
# We'll demonstrate by loading three variables, one each that are on 'c','u', and 'v' points 
# 
# ### A 'c' point variable: ``SSH``
# 
# Let's load tile 3 of the 'c' point sea surface height (``SSH``) variable.  

# In[7]:


data_dir='/Volumes/ECCO_BASE/ECCO_v4r3/nctiles_monthly/SSH/'    
var = 'SSH'
var_type = 'c'
tile_index = 3
ssh_tile_3 = ecco.load_tile_from_netcdf(data_dir, var, var_type, tile_index)


# Since this is the first time we're loading an output variable from the state estimate let's closely examine `ssh_tile_3`.

# In[8]:


ssh_tile_3


# Notice that unlike the grid files, the `ssh_tile_3` `Dataset` object only has one *Data variable*, ``SSH``.  Non-grid ECCO v4 NetCDF tile files only ever store *one* physical variable. Consequently, when loading a non-grid file the *dimensions* of the `Dataset` will be the same as the *dimensions* of the *Data variable*.  Let's take a look at the ``SSH`` `DataArary`:

# In[9]:


ssh_tile_3.SSH


# #### Dimensional coordinates
# As expected, the ``SSH`` `DataArray` uses **i, j** coordinates for its horizontal dimensions.  We also see a **tile** and **time** coordinates.  The ordering of the ``SSH`` dimensions is **time, tile, j, i** with the logic being that **tile** is a kind of space coordinate.
# 
# We find 288 records in the **time** dimensions.  One for each month of the 1992-2015 state estimate.  There is one record in the **tile** dimensions because we have only loaded a single tile file so far.
# 
# #### Non-dimensional coordinates
# 
# Notice the four *Coordinates* that do not have an "\*" in front of their names.  These are so-called [non-dimensional coordinates](http://xarray.pydata.org/en/stable/data-structures.html#coordinates).  Think of non-dimensional coordinates as helpful extra coordinate labels.  Here, the non-dimensional coordinates include two for the **time** dimension: **tim** and **timestep**, and two for the space dimensions: **lon_c, lat_c**.  
# 
# ##### The time non-dimensional coordinates

# In[10]:


ssh_tile_3.time


# The **tim** non-dimensional coordinate provides the calendar date (centered in the month) for the **time** monthly index and the **timestep** non-dimensional coordinate provides the the timestep number corresponding to each month (the ECCO v4 simulation uses 1-hourly timesteps).
# 
# ##### The space non-dimensional coordinates
# **lon_c** and **lat_c** non-dimensional coordinates provide the longitude and latitude for each of the 'c' points of ``SSH``.

# In[11]:


ssh_tile_3.lon_c


# ### A 'u' point variable: ``UVEL``
# 
# Let's load tile 3 of the 'u' point horizontal velocity in the local $x$ direction variable, ``UVEL``.

# In[12]:


data_dir='/Volumes/ECCO_BASE/ECCO_v4r3/nctiles_monthly/UVEL/'    
var = 'UVEL'
var_type = 'u'
tile_index = 3
uvel_tile_3 = ecco.load_tile_from_netcdf(data_dir, var, var_type, tile_index)


# Let's look at `uvel_tile_3`.  This time let's also remove some of the descriptive NetCDF file *Attributes* using a little routine called `minimal_metadata`.  We've already seen these attributes a number of times.

# In[13]:


ecco.minimal_metadata(uvel_tile_3)


# In[14]:


uvel_tile_3


# `uvel_tile_3` has one *Data variable*, ``UVEL``.  Let's take a look at the ``UVEL`` `DataArary`:

# In[15]:


uvel_tile_3.UVEL


# #### Dimensional coordinates
# As expected, ``UVEL`` uses the **i_g, j** coordinates for its horizontal dimensions. Unlike ``SSH``, ``UVEL`` has three-dimensions in space so we find a **k** coordinate.  The ordering of the three-dimensional ECCO v4 output is **time, tile, k, j, i**.
# 
# #### Non-dimensional coordinates
# 
# ``UVEL`` has one new non-dimensional coordinate for **k**: **dep** and two new non-dimensional coordinates for space: **lon_u, lat_u**
# 
# ##### The ``dep`` non-dimensional coordinate
# 
# ``dep`` is the depth of the center of the tracer grid cell in meters:

# In[16]:


uvel_tile_3.dep


# ##### The space non-dimensional coordinates
# **lon_u** and **lat_u** provide the longitude and latitude for each **i_g,j** point.
# 
# ### A 'v' point variable: ``VVEL``
# 
# Let's load tile 3 of the 'v' point horizontal velocity in the local $x$ direction variable, ``VVEL``.

# In[17]:


data_dir='/Volumes/ECCO_BASE/ECCO_v4r3/nctiles_monthly/VVEL/'    
var = 'VVEL'
var_type = 'v'
tile_index = 3
vvel_tile_3 = ecco.load_tile_from_netcdf(data_dir, var, var_type, tile_index)
ecco.minimal_metadata(vvel_tile_3)


# Let's look at the *Coordinates* of ``VVEL``.

# In[18]:


vvel_tile_3.coords


# #### Dimensional coordinates
# As expected, `vvel_tile_3` uses the **i, j_g** coordinates for its horizontal dimensions. 
# 
# #### Non-dimensional coordinates
# 
# ``VVEL`` has two new non-dimensional coordinates for space: **lon_v, lat_v**.  As you might expect, these are the longitude and latitude of this v-point variable.
# 
# ## Conclusion
# 
# Now you know how to load model grid parameters and output variable tiles with coordinates that correspond with their location on the c-grid.
