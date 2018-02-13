
# coding: utf-8

# # Accessing and Subsetting Variables
# 
# ## Objectives
# 
# Introduce several alternative methods for accessing and subsetting the variables stored in the `Datasets` and `DataArrays` that store the ECCO v4 state estimate variables and grid parameters.
# 
# ## Accessing variables inside `Dataset` and `DataArray` objects
# 
# There are two methods for accessing variables stored in `DataArray` and `Dataset` objects, the "dot" method and the "dictionary" method.  The syntax of these methods is as follows:
# 
# 1. The "dot" method: e.g. ,`X.Y`
# 2. The "dictionary" method: e.g., `Y['Y']`
# 
# Both methods work identically to access *Dimensions*, *Coordinates*, and *Data variables*. Accessing *Attribute* variables requires a slightly different approach as we will see.
# 
# ### Accessing *Data variables*, *Coordinates*, and *Data variables*
# 
# Let's demonstrate both methods by accessing the a *Data variable* from a `Dataset`.  First create the `Dataset`.

# In[1]:


import matplotlib.pylab as plt
import numpy as np
import sys
import xarray as xr
from copy import deepcopy 
import ecco_v4_py as ecco
import warnings
warnings.filterwarnings('ignore')
get_ipython().magic(u'pylab inline')
pylab.rcParams['figure.figsize'] = (10, 6)


# In[2]:


data_dir='/Volumes/ECCO_BASE/ECCO_v4r3/nctiles_monthly/SSH/'
var = 'SSH'
var_type = 'c'

ssh_all_tiles = ecco.load_all_tiles_from_netcdf(data_dir, var, var_type)
ecco.minimal_metadata(ssh_all_tiles)


# In[3]:


ssh_all_tiles


# Now we'll use the two methods to access the ``SSH`` `DataArray`,

# In[4]:


ssh_A = ssh_all_tiles.SSH
ssh_B = ssh_all_tiles['SSH']


# In[5]:


print type(ssh_A)
print type(ssh_B)


# We access the `numpy` arrays stored in these `DataArrays` by invoking their `.values`.

# In[6]:


ssh_arr = ssh_A.values
print type(ssh_arr)


# The shape of the `numpy` array can be found by invoking its `.shape`

# In[7]:


ssh_arr.shape


# The order of these four dimensions is consistent with their ordering in the original `DataArray`,

# In[8]:


ssh_A.dims


# ``ssh_A`` and ``ssh_B`` are new variables but they are not **copies** of the original  ``SSH`` `DataArray` object, they both point to the original `numpy` array.  
# 
# We can confirm that ``ssh_A`` and ``ssh_B`` both refer to the same array in memory using the Python `is` command: 

# In[9]:


ssh_A.values is ssh_B.values


# ### Accessing *Attribute* fields
# 
# Accessing *Attribute* variables using the dictionary method requires using the ``attrs`` variable field:

# In[10]:


print ssh_all_tiles.date


# In[11]:


print ssh_all_tiles.attrs['date']


# ## Subsetting variables using the [], .sel, and .isel syntaxes
# 
# So far, a considerable amount of attention has been placed on the *Coordinates* of `Dataset` and `DataArray` objects.  Why?  Labeled coordinates are certainly not necessary for calculations on the basic numerical arrays that store the ECCO state estimate fields.  The reason so much attention has been placed on coordinates is because the `xarray` offers several very useful methods for selecting (or indexing) subsets of data.
# 
# We'll introduce these indexing techniques  with `DataArray` objects first.  To begin let's make a new variable for the ``SSH`` `DataArray`.

# In[12]:


ssh_da = ssh_all_tiles.SSH
type(ssh_da)


# ### Subsetting `numpy` arrays using the **[ ]** syntax
# 
# Subsetting `numpy` arrays is simple with the standard Python **[ ]** syntax.  To demonstrate, subset the first month and second tile of SSH.  
# 
#   > **Note:** *Note `numpy` array indexing starts with 0*

# In[13]:


ssh_jan92_tile2 = ssh_arr[0,1,:,:]
ssh_jan92_tile2


# The resulting subset is itself a `numpy` array,

# In[14]:


type(ssh_jan92_tile2)


# of shape 90x90, as expected.

# In[15]:


ssh_jan92_tile2.shape


# We can always use **[ ]** method to subset our `numpy` arrays.  It is a simple, direct method for accessing our fields.  Just for fun let's plot this SSH subset:

# In[16]:


## origin='lower' is required so that y increases from bottom to top.
plt.imshow(ssh_jan92_tile2, origin='lower')
plt.colorbar()
plt.title('Jan 1992 SSH [m]')
plt.show()


# ### Subsetting `DataArrays` using the **[ ]** syntax
# 
# An interesting and useful alternative to subsetting `numpy` arrays with the **[ ]** method is to subset `DataArray` instead:

# In[17]:


ssh_jan92_tile2_da = ssh_da[0,1,:,:]
ssh_jan92_tile2_da


# The resulting `DataArray` is a subset of the original `DataArray`.  The subset has two fewer dimensions (**tile** and **time** have been eliminated). The horizontal dimensions **j** and **i** are unchanged.
# 
# Even though the **tile** and **time** dimensions have been eliminated, the dimensional and non-dimensional coordinates associated with **time** and **tile** remain.  In fact, these coordinates *tell us when in time and which tile our subset comes from*:
# 
# ``
# Coordinates:
#     time      float64 1.0
#     tim       datetime64[ns] 1992-01-16
#     timestep  float64 732.0
#     tile      int64 2
# ``
# 
#   > **Note:** *The **tile** coordinate is $2$ because coordinates are labels.

# ### Subsetting `DataArrays` using the **.sel( )** syntax
# 
# Another useful method for subsetting `DataArrays` is the **.sel( )** syntax.  The **.sel( )** syntax takes advantage of the fact that coordinates are labels.  We **sel**ect subsets of the `DataArray` by providing a subset of coordinate labels.
# 
# Let's select tile 2 and time (month) 1:

# In[18]:


ssh_jan92_tile2_sel = ssh_da.sel(tile=2, time=1)
ssh_jan92_tile2_sel


# ### Subsetting `DataArrays` using the **.isel( )** syntax
# 
# The last subsetting method is **.isel( )** syntax.  **.isel( )** uses the numerical **index** of coordinates instead of their label.  Subsets are extraccted by providing a set of coordinate indices.
# 
# The equivalent syntax for subsetting tile 2 and month 1 with the **.isel( )** syntax is,

# In[19]:


ssh_jan92_tile2_isel = ssh_da.isel(tile=1, time=0)
ssh_jan92_tile2_isel


# ### More examples of subsetting using the **[ ]**, **.sel( )** and **.isel( )** syntaxes
# 
# In the examples above we only subsetted a single month (Jan 1992) and a single tile (tile 2).  More complex subsetting is possible.  Here are some three examples that yield equivalent, more complex, subsets:
# 
#   > **Note:** Python array indexing goes up to but not including the final number in a range.  Because array indexing starts from 0, array index 41 corresponds to the 42nd element.

# In[20]:


ssh_sub_bracket  = ssh_da[range(0,4), 1, 31:41, 5:22]
ssh_sub_isel     = ssh_da.isel(tile=1, time=[0,1,2,3], i=range(5,22), j=range(31,41))
ssh_sub_sel      = ssh_da.sel(tile=2, time=[1,2,3,4], i=range(6,23), j=range(32,42))

print ssh_sub_bracket.shape
print ssh_sub_isel.shape
print ssh_sub_sel.shape


# In[21]:


print ssh_sub_bracket.coords


# In[22]:


print ssh_sub_isel.coords


# In[23]:


print ssh_sub_sel.coords


# ### Subsetting `Datasets` using the **.sel( )**, and **.isel( )** syntaxes
# 
# Amazingly, we can use the **.sel** and **.isel** methods to simultaneously subset multiple `DataArrays` stored within an single `Dataset`.  Let's make an interesting `Dataset` to subset and then test out the **.sel( )** and **.isel( )** subsetting methods.

# In[24]:


# specify the location of your nctiles_monthly directory
data_dir='/Volumes/ECCO_BASE/ECCO_v4r3/nctiles_monthly/SSH/'    
var = 'SSH'
var_type = 'c'
ssh_all_tiles = ecco.load_all_tiles_from_netcdf(data_dir, var, var_type)
ecco.minimal_metadata(ssh_all_tiles)

# specify the location of your nctiles_grid directory
grid_dir='/Volumes/ECCO_BASE/ECCO_v4r3/nctiles_grid/'
var = 'GRID'
var_type = 'grid'
grid_all_tiles = ecco.load_all_tiles_from_netcdf(grid_dir, var, var_type)
ecco.minimal_metadata(grid_all_tiles)

# Merge these datasets
output_all = xr.merge([ssh_all_tiles, grid_all_tiles])


# Subset tile 2, j = 50 (a single row through the array), and time = 10 (October 1992)

# In[25]:


output_tile2_time10_j50= output_all.sel(tile=2, time=10, j=50)
output_tile2_time10_j50


# All variables that had **tile, time**, or **j** coordinates have been subset while other variables are unchanged.  Let's plot the seafloor depth and sea surface height from west to east along j=50, (see plot at Line 16) which extends across the S. Atlantic, across Africa to Madagascar, and finally into to W. Indian Ocean.

# In[26]:


f, axarr = plt.subplots(2, sharex=True)
(ax1, ax2) = axarr
ax1.plot(output_tile2_time10_j50.lon_c, output_tile2_time10_j50.SSH,color='b')
ax1.set_ylabel('m')
ax1.set_title('SSH (m)')
ax2.plot(output_tile2_time10_j50.lon_c, -output_tile2_time10_j50.Depth,color='k')
ax2.set_xlabel('longitude')
ax2.set_ylabel('m')
ax2.set_title('Seafloor Depth (m)')
plt.show()


# ### Subsetting using **where( )**
# 
# The **where( )** method is quite different than other subsetting methods because  subsetting is done by masking out values with *nans* that do not meet some specified criteria.  
# 
# For more infomation about **where( )** see http://xarray.pydata.org/en/stable/indexing.html#masking-with-where
# 
# Let's demonstrate **where** by masking out all SSH values that do not fall within a box defined between 20S to 60N and 50W to 10E.
# 
# First, we'll extract the ``SSH`` `DataArray`

# In[27]:


ssh_da=output_all.SSH


# Create a matrix that is `True` where latitude is between 20S and 60N and `False` otherwise.

# In[28]:


lat_bounds = np.logical_and(ssh_da.lat_c  > -20, ssh_da.lat_c < 60)


# Create a matrix that is `True` where longitude is between 50W and 10E and `False` otherwise.

# In[29]:


lon_bounds = np.logical_and(ssh_da.lon_c  > -50, ssh_da.lon_c < 10)


# Combine the ``lat_bounds`` and ``lon_bounds`` logical matrices:

# In[30]:


lat_lon_bounds = np.logical_and(lat_bounds, lon_bounds)


# Finally, use **where** to mask out all SSH values that do not fall within our  ``lat_lon_bounds``

# In[31]:


ssh_da_subset_space = ssh_da.where(lat_lon_bounds)


# To visualize the SSH in our box we'll use one of our ECCO v4 custom plotting routines (which will be the subject of another tutorial).  
# 
# Notice the use of **.sel( )** to subset a single time slice (time=1) for plotting.

# In[32]:


ecco.plot_tiles_proj(ssh_da.lon_c, ssh_da.lat_c, 
                     ssh_da_subset_space.sel(time=1), 
                     plot_type='contourf', lon_0=180, 
                     cmin=-.5, cmax=.5, cbar=True);
plt.title('SSH [m]')
plt.show()


# ## Conclusion
# 
# You now know several different methods for accessing and subsetting fields in `Dataset` and `DataArray` objects.    
# 
# To learn a more about indexing/subsetting methods please refer to the `xarray` manual for indexing methods, http://xarray.pydata.org/en/stable/indexing.html.  
