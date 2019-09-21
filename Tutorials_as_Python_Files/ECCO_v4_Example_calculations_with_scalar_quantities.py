#!/usr/bin/env python
# coding: utf-8

# # Example calculations with scalar quantities
# 
# ## Objectives
# 
# To demonstrate basic calculations using scalar fields (e.g., SSH, T, S) from the state estimate including: time series of mean quantities, spatial patterns of mean quantities, spatial patterns of linear trends, and spatial patterns of linear trends over different time periods.
# 
# ## Introduction
# 
# We will demonstrate global calculations with SSH (global mean sea level time series, mean dynamic topography, global mean sea level trend) and a regional calculation with THETA (nino 3.4 index).
# 
# ## Global calculations with SSH
# 
# First, load all 13 tiles for sea surface height and the model grid parameters and merge the two `Datasets`.

# In[1]:


import numpy as np
import sys
import xarray as xr
from copy import deepcopy 
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
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


# In[4]:


## Load the model grid
grid_dir= ECCO_dir + '/nctiles_grid/'

ecco_grid = ecco.load_ecco_grid_nc(grid_dir, k_subset=[0])

## Load one year of 2D daily data, SSH, SST, and SSS 
data_dir= ECCO_dir + '/nctiles_monthly'

ecco_vars = ecco.recursive_load_ecco_var_from_years_nc(data_dir,                                            vars_to_load=['SSH','THETA'],                                            years_to_load = range(2008,2014),                                            k_subset=[0])
## Merge the ecco_grid with the ecco_vars to make the ecco_ds
ecco_ds = xr.merge((ecco_grid , ecco_vars))

# load the ecco_ds into memory (xarray's Dask is not helpful when we only need to load a few fields)
ecco_ds.load()


# In[5]:


print(ecco_ds.time[0].values)
print(ecco_ds.time[-1].values)


# In[6]:


SST = ecco_ds.THETA.isel(k=0)
SST_0=ecco_ds.THETA.squeeze().isel(time=0)


# ## Sea surface height
# 
# ### Global mean sea level
# 
# Global mean sea surface height at time $t$ is defined as follows:
# 
# $$SSH_{\text{global mean}}(t) = \frac{\sum_{i} SSH(i,t) \,\, A(i)}{A_{\text{global ocean}}}$$
# 
# $$A_{\text{global ocean}} = \sum_{i} A(i)$$
# 
# Where $SSH(i,t)$ is dynamic height at model grid cell $i$ and time $t$, $A(i)$ is the area (m^2) of model grid cell $i$
# 
# There are several ways of doing the above calculations.  Since this is the first tutorial with actual calcuations, we'll present a few different approaches for getting to the same answer.
# 
# #### Part 1: $A_{\text{global ocean}}$
# 
# Let's start on the simplest quantity, the global ocean surface area $A_{\text{global ocean}}$.  Our calculation uses *SSH* which is a 'c' point variable.  The surface area of tracer grid cells is provided by the model grid parameter *rA*.  *rA* is a two-dimensional field that is defined over all model grid points, including land.  
# 
# To calculate the total ocean surface area we need to ignore the area contributions from land. 
# 
# We will first construct a 3D mask that is True for model grid cells that are wet and False for model grid cells that are dry cells.  

# In[7]:


# ocean_mask is ceiling of hFacC which is 0 for 100% dry cells and
# 0 > hFacC >= 1 for grid cells that are at least partially wet

# hFacC is the fraction of the thickness (h) of the grid cell which
# is wet.  we'll consider all hFacC > 0 as being a wet grid cell
# and so we use the 'ceiling' function to make all values > 0 equal to 1.

ocean_mask = np.ceil(ecco_ds.hFacC)
ocean_mask = ocean_mask.where(ocean_mask==1, np.nan)


# In[8]:


# the resulting ocean_mask variable is a 2D DataArray because we only loaded 1 vertical level of the model grid
print(type(ocean_mask))
print((ocean_mask.dims))


# In[9]:


plt.figure(figsize=(12,5), dpi= 90)

ecco.plot_tiles(ocean_mask.isel(k=0),layout='latlon', rotate_to_latlon=True)

# select out the model depth at k=1, round the number and convert to string.
z = str((np.round(ecco_ds.Z.values[0])))
plt.suptitle('Wet (1) /dry (0) mask for k=' + str(0) + ',   z=' + z + 'm');


# To calculate $A_{\text{global ocean}}$ we must apply the surface wet/dry mask to $rA$.  

# In[10]:


# Method 1: the array index method, []
#           select land_c at k index 0
total_ocean_area = np.sum(ecco_ds.rA*ocean_mask[0,:])

# these three methods give the same numerical result.  Here are 
# three alternative ways of printing the result
print ('total ocean surface area ( m^2) %d  ' % total_ocean_area.values)
print ('total ocean surface area (km^2) %d ' % (total_ocean_area.values/1.0e6))

# or in scientific notation with 2 decimal points
print ('total ocean surface area (km^2) %.2E' % (total_ocean_area.values/1.0e6))


# This compares favorable with *Global surface area of Earth's Oceans : approx 3.60 x $10^8$ $km^2$* from https://hypertextbook.com/facts/1997/EricCheng.shtml
# 
# ##### Multiplication of DataArrays
# You probably noticed that the multiplication of grid cell area with the land mask was done element by element.  One useful feature of `DataArrays` is that their dimensions are automatically lined up when doing binary operations.  Also, because *rA* and *ocean_mask* are both `DataArrays`, their inner product and their sums are also `DataArrays`.  
# 
# > Note:: *ocean_mask* has a depth (**k**) dimension while *rA* does not (horizontal model grid cell area does not change as a function of depth in ECCOv4).  As a result, when *rA* is multiplied with *ocean_mask* `xarray` **broadcasts** *rA* to all **k** levels.  The resulting matrix  inherits the **k** dimension from *ocean_mask*.
# 
# ##### Another way of summing over `numpy` arrays
# 
# As *rA* and ocean both store `numpy` arrays, you can also calculate the sum of their product by invoking the `.sum()` command inherited in all `numpy arrays`:

# In[11]:


total_ocean_area = (ecco_ds.rA*ocean_mask).isel(k=0).sum()
print ('total ocean surface area (km^2) ' + '%.2E' % (total_ocean_area.values/1e6))


# #### Part2 : $SSH_{\text{global mean}}(t)$
# 
# The global mean *SSH* at each $t$ is given by,
# 
# $$SSH_{\text{global mean}}(t) = \frac{\sum_{i} SSH(i,t) \,\, A(i)}{A_{\text{global ocean}}}$$
# 
# One way of calculating this is to take advantage of `DataArray` coordinate labels and use its *.sum()* functionality to explicitly specify which dimensions to sum over:

# In[12]:


# note no need to multiple RAC by land_c because SSH is nan over land
SSH_global_mean = (ecco_ds.SSH*ecco_ds.rA).sum(dim=['i','j','tile'])/total_ocean_area


# Alternatively we can do the summation over the three non-time dimensions.  The time dimension of SSH is along the first dimension (axis) of the array, axis 0.

# In[13]:


# note no need to multiple RAC by land_c because SSH is nan over land
SSH_global_mean = np.sum(ecco_ds.SSH*ecco_ds.rA,axis=(1,2,3))/total_ocean_area
SSH_global_mean = SSH_global_mean.compute()


# Even though *SSH* has 3 dimensions (time, tile, j, i) and *rA* and *ocean_mask.isel(k=0)* have 2 (j,i), we can mulitply them. With `xarray` the element-by-element multiplication occurs over their common dimension.
# 
# The resulting $SSH_{global-mean}$ `DataArray` has a single dimension, time.
# 
# #### Part 3 : Plotting the global mean sea level time series:
# 
# Before we plot the global mean sea level curve let's remove its time-mean to make it global mean sea level anomaly (the absolute value has no meaning here anyway).

# In[14]:


plt.figure(figsize=(8,4), dpi= 90)

# Method 1: .mean() method of `DataArrays`
SSH_global_mean_anomaly = SSH_global_mean - SSH_global_mean.mean()

# Method 2: numpy's `mean` method 
SSH_global_mean_anomaly = SSH_global_mean - np.mean(SSH_global_mean)

SSH_global_mean_anomaly.plot()
plt.grid()
plt.title('ECCO v4r3 Global Mean Sea Level Anomaly, 2008-2014');
plt.ylabel('m');
plt.xlabel('year');


# In[15]:


np.max(SSH_global_mean.values)


# ### Mean Dynamic Topography
# 
# Mean dynamic topography is calculated as follows,
# 
# $MDT(i) = \frac{\sum_{t} SSH(i,t) - SSH_{\text{global mean}}(t)}{nt} $
# 
# Where $nt$ is the number of time records. 
# 
# For *MDT* we presere the spatial dimensions. Summation and averaging are over the time dimensions (axis 0).

# In[16]:


## Two equivalent methods

# Method 1, specify the axis over which to average
MDT = np.mean(ecco_ds.SSH - SSH_global_mean,axis=0)

# Method 2, specify the coordinate label over which to average
MDT_B = (ecco_ds.SSH - SSH_global_mean).mean(dim=['time'])

# which can be verified using the '.equals()' method to compare Datasets and DataArrays
print(MDT.equals(MDT_B))


# As expected, MDT has preserved its spatial dimensions:

# In[17]:


MDT.dims


# Before plotting the MDT field remove its spatial mean since its spatial mean conveys no dynamically useful information.  

# In[18]:


MDT_no_spatial_mean = MDT - MDT*ecco_ds.rA/total_ocean_area


# In[19]:


MDT_no_spatial_mean.shape


# In[34]:


plt.figure(figsize=(12,5), dpi= 90)

# mask land points to Nan
MDT_no_spatial_mean = MDT_no_spatial_mean.where(ocean_mask[0,:] !=0)

ecco.plot_proj_to_latlon_grid(ecco_ds.XC,                               ecco_ds.YC,                               MDT_no_spatial_mean*ocean_mask,                               user_lon_0=0,                              plot_type = 'pcolormesh',                               show_colorbar=True,                              dx=2,dy=2);

plt.title('ECCO v4r3 Mean Dynamic Topography [m] 2008-2013');


# ### Spatial variations of sea level linear trends  
# 
# To calculate the linear trend for the each model point we will use on the `polyfit` function of `numpy`.  First, define a time variable in years for SSH.

# In[21]:


ecco_ds.SSH.values.shape


# In[22]:


ssh_flat = np.reshape(ecco_ds.SSH.values,[72, 13*90*90])
ssh_flat.shape


# In[23]:


(ecco_ds.time[0]- ecco_ds.time[1])/86400e9


# In[24]:


days_since_first_record = ((ecco_ds.time - ecco_ds.time[0])/(86400e9)).astype(int)
days_since_first_record


# Next, reshape the four dimensional SSH field into two dimensions, time and space (t, i)

# Now set all $SSH$ values that are 'nan' to zero because the polynominal fitting
# routine can't handle nans,

# In[25]:


ssh_flat[np.isnan(ssh_flat)]=0.0
ssh_flat.shape


# Do the polynomial fitting, https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.polyfit.html

# In[26]:


# slope is in m / day
ssh_slope, ssh_intercept = np.polyfit(days_since_first_record, ssh_flat, 1)

print(ssh_slope.shape)

# and reshape the slope result back to 13x90x90    
ssh_slope = np.reshape(ssh_slope, (13, 90,90))

# mask 
ssh_slope_masked = np.where(ocean_mask[0,:] > 0, ssh_slope, np.nan)

# convert from m / day to mm/year
ssh_slope_mm_year = ssh_slope_masked*365*1e3


# In[27]:


plt.figure(figsize=(12,5), dpi= 90)


ecco.plot_proj_to_latlon_grid(ecco_ds.XC,                               ecco_ds.YC,                               ssh_slope_mm_year,                               user_lon_0=-66,                              plot_type = 'pcolormesh',                               show_colorbar=True,                              dx=2,dy=2, cmin=-15, cmax=15)

plt.title('ECCO v4r3 Global Mean Sea Level Trend 2008 - 2013 [mm/yr]');


# In[28]:


((ssh_slope_mm_year*ecco_ds.rA)/(ecco_ds.rA*ocean_mask).sum()).sum()


# ## Regional calculations with THETA

# In[29]:


lat_bounds = np.logical_and(ecco_ds.YC >= -5, ecco_ds.YC <= 5)
lon_bounds = np.logical_and(ecco_ds.XC >= -170, ecco_ds.XC <= -120)

SST = ecco_ds.THETA.isel(k=0)
SST_masked=SST.where(np.logical_and(lat_bounds, lon_bounds))


# In[30]:


plt.figure(figsize=(12,5), dpi= 90)

ecco.plot_proj_to_latlon_grid(ecco_ds.XC,                               ecco_ds.YC,                               SST_masked.isel(time=0),                              user_lon_0 = -66,                              show_colorbar=True);

plt.title('SST in nino 3.4 box: \n %s ' % str(ecco_ds.time[0].values));


# In[31]:


# Create the same mask for the grid cell area
rA_masked=ecco_ds.rA.where(np.logical_and(lat_bounds, lon_bounds));

# Calculate the area-weighted mean in the box
SST_masked_mean=(SST_masked*rA_masked).sum(dim=['tile','j','i'])/np.sum(rA_masked)

# Substract the temporal mean from the area-weighted mean to get a time series
SST_nino_34_anom = SST_masked_mean - np.mean(SST_masked_mean)
SST_nino_34_anom


# ### Load up the Nino 3.4 index values from ESRL

# In[32]:


# https://www.esrl.noaa.gov/psd/gcos_wgsp/Timeseries/Nino34/
# https://www.esrl.noaa.gov/psd/gcos_wgsp/Timeseries/Data/nino34.long.anom.data
# NINA34
# 5N-5S 170W-120W 
# HadISST 
#  Anomaly from 1981-2010
#  units=degC

nino34_2008_2013 =          [24.79 ,  25.07 ,  26.09 ,  26.88 , 27.22 ,  27.24 , 27.19 , 26.83 ,  26.47 ,           26.43 ,  26.29 ,  25.69 ,  25.58 , 26.05 ,  26.54 , 27.52 , 28.04 ,  28.17 ,  27.91 ,           27.49 ,  27.43 ,  27.69 ,  28.15 , 28.40 ,  28.00 , 27.94 , 28.33 ,  28.33 ,  27.72 ,           27.07 ,  26.34 ,  25.55 ,  25.19 , 25.08 ,  25.08 , 24.95 , 24.88 ,  25.50 ,  26.27 ,           27.03 ,  27.34 ,  27.43 ,  27.00 , 26.22 ,  25.99 , 25.81 , 25.56 ,  25.54 ,  25.65 ,           26.15 ,  26.78 ,  27.48 ,  27.69 , 27.82 ,  27.66 , 27.54 , 27.19 ,  26.96 ,  26.98 , 26.45 ,           26.16 ,  26.36 ,  27.12 ,  27.69 , 27.59 ,  27.36 , 26.94 , 26.59 ,  26.66 ,  26.49 , 26.64  , 26.50]


# ### Plot the ECCOv4r3 and ESRL nino 3.4 index

# In[33]:


plt.figure(figsize=(8,5), dpi= 90)
plt.plot(SST_masked.time,SST_nino_34_anom - SST_nino_34_anom.mean(),'ro-')
plt.plot(SST_masked.time, nino34_2008_2013 - np.mean(nino34_2008_2013),'kx-')
plt.title('ECCO v4r3 nino 3.4 SST Anomaly');
plt.legend(('ECCO','ESRL'))
plt.ylabel('deg C');
plt.xlabel('year');
plt.grid()


# Wow, it's almost like SST is an easy measurement for ECCO to fit!

# ## Conclusion
# 
# You should now be familiar with doing some doing calculations using scalar quantities.
# 
# ### Suggested exercises
# 
# 1. Create the SSH time series from 1992-2015
# 2. Create the global mean sea level trend (map) from 1992-2015
# 3. Create the global mean sea level trend (map) for two epochs 1992-2003, 2003-2015
# 4. Compare other nino indices
# 
