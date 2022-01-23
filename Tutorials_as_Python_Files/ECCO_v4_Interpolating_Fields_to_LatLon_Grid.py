#!/usr/bin/env python
# coding: utf-8

# # Interpolating ECCO fields from llc90 to lat-lon Grids + Saving Results to netCDF  
# 
# 
# ## Objectives
# 
# 1. Learn how to interpolate scalar and vector fields from ECCOv4's lat-lon-cap 90 (llc90) model grid to the more common regular latitude-longitude grid.  
# 
# 2. Learn how to save these interpolated fields as netCDF for later analysis 
# 
# ## Introduction
# 
# Recall the orientations of the 13 tiles of the ECCOv4 native llc90 model grid.
# 
# ![llc90 tile layout](../figures/llc90_0.png)
# 
# Tiles 7-12 are rotated 90 degrees counter-clockwise relative to tiles 0-5.
# 
# In this tutorial we demonstrate two methods for mapping scalar and vector fields from the llc90 model grid to "regular" latitude-longitude grids of arbitrary resolution.  
# 
# > **Note:**  *There are many methods on can use to map between the grids (e.g., nearest neighbor, bilinear interpolation, bin-averaging, etc.), each with its own advantages.)*
# 
# ## How to interpolate scalar ECCO fields to a lat-lon grid
# 
# Scalar fields are the most straightforward fields to interpolate.
# 
# ### Preliminaries: Load fields
# 
# First, let's load the all 13 tiles for sea surface height and the model grid parameters.

# In[1]:


import numpy as np
import sys
import xarray as xr
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import warnings
warnings.filterwarnings('ignore')
from pprint import pprint
import importlib


# In[2]:


## Import the ecco_v4_py library into Python
## =========================================

## -- If ecco_v4_py is not installed in your local Python library, 
##    tell Python where to find it.  For example, if your ecco_v4_py
##    files are in /Users/ifenty/ECCOv4-py/ecco_v4_py, then use:

from pathlib import Path
sys.path.append(str(Path('c:/Users/Ian/ECCOv4-py')))
import ecco_v4_py as ecco


# In[3]:


## Set top-level file directory for the ECCO NetCDF files
## =================================================================
# base_dir = homehome/username/'
ECCO_dir = Path('E:/inSync Share/Projects/ECCOv4/Release4/')


# In[4]:


## Load the model grid
grid_dir= ECCO_dir / 'nctiles_grid/'


# In[5]:


## Load the model grid
ecco_grid = ecco.load_ecco_grid_nc(grid_dir, 'ECCO-GRID.nc')                                    


# In[6]:


## Load one year of 2D daily data, SSH, SST, and SSS 
day_mean_dir= ECCO_dir / 'nctiles_monthly/'

ecco_vars = ecco.recursive_load_ecco_var_from_years_nc(day_mean_dir,                                            vars_to_load=['SSH'],                                            years_to_load=2000,less_output=True)

ecco_ds = []

## Merge the ecco_grid with the ecco_vars to make the ecco_ds
ecco_ds = xr.merge((ecco_grid , ecco_vars)).load()

pprint(ecco_ds.data_vars)


# ### Plotting the dataset
# 
# Plotting the ECCOv4 fields was covered in an earlier tutorial.  Before demonstrating interpolation, we will first plot one of our SSH fields. Take a closer look at the arguments of ``plot_proj_to_latlon_grid``, ``dx=2`` and ``dy=2``.

# In[7]:


plt.figure(figsize=(12,6), dpi= 90)

tmp_plt = ecco_ds.SSH.isel(time=1)
tmp_plt = tmp_plt.where(ecco_ds.hFacC.isel(k=0) !=0)

ecco.plot_proj_to_latlon_grid(ecco_ds.XC,                               ecco_ds.YC,                               tmp_plt,                               plot_type = 'pcolormesh',                               dx=2,                              dy=2,                               projection_type = 'robin',                              less_output = True);


# These ``dx`` and ``dy`` arguments tell the plotting function to interpolate the native grid tiles onto a lat-lon grid with spacing of ``dx`` degrees longitude and ``dy`` degrees latitude.  If we reduced ``dx`` and ``dy`` the resulting map would have finer features: Compare with the map when ``dx``=``dy``=0.25 degrees:

# In[8]:


plt.figure(figsize=(12,6), dpi= 90)

tmp_plt = ecco_ds.SSH.isel(time=1)
tmp_plt = tmp_plt.where(ecco_ds.hFacC.isel(k=0) !=0)

ecco.plot_proj_to_latlon_grid(ecco_ds.XC,                               ecco_ds.YC,                               tmp_plt,                               plot_type = 'pcolormesh',                               dx=0.25,                              dy=0.25,                               projection_type = 'robin',                              less_output = True);


# Of course you can interpolate to arbitrarily high resolution lat-lon grids, the model resolution will ultimately determine the smallest resolvable features. 
# 
# Under the hood of ``plot_proj_to_lat_longrid`` is a call to the very useful routine ``resample_to_latlon`` which is the now described in more detail
# 
# ## ``resample_to_latlon``
# 
# ``resample_to_latlon`` takes a field with a cooresponding set of lat lon coordinates (the *source* grid) and interpolates to a new lat-lon *target* grid.  The arrangement of coordinates in the *source* grid is arbitrary. One also provides the $bounds$ of the new lat lon grid.  The example shown above uses -90..90N and -180..180E by default.  In addition, one specifies which interpolation scheme to use (*mapping method*) and the *radius of influence*, the radius around the *target* grid cells within which to search for values from the *source* grid.   

# In[9]:


help(ecco.resample_to_latlon)


# ## Demonstrations of  ``resample_to_latlon``
# 
# ### Global
# First we will map to a global lat-lon grid at 1degree using nearest neighbor

# In[10]:


new_grid_delta_lat = 1
new_grid_delta_lon = 1

new_grid_min_lat = -90+new_grid_delta_lat/2
new_grid_max_lat = 90-new_grid_delta_lat/2

new_grid_min_lon = -180+new_grid_delta_lon/2
new_grid_max_lon = 180-new_grid_delta_lon/2

new_grid_lon, new_grid_lat, field_nearest_1deg =        ecco.resample_to_latlon(ecco_ds.XC,                                 ecco_ds.YC,                                 ecco_ds.SSH.isel(time=0),                                new_grid_min_lat, new_grid_max_lat, new_grid_delta_lat,                                new_grid_min_lon, new_grid_max_lon, new_grid_delta_lon,                                fill_value = np.NaN,                                 mapping_method = 'nearest_neighbor',
                                radius_of_influence = 120000)


# In[11]:


# Dimensions of the new grid:
pprint(new_grid_lat.shape)
pprint(new_grid_lon.shape)


# In[12]:


pprint(new_grid_lon)

# the second dimension of new_grid_lon has the center longitude of the new grid cells
pprint(new_grid_lon[0,0:10])


# In[13]:


# The set of lat points
pprint(new_grid_lat)

# or as a 1D vector
# the first dimension of new_grid_lat has the center latitude of the new grid cells
pprint(new_grid_lat[0:10,0])


# In[14]:


# plot the whole field
plt.figure(figsize=(12,6), dpi= 90)
plt.imshow(field_nearest_1deg,origin='lower')


# Notice that although we specified a fill_value of np.nan, the values on land are zero.  This is because we did not first mask out original field with nans over dry points.  Let's nan out the land points first and interpolate again:

# In[15]:


original_field_with_land_mask= np.where(ecco_ds.maskC.isel(k=0)>0, ecco_ds.SSH.isel(time=0), np.nan)

new_grid_delta_lat = 1
new_grid_delta_lon = 1

new_grid_min_lat = -90+new_grid_delta_lat/2
new_grid_max_lat = 90-new_grid_delta_lat/2

new_grid_min_lon = -180+new_grid_delta_lon/2
new_grid_max_lon = 180-new_grid_delta_lon/2
new_grid_lon, new_grid_lat, field_nearest_1deg =        ecco.resample_to_latlon(ecco_ds.XC,                                 ecco_ds.YC,                                 original_field_with_land_mask,                                new_grid_min_lat, new_grid_max_lat, new_grid_delta_lat,                                new_grid_min_lon, new_grid_max_lon, new_grid_delta_lon,                                fill_value = np.NaN,                                 mapping_method = 'nearest_neighbor',
                                radius_of_influence = 120000)

# plot the whole field, this time land values are nans.
plt.figure(figsize=(12,6), dpi= 90)
plt.imshow(field_nearest_1deg,origin='lower')


# 
# ### Regional
# 
# #### 1 degree, nearest neighbor
# 
# First we'll interpolate only to a subset of the N. Atlantic at 1 degree.

# In[16]:


original_field_with_land_mask= np.where(ecco_ds.maskC.isel(k=0)>0, ecco_ds.SSH.isel(time=0), np.nan)

new_grid_delta_lat = 1
new_grid_delta_lon = 1

new_grid_min_lat = 30+new_grid_delta_lat/2
new_grid_max_lat = 82-new_grid_delta_lat/2

new_grid_min_lon = -90+new_grid_delta_lon/2
new_grid_max_lon = 10-new_grid_delta_lon/2

new_grid_lon, new_grid_lat, field_nearest_1deg =        ecco.resample_to_latlon(ecco_ds.XC,                                 ecco_ds.YC,                                 original_field_with_land_mask,                                new_grid_min_lat, new_grid_max_lat, new_grid_delta_lat,                                new_grid_min_lon, new_grid_max_lon, new_grid_delta_lon,                                fill_value = np.NaN,                                 mapping_method = 'nearest_neighbor',
                                radius_of_influence = 120000)

# plot the whole field, this time land values are nans.
plt.figure(figsize=(6,6), dpi= 90)
plt.imshow(field_nearest_1deg,origin='lower',cmap='jet')


# #### 0.05 degree, nearest neighbor
# 
# Next we'll interpolate the same domain at a much higher resolution, 0.05 degree:

# In[17]:


original_field_with_land_mask= np.where(ecco_ds.maskC.isel(k=0)>0, ecco_ds.SSH.isel(time=0), np.nan)

new_grid_delta_lat = 0.05
new_grid_delta_lon = 0.05

new_grid_min_lat = 30+new_grid_delta_lat/2
new_grid_max_lat = 82-new_grid_delta_lat/2

new_grid_min_lon = -90+new_grid_delta_lon/2
new_grid_max_lon = 10-new_grid_delta_lon/2


new_grid_lon, new_grid_lat, field_nearest_1deg =        ecco.resample_to_latlon(ecco_ds.XC,                                 ecco_ds.YC,                                 original_field_with_land_mask,                                new_grid_min_lat, new_grid_max_lat, new_grid_delta_lat,                                new_grid_min_lon, new_grid_max_lon, new_grid_delta_lon,                                fill_value = np.NaN,                                 mapping_method = 'nearest_neighbor',
                                radius_of_influence = 120000)

# plot the whole field, this time land values are nans.
plt.figure(figsize=(6,6), dpi= 90)
plt.imshow(field_nearest_1deg,origin='lower',cmap='jet');


# The new grid is much finer than the llc90 grid.  Because we are using the nearest neighbor method, you can make out the approximate size and location of the llc90 model grid cells (think about it: nearest neighbor).
# 
# #### 0.05 degree, bin average
# 
# With bin averaging many values from the source grid are 'binned' together and then averaged to determine the value in the new grid.  The numer of grid cells from the source grid depends on the source grid resolution and the radius of influence.  If we were to choose a radius of influence of 120000 m (a little longer than 1 degree longitude at the equator) the interpolated lat-lon map would be smoother than if we were to choose a smaller radius.  To demonstrate:

# bin average with 120000 m radius

# In[18]:


original_field_with_land_mask= np.where(ecco_ds.maskC.isel(k=0)>0, ecco_ds.SSH.isel(time=0), np.nan)

new_grid_lon, new_grid_lat, field_nearest_1deg =        ecco.resample_to_latlon(ecco_ds.XC,                                 ecco_ds.YC,                                 original_field_with_land_mask,                                new_grid_min_lat, new_grid_max_lat, new_grid_delta_lat,                                new_grid_min_lon, new_grid_max_lon, new_grid_delta_lon,                                fill_value = np.NaN,                                 mapping_method = 'bin_average',
                                radius_of_influence = 120000)

# plot the whole field, this time land values are nans.
plt.figure(figsize=(6,6), dpi= 90)
plt.imshow(field_nearest_1deg,origin='lower',cmap='jet');


# bin average with 40000m radius

# In[19]:


original_field_with_land_mask= np.where(ecco_ds.maskC.isel(k=0)>0, ecco_ds.SSH.isel(time=0), np.nan)

new_grid_lon, new_grid_lat, field_nearest_1deg =        ecco.resample_to_latlon(ecco_ds.XC,                                 ecco_ds.YC,                                 original_field_with_land_mask,                                new_grid_min_lat, new_grid_max_lat, new_grid_delta_lat,                                new_grid_min_lon, new_grid_max_lon, new_grid_delta_lon,                                fill_value = np.NaN,                                 mapping_method = 'bin_average',
                                radius_of_influence = 40000)

# plot the whole field, this time land values are nans.
plt.figure(figsize=(6,6), dpi= 90)
plt.imshow(field_nearest_1deg,origin='lower',cmap='jet');


# You may wonder why there are missing (white) values in the resulting map.  It's because there are some lat-lon grid cells whose center is more than 40km away from the center of any grid cell on the source grid. There are ways to get around this problem by specifying spatially-varying radii of influence, but that will have to wait for another tutorial.  
# 
# If you want to explore on your own, explore some of the low-level routines of the pyresample package: https://pyresample.readthedocs.io/en/latest/

# ## Interpolating ECCO vectors fields to lat-lon grids
# 
# Vector fields require a few more steps to interpolate to the lat-lon grid.  At least if what you want are the zonal and meridional vector components.  Why? Because in the llc90 grid, vector fields like ``UVEL`` and ``VVEL``  are defined with positive in the +x and +y directions of the local coordinate system, respectively. A term like ``oceTAUX`` correponds to ocean wind stress in the +x direction and does not correspond to *zonal* wind stress.  To calculate zonal wind stress we need to rotate the vector from the llc90 x-y grid system to the lat-lon grid.
# 
# To demonstrate why rotation is necessary, look at ``oceTAUX`` 

# In[20]:


## Load vector fields
day_mean_dir= ECCO_dir / 'nctiles_monthly/'

ecco_vars = ecco.recursive_load_ecco_var_from_years_nc(day_mean_dir,                                            vars_to_load=['oceTAUX', 'oceTAUY'],                                            years_to_load=2000,less_output=True)

ecco_ds = []

## Merge the ecco_grid with the ecco_vars to make the ecco_ds
ecco_ds = xr.merge((ecco_grid , ecco_vars)).load()

pprint(ecco_ds.data_vars)


# In[21]:


tmp = ecco_ds.oceTAUX.isel(time=0)
tmp_masked = np.where(ecco_ds.maskC.isel(k=0)==1, tmp, np.nan)
ecco.plot_tiles(tmp_masked);


# We can sees the expected positive zonal wind stress in tiles 0-5 because the x-y coordinates of tiles 0-5 happen to approximately line up with the meridians and parallels of the lat-lon grid.  Importantly,  for tiles 7-12 wind stress in the +x direction corresponds to mainly wind stress in the SOUTH direction.  To see the zonal wind stress in tiles 7-12 one has to plot ``oceTAUY`` and recognize that for those tiles positive values correspond with wind stress in the tile's +y direction, which is approximately east-west.  To wit,

# In[22]:


tmp = ecco_ds.oceTAUY.isel(time=0)
tmp_masked = np.where(ecco_ds.maskC.isel(k=0)==1, tmp, np.nan)
ecco.plot_tiles(tmp_masked);


# Great, but if you want have a normal recognizable map of zonal and meridional wind stress, we need to determine the zonal and meridional components of these vectors.  

# ### Vector rotation
# 
# For vector rotation we leverage the very useful XGCM package (https://xgcm.readthedocs.io/en/latest/).
# 
# > **Note:**  *The XGCM documentation contains an MITgcm ECCOv4 Example Page: https://xgcm.readthedocs.io/en/latest/example_eccov4.html.  In that example the dimension 'tile' is called 'face' and the fields were loaded from the binary output of the MITgcm, not the netCDF files that we produce for the official ECCOv4r4 product.  Differences are largely cosmetic.)*
# 
# 
# We use XGCM to map the +x and +y vector components to the grid cell centers from the `u` and `v` grid points of the Arakawa C grid.  The ecco-v4-py routine ``get_llc_grid`` creates the XGCM grid object using a DataSet object containing the following information about the model grid:
# 
# ```
#     i,j,i_g,j_g,k,k_l,k_u,k_p1 .  
# 
# ```
# 
# Our ``ecco_ds`` DataSet does have these fields, they come from the ``ECCO-GRID.nc`` object:

# In[23]:


# dimensions of the ecco_ds DataSet
ecco_ds.dims


# In[24]:


# Make the XGCM object
XGCM_grid = ecco.get_llc_grid(ecco_ds)

# look at the XGCM object.
XGCM_grid

# Depending on how much you want to geek out, you can learn about this fancy XGCM_grid object here:
# https://xgcm.readthedocs.io/en/latest/grid_topology.html


# Once we have the XGCM_grid object, we can use built-in routines of XGCM to interpolate the x and y components of a vector field to the cell centers.

# In[25]:


import xgcm
xfld = ecco_ds.oceTAUX.isel(time=0)
yfld = ecco_ds.oceTAUY.isel(time=0)

velc = XGCM_grid.interp_2d_vector({'X': xfld, 'Y': yfld},boundary='fill')


# velc is a dictionary of the x and y vector components taken to the model grid cell centers. At this point they are not yet rotated!

# In[26]:


velc.keys()


# The magic comes here, with the use of the grid cosine 'cs' and grid  sine 'cs' values of the ``ECCO-GRID`` object:

# In[27]:


# Compute the zonal and meridional vector components of oceTAUX and oceTAUY
oceTAU_E  = velc['X']*ecco_ds['CS'] - velc['Y']*ecco_ds['SN']
oceTAU_N  = velc['X']*ecco_ds['SN'] + velc['Y']*ecco_ds['CS']


# Now we have the zonal and meridional components of the vectors, albeit still on the llc90 grid.

# In[28]:


oceTAU_E_masked = np.where(ecco_ds.maskC.isel(k=0)>0, oceTAU_E, np.nan)
ecco.plot_tiles(oceTAU_E_masked);


# So let's resample ``oceTAU_E`` to a lat-lon grid (that's why you're here, right?) and plot

# In[29]:


new_grid_delta_lat = .5
new_grid_delta_lon = .5

new_grid_min_lat = -90+new_grid_delta_lat/2
new_grid_max_lat = 90-new_grid_delta_lat/2

new_grid_min_lon = -180+new_grid_delta_lon/2
new_grid_max_lon = 180-new_grid_delta_lon/2

new_grid_lon, new_grid_lat, oceTAU_E_masked_latlon =        ecco.resample_to_latlon(ecco_ds.XC,                                 ecco_ds.YC,                                 oceTAU_E_masked,                                new_grid_min_lat, new_grid_max_lat, new_grid_delta_lat,                                new_grid_min_lon, new_grid_max_lon, new_grid_delta_lon,                                fill_value = np.NaN,                                 mapping_method = 'nearest_neighbor',
                                radius_of_influence = 120000)

# plot the whole field, this time land values are nans.
plt.figure(figsize=(12,8), dpi= 90);
plt.imshow(oceTAU_E_masked_latlon,origin='lower',vmin=-.2,vmax=.3,cmap='jet');
plt.title('ECCOv4 Jan 2000 Zonal Wind Stress (N/m^2)');
plt.colorbar(orientation='horizontal');


# Or with extra fanciness:

# In[30]:


plt.figure(figsize=(12,5), dpi= 90)

ecco.plot_proj_to_latlon_grid(ecco_ds.XC, ecco_ds.YC,                               oceTAU_E_masked,                               user_lon_0=180,                              projection_type='PlateCaree',                              plot_type = 'pcolormesh',                               dx=1,dy=1,show_colorbar=True,cmin=-.2, cmax=0.3);


# Compare vs "The Scatterometer Climatology of Ocean Winds (SCOW)",
# 
# http://cioss.coas.oregonstate.edu/scow/zonal_wind_stress.html
# 
# Risien, C.M., and D.B. Chelton, 2008: A Global Climatology of Surface Wind and Wind Stress Fields from Eight Years of QuikSCAT Scatterometer Data. J. Phys. Oceanogr., 38, 2379-2413.

# ![http://cioss.coas.oregonstate.edu/scow/figures/January_SCOW_Zonal_Wind_Stress.png](attachment:image.png)    

# And of course we can't forget about our meridional wind stress:

# In[31]:


oceTAU_N_masked = np.where(ecco_ds.maskC.isel(k=0)>0, oceTAU_N, np.nan)

plt.figure(figsize=(12,6), dpi= 90)
ecco.plot_proj_to_latlon_grid(ecco_ds.XC, ecco_ds.YC,                               oceTAU_N_masked,                               user_lon_0=180,                              projection_type='PlateCaree',                              plot_type = 'pcolormesh',                               dx=1,dy=1,cmin=-.2, cmax=0.3,show_colorbar=True);


# Which we also compare against SCOW:

# ![](http://cioss.coas.oregonstate.edu/scow/figures/January_SCOW_Meridional_Wind_Stress.png)

# ### ``UEVNfromUXVY``
# 
# The ecco-v4-py library include a routine, ``UEVNfromUXVY`` which does the interpolation to the grid cell centers and the rotation in one call:

# In[32]:


xfld = ecco_ds.oceTAUX.isel(time=0)
yfld = ecco_ds.oceTAUY.isel(time=0)

# Compute the zonal and meridional vector components of oceTAUX and oceTAUY
oceTAU_E, oceTAU_N  = ecco.vector_calc.UEVNfromUXVY(xfld, yfld, ecco_ds)


# Before interpolated to the lat-lon grid it is convenient to mask out the land values 

# In[33]:


# mask the rotated vectors
oceTAU_E=oceTAU_E.where(ecco_ds.maskC.isel(k=0))
oceTAU_N=oceTAU_N.where(ecco_ds.maskC.isel(k=0))


# Now plot to verify

# In[34]:


# interpolate to lat-lon
new_grid_lon, new_grid_lat, oceTAU_E_masked_latlon =        ecco.resample_to_latlon(ecco_ds.XC,                                 ecco_ds.YC,                                 oceTAU_E,                                new_grid_min_lat, new_grid_max_lat, new_grid_delta_lat,                                new_grid_min_lon, new_grid_max_lon, new_grid_delta_lon,                                fill_value = np.NaN,                                 mapping_method = 'nearest_neighbor',
                                radius_of_influence = 120000)

# plot the whole field, this time land values are nans.
plt.figure(figsize=(12,8), dpi= 90);
plt.imshow(oceTAU_E_masked_latlon,origin='lower',vmin=-.2,vmax=.3,cmap='jet');
plt.title('ECCOv4 Jan 2000 Zonal Wind Stress (N/m^2)');
plt.colorbar(orientation='horizontal');


# ## Saving interpolated fields to netCDF
# 
# For this demonstration we will rotate the 12 monthly-mean records of oceTAUX and oceTAUY to their zonal and meridional components, interpolate to a lat-lon grids, and then save the output as netCDF format.
# 

# In[35]:


pprint(ecco_ds.oceTAUX.time)


# We will loop through each month, use ``UEVNfromUXVY`` to determine the zonal and meridional components of the ocean wind stress vectors, and then interpolate to a 0.5 degree lat-lon grid.

# In[36]:


oceTAUE = np.zeros((12, 360,720))
oceTAUN = np.zeros((12, 360,720))

new_grid_delta_lat = .5
new_grid_delta_lon = .5

new_grid_min_lat = -90+new_grid_delta_lat/2
new_grid_max_lat = 90-new_grid_delta_lat/2

new_grid_min_lon = -180+new_grid_delta_lon/2
new_grid_max_lon = 180-new_grid_delta_lon/2


for m in range(12):
    cur_oceTAUX = ecco_ds.oceTAUX.isel(time=m)
    cur_oceTAUY = ecco_ds.oceTAUY.isel(time=m)

    # Compute the zonal and meridional vector components of oceTAUX and oceTAUY
    tmp_e, tmp_n  = ecco.vector_calc.UEVNfromUXVY(cur_oceTAUX, cur_oceTAUY, ecco_ds)
    
    # apply landmask
    tmp_e_masked = np.where(ecco_ds.maskC.isel(k=0)>0, tmp_e, np.nan)
    tmp_n_masked = np.where(ecco_ds.maskC.isel(k=0)>0, tmp_n, np.nan)

    # zonal component
    new_grid_lon, new_grid_lat, tmp_e_masked_latlon =        ecco.resample_to_latlon(ecco_ds.XC,                                 ecco_ds.YC,                                 tmp_e_masked,                                new_grid_min_lat, new_grid_max_lat, new_grid_delta_lat,                                new_grid_min_lon, new_grid_max_lon, new_grid_delta_lon,                                fill_value = np.NaN,                                 mapping_method = 'nearest_neighbor',
                                radius_of_influence = 120000)
    
    # meridional component
    new_grid_lon, new_grid_lat, tmp_n_masked_latlon =        ecco.resample_to_latlon(ecco_ds.XC,                                 ecco_ds.YC,                                 tmp_n_masked,                                new_grid_min_lat, new_grid_max_lat, new_grid_delta_lat,                                new_grid_min_lon, new_grid_max_lon, new_grid_delta_lon,                                fill_value = np.NaN,                                 mapping_method = 'nearest_neighbor',
                                radius_of_influence = 120000)

    oceTAUE[m,:] = tmp_e_masked_latlon
    oceTAUN[m,:] = tmp_n_masked_latlon


# In[37]:


# make the new data array structures for the zonal and meridional wind stress fields
oceTAUE_DA = xr.DataArray(oceTAUE,  name = 'oceTAUE', 
                      dims = ['time','latitude','longitude'], 
                      coords = {'latitude': new_grid_lat[:,0],
                                'longitude': new_grid_lon[0,:],
                                'time': ecco_ds.time})

# make the new data array structures for the zonal and meridional wind stress fields
oceTAUN_DA = xr.DataArray(oceTAUN,  name = 'oceTAUN', 
                      dims = ['time','latitude','longitude'], 
                      coords = {'latitude': new_grid_lat[:,0],
                                'longitude': new_grid_lon[0,:],
                                'time': ecco_ds.time})


# Plot the zonal and meridional wind stresses in these new DataArray objects:

# In[38]:


oceTAUE_DA.dims     


# In[39]:


oceTAUE_DA.isel(time=0).plot()


# In[40]:


oceTAUN_DA.isel(time=0).plot()


# ### Saving = easy with xarray
# 

# In[41]:


# meridional component
output_path = Path('C:/Users/Ian/Downloads/oceTAUN_latlon_2000.nc')
oceTAUN_DA.to_netcdf(output_path)

# zonal component
output_path = Path('C:/Users/Ian/Downloads/oceTAUE_latlon_2000.nc')
oceTAUE_DA.to_netcdf(output_path)


# ... done!
