#!/usr/bin/env python
# coding: utf-8

# # Compute meridional heat transport

# This notebook shows how to compute meridional heat transport (MHT) across any particular latitude band. 
# Additionally, we show this for both global and basin specific cases. 
# 
# Oceanographic computations:
# 
# * use [xgcm](https://github.com/xgcm/xgcm) to compute masks and grab values at a particular latitude band
# 
# * use [ecco_v4_py](https://github.com/ECCO-GROUP/ECCOv4-py) to select a specific basin
# 
# * compute meridional heat transport at one or more latitude bands
# 
# Python basics on display:
# 
# * how to [open a dataset using xarray](http://xarray.pydata.org/en/stable/generated/xarray.open_dataset.html) (a one liner!)
# 
# * how to [save a dataset using xarray](http://xarray.pydata.org/en/stable/generated/xarray.Dataset.to_netcdf.html)  (another one liner!)
# 
# * one method for making subplots 
# 
# * some tricks for plotting quantities defined as [dask arrays](http://docs.dask.org/en/latest/array.html)
# 
# Note that each of these tasks can be accomplished more succinctly with [ecco_v4_py](https://github.com/ECCO-GROUP/ECCOv4-py) functions, but are shown explicitly to illustrate these tools. 
# Throughout, we will note the ecco_v4_py (python) and [gcmfaces](https://github.com/ECCO-GROUP/gcmfaces) (MATLAB) functions which can perform these computations.

# In[1]:


import warnings
warnings.filterwarnings('ignore')


# In[2]:


import os
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import cartopy as cart
import sys


# ## Load Model Variables
# 
# Because we're computing transport, we want the files containing 'UVELMASS' and 'VVELMASS' for volumetric transport, and 'ADVx_TH', 'ADVy_TH' and 'DFxE_TH', 'DFyE_TH' for the advective and diffusive components of heat transport, respectively. 

# In[3]:


## Import the ecco_v4_py library into Python
## =========================================

## -- If ecco_v4_py is not installed in your local Python library, 
##    tell Python where to find it.  For example, if your ecco_v4_py
##    files are in /Users/ifenty/ECCOv4-py/ecco_v4_py, then use:

sys.path.append('/Users/ifenty/ECCOv4-py')
import ecco_v4_py as ecco


# In[4]:


## Set top-level file directory for the ECCO NetCDF files
## =================================================================
# base_dir = '/home/username/'
base_dir = '/Users/ifenty/ECCOv4-release'

## define a high-level directory for ECCO fields
ECCO_dir = base_dir + '/Release3_alt'


# In[5]:


## Load the model grid
grid_dir= ECCO_dir + '/nctiles_grid/'

ecco_grid = ecco.load_ecco_grid_nc(grid_dir)

## Load one year of 2D daily data, SSH, SST, and SSS 
data_dir= ECCO_dir + '/nctiles_monthly'

ecco_vars = ecco.recursive_load_ecco_var_from_years_nc(data_dir,                                            vars_to_load=['ADVx_TH', 'ADVy_TH',                                                         'DFxE_TH', 'DFyE_TH'],                                           years_to_load = range(2008,2014))
## Merge the ecco_grid with the ecco_vars to make the ecco_ds
ecco_ds = xr.merge((ecco_grid , ecco_vars))


# ## Grab latitude band: 26$^\circ$N array as an example
# 
# Here we want to grab the transport values which along the band closest represented in the model to 26$^\circ$N. 
# In a latitude longitude grid this could simply be done by, e.g. `U.sel(lat=26)`. 
# However, the LLC grid is slightly more complicated. 
# Luckily, the functionality enabled by the [xgcm Grid object](https://xgcm.readthedocs.io/en/latest/api.html#grid) makes this relatively easy. 

# Note that this subsection can be performed with the with the ecco_v4_py modules [vector_calc](https://github.com/ECCO-GROUP/ECCOv4-py/blob/master/ecco_v4_py/vector_calc.py) and [scalar_calc](https://github.com/ECCO-GROUP/ECCOv4-py/blob/master/ecco_v4_py/scalar_calc.py) as follows:
# 
# ```
# from ecco_v4_py import vector_calc, scalar_calc
# grid = ecco_v4_py.get_llc_grid(ds)
# rapid_maskW, rapid_maskS = vector_calc.get_latitude_masks(lat_val=26,yc=ds.YC,grid=grid)
# rapid_maskC = scalar_calc.get_latitude_mask(lat_val=26,yc=ds.YC,grid=grid)
# ```
# 
# One can also use the gcmfaces function [gcmfaces_calc/gcmfaces_lines_zonal.m](https://github.com/ECCO-GROUP/gcmfaces/blob/master/gcmfaces_calc/gcmfaces_lines_zonal.m).

# In[6]:


# Get array of 1's at and north of latitude
lat = 26
ones = xr.ones_like(ecco_ds.YC)
dome_maskC = ones.where(ecco_ds.YC>=lat,0)


# In[7]:


plt.figure(figsize=(12,6))
fig,ax,p,cb=ecco.plot_proj_to_latlon_grid(ecco_ds.XC,ecco_ds.YC,dome_maskC,
                              projection_type='robin',cmap='viridis',user_lon_0=0,show_colorbar=True)
ax.add_feature(cart.feature.LAND,facecolor='0.7',zorder=2)
plt.show()


# Again, if this were a lat/lon grid we could simply take a finite difference in the meridional direction. 
# The only grid cells with 1's remaining would be at the southern edge of grid cells at approximately 26$^\circ$N.
# 
# However, recall that the LLC grid has a different picture.

# In[8]:


# Mask masks and add them to ecco_ds

maskC = ecco_ds.hFacC.where(ecco_ds.hFacC>0, np.nan)
maskW = ecco_ds.hFacW.where(ecco_ds.hFacW>0, np.nan)
maskS = ecco_ds.hFacS.where(ecco_ds.hFacS>0, np.nan)

maskW.name = 'maskW'
maskS.name = 'maskS'
maskC.name = 'maskC'
ecco_ds = xr.merge((ecco_ds, maskW, maskC, maskS))


# In[9]:


plt.figure(figsize=(12,6))
ecco.plot_tiles(dome_maskC+maskC.isel(k=0), cmap='viridis')
plt.show()


# Recall that for tiles 7-12, the y-dimension actually runs East-West. 
# Therefore, we want to compute a finite difference in the x-dimension on these tiles to get the latitude band at 26$^\circ$N. 
# For tiles 1-5, we clearly want to difference in the y-dimension. 
# Things get more complicated on tile 6.

# Here we make the [xgcm Grid object](https://xgcm.readthedocs.io/en/latest/api.html#grid) which allows us to compute finite differences in simple one liners. 
# This object understands how each of the tiles on the LLC grid connect, because we provide that information. 
# To see under the hood, checkout the utility function [get_llc_grid](https://github.com/ECCO-GROUP/ECCOv4-py/blob/master/ecco_v4_py/ecco_utils.py) where these connections are defined. 

# In[10]:


ecco_ds


# In[11]:


grid = ecco.get_llc_grid(ecco_ds)


# In[12]:


lat_maskW = grid.diff(dome_maskC,'X',boundary='fill')
lat_maskS = grid.diff(dome_maskC,'Y',boundary='fill')


# In[13]:


plt.figure(figsize=(12,6))
ecco.plot_tiles(lat_maskW+maskW.isel(k=0), cmap='viridis')
plt.show()


# In[14]:


plt.figure(figsize=(12,6))
ecco.plot_tiles(lat_maskS+maskS.isel(k=0), cmap='viridis',show_colorbar=True)
plt.show()


# ## Select the Atlantic ocean basin for RAPID-MOCHA MHT

# Now that we have 26$^\circ$N we just need to select the Atlantic. 
# This can be done with the [ecco_v4_py.get_basin](https://github.com/ECCO-GROUP/ECCOv4-py/blob/master/ecco_v4_py/get_basin.py) module, specifically `ecco_v4_py.get_basin.get_basin_mask`.
# Note that this function requires a mask as an input, and then returns the values within a particular ocean basin.
# Therefore, provide the function with `ds['maskC']` for a mask at tracer points, `ds['maskW']` for west grid cell edges, etc.
# 
# Note: this mirrors the gcmfaces functionality [ecco_v4/v4_basin.m](https://github.com/ECCO-GROUP/gcmfaces/blob/master/ecco_v4/v4_basin.m). 

# Here we just want the Atlantic ocean, but lets look at what options we have ...

# In[15]:


print(ecco.get_available_basin_names())


# Notice that, for instance, 'mexico' exists separately from the Atlantic ('atl'). 
# This is to provide as fine granularity as possible (and sensible). 
# To grab the broader Atlantic ocean basin, i.e. the one people usually refer to, use the option 'atlExt'. 
# Similar options exist for the Pacific and Indian ocean basins.

# In[16]:


atl_maskW = ecco.get_basin_mask(basin_name='atlExt',mask=maskW.isel(k=0))
atl_maskS = ecco.get_basin_mask(basin_name='atlExt',mask=maskS.isel(k=0))


# Notice that we pass the routine a 2D mask by selecting the first depth level. This is simply to make things run faster.

# In[17]:


plt.figure(figsize=(12,6))
fig,ax,p,cb=ecco.plot_proj_to_latlon_grid(ecco_ds.XC,ecco_ds.YC,atl_maskW,
                              projection_type='robin',cmap='viridis',user_lon_0=0,show_colorbar=True)
ax.add_feature(cart.feature.LAND,facecolor='0.7',zorder=2)
plt.show()


# ## MHT at the approximate RAPID array latitude
# 
# This can be done with the [ecco_v4_py.calc_meridional_trsp](https://github.com/ECCO-GROUP/ECCOv4-py/blob/master/ecco_v4_py/calc_meridional_trsp.py) module for heat, salt, and volume transport as follows:
# 
# ```
# mvt = ecco_v4_py.calc_meridional_vol_trsp(ecco_ds,lat_vals=26,basin_name='atlExt')
# mht = ecco_v4_py.calc_meridional_heat_trsp(ecco_ds,lat_vals=26,basin_name='atlExt')
# mst = ecco_v4_py.calc_meridional_salt_trsp(ecco_ds,lat_vals=26,basin_name='atlExt')
# ```
# 
# Additionally, one could compute the overturning streamfunction at this latitude band with `ecco_v4_py.calc_meridional_stf`. 
# The inputs are the same as the functions above, see the module [ecco_v4_py.calc_stf](https://github.com/ECCO-GROUP/ECCOv4-py/blob/master/ecco_v4_py/calc_stf.py).
# 
# In MATLAB, one can use the functions: 
# 
# * compute meridional transports: [gcmfaces_calc/calc_MeridionalTransport.m](https://github.com/ECCO-GROUP/gcmfaces/blob/master/gcmfaces_calc/calc_MeridionalTransport.m)
# 
# * compute the overturning streamfunction: [gcmfaces_calc/calc_overturn.m](https://github.com/ECCO-GROUP/gcmfaces/blob/master/gcmfaces_calc/calc_overturn.m)
# 
# 

# ### A note about computational performance
# 
# When we originally open the dataset with all of the variables, we don't actually load anything into memory. In fact, nothing actually happens until "the last minute". 
# For example, the data are only loaded once we do any computation like compute a mean value, plot something, or explicitly provide a `load` command for either the entire dataset or an individual DataArray within the dataset. 
# This 'lazy execution' is enabled by the data structure underlying the xarray Datasets and DataArrays, the [dask array](https://docs.dask.org/en/latest/array.html). 
# 
# In short, the when the data are opened, dask builds an execution task graph which it saves up to execute at the last minute. 
# Dask also allows for parallelism, and by default runs in parallel across [threads for a single core architecture](https://docs.dask.org/en/latest/scheduling.html#local-threads). 
# For now, this is what we will show. 
# 
# Some quick tips are: 
# 
# * Don't load the full 25 years of ECCOv4r3 output unless you're on a machine with plenty of memory. I am doing this in the cell below because I'm on a Skylake node with 192GB. Proceed with caution before copying and pasting.
# 
# * If you're in this situation where you can't load all months into memory, it's a good idea to load 
#   a final result before plotting, in case you need to plot it many times in a row, see below...

# In[18]:


get_ipython().run_cell_magic('time', '', "ecco_ds['ADVx_TH'].load();\necco_ds['DFxE_TH'].load();\necco_ds['ADVy_TH'].load();\necco_ds['DFyE_TH'].load();\necco_ds['dxG'].load();\necco_ds['dyG'].load();\necco_ds['drF'].load();")


# In[19]:


get_ipython().run_cell_magic('time', '', "trsp_x = ((ecco_ds['ADVx_TH'] + ecco_ds['DFxE_TH']) *  atl_maskW).sum('k')\ntrsp_y = ((ecco_ds['ADVy_TH'] + ecco_ds['DFyE_TH']) *  atl_maskS).sum('k')")


# In[20]:


get_ipython().run_cell_magic('time', '', "trsp_x = trsp_x * lat_maskW\ntrsp_y = trsp_y * lat_maskS\n\n# Sum horizontally \ntrsp_x = trsp_x.sum(dim=['i_g','j','tile'])\ntrsp_y = trsp_y.sum(dim=['i','j_g','tile'])\n\n# Add together to get transport at depth level \ntrsp_at_depth = trsp_x + trsp_y\n\n# Sum over full depth and convert to PW\nmht = (trsp_at_depth *10**-15 * 1000 * 4000)\nmht.attrs['units']='PW'")


# #### Now that we have computed MHT, lets load the result for iterative plotting
# 
# For some reason when dask arrays are plotted, they are computed on the spot but don't stay in memory.
# This takes a bit to get the hang of, but keep in mind that this allows us to scale the same code on distributed architecture, so we could use these same routines for high resolution output. This seems worthwhile!
# 
# Note that we probably don't need this load statement if we have already loaded the underlying datasets.

# In[21]:


plt.figure(figsize=(12,6))
plt.plot(mht['time'],mht,'r')
plt.grid()
plt.title('Monthly Meridional Heat Transport at 26N')
plt.ylabel(f'[{mht.attrs["units"]}]')
plt.show()


# ## Now compare global and Atlantic MHT at many latitudes

# In[22]:


# pick a temp directory for yourself
nc_save_dir = '/Users/ifenty/tmp/'

if not os.path.isdir(nc_save_dir):
    os.makedirs(nc_save_dir)
    
nc_file = f'{nc_save_dir}/eccov4r3_mht.nc'


# In[23]:


ecco_ds


# In[24]:


del(ecco)
import ecco_v4_py as ecco


# In[25]:


list(ecco_ds.keys())


# In[26]:


global_lats = np.arange(-60,60,1)
mht = ecco.calc_meridional_heat_trsp(ecco_ds, lat_vals=global_lats)
mht = mht.rename({'heat_trsp':'global_heat_trsp'})
mht = mht.rename({'heat_trsp_z':'global_heat_trsp_z'})
print(' --- Done with global --- ')


# In[27]:


plt.figure(figsize=(12,6))
plt.plot(mht['lat'], mht['global_heat_trsp'].mean('time'))


# In[28]:


basin_lats = np.arange(-35,60,1)
atl = ecco.calc_meridional_heat_trsp(ecco_ds,lat_vals=basin_lats,basin_name='atlExt')
atl = atl.rename({'heat_trsp':'atl_heat_trsp'})
atl = atl.rename({'heat_trsp_z':'atl_heat_trsp_z'})
print(' --- Done with Atlantic --- ')


# In[29]:


atl


# In[30]:


plt.figure(figsize=(12,6))
plt.plot(atl['lat'], atl['atl_heat_trsp'].mean('time'))


# In[31]:


plt.figure(figsize=(12,6))
plt.plot(mht['lat'], mht['global_heat_trsp'].mean('time'))
plt.plot(atl['lat'], atl['atl_heat_trsp'].mean('time'))
plt.legend(('Global','Atlantic'))
plt.grid(linestyle='--')
plt.title(f'Meridional Heat Transport [{mht["global_heat_trsp"].attrs["units"]}]')
plt.ylabel(f'[{mht["global_heat_trsp"].attrs["units"]}]')
plt.xlabel('Latitude')
plt.show()


# ## MHT as a function of depth
# 

# In[43]:


def lat_depth_plot(mht,fld,label):
    fig = plt.figure(figsize=(12,6))
    
    # Set up depth coordinate
    depth = -mht['Z']
    stretch_depth = 200
    
    # Set up colormap and colorbar
    cmap = 'RdBu_r'
    fld_mean = mht[fld].mean('time')
    abs_max = np.max(np.abs([fld_mean.min(),fld_mean.max()]))
    cmin = -abs_max*.15
    cmax = -cmin
    
    # First top 500m
    ax1 = plt.subplot(2,1,1)
    p1 = ax1.pcolormesh(mht['lat'],depth,fld_mean,cmap=cmap,vmin=cmin,vmax=cmax)

    # Handle y-axis
    ax1.invert_yaxis()
    plt.ylim([stretch_depth, 0])
    ax1.yaxis.axes.set_yticks(np.arange(stretch_depth,0,-50))
    plt.ylabel(f'Depth [{mht["Z"].attrs["units"]}]')

    # Remove top plot xtick label
    ax1.xaxis.axes.set_xticklabels([])

    # Now the rest ...
    ax2 = plt.subplot(2,1,2)
    p2 = ax2.pcolormesh(mht['lat'],depth,fld_mean,cmap=cmap,vmin=cmin,vmax=cmax)

    # Handle y-axis
    ax2.invert_yaxis()
    plt.ylim([depth.max(), stretch_depth])
    yticks = np.flip(np.arange(6000,2*stretch_depth,-1000))
    ax2.yaxis.axes.set_yticks(yticks)
    plt.ylabel(f'Depth [{mht["Z"].attrs["units"]}]')
               
    # Label  axis
    plt.xlabel('Latitude')

    # Reduce space between subplots
    fig.subplots_adjust(hspace=0.05)

    # Make a single title
    fig.suptitle(f'{label} time mean meridional heat transport',verticalalignment='top',fontsize=24)

    # Make an overyling colorbar
    fig.subplots_adjust(right=0.83)
    cbar_ax = fig.add_axes([0.87, 0.1, 0.025, 0.8])
    fig.colorbar(p2,cax=cbar_ax)
    cbar_ax.set_ylabel(f'[{mht[fld].attrs["units"]}]')
                       
    plt.show()


# In[44]:


lat_depth_plot(mht,'global_heat_trsp_z','Global')
lat_depth_plot(atl,'atl_heat_trsp_z','Atlantic')


# ## Exercise: reproduce figure from (Ganachaud and Wunsch, 2000)

# In[34]:


from IPython.display import Image
Image('../figures/buckley_mht.png')


# Figure from (Buckley and Marshall, 2016), which is adapted from (Ganachaud and Wunsch, 2000).
# 
# Note: to do this, you may need to pair latitude masks with those defined through the `get_section_masks` module. An example of this functionality is shown in the next tutorial. 

# ## References
# 
# Buckley, M. W. and Marshall, J. ( 2016), Observations, inferences, and mechanisms of the Atlantic Meridional Overturning Circulation: A review, Rev. Geophys., 54, doi:10.1002/2015RG000493.
# 
# Ganachaud, A., & Wunsch, C. (2000). Improved estimates of global ocean circulation, heat transport and mixing from hydrographic data. Nature, 408(6811), 453-7. doi:http://dx.doi.org.ezproxy.lib.utexas.edu/10.1038/35044048
