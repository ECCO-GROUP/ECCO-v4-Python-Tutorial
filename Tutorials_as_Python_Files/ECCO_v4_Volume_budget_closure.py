#!/usr/bin/env python
# coding: utf-8

# # ECCOv4 Global Volume Budget Closure
# 
# Here we demonstrate the closure of volume budgets in ECCOv4 configurations.  This notebook is draws heavily from  `evaluating_budgets_in_eccov4r3.pdf` which explains the procedure with Matlab code examples by Christopher G. Piecuch.  See ECCO Version 4 release documents: */doc/evaluating_budgets_in_eccov4r3.pdf*).
# 
# 
# ## Objectives
# 
# Illustrate how volume budgets are closed globally.  
# 
# 
# ## Introduction
# 
# 
# ECCOv4 uses the $z^*$ coordinate system in which the depth of the vertical coordinate, $z^*$ varies with time as:
# 
# \begin{equation}
# z^* = \frac{z - \eta(x,y,t)}{H(x,y) + \eta(x,y,t)} H(x,y)
# \end{equation}
# 
# With $H$ being the model depth, $\eta$ being the model sea level anomaly, and $z$ being depth.
# 
# If the vertical coordinate didn't change through time then volume fluxes across the 'u' and 'v' grid cell faces of a tracer cell could be calculated by multiplying the velocities at the face with the face area:
# 
# volume flux across 'u' face in the +x direction = $\mathit{UVEL}(x,y,k) \times \mathit{drF}(k) \times \mathit{dyG}(x,y) \times \mathit{hFacW}(x,y,k)$
# 
# volume flux across 'v' face in the +y direction = $\mathit{VVEL}(x,y,k) \times \mathit{drF}(k) \times \mathit{dxG}(x,y) \times \mathit{hFacS}(x,y,k)$
# 
# With ``dyG`` and ``dxG`` being the lengths of the 'u' and 'v' faces, ``drF`` being the grid cell height and ``hFacW`` and ``hFacS`` being the vertical fractions of the 'u' and 'v' grid cell faces that are open water (ECCOv4 uses partial cells to better represent bathymetry which can allows 0 < hfac $\le$ 1).
# 
# However, because the vertical coordiate varies with time in the $z^*$ system, the grid cell height ``drF`` varies with time as ``drF``$\times s^*(t)$, with
# 
# \begin{equation}
# s^*(x,y,k,t) = 1 + \frac{\eta(x,y,t)}{H}
# \end{equation}
# 
# with $s^* > 1$ when $\eta > 0$
# 
# Thus, to calculate the volume fluxes grid cell through horizontal faces we must account for the time-varying grid cell face areas:
# 
# volume flux across 'u' in the +x direction face with $z^*$ coordinates = $\mathit{UVEL}(x,y,k) \times \mathit{drF}(k) \times \mathit{dyG}(x,y) \times \mathit{hFacW}(x,y,k) \times s^*(x,y,k,t)$
# 
# volume flux across 'v' in the +y direction face with $z^*$ coordinates = $\mathit{VVEL}(x,y,k) \times \mathit{drF}(k) \times \mathit{dxG}(x,y) \times \mathit{hFacS}(x,y,k) \times s^*(x,y,k,t)$
# 
# 
# To make budget calculations easier we provide the scaled velocities quantities ``UVELMASS`` and ``VVELMASS``,
# 
# \begin{align}
# \mathit{UVELMASS}(x,y,k) = \mathit{UVEL}(x,y,k) \times \mathit{hFacW}(x,y,k) \times s^{*}(x,y,k,t)
# \end{align}
# and 
# \begin{align}
# \mathit{VVELMASS}(x,y,k) = \mathit{VVEL}(x,y,k) \times \mathit{hFacS}(x,y,k) \times s^{*}(x,y,k,t)
# \end{align}
# 
# It is worth noting that the word **mass** in ``UVELMASS`` and ``VVELMASS`` is confusing since there is no mass involved here.  Think of these terms as simply being ``UVEL`` and ``VVEL`` multiplied by the fraction of the grid cell height that is open grid cell face across which the volume transport occurs.  Partial cell bathymetry can make this fraction (hFacW, hFacS) less than one, and the $s^*$ scaling factor further adjusts this fraction higher or lower through time.
# 
# Fully closing the budget requires the vertical volume fluxes across the top and bottom 'w' faces of the grid cell and surface freshwater fluxes.  Regarding vertical volume fluxes, there are no $s^*$ or ``hFac`` equivalent scaling factors that modify our top and bottom grid cell areas.  Therefore, vertical volume fluxes through 'w' faces are simply:
# 
# volume flux across 'w' face in the +z direction = $\mathit{WVEL}(x,y,k) \times \mathit{rA}(x,y)$
# 
# > **Note:** Inexplicably, the term ``WVEL`` is provided with the silly name ``WVELMASS``.  Sometimes it's difficult to ignore other people's poor life choices, but please try to do so here.  Ignore the confusing name, ``WVELMASS`` is identical to ``WVEL``.
# 
# In the $z^*$ coordinate system the depth of the surface grid cell is always $z^* = 0$.  In the MITgcm, ``WELMASS`` at the top of the surface grid cell is the liquid volume flux out of the ocean surface and is proportional to the vertical ocean mass flux, ``oceFWflx``
# 
# ## ``ETAN`` in a Boussinesq Model
# 
# ECCOv4 uses a volume-conserving Boussinesq formulation of the MITgcm.  Because volume is conserved in Boussinesq formulations, seawater density changes do not change model sea level anomaly, ``ETAN``.  The following demonstration of``ETAN`` budget closure considers volumetric fluxes but does not take into consideration expansion/contraction due to changes in density.  Furthermore, model ``ETAN`` changes with the exchange of water between ocean and sea-ice while  in reality the Archimedes principle holds that sea level should not change following the growth or melting of sea-ice because floating sea-ice displaces a volume of seawater equal to its weight.  Thus, ``ETAN`` is not comparable to observed sea level. 
# 
# We correct ``ETAN`` to make a sea surface height field that is comaprable with observations by making three corrections: with a) the "Greatbatch correction", a time varying, globally-uniform correction to ocean volume due to changes in global mean density, b) the inverted barometer (IB) correction (see ``SSHIBC``) and c) the 'sea ice load' correction to account for the displacement of seawater due to submerged sea-ice and snow (see ``sIceLoad``).   A demonstration of these corrections is outside the scope of this tutorial.   Here we focus on closing the model budget keeping in mind that we are neglecting sea level changes from changes in global mean density and the fact that ``ETAN`` does not account for volume displacement due to submerged sea-ice.
# 
# Greatbatch, 1994. J. of Geophys. Res. Oceans, https://doi.org/10.1029/94JC00847

# ## Evaluating the model sea level anomaly ``ETAN`` volume budget

# We will evalute 
# 
# \begin{align}
# \underbrace{\frac{\partial \eta}{\partial t}}_{G_\text{total tendency}} = \underbrace{\int_{-H}^0 \left( -\nabla_{z^*}(s^*\,{\bf v}- \frac{\partial w}{\partial z^*} \right) dz^*}_{G_{\text{volumetric divergence}}} + \underbrace{F}_{{G_{surface fluxes}}}
# \end{align}
# 
# The total tendency of $\eta$, $G_{\text{total tendency}}$ is the sum of the $\eta$ tendencies from volumetric divergence, $G_{\text{volumetric divergence}}$, and volumetric surface fluxes, $G_{\text{surface fluxes}}$. 
# 
# In discrete form, using indexes that start from k=0 (surface tracer cell) and running to k=nk-1 (bottom tracer cell)
# 
# \begin{align}
# \frac{\eta(i,j)}{\partial t} = \sum_{k=nk-1}^{0} \underbrace{\left[\mathit{UVELMASS}(i_g,j,k)-\mathit{UVELMASS}(i_g+1,j,k) \right]\, \mathit{dyG}(i_g,j) \, \mathit{drF}(k) }_{\text{volumetric flux in minus out in x direction}}  + 
# \\
# \sum_{k=nk-1}^{0} \underbrace{\left[\mathit{VVELMASS}(i,j_g,k)-\mathit{VVELMASS}(i,j_g+1,k) \right] \, \mathit{dxG}(i,j_g) \, \mathit{drF}(k)}_{\text{volumetric flux in minus out in y direction}} +
# \\
# \sum_{k_l=nk}^{1} \underbrace{\mathit{WVELMASS}(i,j,k_l)\, \mathit{drA}(i,j)}_{\text{volumetric flux through grid cell bottom surface}} + 
# \\
# \underbrace{\mathit{oceFWflx}(i,j)/ \mathit{rhoConst}}_{\text{volumetric flux through  the top surface of the uppermost tracer cell}}
# \end{align}
# 
# In the above we intentionally sum ``WVELMASS`` fluxes from the BOTTOM surface of the lowermost grid cell (at $k_l = 50$) to the BOTTOM face of the uppermost grid cell ($k_l = 1$) so that we can explicitly include the surface volume flux (forcing) term, $\mathit{oceFWflx}(i,j)/\mathit{rhoConst}$.
# 
# We will calculate $\partial \eta / \partial t$ by differencing instantaneous monthly snapshots of $\eta$ as
# 
# $$\frac{\partial \eta}{\partial t} = \frac{\eta(i,j,t+1) - \eta(i,j,t)}{\Delta t}$$
# 
# The ``UVELMASS, VVELMASS, WVELMASS`` and ``oceFWflx`` terms must be time-average quantities between the monthly $\eta$ snapshots.
# 
# ### Prepare environment and loading the relevant model variables

# In[1]:


import numpy as np
import sys
import xarray as xr
from copy import deepcopy 
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import warnings
import cmocean
warnings.filterwarnings('ignore')


# In[62]:


# Density kg/m^3
rhoconst = 1029

## needed to convert surface mass fluxes to volume fluxes

# lat/lon resolution in degrees to interpolate the model 
# fields for the purposes of plotting
map_dx = .2
map_dy = .2


# In[3]:


## Import the ecco_v4_py library into Python
## =========================================

## -- If ecco_v4_py is not installed in your local Python library, 
##    tell Python where to find it.  For example, if your ecco_v4_py
##    files are in /Users/ifenty/ECCOv4-py/ecco_v4_py, then use:

sys.path.append('/home/ifenty/ECCOv4-py')
import ecco_v4_py as ecco


# ### Load ecco_grid

# In[4]:


## Set top-level file directory for the ECCO NetCDF files
## =================================================================
# base_dir = '/home/username/'
base_dir = '/home/ifenty/ECCOv4-release'

ecco_version = 'v4r3'
    
## define a high-level directory for ECCO fields
ECCO_dir = base_dir + '/Release3_alt'


# In[5]:


## Load the model" grid
grid_dir= ECCO_dir + '/nctiles_grid/'
ecco_grid = ecco.load_ecco_grid_nc(grid_dir, 'ECCOv4r3_grid.nc')


# ### Load 2D MONTHLY $\eta$ snapshots

# In[6]:


data_dir= ECCO_dir + '/nctiles_monthly_snapshots'

year_start = 1993
year_end = 2017


# load one extra year worth of snapshots
ecco_monthly_snaps = ecco.recursive_load_ecco_var_from_years_nc(data_dir,                          vars_to_load=['ETAN'],                         years_to_load=range(year_start, year_end+1)).load()

num_months = len(ecco_monthly_snaps.time.values)
# drop the last 11 months so that we have one snapshot at the beginning and end of each month within the 
# range t1993/1/1 to 2015/1/1)

ecco_monthly_snaps = ecco_monthly_snaps.isel(time=np.arange(0, num_months-11))


# In[7]:


# 1993-01 (beginning of first month) to 2015-01-01 (end of last month, 2014-12)
print(ecco_monthly_snaps.ETAN.time.isel(time=[0, -1]).values)


# In[8]:


# find the record of the last ETAN snapshot
last_record_date = 
    ecco.extract_yyyy_mm_dd_hh_mm_ss_from_datetime64(ecco_monthly_snaps.time[-1].values)
print(last_record_date)
last_record_year = last_record_date[0]


# ### Load MONTHLY mean data 

# In[9]:


data_dir= ECCO_dir + '/nctiles_monthly'

year_end = last_record_year
ecco_monthly_mean = ecco.recursive_load_ecco_var_from_years_nc(data_dir,                         vars_to_load=['oceFWflx',
                                      'UVELMASS',
                                      'VVELMASS', 
                                      'WVELMASS'],\
                        years_to_load=range(year_start, year_end)).load()


# In[10]:


# first and last monthly-mean records
print(ecco_monthly_mean.time.isel(time=[0, -1]).values)


# In[11]:


# each monthly mean record is bookended by a snapshot.  
#we should have one more snapshot than monthly mean record
print('number of monthly mean records: ', len(ecco_monthly_mean.time))
print('number of monthly snapshot records: ', len(ecco_monthly_snaps.time))


# ### Create the xgcm 'grid' object
# 
# the xgcm grid object makes it easy to make flux divergence calculations across different tiles of the lat-lon-cap grid.

# In[12]:


ecco_xgcm_grid = ecco.get_llc_grid(ecco_grid)
ecco_xgcm_grid


# ## Calculate LHS: $\eta$ time tendency: $G_{\text{total tendency}}$
# 
# We calculate the monthly-averaged time tendency of ``ETAN`` by differencing monthly ``ETAN`` snapshots.
# Subtract the numpy arrays $\eta(t+1)$ - $\eta(t)$.  This operation gives us $\Delta$ ``ETAN`` $/ \Delta$ t (month) records.

# In[13]:


num_months = len(ecco_monthly_snaps.time)
G_total_tendency_month =     ecco_monthly_snaps.ETAN.isel(time=range(1,num_months)).values -     ecco_monthly_snaps.ETAN.isel(time=range(0,num_months-1)).values

# The result is a numpy array of 264 months
print('shape of G_total_tendency_month: ', G_total_tendency_month.shape)


# In[14]:


ecco_monthly_mean.oceFWflx.shape


# In[15]:


# Convert this numpy array to an xarray Dataarray to take advantage 
# of xarray time indexing.
# The easiest way is to copy an existing DataArray that has the 
# dimensions and time indexes that we want,
# replace its values, and change its name.  

tmp = ecco_monthly_mean.oceFWflx.copy(deep=True)
tmp.values = G_total_tendency_month
tmp.name = 'G_total_tendency_month'
G_total_tendency_month = tmp

# the nice thing is that now the time values of G_total_tendency_month now line 
# up with the time values of the time-mean fields (middle of the month)
print('\ntime of first array of G_total_tendency_month');
print(G_total_tendency_month.time[0].values)


# Now convert $\Delta$ ``ETAN`` $/ \Delta$ t (month) to $\Delta$ ``ETAN`` $/ \Delta$ t (seconds) by dividing by the number of seconds in each month. To find the number of seconds in each month, subtract the model time step number (which is hourly) from the beginning and end of each month:

# In[16]:


if ecco_version == 'v4r4':
    hrs_per_month = ecco_monthly_snaps.timestep[1:].values -     ecco_monthly_snaps.timestep[0:-1].values
elif ecco_version == 'v4r3':
    hrs_per_month = ecco_monthly_snaps.iter[1:].values -     ecco_monthly_snaps.iter[0:-1].values

# convert hours per month to seconds per month:
secs_per_month = hrs_per_month * 3600

# Make a DataArray with the number of seconds in each month, 
#time indexed to the times in dETAN_dT_perMonth (middle of each month)
secs_per_month = xr.DataArray(secs_per_month,                               coords={'time': G_total_tendency_month.time.values},                               dims='time')

# show number of seconds in the first two months:
print('# of seconds in Jan and Feb 1993 ', secs_per_month[0:2].values)

# sanity check: show number of days in the first two months:
print('# of days in Jan and Feb 1993 ', secs_per_month[0:2].values/3600/24)


# Convert the `ns_in_month` from timedelta64 object to float so we can use it to use it for a mathematical operation: converting G_total_tendency_month to G_total_tendency.  Also, convert from ns to seconds.

# In[17]:


# convert dETAN_dT_perMonth to perSeconds
G_total_tendency = G_total_tendency_month / secs_per_month


# ## Plot the time-mean $\partial \eta / \partial t$, total $\Delta \eta$, and one example $\partial \eta / \partial t$ field
# 
# ### Time-mean $\partial \eta / \partial t$
# 
# To calculate the time averaged $G_{total_tendency}$ one might be tempted to do the following
# 
# ``> G_total_tendency_mean = G_total_tendency.mean('time')``
# 
# But that would be folly because the 'mean' function does not know that the number of days in each month is different!  The result would downweight Februarys and upweight Julys.  We have to weight the tendency records by the length of each month.  A clever way of doing that is provided in the xarray documents: https://xarray-test.readthedocs.io/en/latest/examples/monthly-means.html
# 
# In our case we know the length of each month, we just calculated it above in ``secs_per_month``.  We will weight each month by the relative # of seconds in each month and sum to get a weighted average.

# In[18]:


# the weights are just the # of seconds per month divided by total seconds
month_length_weights = secs_per_month / secs_per_month.sum()


# The time mean of the ``ETAN`` tendency, $\overline{G_{\text{total tendency}}}$, is given by 
# 
# $\overline{G_{\text{total tendency}}} = \sum_{i=1}^{nm} w_i G_{\text{total tendency}}$
# 
# with $\sum_{i=1}^{nm} w_i = 1$ and  nm=number of months

# In[19]:


# the weights sum to 1
print(month_length_weights.sum())


# In[20]:


# the weighted mean weights by the length of each month (in seconds)
G_total_tendency_mean = (G_total_tendency*month_length_weights).sum('time')


# In[21]:


plt.figure(figsize=(20,8));
ecco.plot_proj_to_latlon_grid(ecco_grid.XC, ecco_grid.YC,                              G_total_tendency_mean,show_colorbar=True,                              cmin=-1e-9, cmax=1e-9,                               cmap=cmocean.cm.balance, user_lon_0=-67,                              dx=map_dx,dy=map_dy);
plt.title('Average $\partial \eta / \partial t$ [m/s]', fontsize=20);


# ### Total $\Delta \eta$
# 
# 
# The time average eta tendency is small, about 1 billionth of a meter per second.  The ECCO period is coming up to a billion seconds though...  How much did ``ETAN`` change over the analysis period?

# In[22]:


# the number of seconds in the entire period 
seconds_in_entire_period =     float(ecco_monthly_snaps.time[-1] - ecco_monthly_snaps.time[0])/1e9
print ('seconds in analysis period: ', seconds_in_entire_period)

# which is also the sum of the number of seconds in each month
print('sum of seconds in each month ', secs_per_month.sum().values)


# In[23]:


ETAN_delta = G_total_tendency_mean*seconds_in_entire_period


# In[24]:


plt.figure(figsize=(20,8));
ecco.plot_proj_to_latlon_grid(ecco_grid.XC, ecco_grid.YC,                               ETAN_delta,show_colorbar=True,                              cmin=-1, cmax=1,                               cmap=cmocean.cm.balance, user_lon_0=-67,                               dx=map_dx,dy=map_dy);
plt.title('Predicted $\Delta \eta$ [m] from $\partial \eta / \partial t$',           fontsize=20);


# We can sanity check the total ``ETAN`` change that we found by multipling the time-mean ``ETAN`` tendency with the number of seconds in the simulation by comparing that with the difference in ``ETAN`` between the end of the last month and start of the first month. 

# In[25]:


ETAN_delta_method_2 = ecco_monthly_snaps.ETAN.isel(time=-1).values -     ecco_monthly_snaps.ETAN.isel(time=0).values


# In[26]:


plt.figure(figsize=(20,8));
ecco.plot_proj_to_latlon_grid(ecco_grid.XC, ecco_grid.YC,                               ETAN_delta_method_2,
                              show_colorbar=True,\
                              cmin=-1, cmax=1, \
                              cmap=cmocean.cm.balance, user_lon_0=-67,\
                              dx=map_dx,dy=map_dy);

plt.title('Actual $\Delta \eta$ [m]', fontsize=20);


# In[27]:


plt.figure(figsize=(20,8));
ecco.plot_proj_to_latlon_grid(ecco_grid.XC, ecco_grid.YC,                               ETAN_delta_method_2-ETAN_delta,                              show_colorbar=True,                              cmin=-1e-6, cmax=1e-6,                               cmap=cmocean.cm.balance, user_lon_0=-67,                              dx=map_dx,dy=map_dy);
plt.title('Difference between actual and predicted $\Delta \eta$ [m]',           fontsize=20);


# That's a big woo, these are the same to within 10^-6 meters!

# ### Example $\partial \eta / \partial t$ field

# In[28]:


plt.figure(figsize=(20,8));

# get an array of YYYY, MM, DD, HH, MM, SS for 
#dETAN_dT_perSec at time index 100
tmp = ecco.extract_yyyy_mm_dd_hh_mm_ss_from_datetime64(G_total_tendency.time[100].values)
print(tmp)
ecco.plot_proj_to_latlon_grid(ecco_grid.XC, ecco_grid.YC,                               G_total_tendency.isel(time=100),                              show_colorbar=True,                              cmin=-1e-7, cmax=1e-7,                              cmap=cmocean.cm.balance, user_lon_0=-67,                              dx=map_dx,dy=map_dy);

plt.title('$\partial \eta / \partial t$ [m/s] during ' + 
          str(tmp[0]) +'/' + str(tmp[1]), fontsize=20);


# For any given month the time rate of change of ``ETAN`` is two orders of magnitude smaller than the 1993-2015 mean. In the above we are looking at May 2001.  We see positive ``ETAN`` tendency due sea ice melting in the northern hemisphere (e.g., Baffin Bay, Greenland Sea, and Chukchi Sea).

# ## Calculate RHS: $\eta$ tendency due to surface fluxes, $G_\text{surface fluxes}$
# 
# Surface mass fluxes are given in `oceFWflx`.  Convert surface mass flux to a vertical velocity by dividing by the reference density ``rhoConst``= 1029 kg m-3

# In[29]:


# tendency of eta implied by surface volume fluxes (m/s)
G_surf_fluxes = ecco_monthly_mean.oceFWflx/rhoconst


# ## Plot the time-mean, total, and one month average of $G_{\text{surface fluxes}}$
# 
# ### Time-mean $G_{\text{surface fluxes}}$
# 
# We calculate the time-mean surface flux $\eta$ tendency using the same weights as the total $\eta$ tendency.

# In[30]:


G_surf_fluxes_mean= (G_surf_fluxes*month_length_weights).sum('time')
G_surf_fluxes_mean.shape


# In[31]:


plt.figure(figsize=(20,8));
ecco.plot_proj_to_latlon_grid(ecco_grid.XC, ecco_grid.YC,                               G_surf_fluxes_mean,
                              show_colorbar=True,\
                              cmin=-1e-7, cmax=1e-7, \
                              cmap=cmocean.cm.balance, user_lon_0=-67,
                              dx=map_dx,dy=map_dy)
plt.title('Average $\partial \eta / \partial t$ [m/s] implied by oceFWflx surface mass fluxes\n Negative = Water out',
          fontsize=20);


# ### Total $\Delta \eta$ due to surface fluxes
# 
# If there were no other terms on the RHS to balance surface fluxes, the total change in ``ETAN`` between 1993 and 2015 would be order of h10s of meters almost everwhere.  

# In[32]:


ETAN_delta_surf_fluxes = G_surf_fluxes_mean*seconds_in_entire_period


# In[33]:


plt.figure(figsize=(20,8));
ecco.plot_proj_to_latlon_grid(ecco_grid.XC, ecco_grid.YC,                               ETAN_delta_surf_fluxes,show_colorbar=True,                              cmin=-100, cmax=100,                               cmap=cmocean.cm.balance, user_lon_0=-67,                               dx=map_dx,dy=map_dy);
plt.title('$\Delta \eta$ [m] implied by oceFWflx surface mass fluxes', 
          fontsize=20);


# ### Example $\partial \eta / \partial t$ impliedy by surface fluxes

# In[34]:


plt.figure(figsize=(20,8));

# get an array of YYYY, MM, DD, HH, MM, SS for 
# dETAN_dT_perSec at time index 100
tmp = ecco.extract_yyyy_mm_dd_hh_mm_ss_from_datetime64(G_surf_fluxes.time[100].values)
print(tmp)
ecco.plot_proj_to_latlon_grid(ecco_grid.XC, ecco_grid.YC,                               G_surf_fluxes.isel(time=100),show_colorbar=True,                             cmin=-1e-7, cmax=1e-7,                               cmap=cmocean.cm.balance, user_lon_0=-67,                              dx=map_dx,dy=map_dy);
plt.title('$\partial \eta / \partial t$ [m/s] implied by oceFWflx surface mass fluxes ' + 
          str(tmp[1]) +'/' + str(tmp[0]), fontsize=20);


# For any given month the time rate of change of ``ETAN`` is almost the same as its 22 year mean.  Differences are largest in the high latitudes where sea-ice melt and growth during any particular month induce large changes in ``ETAN``.

# ## Calculate RHS: $\eta$ tendency due to volumetric flux divergence, $G_\text{volumetric fluxes}$
# 
# First we will look at vertical volumetric flux divergence, then horizontal volumetric flux divergence.
# 
# ### Vertical volumetric flux divergence

# It turns out that `WVELMASS` at k_l=0 (the top face of the top tracer cell) is proportional to the ocean surface mass flux `oceFWflx`.  The vertical velocity of the ocean surface is equal to the rate at which water is being added or removed across the top surface of the uppermost grid cell.  This is demonstrated by differencing the velocity at the top 'w' face of the uppermost tracr cell `WVELMASS` (k_l = 0) and the velocity equivalent of transporting the surface mass flux term `oceFWFlx` through this same face.
# 
# First, find the time-mean vertical velocity at the liquid ocean surface

# In[35]:


WVELMASS_surf_mean =     (ecco_monthly_mean.WVELMASS.isel(k_l=0)*month_length_weights).sum('time')


# Next, find the time-mean vertical velocity implied by the `oceFWflx` at k_l=0:

# In[36]:


WVEL_from_oceFWflx_mean =     -(ecco_monthly_mean.oceFWflx*month_length_weights).sum('time')/rhoconst


# In[37]:


plt.figure(figsize=(15,15))
#plt.sca(axs[0,0])
F=ecco.plot_proj_to_latlon_grid(ecco_grid.XC, ecco_grid.YC,                                 WVELMASS_surf_mean,                                show_colorbar=True,                                cmin=-1e-7, cmax=1e-7,                                 cmap=cmocean.cm.balance, user_lon_0=-67,                                dx=2,dy=2, subplot_grid=[3,1,1]);
F[1].set_title('A: surface velocity from WVEL [m/s]')

F=ecco.plot_proj_to_latlon_grid(ecco_grid.XC, ecco_grid.YC,                                 WVEL_from_oceFWflx_mean,                                show_colorbar=True,                                cmin=-1e-7, cmax=1e-7,                                 cmap=cmocean.cm.balance, user_lon_0=-67,                                dx=2,dy=2, subplot_grid=[3,1,2])
F[1].set_title("B: surface velocity implied by oceFWflx [m/s]")

F=ecco.plot_proj_to_latlon_grid(ecco_grid.XC, ecco_grid.YC,                                 WVELMASS_surf_mean-WVEL_from_oceFWflx_mean,                                show_colorbar=True,                                cmin=-1e-7, cmax=1e-7,                                 cmap=cmocean.cm.balance, user_lon_0=-67,                                dx=2,dy=2, subplot_grid=[3,1,3])
F[1].set_title("difference between A and B");


# ``WVELMASS`` at the surface evidently is the same as the surface velocity implied by the surface mass flux ``oceFWflx``.  Thus, we do not actually need ``oceFWflx`` to close the volume budget.  However, to keep the surface forcing term explicitly represented, we will keep ``oceFWflx`` and instead zero out the values of ``WVELMASS`` at the surface so as to avoid double counting. 
# 
# Calculate the vertical volume fluxes at each level: (velocity x grid cell area) [m3 s-1]

# In[38]:


vol_transport_z = ecco_monthly_mean.WVELMASS * ecco_grid.rA


# Set the volume transport at the surface level to be zero because we already counted the fluxes out of the domain with ``oceFWflx``.

# In[39]:


vol_transport_z.isel(k_l=0).values[:] = 0


# Each grid cell has a top and bottom surface and therefore ``WEVELMASS ``should have 51 vertical levels (one more than the number of tracer cells).  For some reason we only have 50 vertical levels, with the bottom of the 50th tracer cell missing.  To calculate vertical flux divergence we need to add this 51st ``WVELMASS`` which is everywhere  zero (no volume flux from the seafloor). The xgcm library handles this in its `diff` routine by specifying the boundary='fill' and fill_value = 0.  

# In[40]:


# volume flux divergence into each grid cell, m^3 / s 
vol_vert_divergence = ecco_xgcm_grid.diff(vol_transport_z, 'Z',                                           boundary='fill', fill_value=0)

# change in eta per unit time due to volumetric vertical convergence 
# at each depth level: m/s
G_vertical_flux_divergence = vol_vert_divergence / ecco_grid.rA;


# In[41]:


# change in eta per unit time due to vertical integral of 
# volumetric horizonal convergence: m/s
G_vertical_flux_divergence_depth_integrated = G_vertical_flux_divergence.sum('k')


# In[42]:


# Calculate the time-mean surface flux $\eta$ tendency using 
# the same weights as the total $\eta$ tendency.
G_vertical_flux_divergence_depth_integrated_time_mean =     (G_vertical_flux_divergence_depth_integrated * month_length_weights).sum('time')


# ### Plot the time-mean $G_{\text{vertical flux divergence}}$

# In[43]:


plt.figure(figsize=(20,8));
ecco.plot_proj_to_latlon_grid(ecco_grid.XC, ecco_grid.YC,                              G_vertical_flux_divergence_depth_integrated_time_mean,
                              show_colorbar=True,\
                              cmap=cmocean.cm.balance, \
                              cmin=-1e-9, cmax=1e-9,
                              dx=map_dx,dy=map_dy);
plt.title('Average $\partial \eta / \partial t$ [m/s] due to vertical flux divergence\n Negative = Water out',
          fontsize=20);


# These values are everywhere essentially zero (numerical noise). On average, the only vertical flux divergence in the column is across the ocean surface. Below the surface, the sum of vertical flux divergence in all tracer cells in the column must be zero because any divergence in any one particular cell is exactly offset by convergence in another cell. Net convergence into the column manifests as a positive vertical velocity at the surface which is equal to oceFWflux in the time-mean. **Thanks to Hong Zhang for comments that improved this explanation**

# ### Horizontal Volume Flux Divergence

# In[44]:


# Volumetric transports in x and y(m^3/s)
vol_transport_x = ecco_monthly_mean.UVELMASS * ecco_grid.dyG * ecco_grid.drF
vol_transport_y = ecco_monthly_mean.VVELMASS * ecco_grid.dxG * ecco_grid.drF


# In[45]:


# Difference of horizontal transports in x and y directions
vol_flux_diff = ecco_xgcm_grid.diff_2d_vector({'X': vol_transport_x,                                                'Y': vol_transport_y},                                              boundary='fill')

# volume flux divergence into each grid cell, m^3 / s 
vol_horiz_divergence = (vol_flux_diff['X'] + vol_flux_diff['Y'])

# change in eta per unit time due to volumetric horizonal
# convergence at each depth level: m/s
# a positive DIVERGENCE leads to negative eta tendency
G_vol_horiz_divergence = -vol_horiz_divergence / ecco_grid.rA

# change in eta in each grid cell per unit time due to horiz. divergence: m/s
G_vol_horiz_divergence_depth_integrated = G_vol_horiz_divergence.sum('k')


# In[46]:


# calculate time-mean using the month length weights
G_vol_horiz_divergence_depth_integrated_mean =     (G_vol_horiz_divergence_depth_integrated * month_length_weights).sum('time')


# ### Plot the time-mean $G_{\text{horizontal flux divergence}}$

# In[47]:


plt.figure(figsize=(20,8));
ecco.plot_proj_to_latlon_grid(ecco_grid.XC, ecco_grid.YC,                               G_vol_horiz_divergence_depth_integrated_mean,
                              show_colorbar=True,\
                              cmin=-1e-7, cmax=1e-7, \
                              cmap=cmocean.cm.balance, user_lon_0=-67,\
                              dx=map_dx,dy=map_dy);
plt.title('Average $\partial \eta / \partial t$ [m/s] implied by horizontal flux divergence\n Positive = Water out',
          fontsize=20);


# ### Plot one example $G_{\text{horizontal flux divergence}}$

# In[48]:


plt.figure(figsize=(20,8));
ecco.plot_proj_to_latlon_grid(ecco_grid.XC, ecco_grid.YC,                               G_vol_horiz_divergence_depth_integrated.isel(time=100),
                              show_colorbar=True,\
                              cmin=-2e-7, cmax=2e-7, \
                              cmap=cmocean.cm.balance, user_lon_0=-67,\
                              dx=map_dx,dy=map_dy);
plt.title('$\partial \eta / \partial t$ [m/s] implied by horizontal flux divergence for May 2001\nPositive = Water out',
          fontsize=20);


# ### Plot $G_{\text{horizontal flux divergence}}$ at 100m depth

# In[49]:


plt.figure(figsize=(20,8));
ecco.plot_proj_to_latlon_grid(ecco_grid.XC, ecco_grid.YC,                               G_vol_horiz_divergence.isel(k=10,time=100),
                              show_colorbar=True,\
                              cmin=-2e-6, cmax=2e-6, \
                              cmap=cmocean.cm.balance, user_lon_0=-67,\
                              dx=map_dx,dy=map_dy);
plt.title('$\partial \eta / \partial t$ [m/s] at 100 m implied by horizontal flux divergence for May 2001\nPositive = Water out',
          fontsize=20);


# ## Comparison of LHS and RHS

# ### Time mean difference of LHR and RHS

# Plot the time-mean difference between the LHS and RHS of the volume budget equation.  

# In[50]:


#LHS ETA TENDENCY
a = G_total_tendency_mean
# RHS ETA TENDENCY FROM VOL DIVERGENCE AND SURFACE FLUXES
b = G_vol_horiz_divergence_depth_integrated_mean
c = G_surf_fluxes_mean

delta = a - b -c 


# In[51]:


plt.figure(figsize=(20,8));
ecco.plot_proj_to_latlon_grid(ecco_grid.XC, ecco_grid.YC, delta,                               show_colorbar=True,                              cmin=-1e-9, cmax=1e-9,                               cmap=cmocean.cm.balance,                               dx=map_dx,dy=map_dy);
plt.title('Residual $\partial \eta / \partial t$ [m/s]: LHS - RHS ', 
          fontsize=20);


# The residual of the time-mean surface velocity tendency terms is essentially zero.  We can look at the distribution of residuals to get a little more confidence.

# ### Histogram of residuals

# In[52]:


tmp = np.abs( a-b-c).values.ravel();
plt.figure(figsize=(10,3));

plt.hist(tmp[np.nonzero(tmp > 0)],np.linspace(0, .5e-10,1000));
plt.grid()


# Almost all residuals < $10^{-11}$ m/s.  We can close the ETAN budget using UVELMASS, VVELMASS, ``WVELMASS`` and ``oceFWflx``.  
# 
# > Note: As stated earlier, we don't actually need ``oceFWflx`` because the surface velocity of ``WVELMASS`` is proportional to ``oceFWflx`` but we kept it so that the surface forcing term is explicit.
# 
# ## ETAN budget closure through time
# 
# ### Global average ETAN budget closure
# 
# Another way of demonstrating volume budget closure is to show the global spatially-averaged ETAN tendency terms through time

# In[53]:


# LHR and RHS through time
a = G_total_tendency
b = G_vol_horiz_divergence_depth_integrated
c = G_surf_fluxes
# residuals
d = a-b-c

area_masked = ecco_grid.rA.where(ecco_grid.hFacC.isel(k=0)> 0)

# take area-weighted mean of these terms
tmp_a=(a*area_masked).sum(dim=('i','j','tile'))/area_masked.sum()
tmp_b=(b*area_masked).sum(dim=('i','j','tile'))/area_masked.sum()
tmp_c=(c*area_masked).sum(dim=('i','j','tile'))/area_masked.sum()
tmp_d=(d*area_masked).sum(dim=('i','j','tile'))/area_masked.sum()

# result is four time series
tmp_a.dims


# In[54]:


fig, axs = plt.subplots(2, 2, figsize=(12,8))

plt.sca(axs[0,0])
tmp_a.plot(color='orange')
axs[0,0].set_title("d(eta)/d(t) total [m/s]")
plt.grid()
plt.ylim([-1e-8, 1e-8]);

plt.sca(axs[0,1])
tmp_b.plot(color='g')
axs[0,1].set_title("d(eta)/d(t) vol divergence [m/s]")
plt.grid()
plt.ylim([-1e-8, 1e-8]);

plt.sca(axs[1,0])
tmp_c.plot(color='r')
axs[1,0].set_title("d(eta)/d(t) surf fluxes [m/s]")
plt.grid()
plt.ylim([-1e-8, 1e-8]);

plt.sca(axs[1,1])
tmp_d.plot(color='b')
axs[1,1].set_title("d(eta)/d(t) LHS - RHS [m/s]")
plt.grid()
plt.subplots_adjust(hspace = .5, wspace=.2)
plt.ylim([-1e-8, 1e-8])
plt.suptitle('Global Volume Budget',fontsize=20);


# When averaged over the entire ocean surface the volumetric divergence has no net impact on $\partial \eta / \partial t$.  This makes sence because horizontal flux divergence can only redistributes volume.  Globally, $\eta$ can only change via net surface fluxes.
# 
# ### Local ETAN budget closure
# 
# Locally we expect that volume divergence can impact $\eta$.  This is demonstrated for a single point the Arctic.

# In[55]:


# Recall, from above...
#a = G_total_tendency
#b = G_vol_horiz_divergence_depth_integrated
#c = G_surf_fluxes
#d = a-b-c

tmp_aa = a.isel(tile=6,j=40,i=29)
tmp_bb = b.isel(tile=6,j=40,i=29)
tmp_cc = c.isel(tile=6,j=40,i=29)
tmp_dd = d.isel(tile=6,j=40,i=29)

fig, axs = plt.subplots(2, 2, figsize=(12,8))
plt.sca(axs[0,0])
tmp_aa.plot(color='orange')
axs[0,0].set_title("d(eta)/d(t) total [m/s]")
plt.ticklabel_format(axis='y', style='sci', useMathText=True)
plt.ylim([-5e-7, 5e-7]);
plt.grid()

plt.sca(axs[0,1])
tmp_bb.plot(color='g')
axs[0,1].set_title("d(eta)/d(t) vol divergence [m/s]")
plt.ticklabel_format(axis='y', style='sci', useMathText=True)
plt.ylim([-5e-7, 5e-7]);
plt.grid()

plt.sca(axs[1,0])
tmp_cc.plot(color='r')
axs[1,0].set_title("d(eta)/d(t) surf fluxes [m/s]")
plt.ticklabel_format(axis='y', style='sci', useMathText=True)
plt.ylim([-5e-7, 5e-7]);
plt.grid()

plt.sca(axs[1,1])
tmp_dd.plot(color='b')
axs[1,1].set_title("d(eta)/d(t) LHS - RHS [m/s]")
plt.ticklabel_format(axis='y', style='sci', useMathText=True)
plt.ylim([-5e-7, 5e-7]);
plt.grid()

plt.subplots_adjust(hspace = .5, wspace=.2)
plt.suptitle('Volume Budget for one Arctic Ocean point',
                 fontsize=20);


# Indeed, the volume divergence term does contribute to $\eta$ variations at this one point.  
# 
# The seasonal cycle of surface fluxes from sea-ice growth and melt can be seen in the surface fluxes term (plotted below for just the year 1995)

# In[56]:


tmp_cc.sel(time='1995').plot(color='r')
plt.ticklabel_format(axis='y', style='sci', useMathText=True)
plt.ylim([-5e-7, 5e-7]);
plt.grid()
plt.title('d(eta)/d(t) from surf. fluxes in 1995 [m/s]');


# ## Predicted vs. actual $\eta$ 
# 
# As we have shown, in our Boussinesq model the only term that can change global mean model sea level anomaly $\eta$ is net surface freshwater flux.  Let us compare the time-evolution of $\eta$ implied by surface freshwater fluxes and the actual $\eta$ from the model output.
# 
# The predicted $\eta$ time series is calculated by time integrating **$G_{surface fluxes}$**.  
# This time series is compared against the actual $\eta$ time series anomaly relative to the $\eta(t=0)$.

# In[57]:


area_masked = ecco_grid.rA.where(ecco_grid.maskC.isel(k=0) == 1)

dETA_per_month_predicted_from_surf_fluxes =     ((G_surf_fluxes * area_masked).sum(dim=('i','j','tile')) / 
     area_masked.sum())*secs_per_month

ETA_predicted_by_surf_fluxes =     np.cumsum(dETA_per_month_predicted_from_surf_fluxes.values)

ETA_from_ETAN = 
    (ecco_monthly_snaps.ETAN * area_masked).sum(dim=('i','j','tile')) /
    area_masked.sum()

# plotting
plt.figure(figsize=(14,5));

plt.plot(dETA_per_month_predicted_from_surf_fluxes.time,          ETA_predicted_by_surf_fluxes,'b.')
plt.plot(ETA_from_ETAN.time.values, ETA_from_ETAN-ETA_from_ETAN[0],'r-')
plt.grid()
plt.ylabel('global mean $\eta$');
plt.legend(('predicted', 'actual'));
plt.title('$\eta(t)$ as predicted from net surface fluxes and model ETAN [m]', 
          fontsize=20);


# The first predicted $\eta$ occurs at the end of the first month (one month time integral of $\partial \eta / \partial t$.  The first *actual* $\eta$ is set to be zero.
# 
# The above plot is another way of confirming that in our Boussinesq model the only term that can change global mean model sea level anomaly $\eta$ is net surface freshwater flux.  To account for changes in global mean density we must apply the Greatbatch correction, inverse-barometer correction, and a correction term to account for the fact that sea-ice does not 'float' on top of the ocean but in fact displaces seawater upwards.  All of these corrections are made for the term ``SSH``, dynamic sea surface height anomaly (not shown here).

# One can compare the sea level rise from mass fluxes in ECCO vs those estimated from GRACE and other datasets published  the WCRP Global Sea Level Budget Group: Global sea-level budget 1993-present, Earth Syst. Sci. Data, 10, 1551-1590, https://doi.org/10.5194/essd-10-1551-2018, 2018. 
# 
# Available here.  See Figure 16.
# https://www.earth-syst-sci-data.net/10/1551/2018/essd-10-1551-2018.pdf

# In[58]:


# Annual mean SL calculation must account for different lengths of each month.

# step 1. weight predicted eta by seconds in each month
tmp1=ETA_predicted_by_surf_fluxes*secs_per_month
# step 2, group records by year and sum across each year
tmp2=tmp1.groupby('time.year').sum()
# step 3, group secs per month by year and sum across each year
secs_per_year = secs_per_month.groupby('time.year').sum()

# step 3, divide time-weighted ETA by seconds per year
annual_mean_GMSL_due_to_mass_fluxes =tmp2/secs_per_year
num_years = len(annual_mean_GMSL_due_to_mass_fluxes.year.values)

plt.figure(figsize=(8,5));
# the -0.13 is to make the starting value comparable with WCRP fig 16.
plt.bar(annual_mean_GMSL_due_to_mass_fluxes.year.values[12:num_years],        1000*(annual_mean_GMSL_due_to_mass_fluxes.values[12:num_years]-.017),         color='k')
plt.grid()
plt.xticks(np.arange(2005, 
                     annual_mean_GMSL_due_to_mass_fluxes.year.values[-1]+1,step=1));
plt.title('Sea level (mm) caused by mass fluxes');


# ## Time-mean $\partial \eta / \partial t$ due to net freshwater fluxes
# 
# We can easily calculate the time mean rate of global sea level rise due to net freshwater flux.

# In[59]:


x=(G_surf_fluxes_mean * area_masked).sum() / area_masked.sum()
total_GMSLR_due_to_mass_fluxes = x*seconds_in_entire_period
print('total global mean sea level rise due to mass fluxes m: ',       total_GMSLR_due_to_mass_fluxes.values)


# Dividing this by the total number of years in the analysis period gives a average rate per year:

# In[60]:


total_number_of_years = len(secs_per_month)/12
total_number_of_years


# In[61]:


mean_rate_of_GMSLR_due_to_mass_fluxes =     total_GMSLR_due_to_mass_fluxes/total_number_of_years

print('mean rate of GMSLR due to mass fluxes [mm/yr] ',        1000*mean_rate_of_GMSLR_due_to_mass_fluxes.values)


# Compare with other estimates of GMSLR due to mass fluxes.
