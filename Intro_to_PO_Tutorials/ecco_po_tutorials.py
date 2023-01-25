def cmap_zerocent_scale(plot,scale_factor):
    """
    Center colormap at zero and scale relative to existing |maximum| value, 
    given plot object and scale_factor, a number of type float.
    Returns new colormap limits as new_clim.
    """
    
    import numpy as np
    
    curr_clim = plot.get_clim()
    new_clim = (scale_factor*np.max(np.abs(curr_clim)))*np.array([-1,1])
    plot.set_clim(new_clim)
    return new_clim

#----------------------------------------------------------------------------------------


def plot_mask(*args,ax=None,color):
    """
    Plot mask, given input parameters:
    - X, Y: (optional) coordinates as 1-D or 2-D NumPy arrays or xarray DataArrays
    - mask: 2-D array of boolean values (True/False or 1/0), NumPy or xarray
    - axes: axes to plot on, defaults to current axes
    - color: a string indicating a color in Matplotlib, or a 3-element tuple or NumPy array indicating RGB color values
    
    Returns plot_obj, the plot object of the mask 
    """
    
    import numpy as np

    if len(args) == 1:
        mask = args[0]
    else:
        X = args[0]
        Y = args[1]
        mask = args[2]
    # set alpha values to 1 where mask is plotted, 0 otherwise
    if str(type(mask))[0:5] == 'xarray':
        mask = mask.values
    # get color for mask
    if isinstance(color,str):
        import matplotlib.colors as mcolors
        color_rgb = mcolors.to_rgb(color)
    elif (isinstance(color,tuple)) and (len(color) == 3):
        color_rgb = np.asarray(color)
    elif (isinstance(color,np.ndarray)) and (len(color) == 3):
        color_rgb = color
    else:
        raise TypeError("input parameter 'color' has incorrect type or number of elements")
    # create a colormap using a 2x4 array with two RGBA entries
    # the RGB entries are the same in each row
    # in the 1st row alpha=0, in the 2nd row alpha=1
    cmap_array = np.hstack((np.tile(color_rgb,(2,1)),np.array([[0],[1]])))
    from matplotlib.colors import ListedColormap
    colormap = ListedColormap(cmap_array)
    # get axis limits of existing plot
    if ax is None:
        ax = plt.gca()
    existing_xlim = ax.get_xlim()
    existing_ylim = ax.get_ylim()
    # plot mask using colormap just created, with alpha=1 where mask=1 or True
    if len(args) == 1:
        plot_obj = ax.pcolormesh(mask,cmap=colormap,vmin=0.,vmax=1.,zorder=50)
    else:
        plot_obj = ax.pcolormesh(X,Y,mask,cmap=colormap,vmin=0.,vmax=1.,zorder=50)
    # set axis limits of mask to axis limits of existing plot
    ax.set_xlim(existing_xlim)
    ax.set_ylim(existing_ylim)
    
    return plot_obj

#----------------------------------------------------------------------------------------


def mean_weighted_binned(value_field,weighting,bin_field,bin_bounds):
    """
    Compute normalized difference in bins, given:
    - value_field: field of values to average
    - weighting: weighting of individual grid cells
    - bin_field: field to use in binning
    - bin_bounds: bound values of bins to use, as numpy array of size (N,2)
    """
    
    import numpy as np
    
    mean_weighted_inbins = np.empty((bin_bounds.shape[0]))
    mean_weighted_inbins.fill(np.nan)
    bin_field_broadcast_flat = (np.ones(value_field.shape)*bin_field).flatten()
    weighting_flat = weighting.flatten()
    value_field_flat = value_field.flatten()
    bin_idx_sorted = np.argsort(bin_field_broadcast_flat)   # flatten and sort bin field values
    sort_idx = (bin_field_broadcast_flat[bin_idx_sorted] >= bin_bounds[0,0]).nonzero()[0][0]
    curr_idx = bin_idx_sorted[sort_idx]
    curr_bin_val = bin_field_broadcast_flat[curr_idx]
    curr_bin_num = ((bin_bounds[:,0] <= curr_bin_val) & (bin_bounds[:,1] > curr_bin_val)).nonzero()[0][0]
    value_sum_inbin = 0.
    weighting_sum_inbin = 0.
    while curr_bin_val < bin_bounds[-1,1]:
        if np.logical_and(~np.isnan(value_field_flat[curr_idx]),~np.isinf(value_field_flat[curr_idx])):
            value_sum_inbin += weighting_flat[curr_idx]*value_field_flat[curr_idx]
            weighting_sum_inbin += weighting_flat[curr_idx]
        sort_idx += 1
        if sort_idx >= len(bin_idx_sorted):
            if weighting_sum_inbin > 0:
                mean_weighted_inbins[curr_bin_num] = value_sum_inbin/weighting_sum_inbin
            break
        curr_idx = bin_idx_sorted[sort_idx]
        curr_bin_val = bin_field_broadcast_flat[curr_idx]
        if curr_bin_val >= bin_bounds[curr_bin_num,1]:
            if weighting_sum_inbin > 0:
                mean_weighted_inbins[curr_bin_num] = value_sum_inbin/weighting_sum_inbin
            curr_bin_num = ((bin_bounds[:,0] <= curr_bin_val) & (bin_bounds[:,1] > curr_bin_val)).nonzero()[0][0]
            value_sum_inbin = 0.
            weighting_sum_inbin = 0.
    return mean_weighted_inbins

#----------------------------------------------------------------------------------------


def geos_vel_compute(dens_press_filename,grid_filename="~/Downloads/ECCO_V4r4_PODAAC/ECCO_L4_GEOMETRY_LLC0090GRID_V4R4/GRID_GEOMETRY_ECCO_V4r4_native_llc0090.nc",fc_filename="llc_13tile_fc.txt"):
    """
    This routine computes geostrophic velocities from an input netCDF file containing ECCO v4r4 density and pressure anomalies on the native llc90 grid.
    
    Parameters
    ----------
    dens_press_filename: the name (including path if not in current directory) of the netCDF file containing the density and pressure anomalies
    
    grid_filename: the name (including path if not in current directory) of the netCDF file containing the latitudes ('YC') and grid parameters ('dxC','dyC') needed to compute horizontal derivatives
    
    fc_filename: text file containing a Python dictionary specifying the tile/face connections in the ECCO llc90 grid
    """
    
    import numpy as np
    import xarray as xr
    import json
    import xgcm
    from os.path import expanduser,join
    import sys
    user_home_dir = expanduser('~')
    sys.path.append(join(user_home_dir,'ECCOv4-py'))
    import ecco_v4_py as ecco
    
    # load file into workspace
    ds_denspress = xr.open_dataset(dens_press_filename, data_vars='minimal',\
                                     coords='minimal', compat='override')
    
    densanom = ds_denspress.RHOAnoma
    
    rhoConst = 1029.
    dens = rhoConst + densanom
    
    pressanom = ds_denspress.PHIHYDcR
    
    
    # load grid parameters file
    ds_grid = xr.open_dataset(grid_filename)
    
    
    # create xgcm Grid object
    xgcm_grid = ecco.get_llc_grid(ds_grid)
    
    # compute derivatives of pressure in x and y
    d_press_dx = (xgcm_grid.diff(rhoConst*pressanom,axis="X",boundary='extend'))/ds_grid.dxC
    d_press_dy = (xgcm_grid.diff(rhoConst*pressanom,axis="Y",boundary='extend'))/ds_grid.dyC
    
    # interpolate (vector) gradient values to center of grid cells
    press_grads_interp = xgcm_grid.interp_2d_vector({'X':d_press_dx,'Y':d_press_dy},boundary='extend')
    dp_dx = press_grads_interp['X']
    dp_dy = press_grads_interp['Y']
    
    # compute f from latitude of grid cell centers
    lat = ds_grid.YC
    Omega = (2*np.pi)/86164
    lat_rad = (np.pi/180)*lat    # convert latitude from degrees to radians
    f = 2*Omega*np.sin(lat_rad)
    
    # compute geostrophic velocities
    v_g = dp_dx/(f*dens)
    u_g = -dp_dy/(f*dens)
    
    # assign attributes to DataArrays (names) and units
    u_g.name = 'u_g'
    u_g.attrs.update({'long_name': 'Geostrophic velocity in model-x direction',\
                      'units': 'm s-1'})
    v_g.name = 'v_g'
    v_g.attrs.update({'long_name': 'Geostrophic velocity in model-y direction',\
                      'units': 'm s-1'})
    
    # create xarray Dataset containing geostrophic velocities
    ds_geos_vel = u_g.to_dataset(name='u_g',promote_attrs=True)
    ds_geos_vel['v_g'] = v_g
    
    return ds_geos_vel

#----------------------------------------------------------------------------------------

def depth_two_subplots(horiz_coords,depth_coords,data,k_split,cmap,mask=None,fig=None,axs=None):
    """
    Make 2 subplots with depth on y-axis, for shallow and deeper depths, given parameters:
    horiz_coords: horizontal coordinate, xarray DataArray
    depth_coords: depth_coordinate, xarray DataArray
    data: 2-D xarray DataArray
    k_split: k coordinate to split the plot at, integer
    cmap: string specifying colormap to use
    mask: 2-D boolean (land) mask, optional
    fig: figure object, optional, default is new figure is created
    axs: axes with two subplots oriented vertically, default is they are created
    """

    import numpy as np
    
    if (fig is None) and (axs is None):
        fig,axs = plt.subplots(2,1,figsize=(10,10))    # subplots for different depths
    elif (fig is None) or (axs is None):
        print("Warning: Only one of fig or axs has been supplied, not both")
    curr_ax = axs[0]
    curr_plot_0 = curr_ax.pcolormesh(horiz_coords,depth_coords[:k_split],\
                                      data.isel(k=np.arange(k_split)),cmap=cmap)                                      
    if mask is not None:
        curr_mask = mask.isel(k=np.arange(k_split))
        plot_mask(horiz_coords,depth_coords[:k_split],curr_mask,ax=curr_ax,color=np.zeros(3,))
    curr_ax.set_ylim(curr_ax.get_ylim()[::-1])
    curr_ax.set_ylabel('Depth [m]')    
    clim_0 = curr_plot_0.get_clim()
    
    curr_ax = axs[1]
    curr_plot_1 = curr_ax.pcolormesh(horiz_coords,depth_coords[k_split:],\
                                      data.isel(k=np.arange(25,len(depth_coords))),cmap=cmap)                                      
    if mask is not None:
        curr_mask = mask.isel(k=np.arange(25,len(depth_coords)))
        plot_mask(horiz_coords,depth_coords[k_split:],curr_mask,ax=curr_ax,color=np.zeros(3,))
    curr_ax.set_ylim(curr_ax.get_ylim()[::-1])
    curr_ax.set_ylabel('Depth [m]')
    clim_1 = curr_plot_1.get_clim()
    
    # create shared colorbar for 2 subplots
    new_clim = [np.fmin(clim_0[0],clim_1[0]),np.fmax(clim_0[1],clim_1[1])]
    curr_plot_0.set_clim(new_clim)
    curr_plot_1.set_clim(new_clim)
    fig.colorbar(curr_plot_1,ax=axs[:])
    return fig,axs

#----------------------------------------------------------------------------------------

# function to scale colormap, coordinating among multiple plots
def cmap_zerocent_scale_multiplots(plot_objs,scale_factor):
    """
    Center colormap at zero and scale relative to existing |maximum| value across multiple plots, 
    given plot objects (as list) and scale_factor, a number of type float.
    Returns new colormap limits as new_clim.
    """

    import numpy as np
    
    clim_plots = np.empty((len(plot_objs),2))
    for count,curr_obj in enumerate(plot_objs):
        clim_plots[count,:] = curr_obj.get_clim()
    new_clim = (scale_factor*np.nanmax(np.abs(clim_plots)))*np.array([-1,1])
    for curr_obj in plot_objs:
        curr_obj.set_clim(new_clim)
    return new_clim

#----------------------------------------------------------------------------------------

def mod_360_range(n,low_bound):
    "Compute n mod 360; output is in the range [n,n+360)"
    out = ((n - low_bound) % 360) + low_bound
    return out

#----------------------------------------------------------------------------------------

def llc_grid_idx_along_lat(lat_transect,lon_bnds):
    """
    Finds grid indices along a given line of latitude.
    Input parameters:
    - lat_transect: line of latitude to plot along, float
    - lon_bnds: 2 elements specifying western and eastern longitude bounds of transect, list or numpy array
                The difference between the bounds must be >0 and <=360
    Outputs:
    - idx_along_lat: grid indices along line of latitude, dict containing 'tile','j','i' as keys
    - XC_along_lat: longitude of grid cell centers, numpy array
    - XG_along_lat: longitude of western grid cell edges, numpy array
    """
    
    import numpy as np
    
    # identify indices of grid cells along latitude transect
    mask_along_lat = np.logical_and((ds_grid.YC_bnds > lat_transect).sum("nb") > 0,\
                                    (ds_grid.YC_bnds > lat_transect).sum("nb") < 4)
    mask_in_lon_bnds = np.logical_and((mod_360_range(ds_grid.XC,lon_bnds[1]) > mod_360_range(lon_bnds[0],lon_bnds[1])),\
                                      (mod_360_range(ds_grid.XC,lon_bnds[0]) < mod_360_range(lon_bnds[1],lon_bnds[0] + 1.e-5)))
    mask_along_lat = np.logical_and(mask_along_lat,mask_in_lon_bnds)
    idx_along_lat = mask_along_lat.values.nonzero()
    # create longitude arrays along transect
    XC_along_lat = np.zeros((len(idx_along_lat[0]),))
    XG_along_lat = np.zeros((len(idx_along_lat[0]),))
    # for loop through indices along transect, using zip to iterate through three spatial indices
    for count,(tile_idx,j_idx,i_idx) in enumerate(zip(idx_along_lat[0],idx_along_lat[1],idx_along_lat[2])):
        curr_XC = mod_360_range(ds_grid.XC[tile_idx,j_idx,i_idx],lon_bnds[0])
        XC_along_lat[count] = curr_XC
        XG_along_lat[count] = np.min(mod_360_range(ds_grid.XC_bnds[tile_idx,j_idx,i_idx,:],curr_XC - 180))
    # sort grid cells in order of increasing longitude
    idx_sorted = np.argsort(XC_along_lat)
    XC_along_lat = XC_along_lat[idx_sorted]
    XG_along_lat = XG_along_lat[idx_sorted]
    tile_idx_sorted = idx_along_lat[0][idx_sorted]
    j_idx_sorted = idx_along_lat[1][idx_sorted]
    i_idx_sorted = idx_along_lat[2][idx_sorted]
    idx_along_lat = {'tile':tile_idx_sorted,'j':j_idx_sorted,'i':i_idx_sorted}
    return idx_along_lat,XC_along_lat,XG_along_lat

#----------------------------------------------------------------------------------------

def data_along_lat(data_in,idx_along_lat,XC_along_lat,XG_along_lat):
    """
    Creates xarray DataArray along latitude transect, given grid indices
    Input parameters:
    - data_in: data to plot along latitude, xarray DataArray or numpy array
            must have 4 or 5 dimensions, 'time','k','tile','j','i', 'time' is optional
    - XC_along_lat: longitudes of grid cell centers, numpy array
    - XG_along_lat: longitudes of western grid cell edges, numpy array
    Outputs:
    - data_xrarray: data along latitude transect, xarray DataArray
    """
    
    import numpy as np
    
    # retrieve grid indices along transect
    tile_idx_along = idx_along_lat['tile']
    j_idx_along = idx_along_lat['j']
    i_idx_along = idx_along_lat['i']
    
    # construct DataArray along transect
    data_along_lat = np.empty(data_in.shape[:-3] + (len(idx_along_lat['i']),))
    data_along_lat.fill(np.nan)
    for count,(tile_idx,j_idx,i_idx) in enumerate(zip(tile_idx_along,j_idx_along,i_idx_along)):
        if len(data_in.shape) == 4:
            data_along_lat[:,count] = data_in[:,tile_idx,j_idx,i_idx]
        elif len(data_in.shape) == 5:
            data_along_lat[:,:,count] = data_in[:,:,tile_idx,j_idx,i_idx]
    if len(data_in.shape) == 4:    # no time dimension
        data_xrarray = xr.DataArray(
                    data=data_along_lat,
                    dims=["k","lon"],
                    coords=dict(
                        Z=(["k"],ds_grid.Z.data),
                        lon=(["lon"],XC_along_lat),
                        lonW=(["lon"],XG_along_lat),
                        lat=lat_transect,
                    ),
        )
    elif len(data_in.shape) == 5:    # include time dimension
        try:
            time_coord = data_in.time.data
        except:
            time_coord = np.nan
        data_xrarray = xr.DataArray(
                    data=data_along_lat,
                    dims=["time","k","lon"],
                    coords=dict(
                        time=(["time"],time_coord),
                        Z=(["k"],ds_grid.Z.data),
                        lon=(["lon"],XC_along_lat),
                        lonW=(["lon"],XG_along_lat),
                        lat=lat_transect,
                    ),
        )
    
    return data_xrarray    
    
#----------------------------------------------------------------------------------------    

def lon_depth_along_lat(lat_transect,lon_bnds,data_in):
    """
    Function to extract data along line of latitude.
    Input parameters:
    - lat_transect: line of latitude to plot along, float
    - lon_bnds: 2 elements specifying western and eastern longitude bounds of transect, list or numpy array
    - data_in: data, numpy array or xarray DataArray
    Outputs:
    - XC_along_lat: longitude of grid cell centers, numpy array
    - XG_along_lat: longitude of western grid cell edges, numpy array
    - data_xrarray: data along latitude transect, xarray DataArray
    """
    
    # identify indices of grid cells along latitude transect
    idx_along_lat,XC_along_lat,XG_along_lat = llc_grid_idx_along_lat(lat_transect,lon_bnds)
    
    # construct DataArray along transect
    data_xrarray = data_along_lat(data_in,idx_along_lat,XC_along_lat,XG_along_lat)
    
    return XC_along_lat,XG_along_lat,data_xrarray

#----------------------------------------------------------------------------------------    

def plot_mask_ecco_tiles(mask,color):
    """
    Plot mask in global ECCO tiles plot on current axes, 
    given 2-D mask (xarray DataArray) and color (a string, RGB tuple or 3-element NumPy array).
    """
    
    import numpy as np
        
    # loop through tiles to add mask
    tile_order = np.array([-1,-1,-1,6, \
                             2,5,7,10,  \
                             1,4,8,11, \
                             0,3,9,12])
    for idx,curr_ax in enumerate(curr_fig.get_axes()):
        if len(curr_ax.get_images()) > 0:
            # plot land mask
            array_plot = mask.isel(tile=tile_order[idx]).squeeze()
            if tile_order[idx] == 6:
                array_plot = np.rot90(array_plot,2)
            elif tile_order[idx] > 6:
                array_plot = np.rot90(array_plot)
            plot_mask(array_plot,ax=curr_ax,color=color)