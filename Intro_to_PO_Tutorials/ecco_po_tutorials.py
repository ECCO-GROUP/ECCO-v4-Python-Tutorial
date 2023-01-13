def cmap_zerocent_scale(plot,scale_factor):
    """
    Center colormap at zero and scale relative to existing |maximum| value, 
    given plot object and scale_factor, a number of type float.
    Returns new colormap limits as new_clim.
    """
    curr_clim = plot.get_clim()
    new_clim = (scale_factor*np.max(np.abs(curr_clim)))*np.array([-1,1])
    plot.set_clim(new_clim)
    return new_clim



def plot_mask(*args,ax=None,color):
    """
    Plot mask, given input parameters:
    - X, Y: (optional) coordinates as 1-D or 2-D NumPy arrays or xarray DataArrays
    - mask: 2-D array of boolean values (True/False or 1/0), NumPy or xarray
    - axes: axes to plot on, defaults to current axes
    - color: a string indicating a color in Matplotlib, or a 3-element tuple or NumPy array indicating RGB color values
    
    Returns plot_obj, the plot object of the mask 
    """
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
