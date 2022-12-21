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
    
    # load file into workspace
    ds_denspress = xr.open_dataset(dens_press_filename, data_vars='minimal',\
                                     coords='minimal', compat='override')
    
    densanom = ds_denspress.RHOAnoma
    
    rhoConst = 1029.
    dens = rhoConst + densanom
    
    pressanom = ds_denspress.PHIHYDcR
    
    
    # load grid parameters file
    ds_grid = xr.open_dataset(grid_filename)
    
    
    # load face_connections dictionary
    
    # read in as string
    with open(fc_filename) as fc:
        data = fc.read()
    # convert string to dictionary
    face_connections = json.loads(data)
    
    
    # create xgcm Grid object
    grid = xgcm.Grid(ds_grid,periodic=False,face_connections=face_connections)
    
    # compute derivatives of pressure in x and y
    d_press_dx = (grid.diff(rhoConst*pressanom,axis="X",boundary='extend'))/ds_grid.dxC
    d_press_dy = (grid.diff(rhoConst*pressanom,axis="Y",boundary='extend'))/ds_grid.dyC
    
    # interpolate (vector) gradient values to center of grid cells
    press_grads_interp = grid.interp_2d_vector({'X':d_press_dx,'Y':d_press_dy},boundary='extend')
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
