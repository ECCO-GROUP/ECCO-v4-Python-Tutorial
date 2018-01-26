
An Improved Method for Loading ECCOv4 netCDF files
==================================================

Objectives:
-----------

To introduce a method for loading data from the ECCO v4 netCDF tile
files that returns ``Dataset`` and ``DataArray`` objects that properly
distinguish *where* on the Arakawa-C grid the variables are situated.

This custom routine written for the *ecco\_v4\_py* package is:
``load_tile_from_netcdf``.

Introduction
------------

As we showed in the first tutorial, we can use the
``xr.open_dataset(data_dir + fname)`` routine from xarray to load a
netCDF tile file into Python. This routine parses the netCDF file and
extracts the dimensions, coordinates, variables, and metadata
information and uses those to construct the ``Dataset`` object.

In the last tutorial we read in a single GRID netCDF tile file and
examined its contents. We found that it listed three dimensions, i1, i2,
and i3, for its Data variables. Let's load it up again and take a closer
look. This time we'll name the ``Dataset`` object as ``grid_default``
since we are loading it with the default method from ``xarray``.

.. code:: ipython2

    import matplotlib.pylab as plt
    import numpy as np
    import sys
    import xarray as xr
    from copy import deepcopy 
    sys.path.append('/Users/ifenty/git_repo/ECCOv4-py')
    import ecco_v4_py as ecco

.. code:: ipython2

    data_dir='/Volumes/ECCO_BASE/ECCO_v4r3/nctiles_grid/'    
    fname = 'GRID.0003.nc'
    grid_default = xr.open_dataset(data_dir + fname)

.. code:: ipython2

    grid_default




.. parsed-literal::

    <xarray.Dataset>
    Dimensions:  (i1: 50, i2: 90, i3: 90)
    Coordinates:
      * i1       (i1) float64 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 ...
      * i2       (i2) float64 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 ...
      * i3       (i3) float64 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 ...
    Data variables:
        hFacC    (i1, i2, i3) float64 ...
        hFacW    (i1, i2, i3) float64 ...
        hFacS    (i1, i2, i3) float64 ...
        XC       (i2, i3) float64 ...
        YC       (i2, i3) float64 ...
        XG       (i2, i3) float64 ...
        YG       (i2, i3) float64 ...
        RAC      (i2, i3) float64 ...
        RAZ      (i2, i3) float64 ...
        DXC      (i2, i3) float64 ...
        DYC      (i2, i3) float64 ...
        DXG      (i2, i3) float64 ...
        DYG      (i2, i3) float64 ...
        Depth    (i2, i3) float64 ...
        AngleCS  (i2, i3) float64 ...
        AngleSN  (i2, i3) float64 ...
        RC       (i1) float64 ...
        RF       (i1) float64 ...
        DRC      (i1) float64 ...
        DRF      (i1) float64 ...
    Attributes:
        description:    C-grid parameters (see MITgcm documentation for details)....
        A:              :Format      = native grid (nctiles w. 13 tiles)
        B:              :source      = ECCO consortium (http://ecco-group.org/)
        C:              :institution = JPL/UT/MIT/AER
        D:              :history     = files revision history :
        E:                                 04/20/2017: fill in geometry info for ...
        F:                                 11/06/2016: third release of ECCO v4 (...
        G:                             estimates revision history (from second re...
        H:                                 employs bi-harmonic viscosity (enhance...
        I:                                 sea-ice parameters, updated or novel o...
        J:                                 GRACE OBP, Aquarius SSS, global mean s...
        K:                                 time-series, extended and/or expanded ...
        L:                                 revised weights including data and con...
        M:                                 to account for grid-size variation and...
        N:                                 separate time-mean and time-variable d...
        O:                                 and controls, sea-ice costs, and initi...
        P:                                 additional controls.\n 
        Q:              :references  = Forget, G., J.-M. Campin, P. Heimbach, C. ...
        R:                              and C. Wunsch, 2015: ECCO version 4: an i...
        S:                              non-linear inverse modeling and global oc...
        T:                              Geoscientific Model Development, 8, 3071-...
        U:                             Forget, G., J.-M. Campin, P. Heimbach, C. ...
        V:                              ECCO version 4: Second Release, 2016, htt...
        W:              file created using gcmfaces_IO/write2nctiles.m
        date:           21-Apr-2017
        Conventions:    CF-1.6
        _FillValue:     nan
        missing_value:  nan



We see that all of the Data variables in ``grid_default`` use one of
three dimensions, **i1**,\ **i2**,\ **i3**. As we saw before, some
variables are 3D (e.g., hFacC), others are 2D (e.g., XC), and others are
1D (e.g., RF).

Now, while the default format of this Dataset object is already quite
useful, it falls short of taking full advantage of the 'coordinate'
feature afforded by the Dataset object.

Model variables on Arakawa-C grids are staggered in space. On the
horizontal plane, model variables can be situated at one of four
different classes of point.

figure

The four horizontal points of the Arakawa-C grid
------------------------------------------------

*c* points
~~~~~~~~~~

Variables that do not have a horizontal velocity component (e.g., T, S,
SSH, OBP, sea ice concentration, vertical velocity) are situated at
:math:`c` points in the horizontal plane. :math:`c` points are at the
center of the *tracer* grid cell in the horizontal plane.

Let us define the :math:`(i,j)` coordinate system for the indices of
:math:`c` points.

:math:`c(0,0)` is the -x most and -y most tracer grid cell.

Moving in the +\ :math:`y` direction, the next grid cell is
:math:`c(0,1)` Moving in the +\ :math:`x` direction, the next grid cell
is :math:`c(1,0)`

*u* points
~~~~~~~~~~

Variables that are explictly related to horizontal velocity or
horizontal fluxes in the model's :math:`x` direction are situated at
:math:`u` points in the horizontal plane. Examples include horizontal
velocity in the :math:`x` direction (:math:`UVEL`) and horizontal
advective flux of snow in the :math:`x` direction (:math:`ADVxSNOW`).

In the :math:`x` direction they are situated on the edges or faces of
the tracer grid cell. In the :math:`y` direction they are at the same
locations as the :math:`c` points.

Let us define the :math:`(i_g, j)` coordinate system for :math:`u`
points. We use :math:`i_g` as the coordinate in the :math:`x` direction
because :math:`u` points are situated along the tracer grid cell
ed\ ***G***\ es. We use the :math:`j` for its :math:`y` coordinate
because it is the same as the :math:`y` coordinate of the :math:`c`
points.

:math:`u(0,0)` is the first :math:`u` point situated in the :math:`-x`
direction of $c(0,0).

Moving in the +\ :math:`y` direction, the next grid cell is
:math:`u(0,1)` Moving in the +\ :math:`x` direction, the next grid cell
is :math:`u(1,0)`

*v* points
~~~~~~~~~~

Variables that are explictly related to horizontal velocity or
horizontal fluxes in the model's :math:`y` direction are situated at
:math:`v` points in the horizontal plane. Examples include horizontal
velocity in the :math:`y` direction (:math:`VVEL`) and horizontal
advective flux of snow in the :math:`y` direction (:math:`ADVySNOW`).

In the :math:`x` direction they are at the same locations as the
:math:`c` points. In the :math:`y` direction they are situated on the
edges (or faces) of the tracer grid cell.

Let us define the :math:`(i, j_g)` coordinate system for :math:`v`
points. We use the :math:`i` for its :math:`x` coordinate because it is
the same as the :math:`x` coordinate of the :math:`c` points. We use
:math:`j_g` as the coordinate in the :math:`y` direction because
:math:`v` points are situated along the tracer grid cell
ed\ ***G***\ es.

:math:`v(0,0)` is the first :math:`v` point situated in the :math:`-y`
direction of $c(0,0).

Moving in the +\ :math:`y` direction, the next grid cell is
:math:`v(0,1)` Moving in the +\ :math:`x` direction, the next grid cell
is :math:`v(1,0)`

*g* points
~~~~~~~~~~

Variables that are explictly related to horizontal velocities in the
model in both the :math:`x` and :math:`y` direction are situated at
:math:`g` points in the horizontal plane. Vorticity is an example.

:math:`g` points are situated along the edges of the grid cell in both
:math:`x` and :math:`y`. In other words, they are at the **corners** of
tracer grid cells.

Let us define the :math:`(i_g, j_g)` coordinate system for :math:`g`
points following the same reasoning as described above: in both the
:math:`x` and :math:`y` directions, :math:`g` points are on the
ed\ ***G***\ es of tracer grid cells.

:math:`g(0,0)` is the first :math:`g` point situated in the :math:`-y`
and :math:`-x` directions of $c(0,0).

Moving in the +\ :math:`y` direction, the next grid cell is
:math:`g(0,1)` Moving in the +\ :math:`x` direction, the next grid cell
is :math:`g(1,0)`

The two vertical points of the Arakawa-C gird
---------------------------------------------

There are two coordinates in the vertical :math:`z` dimension:

*w* points
~~~~~~~~~~

Variables related to vertical velocity or vertical fluxes are situated
at :math:`w` in the vertical direction. These variables are situated on
the upper and lower faces of the tracer grid cell.

Let us define the :math:`k_g` coordinate system for :math:`w` points by
following the same reasoning as we used above: :math:`w` points fall
along the the ed\ ***G***\ es of tracer grid cells in the :math:`z`
direction.

Indexing begins at the sea surface, k\_g=0.

*k* points
~~~~~~~~~~

Let us define the :math:`k` coordinate system for variables situated in
the middle of a tracer grid cell in the vertical :math:`z` direction.
Examples include all tracers.

Indexing begins in the uppermost grid cell surface, k=0.

| The default coordinate names in the GRID netcdf tile files are not
  sufficiently descriptive to distinguish between the four horizontal
  coordinates ('i','i\_g','j','j\_g') and the three vertical
  coordinates.
| Therefore, we provide a routine especially for reading in ECCOv4 llc90
  GRID netcdf files and giving the default i1, i2, and i3 coordinates
  more descriptive names.

Load the file containing the grid parameter information for one tile.
---------------------------------------------------------------------

To load ECCO v4's netcdf files we will use the *open\_dataset* command
from the `xarray <http://xarray.pydata.org/en/stable/index.html>`__
Python package. *open\_dataset* creates a **Dataset** object and loads
the contents of the netcdf file, including its metadata, into a data
structure.

Let's open the grid file for *tile 3* (North East Atlantic Ocean), of
the 13 ECCO v4 llc90 grid files.

Change ``data_dir`` to match the location of your ``nctiles_grid``
directory.

.. code:: ipython2

    data_dir='/Volumes/ECCO_BASE/ECCO_v4r3/nctiles_grid/'    
    var = 'GRID'
    var_type = 'grid'
    tile_index = 3
    grid_tile_3 = ecco.load_tile_from_netcdf(data_dir, 
                                             var, 
                                             var_type, 
                                             tile_index)


.. parsed-literal::

    loading /Volumes/ECCO_BASE/ECCO_v4r3/nctiles_grid/GRID.0003.nc


.. code:: ipython2

    grid_tile_3




.. parsed-literal::

    <xarray.Dataset>
    Dimensions:  (i: 90, i_g: 90, j: 90, j_g: 90, k: 50, k_l: 50, k_u: 50)
    Coordinates:
        tile     int64 3
      * k        (k) float64 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 ...
      * i        (i) float64 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 ...
      * j        (j) float64 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 ...
      * i_g      (i_g) float64 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 ...
      * j_g      (j_g) float64 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 ...
      * k_u      (k_u) float64 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 ...
      * k_l      (k_l) int64 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 ...
    Data variables:
        XC       (j, i) float64 ...
        YC       (j, i) float64 ...
        RAC      (j, i) float64 ...
        Depth    (j, i) float64 ...
        AngleCS  (j, i) float64 ...
        AngleSN  (j, i) float64 ...
        hFacC    (k, j, i) float64 ...
        land_c   (k, j, i) float64 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 ...
        XG       (j_g, i_g) float64 ...
        YG       (j_g, i_g) float64 ...
        RAZ      (j_g, i_g) float64 ...
        DXC      (j, i_g) float64 ...
        DYG      (j, i_g) float64 ...
        hFacW    (k, j, i_g) float64 ...
        land_u   (k, j, i_g) float64 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 ...
        DYC      (j_g, i) float64 ...
        DXG      (j_g, i) float64 ...
        hFacS    (k, j_g, i) float64 ...
        land_v   (k, j_g, i) float64 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 ...
        RF       (k_u) float64 ...
        DRC      (k_u) float64 ...
        RC       (k) float64 ...
        DRF      (k) float64 ...
        RB       (k_l) float64 -10.0 -20.0 -30.0 -40.0 -50.0 -60.0 -70.0 -80.01 ...
    Attributes:
        description:    C-grid parameters (see MITgcm documentation for details)....
        A:              :Format      = native grid (nctiles w. 13 tiles)
        B:              :source      = ECCO consortium (http://ecco-group.org/)
        C:              :institution = JPL/UT/MIT/AER
        D:              :history     = files revision history :
        E:                                 04/20/2017: fill in geometry info for ...
        F:                                 11/06/2016: third release of ECCO v4 (...
        G:                             estimates revision history (from second re...
        H:                                 employs bi-harmonic viscosity (enhance...
        I:                                 sea-ice parameters, updated or novel o...
        J:                                 GRACE OBP, Aquarius SSS, global mean s...
        K:                                 time-series, extended and/or expanded ...
        L:                                 revised weights including data and con...
        M:                                 to account for grid-size variation and...
        N:                                 separate time-mean and time-variable d...
        O:                                 and controls, sea-ice costs, and initi...
        P:                                 additional controls.\n 
        Q:              :references  = Forget, G., J.-M. Campin, P. Heimbach, C. ...
        R:                              and C. Wunsch, 2015: ECCO version 4: an i...
        S:                              non-linear inverse modeling and global oc...
        T:                              Geoscientific Model Development, 8, 3071-...
        U:                             Forget, G., J.-M. Campin, P. Heimbach, C. ...
        V:                              ECCO version 4: Second Release, 2016, htt...
        W:              file created using gcmfaces_IO/write2nctiles.m
        date:           21-Apr-2017
        Conventions:    CF-1.6
        _FillValue:     nan
        missing_value:  nan


