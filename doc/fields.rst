############################################################
ECCO v4 state estimate ocean, sea-ice, and atmosphere fields
############################################################

The complete state estimate consists of a set of ocean, sea-ice, air-sea flux, and atmosphere state variables that are the output from a free-running ocean and sea-ice general circulation model. 

.. _in-layout:

*******************
Geographical layout
*******************

Ocean, sea-ice, air-sea flux, and atmosphere fields are provided in two spatial layouts:

- 13-tile *native* lat-lon-cap 90 (llc90) grid
- 0.5° x 0.5° latitude and longitude grid

13-tile *native* lat-lon-cap 90 grid
====================================

The lat-lon-cap (llc) is the decomposition of the spherical Earth into a Cartesian curvilinear coordinate system .  It is a topologically non-trivial cubed-sphere rendering in the northern hemisphere and a dipolar grid in the southern hemisphere.  Between 70°S and ~57°N, model grid cells are approximately oriented to lines of latitude and longitude.  A special Arctic "cap" is situated north of ~57°N.  

The Cartesian curvilinear coordinate system is divided into 13 tiles, each consisting of 90x90 grid cells in the horizontal and 50 vertical levels.  Horizontal model grid resolution varies spatially from 22km to 110km, with the highest resolutions at high latitudes and lowest resolution in mid latitudes. Vertical grid spacing increases with depth from 10m to 456.5m.  The bottom of the deepest model grid cell is 6145m below the surface.

The Cartesian (x,y) coordinates of llc tiles do not coorespond to longitude and latitude.  Horizontal velocities are defined relative to the **local orientation** of x and y in the tile.  Velocities in the positive *x* direction are defined as positive *u*.  Velocities in the positive *y* direction are defined as positive *v*.

.. figure:: ../figures/llc90_0.png
    :align: center
    :alt: alternate text
    :figclass: align-center


Available fields on the llc90 grid
----------------------------------

`Monthly-averaged fields <https://raw.githubusercontent.com/ECCO-GROUP/ECCO-v4-Python-Tutorial/master/varlist/v4r4_nctiles_monthly_varlist.txt>`_

`Daily-averaged fields <https://raw.githubusercontent.com/ECCO-GROUP/ECCO-v4-Python-Tutorial/master/varlist/v4r4_nctiles_daily_varlist.txt>`_

`Daily snapshot fields <https://raw.githubusercontent.com/ECCO-GROUP/ECCO-v4-Python-Tutorial/master/varlist/v4r4_nctiles_snapshots_varlist.txt>`_


*interpolated* 0.5° x 0.5° latitude-longitude grid
==================================================

Many fields from the *native* lat-lon-cap model output have been interpolated to a more user-friendly 0.5° latitude-longitude grid.  Note that the interpolated fields generally can not be used to close budgets, and so flux fields that are primarily used for budget calculations have not been interpolated.

Available fields on the 0.5° x 0.5° latitude-longitude grid
-----------------------------------------------------------

`0.5° x 0.5° monthly-averaged fields <https://raw.githubusercontent.com/ECCO-GROUP/ECCO-v4-Python-Tutorial/master/varlist/v4r4_latlon_monthly_varlist.txt>`_

`0.5° x 0.5° daily-averaged fields <https://raw.githubusercontent.com/ECCO-GROUP/ECCO-v4-Python-Tutorial/master/varlist/v4r4_latlon_daily_varlist.txt>`_


Miscellaneous fields and data
==================================================

A few time series that do not have spatial dimensions (e.g. averages over the global ocean), as well as the grid parameter fields, are listed `here <https://raw.githubusercontent.com/ECCO-GROUP/ECCO-v4-Python-Tutorial/master/varlist/v4r4_tseries_grid_varlist.txt>`_.

See the "Using Python to Download ECCO Datasets" tutorial for information on how to download the output.


*******************************************
Temporal frequency of state estimate fields
*******************************************

All ocean, sea-ice, surface atmosphere, and air-sea flux fields are provided as monthly and daily averages on the native LLC90 grid.  In addition, certain fields are provided as daily snapshots (at 0Z UTC time) to support budget closure calculations.

*************
Custom output
*************

Because the state estimate fields are the output from a free-running ocean model, users can re-run the model to generate custom output on the native lat-lon-cap model grid.  Instructions for doing so are provided `here <https://www.ecco-group.org/docs/v4r4_reproduction_howto.pdf>`_ .


