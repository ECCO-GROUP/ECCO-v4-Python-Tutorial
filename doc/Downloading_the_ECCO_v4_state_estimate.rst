######################################
Downloading the ECCO v4 State Estimate
######################################


.. _in-ftp-site:
The ECCO v4 r3 ftp site
=======================

The ECCO v4 r3 state estimate is hosted on ftp://ecco.jpl.nasa.gov/Version4/Release3

The directory layout of the ftp site is described in the following document,
ftp://ecco.jpl.nasa.gov/Version4/Release3/doc/v4r3_overview.pdf

.. _in-grid:
lat-lon-cap 90 fields
=====================

geometric model grid parameters
-------------------------------

Calculations involving the state estimate variables often require the geometric model grid parameters.  These parameters are packaged together as 13 grid NetCDF files, one for each llc tile, in the *nctiles_grid/* directory.

.. _in-monthly:
monthly-averaged ocean and sea-ice variables
--------------------------------------------
Monthly-averaged ocean and sea-ice fields are in provided in subdirectories of *nctiles_monthly/*. Each subdirectory corresponds to a single variable and contains 13 NetCDF files, one for each different llc tile.

.. _in-daily:
daily-averaged ocean and sea-ice variables
--------------------------------------------
Daily-averaged ocean and sea-ice fields are in the subdirectories of *nctiles_daily/*. Each subdirectory corresponds to a single variable and contains 13 NetCDF files, one for each llc tile.


6-hourly atmosphere variables
-----------------------------
6-hourly atmospheric fields that can be used as atmospheric boundary conditions for the model are provided in *input_forcing/*. Each atmospheric state variable is divided by year.  These files are *not* divided into 13 tiles but are instead provided in the special *native* lat-lon-cap flat binary format required by the model.  Tools for reading and plotting files in llc binary format are provided in this tutorial.

For reference, these are the atmospheric state fields provided in *input_forcing*

::
  eccov4r3_dlw_YYYY                 downward longwave radiation
  eccov4r3_dsw_YYYY                 downward shortwave radiation
  eccov4r3_rain_YYYY                precipitation
  eccov4r3_spfh2m_YYYY              2m near-surface atmospheric specific humidity
  eccov4r3_tmp2m_degC_YYYY          2m near-surface air temperature
  eccov4r3_ustr_YYYY                zonal wind stress
  eccov4r3_vstr_YYYY                meridional wind stress
  eccov4r3_wspeed_YYYY              near-surface wind speed


monthly-averaged *interpolated* 1° x 1° latitude-longitude variables
====================================================================

Monthly-averaged ocean, sea ice, and air-sea flux terms are in the subdirectory *interp_monthly/*. Each NetCDF file in *interp_monthly* corresponds to a single variable.


Downloading the State Estimate
==============================

There are many ways of downloading files from the ecco v4 r3 ftp site.   

On windows machine, I typically use an ftp client.

On linux machines, I use *wget*.  To download the model grid, and monthly-averaged sea surface height, ocean bottom pressure, temperature, and salinity fields use the following:

.. code-block:: bash

    wget -r ftp://ecco.jpl.nasa.gov/Version4/Release3/nctiles_grid/
    wget -r ftp://ecco.jpl.nasa.gov/Version4/Release3/nctiles_monthly/SSH/
    wget -r ftp://ecco.jpl.nasa.gov/Version4/Release3/nctiles_monthly/OBP/
    wget -r ftp://ecco.jpl.nasa.gov/Version4/Release3/nctiles_monthly/THETA/
    wget -r ftp://ecco.jpl.nasa.gov/Version4/Release3/nctiles_monthly/SALT/


On osx machines, use the *curl* command:

.. code-block:: bash

    curl ftp://ecco.jpl.nasa.gov/Version4/Release3/nctiles_grid/
    curl ftp://ecco.jpl.nasa.gov/Version4/Release3/nctiles_monthly/SSH/
    curl ftp://ecco.jpl.nasa.gov/Version4/Release3/nctiles_monthly/OBP/
    curl ftp://ecco.jpl.nasa.gov/Version4/Release3/nctiles_monthly/THETA/
    curl ftp://ecco.jpl.nasa.gov/Version4/Release3/nctiles_monthly/SALT/


Take note of the location of your files.  You'll need to specify their path location to load them in the tutorial.
