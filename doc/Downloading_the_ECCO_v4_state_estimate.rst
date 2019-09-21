#####################################
How to get the ECCO v4 State Estimate
#####################################

.. _in-ftp-site:

ECCO v4 r3 on PO.DAAC drive and UT Austin FTP mirror
====================================================

ECCOv4 output fields are provided as NetCDF files.  We have been working hard at improving  these NetCDF files so that they contain more useful metadata and variable descriptions, fewer extraneous fields, and have more consistent naming.

The Python examples in this tutorial are compatible with the ECCOv4 NetCDF grid files provided in the 'Release3_alt' directories:

The ECCO v4 state estimate is now available in two places.

1. PO.DAAC drive: https://ecco.jpl.nasa.gov/drive/files/Version4/ (po.daac drive, recommended)

2. UT Austin FTP mirror: https://web.corral.tacc.utexas.edu/OceanProjects/ECCO/ECCOv4/


IMPORTANT: The NetCDF files in the ECCO Version 4 'Release 1', 'Release 2', and 'Release3' directories are not fully compatible with the tutorial examples and the ecco-v4-py library.  We may update these earlier releases to be compatible with the ecco-v4-py library in the future.  As of now, the 'Release3_alt' has the most up-to-date ECCO solution and file format.  These files can be found here:

1. PO.DAAC drive: https://ecco.jpl.nasa.gov/drive/files/Version4/Release3_alt (po.daac drive, recommended)
2. UT Austin FTP mirror: https://web.corral.tacc.utexas.edu/OceanProjects/ECCO/ECCOv4/Release3_alt

Please see the ECCO website, ecco.jpl.nasa.gov, for udpates.

.. _in-grid:

fields on the 13-tile lat-lon-cap (llc) *native* model grid
===========================================================

geometric model grid parameters
-------------------------------

Calculations involving the state estimate variables often require the geometric model grid parameters.  These parameters are packaged together as 13 grid NetCDF files, one for each llc tile, in the *nctiles_grid/* directory.

.. _in-monthly:

monthly-averaged ocean and sea-ice variables
--------------------------------------------

Monthly-averaged ocean and sea-ice fields are in provided in subdirectories of *nctiles_monthly/*. Each subdirectory corresponds to a single variable and contains 13 NetCDF files, one for each different llc tile.

.. _in-daily:

daily-averaged ocean and sea-ice variables
------------------------------------------

Daily-averaged ocean and sea-ice fields are in the subdirectories of *nctiles_daily/*. Each subdirectory corresponds to a single variable and contains 13 NetCDF files, one for each llc tile.


6-hourly atmosphere variables
-----------------------------

6-hourly atmospheric fields that can be used as atmospheric boundary conditions for the model are provided in *input_forcing/*. Each atmospheric state variable is divided by year.  These files are *not* divided into 13 tiles but are instead provided in the special *native* lat-lon-cap flat binary format required by the model.  Tools for reading and plotting files in llc binary format are provided in this tutorial.

For reference, these are the atmospheric state fields provided in *input_forcing*

.. code-block:: bash

  eccov4r3_dlw_YYYY                 downward longwave radiation
  eccov4r3_dsw_YYYY                 downward shortwave radiation
  eccov4r3_rain_YYYY                precipitation
  eccov4r3_spfh2m_YYYY              2m near-surface atmospheric specific humidity
  eccov4r3_tmp2m_degC_YYYY          2m near-surface air temperature
  eccov4r3_ustr_YYYY                zonal wind stress
  eccov4r3_vstr_YYYY                meridional wind stress
  eccov4r3_wspeed_YYYY              near-surface wind speed


monthly-averaged *interpolated* 0.5° x 0.5° latitude-longitude variables
========================================================================

Monthly-averaged ocean, sea ice, and air-sea flux terms are in the subdirectory *interp_monthly/*. Each NetCDF file in *interp_monthly* corresponds to a single variable.


Downloading the State Estimate
==============================

ECCO v4 solutions are now hosted on the PO.DAAC drive.  This service is very useful because one can mount the ECCO file directory on PO.DAAC drive to your local machine.  

You can even use wget to download files through PO.DAAC drive.  See the ECCO website for more details:
https://ecco.jpl.nasa.gov/products/latest/

Take note of the location of your files.  You'll need to specify their path location to load them in the tutorial.
