**********************************
Getting the ECCO v4 State Estimate
**********************************

ECCO Version 4 Release 3 is the most recent edition of the
global ocean state estimate and estimation system described by Forget et al. (2015b, 2016).  

A brief synopsis describing Release 3 can be found here:  ftp://ecco.jpl.nasa.gov/Version4/Release3/doc/v4r3_summary.pdf

A high-level analysis of the state estimate can be found here:
ftp://ecco.jpl.nasa.gov/Version4/Release3/doc/v4r3_depiction.pdf


.. _in-ftp-site:
The ECCO v4 r3 ftp site
=======================

The ECCO v4 r3 state estimate is hosted on ftp://ecco.jpl.nasa.gov/Version4/Release3

The directory layout of the ftp site is described in the following document,
ftp://ecco.jpl.nasa.gov/Version4/Release3/doc/v4r3_overview.pdf

.. _in-layout:
The geographical layout of the fields
=====================================

Each variable of the state estimate and the model grid parameters are provided as a set of 13 "tiles".  These 13 tiles map the spherical Earth into a Cartesian curvilinear coordinate system - a requirement of the model.  Each tile is comprised of 90x90 horizontal grid cells.  The model has 50 vertical levels.  Between 70°S and ~57°N and the model grid is approximately latitude-longitude in geometry.  A special Arctic "Cap" is found north of ~57°N (Forget et al., 2015). Together, the 13 tiles of 90x90 horizontal grid cells in the latitude-longitude geometry plus the Arctic "Cap" give us the so-called "LLC90" grid.

The horizontal resolution of model grid cells varies spatially from 22km to 110km, with the highest resolution in high latitudes and lowest resolution in mid latitudes. The deepest ocean bottom is set to 6145m below the surface, with the vertical grid spacing increasing 


.. _in-grid:
Model Grid Fields
=================
Analysis of the state estimate requires knowledge of the geometric parameters that describe the model's LLC90 (Lat-Lon-Cap 90) grid.  
These parameters are found in the *nctiles_grid/* directory.

.. _in-monthly:
Monthly-Averaged Model fields
=============================
The nominal output of the state estimate is in the form of monthly-averaged fields (found in the subdirectory *nctiles_monthly/*). Each subdirectory inside *nctiles_monthly* contains netCDF files for a particular variable, as indicated by the name of the subdirectory, split into 13 tile files as described above.  Some of the most commonly used fields, like velocity components, potential temperature, salinity, sea surface height, and ocean bottom pressure are UVEL, VVEL, THETA, SALT, SSH, and OBP, respectively

.. _in-daily:
Daily-Averaged Model fields
===========================
Daily averages are also provided for the following variables in the directory nctiles_daily to facilitate studies of the ocean’s high-frequency variations: SSH, OBP, sea surface temperature (SST), sea surface salinity (SSS), sea-ice concentration (SIarea), mean sea-ice thickness (SIheff), mean snow thickness (SIhsnow), and sea-ice and snow loading [kg/m^2] (sIceLoad).
surface to 457m near the ocean bottom.


Downloading the State Estimate
==============================

There are many ways of downloading files from an ftp site.   The simplest is with the `wget command`.  To download the model grid, and  monthly-averaged sea surface height, ocean bottom pressure, temperature, and salinity fields use the following commands.

.. code-block:: bash

    wget -r ftp://ecco.jpl.nasa.gov/Version4/Release3/nctiles_grid/
    wget -r ftp://ecco.jpl.nasa.gov/Version4/Release3/nctiles_monthly/SSH/
    wget -r ftp://ecco.jpl.nasa.gov/Version4/Release3/nctiles_monthly/OBP/
    wget -r ftp://ecco.jpl.nasa.gov/Version4/Release3/nctiles_monthly/THETA/
    wget -r ftp://ecco.jpl.nasa.gov/Version4/Release3/nctiles_monthly/SALT/

Take note of the location of your files.  You'll need to specify their path location to load them in the tutorial.
