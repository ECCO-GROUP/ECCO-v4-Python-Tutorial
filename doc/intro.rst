#########################################
The ECCO Ocean and Sea-Ice State Estimate
#########################################

***************************************************
What is the ECCO Central Production State Estimate?
***************************************************

The Estimating the Circulation and Climate of the Ocean (ECCO) Central Production state estimate is a reconstruction of the three-dimensional time-varying ocean and sea-ice state.  Currently in Version 4 Release 3, the state estimate covers the period Jan 1, 1992 to Dec 31, 2015.  The state estimate is provided on an approximately 1-degree horizontal grid (cells range from ca. 20 to 110 km in length) and 50 vertical levels of varying thickness.

The ECCO CP state estimate has two defining features: (1) it reproduces a large set of remote sensing and in-situ observational data within their prior quantified uncertainties and (2) the dynamical evolution of its ocean circulation, hydrography, and sea-ice through time perfectly satisfies the laws of physics and thermodynamics.  The state estimate is the solution of a free-running ocean and sea-ice general circulation model and consequently can be used to assess budgets of quantities like heat, salt and vorticity.

Relation to other ocean reanalyses
==================================

ECCO state estimates share many similarities with conventional ocean reanalyses but differ in several key respects.  Both are ocean reconstructions that use observational data to fit an ocean model to the data so that the model agrees with the data in a statistical sense.  Ocean reanalyses are constructed by directly adjusting the ocean model's state to reduce its misfit to the data. Furthermore, information contained in the data is only explored forward in time. In contrast, ECCO state estimates are constructed by identifying a set of ocean model initial conditions, parameters, and atmospheric boundary conditions such that a free-running simulation of the ocean model reproduces the observations as a result of the governing equations of motion. These equations also provide a means for propagating information contained in the data back in time ("upstream" of when/where observations have been made).  Therefore, while both ocean reanalyses and ECCO state estimates reproduce observations of ocean variability, only ECCO state estimates provide an explanation for the underlying physical causes and mechanisms responsible for bringing them into existence (e.g., Stammer et al. 2017).

Conservation properties of ECCO state estimates
===============================================

By design, ECCO state estimates perfectly satisfy the laws of physics and thermodynamics and therefore conserve heat, salt, volume, and momentum (Wunsch and Heimbach, 2007, 2013).  Indeed, it is because of these conservation properties that ECCO state estimates are used to investigate the origins of ocean heat, salt, mass, sea-ice, and regional sea level variability (e.g., Wunsch et al., 2007; Köhl and Stammer, 2008; Piecuch and Ponte 2011; Fenty and Heimbach, 2013b; Wunsch and Heimbach 2013, 2014; Fukumori and Wang 2013; Buckley et al. 2014, 2015; Forget and Ponte, 2015, Fenty et al., 2017; Piecuch et al. 2017).  

*******************************************************
How is the ECCO Central Production State Estimate Made?
*******************************************************

The ECCO ocean reanalysis system is a mature, state-of-the-art data tool for synthesizing a diverse Earth System observations, including satellite and in-situ data, into a complete and physically-consistent description of Earth’s global ocean and sea-ice state.  Ocean observations, even those from satellites with global coverage, are still sparse in both space and time, relative to the inherent scales of ocean variability.  The ECCO reanalysis system is able to reconstruct the global ocean and sea-ice state by synthesizing hundreds of millions of sparse and heterogeneous ocean and sea-ice satellite and in-situ data with an ocean and sea-ice general circulation model.  Through iterative high-dimension nonlinear optimization using the adjoint of the ocean and sea-ice model the ECCO reanalysis system identifies a particular solution to the system of equations describing ocean and sea-ice dynamics and thermodynamics that reproduces a set of constraining observations in a least-squares sense.

By simultaneously integrating numerous diverse and heterogeneous data streams into the dynamically-consistent framework of the physical model we make optimal use of the data [Munk and Wunsch 1982].  Users of the ECCO reanalysis are not only provided a comprehensive description of the Earth’s changing ocean and sea-ice states but also information about the underlying physical processes responsible for driving those changes.


*************************
What fields are provided? 
*************************

The complete state estimate consists of a 53 ocean, 13 sea-ice variables, and 8 atmosphere state variables that are the output from a free-running ocean and sea-ice general circulation model.  These variables have been packaged as NetCDF files.

Temporal Frequency
==================

All three-dimensional ocean and sea-ice fields are provided as monthly averages and a select number of two-dimensional ocean and sea-ice fields are provided as daily averages.  Atmospheric state fields are provided as 6-hourly records.  In addition, potential temperature (theta), salinity, and free surface height anomaly at the ocean/sea-ice interface (etan) are provided as monthly snapshots to support budget closure calculations.

Spatial Layout 
==============

Ocean, sea-ice, and atmosphere fields are provided in two spatial layouts:
1. The 13-tile *native* lat-lon-cap 90 (llc90) grid
1. A 1° x 1° latitude and longitude grid

13-tile *native* lat-lon-cap 90 grid
------------------------------------

The lat-lon-cap (llc) is the decomposition of the spherical Earth into a Cartesian curvilinear coordinate system required by the ocean model.  It is a topologically non-trivial cubed-sphere rendering in the northern hemisphere and a dipolar grid in the southern hemisphere.
This decomposition consists of 13 tiles, each with 90x90 grid cells in the horizontal and 50 vertical levels.  The effective resolution of the llc90 grid is 1°.

.. figure:: ../figures/llc90.png
    :width: 200px
    :align: center
    :height: 100px
    :alt: alternate text
    :figclass: align-center


llc90 monthly-averaged ocean and sea-ice fields: ftp://ecco.jpl.nasa.gov/Version4/Release3/nctiles_monthly/README

llc90 monthly-snapshot ocean and sea-ice fields: ftp://ecco.jpl.nasa.gov/Version4/Release3/nctiles_monthly_snapshots/README

llc90 daily-averaged ocean and sea-ice fields: ftp://ecco.jpl.nasa.gov/Version4/Release3/nctiles_daily/README

llc90 6-hourly atmosphere fields: ftp://ecco.jpl.nasa.gov/Version4/Release3/input_forcing/README


*interpolated* 1° x 1° latitude-longitude
-----------------------------------------

Select monthly-average fields from the *native* lat-lon-cap model output have been interpolated to a more user-friendly 1° latitude-longitude grid.  

1° x 1° monthly-averaged ocean, sea-ice, and atmosphere fields: 
ftp://ecco.jpl.nasa.gov/Version4/Release3/interp_monthly/README


