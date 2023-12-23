.. ECCO Version 4 Tutorial documentation master file, created by
   sphinx-quickstart on Mon Jan 22 20:58:25 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the ECCO Version 4 Tutorial
======================================

This website contains a set of tutorials about how to use the ECCO Central Production Version 4 (ECCO v4) global ocean and sea-ice state estimate.  The tutorials were written in Python and make use of the `ecco_v4_py` Python library, a library written specifically for loading, plotting, and analyzing ECCO v4 state estimate fields.


Additional Resources
--------------------

The ECCO v4 state estimate is the output of a free-running simulation of a global ca. 1-degree configuration of the MITgcm.  Prior to public release, the model output files model are assembled into NetCDF files.  If you would like to work directly with the flat binary "MDS" files provided by the model then take a look at the `xmitgcm`_ Python package.  The `xgcm`_ Python package provides tools for operating on model output fields loaded with `xmitgcm`_.   If you wish to analyze the MITgcm model output using Matlab then we recommend the `gcmfaces`_ toolbox.  

The `ecco_v4_py`_ package used in this tutorial was inspired by the `xmitgcm`_ package and `gcmfaces`_ toolboxes.

.. _ecco_v4_py : https://github.com/ECCO-GROUP/ECCOv4-py
.. _xmitgcm : http://xmitgcm.readthedocs.io/en/latest/index.html
.. _xgcm : https://github.com/xgcm/xgcm
.. _gcmfaces : https://github.com/gaelforget/gcmfaces

.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   intro
   fields
   Installing_Python_and_Python_Packages
   Downloading_ECCO_Datasets_from_PODAAC_Python.ipynb
   Downloading_Subsets_of_ECCO_Datasets.ipynb
   Tutorial_wget_Command_Line_HTTPS_Downloading_ECCO_Datasets_from_PODAAC
   Tutorial_Introduction

.. toctree::
   :maxdepth: 2
   :caption: ECCO Data Structures
   
   ECCO_v4_data_structure_basics
   ECCO_v4_Coordinates_and_Dimensions_of_ECCOv4_NetCDF_files

.. toctree::
   :maxdepth: 2
   :caption: Input/Output, Data Structure Manipulation

   ECCO_v4_Loading_the_ECCOv4_native_model_grid_parameters.ipynb
   ECCO_v4_Loading_the_ECCOv4_state_estimate_fields_on_the_native_model_grid.ipynb
   ECCO_v4_Loading_LLC_compact_binary_files.ipynb
   ECCO_v4_Combining_Multiple_Datasets.ipynb
   ECCO_v4_Saving_Datasets_and_DataArrays_to_NetCDF.ipynb

.. toctree::
   :maxdepth: 2
   :caption: Operating on Data Variables

   ECCO_v4_Accessing_and_Subsetting_Variables
   ECCO_v4_Operating_on_Numpy_Arrays

.. toctree::
   :maxdepth: 2
   :caption: Plotting & Interpolation

   ECCO_v4_Plotting_Tiles.ipynb
   ECCO_v4_Interpolating_Fields_to_LatLon_Grid.ipynb

.. toctree::
   :maxdepth: 2
   :caption: Scalar Calculations

   ECCO_v4_Example_calculations_with_scalar_quantities.ipynb

.. toctree::
   :maxdepth: 2
   :caption: Intro to PO Tutorials
   
   Intro_to_PO_start
   Geostrophic_balance.ipynb
   Thermal_wind.ipynb
   Steric_height.ipynb

.. toctree::
   :maxdepth: 2
   :caption: More Advanced Calculations

   ECCO_v4_Example_MHT.ipynb
   ECCO_v4_Example_OSNAP.ipynb
   ECCO_v4_Volume_budget_closure.ipynb 
   ECCO_v4_Heat_budget_closure.ipynb 
   ECCO_v4_Salt_and_salinity_budget.ipynb
   ECCO_v4_Calculating_the_ECCOv4_ocean_thermal_forcing.ipynb

.. toctree::
   :maxdepth: 2
   :caption: Support 

   support

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
