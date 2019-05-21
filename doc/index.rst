.. ECCO Version 4 Tutorial documentation master file, created by
   sphinx-quickstart on Mon Jan 22 20:58:25 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the ECCO Version 4 Tutorial
======================================

This website contains a set of tutorials about how to use the ECCO Central Production Version 4 (ECCO v4) global ocean and sea-ice state estimate.  The tutorials were written in Python and make use of the `ecco_v4_py`_ Python library, a library written specifically for loading, plotting, and analyzing ECCO v4 state estimate fields.


Additional Resources
--------------------

The ECCO v4 state estimate is the output of a free-running simulation of a global ca. 1-degree configuration of the MITgcm.  Prior to public release, the model output files model are assembled into NetCDF files.  If you would like to work directly with the flat binary "MDS" files provided by the model then take a look at the `xmitgcm`_ Python package.  The `xgcm`_ Python package provides tools for operating on model output fields loaded with `xmitgcm`_.   If you wish to analyze the MITgcm model output using Matlab then we strongly recommend the extensive set of tools provided by the `gcmfaces`_ toolbox.  

The `ecco_v4_py`_ package used in this tutorial was inspired by both the `xmitgcm`_ package and `gcmfaces`_ toolbox.

.. _ecco_v4_py : https://pypi.python.org/pypi/ecco_v4_py/
.. _xmitgcm : http://xmitgcm.readthedocs.io/en/latest/index.html
.. _xgcm : https://github.com/xgcm/xgcm
.. _gcmfaces : https://github.com/gaelforget/gcmfaces

.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   intro
   fields
   Installing_Python_and_Python_Packages
   Downloading_the_ECCO_v4_state_estimate
   Tutorial_Introduction

.. toctree::
   :maxdepth: 2
   :caption: ECCO Data Structures
   
   ECCO_v4_data_structure_basics
   ECCO_v4_Coordinates_and_Dimensions_of_ECCOv4_NetCDF_files

.. toctree::
   :maxdepth: 2
   :caption: Input/Output, Data Structure Manipulation

   ECCO_v4_Loading_ECCO_NetCDF_state_estimate_and_model_grid_fields.ipynb
   ECCO_v4_Combining_Multiple_Datasets
   ECCO_v4_Saving_Datasets_and_DataArrays_to_NetCDF 

.. toctree::
   :maxdepth: 2
   :caption: Operating on Data Variables

   ECCO_v4_Accessing_and_Subsetting_Variables
   ECCO_v4_Operating_on_Numpy_Arrays

.. toctree::
   :maxdepth: 2
   :caption: Plotting

   ECCO_v4_Plotting_Tiles

.. toctree::
   :maxdepth: 2
   :caption: Scalar Calculations

   ECCO_v4_Example_calculations_with_scalar_quantities

.. toctree::
   :maxdepth: 2
   :caption: More Advanced Calculations

   compute_mht
 

.. toctree::
   :maxdepth: 2
   :caption: UNDER CONSTRUCTION: Budget and Vorticity Calculations

   VectorCalculus_ECCO_barotropicVorticity
   ecco_budgets

.. toctree::
   :maxdepth: 2
   :caption: Support 

   support

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
