**************************
Python and Python Packages
**************************

Our Python tutorial is compatible with Python 3.  It relies on several packages including **ecco_v4-py** which include codes to facilitate loading, plotting, and performing calculations on ECCOv4 state estimate fields.  

.. _in-python:

Why Python?
-----------

`Python <https://www.python.org/>`_ is an easy to learn open source programming language.  In addition to the standard language library, there are thousands of free third-party modules (code libraries) available on code repositories such as `Python Package Index <https://pypi.org/>`_ (PyPI), `_Conda <https://anaconda.org/anaconda/repo>`_ and `Conda Forge <https://conda-forge.org/feedstock-outputs/>`_.  Unlike commerical numerical computing environments like Matlab and IDL, Python is free for everyone to use.  In addition, Python code can be run on multiple platforms such as Windows, Linux, and OS X.

Here are some links to help you learn more about Python.

- `Python 3.x Documentation <https://docs.python.org/3/>`_
- `Python 3 Tutorial <https://docs.python.org/3/tutorial/>`_ 
- `Scientific Python Lectures <http://www.scipy-lectures.org/>`_ 
- `Using the NumPy module for Matlab Users <http://scipy.github.io/old-wiki/pages/NumPy_for_Matlab_Users>`_ 
- `Learning Python with Anaconda <https://www.datacamp.com/learn-python-with-anaconda>`_ 


.. _in-Installing:

Installing Python
-----------------------------------------------

There are several ways of installing Python on your machine. You can install compiled binaries directly from the  `Python website <https://www.python.org/downloads/release/python-2714/>`_, or one can install via a package manager such as Anaconda or Miniconda. I personally find the Anaconda or Miniconda route to be simplest. 

Anaconda
^^^^^^^^
For scientific computing, the `Anaconda`_ Python distribution is quite convenient because it comes with a `large collection`_ of useful modules, a good open source IDE, `Spyder`_., and the ability to open and execute `Jupyter Notebooks`_

The latest installers for the Anaconda Distribution can be found on the `Anaconda website`_

.. _Anaconda : https://www.anaconda.com/
.. _Anaconda website: https://www.anaconda.com/download/
.. _pip : https://pypi.python.org/pypi/pip
.. _large collection : https://docs.anaconda.com/anaconda/packages/pkg-docs
.. _Spyder : https://pythonhosted.org/spyder/index.html
.. _P2v3 : https://www.digitalocean.com/community/tutorials/python-2-vs-python-3-practical-considerations-2
.. _Jupyter Notebooks : https://jupyter.org/


.. _in-libraries:



Downloading the *ecco_v4_py* Python Package
-------------------------------------------

The *ecco_v4_py* package is a library of routines that are helpful for analyzing the ECCO the Version 4 state estimate.  The latest version can always be found on our `github repository`_ 


Below are three **options** or installing the *ecco_v4_py* Python package.

.. attention::

    Use only one of the options below!


Option 1: Clone into the repository using git (recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Cloning into the *ecco_v4_py* repository using `git` 
is recommended because 

a) you can easily see and modify the ecco_v4_py source code
b) you can improve the source code and share your improvements with the community.

To use `git` to clone into the project simply run the following commands
(in the example below the Python files will go into ~/ECCOv4-py/)

.. code-block:: bash

    > mkdir ~/ECCOv4-py
    > cd ~/ECCOv4-py
    > git clone https://github.com/ECCO-GROUP/ECCOv4-py.git


Option 2: Download the repository using git (not recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This method gets you the source code but if you make changes it is harder to share those changes with the community. Use this method if you don't have access to git.

.. code-block:: bash
	
    > mkdir ~/ECCOv4-py
    > cd ~/ECCOv4-py
    > wget https://github.com/ECCO-GROUP/ECCOv4-py/archive/master.zip
    > unzip master.zip
    > rm master.zip



Option 3: Use the *pip* Python package tool (not at all recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you use *pip* to install the *ecco_v4_py* package the source code will be installed in your Python library directory from https://pypi.org/project/ecco-v4-py/.  This method is OK if you don't plan to look at or modify the library code, but what fun would that be?  Also, YMMV with respect to installing dependencies.

.. code-block:: bash
	
    pip install ecco_v4_py
    
    
       

Installing Required Python Packages
-----------------------------------

The following additional packages must be installed 

.. code-block:: bash

  - aiohttp
  - codecov
  - cartopy>=0.18.0
  - cmocean
  - dask
  - docrep
  - fsspec
  - future
  - geos
  - matplotlib
  - netcdf4
  - numpy
  - pathlib
  - proj
  - pytest
  - python=3.8
  - pytest-cov
  - pyresample
  - scipy
  - xarray
  - xgcm>=0.5.0
  - xmitgcm>=0.5.1

Below are several **options** or installing these packages. PICK ONE!


Option 1: *Conda*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. attention::
    Conda is recommended because it automatically installs the GEOS (Geometry Engine, Open Source) library which is needed to make projection plots.

To simply install all required dependencies using conda, create a new conda environment using the ECCOv4-py environment YAML.


.. code-block:: bash

   cd your_ECCOv4_py_directory/
   conda env create --name ECCOv4_py --file ci/environment-py38.yml

To activate the conda environment with the required dependencies call:

.. code-block:: bash

   conda activate ECCOv4_py
   
   
Alternatively, one can install the dependencies one at a time into an existing conda environment:

.. code-block:: bash

  conda activate yourExistingCondaEnvironment
  
  conda install -c conda-forge aiohttp
  conda install -c conda-forge codecov
  conda install -c conda-forge cartopy 
  conda install -c conda-forge cmocean
  conda install -c conda-forge dask
  conda install -c conda-forge docrep
  conda install -c conda-forge fsspec
  conda install -c conda-forge future
  conda install -c conda-forge geos
  conda install -c conda-forge matplotlib
  conda install -c conda-forge netcdf4
  conda install -c conda-forge numpy
  conda install -c conda-forge pathlib
  conda install -c conda-forge proj
  conda install -c conda-forge pytest
  conda install -c conda-forge pytest-cov
  conda install -c conda-forge pyresample
  conda install -c conda-forge scipy
  conda install -c conda-forge xarray
  conda install -c conda-forge xgcm 
  conda install -c conda-forge xmitgcm


Option 2: *pip* alone (not recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. DANGER::
    The Python module Cartopy requires the GEOS (Geometry Engine, Open Source) library.  Instructions for installing this library can be found on the `geos website`_.   Some users have reported difficulties installing GEOS libraries on their platforms.  For that reason, we recommend using Conda (Option 1).   

To install geos, following instructions here:
https://pygeos.readthedocs.io/en/latest/installation.html


.. code-block:: bash

  pip install geos
  pip install aiohttp
  pip install codecov
  pip install cmocean
  pip install dask
  pip install docrep
  pip install fsspec
  pip install future
  pip install matplotlib
  pip install netcdf4
  pip install numpy
  pip install pathlib
  pip install pytest
  pip install pytest-cov
  pip install pyresample
  pip install scipy
  pip install xarray
  pip install xgcm
  pip install xmitgcm
  pip install proj
  pip install cartopy
  
 



Using the *ecco_v4_py* in your programs
------------------------------------------------------

Assuming you downloaded the *ecco_v4_py* routines to ``/home/username/ECCOv4-py`` then simply add these three lines to the top of your Python programs (or Jupyter Notebooks)

.. code-block:: python

    import sys
    sys.path.append('/home/username/ECCOv4-py')
    import ecco_v4_py as ecco


If you you installed the package using pip then the *ecco_v4_py* library will be automatically installed and will be ready to import into your Python program via the following commands:  

.. code-block:: python

    import ecco_v4_py as ecco

.. _geos website: https://trac.osgeo.org/geos

.. _github repository: https://github.com/ECCO-GROUP/ECCOv4-py/tree/master/ecco_v4_py
