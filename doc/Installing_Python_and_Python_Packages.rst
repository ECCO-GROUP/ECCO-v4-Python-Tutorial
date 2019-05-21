**************************
Python and Python Packages
**************************

Our Python tutorial is compatible with Python 2.7 and 3.  It relies on several packages including **ecco_v4-py** which include codes to facilitate loading, plotting, and performing calculations on ECCOv4 state estimate fields.  

.. _in-python:

Why Python?
-----------

`Python <https://www.python.org/>`_ is an easy to learn, open source programming language.  In addition to the standard language library, there are thousands of free third-party modules (code libraries) available on the `Python Package Index <https://pypi.org/>`_ (PyPI).  Unlike commerical numerical computing environments like Matlab and IDL, Python is free for everyone to use.  In addition, Python code can be run on just any platform whether Windows, Linux, or OS X.

Here are some links to help you learn more about Python.

- `Python 3.x Documentation <https://docs.python.org/3/>`_
- `Python 3 Tutorial <https://docs.python.org/3/tutorial/>`_ 
- `Python 2.7 Documentation <https://docs.python.org/2.7/>`_   (Note Python 2.7 will soon be unsupported)
- `Python 2.7 Tutorial <https://docs.python.org/2.7/tutorial/index.html>`_ 
- `Scientific Python Lectures <http://www.scipy-lectures.org/>`_ 
- `Using the NumPy module for Matlab Users <http://scipy.github.io/old-wiki/pages/NumPy_for_Matlab_Users>`_ 
- `Learning Python with Anaconda <https://www.datacamp.com/learn-python-with-anaconda>`_ 


.. _in-Installing:

Installing Python and the Anaconda Distribution
-----------------------------------------------

Python
^^^^^^
The latest installers for Python for many platforms can be found on the `Python website <https://www.python.org/downloads/release/python-2714/>`_.


Anaconda
^^^^^^^^
Python code can be written in any text editor and run from the command line.  Third-party modules can be manually installed from the command line using the `pip`_ package manager.  

Python code can also be written and executed in an interactive environment (integrated development environment, IDE) similar to the Matlab console.  For scientific computing, the `Anaconda`_ Python distribution is quite convenient because it comes with a `large collection`_ of useful modules, a good open source IDE, `Spyder`_., and the ability to open and execute `Jupyter Notebooks`_

The latest installers for the Anaconda Distribution can be found on the `Anaconda website`_

.. _Anaconda : https://www.anaconda.com/
.. _Anaconda website: https://www.anaconda.com/download/
.. _pip : https://pypi.python.org/pypi/pip
.. _large collection : https://docs.anaconda.com/anaconda/packages/pkg-docs
.. _Spyder : https://pythonhosted.org/spyder/index.html
.. _P2v3 : https://www.digitalocean.com/community/tutorials/python-2-vs-python-3-practical-considerations-2
.. _Jupyter Notebooks : https://jupyter.org/


.. _in-libraries:

Installing Required Python Packages
-----------------------------------

After installing Anaconda the following packages must be installed: 
*netcdf4*, *cartopy*, *pyresample*, *xarray*, *xmitgcm*, *xgcm*
  

Below are two **options** or installing these packages. PICK ONE!


Option 1: Conda (recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. attention::
    Conda is recommended because it automatically installs the GEOS (Geometry Engine, Open Source) library which is needed to make projection plots.


.. code-block:: bash

    conda install netcdf4
    conda install -c scitools cartopy
    conda install -c conda-forge pyresample

.. code-block:: bash

    pip install xarray
    pip install xmitgcm
    pip install xgcm


Option 2: *pip* (not recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. DANGER::
    The Python module Cartopy requires the GEOS (Geometry Engine, Open Source) library.  Instructions for installing this library can be found on the `geos website`_.   Some users have reported difficulties  installing GEOS libraries on their platforms.  For that reason, we recommend using Conda (Option 1).  


.. code-block:: bash

    pip install netcdf4
    pip install pyresample
    pip install cartopy
    pip install xarray
    pip install xmitgcm
    pip install xgcm



Downloading the *ecco_v4_py* Python Package
-------------------------------------------

The *ecco_v4_py* package is a library of routines that are helpful for analyzing the ECCO v4 state estimate.  It is stored on the `github repository`_ 


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
This method gets you the source code but if you make changes it is harder to share those changes with the community.

.. code-block:: bash
	
    > mkdir ~/ECCOv4-py
    > cd ~/ECCOv4-py
    > wget https://github.com/ECCO-GROUP/ECCOv4-py/archive/master.zip
    > unzip master.zip
    > rm master.zip

Of course you may want to use this method if you don't have access to git.

Option 3: Use the *pip* Python package tool (not recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you use *pip* to install the *ecco_v4_py* package the source code will be installed in your Python library directory from https://pypi.org/project/ecco-v4-py/.  This method is OK if you don't plan to look at or modify the library code.   

.. code-block:: bash
	
    pip install ecco_v4_py


Using the *ecco_v4_py* Python Package in your programs
------------------------------------------------------

If you use Options 1 or 2 to download the *ecco_v4_py* source code then you must tell Python the location of the files before Python can it.  This is easy, you just you just have to remember to do it at the top of all of your programs!  

Assuming you downloaded the *ecco_v4_py* routines to ``/home/username/ECCOv4-py`` then simply add these three lines to the top of your Python programs (or Jupyter Notebooks)

.. code-block:: python

    import sys
	sys.path.append('/home/username/ECCOv4-py')
	import ecco_v4_py as ecco


If you used Method 3 (pip install) then the *ecco_v4_py* library will be automatically installed and will be ready to import into your Python program via the following commands:  

.. code-block:: python

    import ecco_v4_py as ecco

.. _geos website: https://trac.osgeo.org/geos

.. _github repository: https://github.com/ECCO-GROUP/ECCOv4-py/tree/master/ecco_v4_py