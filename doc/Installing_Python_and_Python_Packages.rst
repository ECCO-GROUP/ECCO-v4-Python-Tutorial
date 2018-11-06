**************************
Python and Python Packages
**************************

Our tutorials are written for Python 2.7 and use several modules including **ecco_v4-py**, a module specifically written to facilitate loading, plotting, and performing *unary* and *binary* operations on the state estimate fields.  

.. _in-python:

Why Python?
-----------

`Python <https://www.python.org/>`_ is an easy to learn, open source programming language.  In addition to the standard language library, there are thousands of free third-party modules (code libraries) available on the `Python Package Index <https://pypi.org/>`_ (PyPI).  Unlike commerical numerical computing environments like Matlab and IDL, Python is free for everyone to use.  In addition, Python code can be run on just any platform whether Windows, Linux, or OS X.

Here are some links to help you learn more about Python.

- `Python 2.7 Documentation <https://docs.python.org/2.7/>`_ 
- `Python 2.7 Tutorial <https://docs.python.org/2.7/tutorial/index.html>`_ 
- `Scientific Python Lectures <http://www.scipy-lectures.org/>`_ 
- `Using the NumPy module for Matlab Users <http://scipy.github.io/old-wiki/pages/NumPy_for_Matlab_Users>`_ 
- `Learning Python with Anaconda <https://www.datacamp.com/learn-python-with-anaconda>`_ 


.. _in-Installing:

Installing Python and the Anaconda Distribution
-----------------------------------------------

Python
^^^^^^
The latest installers for Python 2.7 for many platforms can be found on the `Python website <https://www.python.org/downloads/release/python-2714/>`_.

.. note::  The code for this tutorial has been developed for Python 2.7 and has not yet been tested with `Python 3 <https://www.digitalocean.com/community/tutorials/python-2-vs-python-3-practical-considerations-2>`_.  

Anaconda
^^^^^^^^
Python code can be written in any text editor and run from the command line.  Third-party modules can be manually installed from the command line using the `pip`_ package manager.  

Alternatively, Python code can be written and executed in an interactive environment (integrated development environment, IDE) similar to what one may be familiar with in Matlab.  For scientific computing, the `Anaconda`_ Python distribution is convenient because it comes with a `large collection`_ of useful modules and a good open source IDE, `Spyder`_.

The latest installers for the Anaconda Distribution can be found on the `Anaconda website`_

.. note::  For this tutorial be sure to install the Anaconda Distribution for Python 2.7.  

.. _Anaconda : https://www.anaconda.com/
.. _Anaconda website: https://www.anaconda.com/download/
.. _pip : https://pypi.python.org/pypi/pip
.. _large collection : https://docs.anaconda.com/anaconda/packages/pkg-docs
.. _Spyder : https://pythonhosted.org/spyder/index.html
.. _P2v3 : https://www.digitalocean.com/community/tutorials/python-2-vs-python-3-practical-considerations-2

.. _in-libraries:

Installing Required Python Packages
-----------------------------------

After installing Anaconda the following packages must be installed: *xarray*, *netcdf4*, *pyresample*, *basemap*.  

.. code-block:: bash

    pip install xarray
    pip install netcdf4
    pip install pyresample
    
.. note:: *pyresample* can be installed on a mac using
    conda install -c conda-forge pyresample    
    
For *basemap*, see the installation instructions on the project page, https://github.com/matplotlib/basemapbrew 

On a mac you should consider installing the geos package using homebrew as,

.. code-block:: bash
    brew install geos 


Installing the *ecco_v4_py* Python Package
------------------------------------------

The *ecco_v4_py* package is a library of routines that are helpful for analyzing the ECCO v4 state estimate.  It is stored on the `github repository`_ 

There are three methods for installing the *ecco_v4_py* Python package, 

1. [RECOMMENDED] use *pip*, a tool for installing Python packages.

.. code-block:: bash
	
    pip install ecco_v4_py

2. Download the latest version from the git repository, https://github.com/ECCO-GROUP/ECCOv4-py/archive/master.zip

3. Use `git` to clone the project:

.. code-block:: bash
	
    git clone https://github.com/ECCO-GROUP/ECCOv4-py.git


If you use Method 1 (pip install) then the *ecco_v4_py* library will be automatically installed and will be ready to import into your Python program.  

If you use Methods 2 or 3 you'll need to take note of the location of your *ecco_v4_py* directory and add it to the Python system path in the header of your routines to allow the library to be imported:  

.. code-block:: python

	sys.path.append('/PATH/TO/YOUR/COPY/OF/ECCOv4-py/ecco_v4_py')
	import ecco_v4_py as ecco


.. _github repository: https://github.com/ECCO-GROUP/ECCOv4-py/tree/master/ecco_v4_py
