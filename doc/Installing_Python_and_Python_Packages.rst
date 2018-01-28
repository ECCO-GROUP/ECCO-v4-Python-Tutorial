*************************************
Installing Python and Python Packages
*************************************

Our tutorials are written for Python 2.7 and use several modules including **ECCOv4-py**, a module specifically written to facilitate loading, plotting, and performing basic unary and binary operations on the physical variables of the state estimate.  

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
Python code can be written in any text editor and run from the command line.  Third-party modules can be manually installed from then command line using the `pip`_ package manager.  

Alternatively, Python code can be written and executed in an interactive environment (integrated development environment, IDE) similar to what one may be familiar with in Matlab.  For scientific computing, the `Anaconda`_ Python distribution is very convienent because it comes with a `large collection`_ of useful modules and a good open source IDE, `Spyder`_.

The latest installers for the Anaconda Distribution can be found on the `Anaconda website`_

.. note::  For this tutorial be sure to install the Anaconda Distribution for Python 2.7.  

.. _Anaconda website: https://www.anaconda.com/download/
.. _pip : https://pypi.python.org/pypi/pip
.. _large collection : https://docs.anaconda.com/anaconda/packages/pkg-docs
.. _Spyder : https://pythonhosted.org/spyder/index.html
.. _P2v3 : https://www.digitalocean.com/community/tutorials/python-2-vs-python-3-practical-considerations-2

.. _in-libraries:

Installing Required Python Packages
-----------------------------------

The following packages must be installed: *xarray*, *netcdf4*, *pyresample*, *basemap*.  After installing Anaconda, you can install these packages with the following command:

.. code-block:: bash

    conda install -c conda-forge xarray
    conda install -c anaconda netcdf4
    conda install -c conda-forge pyresample
    conda install -c anaconda basemap 


The development version of the *ecco_v4_py* can be found on the `github project page`_ 

From there you have two options:

1. Download https://github.com/ECCO-GROUP/ECCOv4-py/archive/master.zip

2. Use `git` to clone the project:

.. code-block:: bash
	
    git clone git@github.com:ECCO-GROUP/ECCOv4-py.git

.. _github project page: https://github.com/ECCO-GROUP/ECCOv4-py/tree/master/ecco_v4_py

Take note of the location of your *ecco_v4_py* directory.  You'll need to specify its location during the *import* header of your Python programs.  
