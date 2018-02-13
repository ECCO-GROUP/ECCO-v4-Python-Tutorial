**************************
Python and Python Packages
**************************

Our tutorials are written for Python 2.7 and use several modules including **ECCOv4-py**, a module specifically written to facilitate loading, plotting, and performing basic *unary* and *binary* operations on the state estimate fields.  

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

Of particular interest to this tutorial are the following Python packages (libraries)

NumPy
^^^^^
From the `NumPy website <http://www.numpy.org/>`_

"NumPy is the fundamental package for scientific computing with Python. It contains among other things:

a powerful N-dimensional array object
sophisticated (broadcasting) functions
tools for integrating C/C++ and Fortran code
useful linear algebra, Fourier transform, and random number capabilities"

Matplotlib
^^^^^^^^^^
From the `Matplotlib website <https://matplotlib.org/>`_

"Matplotlib is a Python 2D plotting library which produces publication quality figures in a variety of hardcopy formats and interactive environments across platforms. Matplotlib can be used in Python scripts, the Python and IPython shell, the jupyter notebook, web application servers, and four graphical user interface toolkits.

Matplotlib tries to make easy things easy and hard things possible. You can generate plots, histograms, power spectra, bar charts, errorcharts, scatterplots, etc., with just a few lines of code. For examples, see the sample plots and thumbnail gallery.

For simple plotting the pyplot module provides a MATLAB-like interface, particularly when combined with IPython. For the power user, you have full control of line styles, font properties, axes properties, etc, via an object oriented interface or via a set of functions familiar to MATLAB users."


xarray
^^^^^^
From the `xarray website <http://xarray.pydata.org/en/stable/why-xarray.html>`_

"Adding dimensions names and coordinate indexes to numpy's ndarray_ makes many
powerful array operations possible:

-  Apply operations over dimensions by name: ``x.sum('time')``.
-  Select values by label instead of integer location:
   ``x.loc['2014-01-01']`` or ``x.sel(time='2014-01-01')``.
-  Mathematical operations (e.g., ``x - y``) vectorize across multiple
   dimensions (array broadcasting) based on dimension names, not shape.
-  Flexible split-apply-combine operations with groupby:
   ``x.groupby('time.dayofyear').mean()``.
-  Database like alignment based on coordinate labels that smoothly
   handles missing values: ``x, y = xr.align(x, y, join='outer')``.
-  Keep track of arbitrary metadata in the form of a Python dictionary:
   ``x.attrs``.


The N-dimensional nature of xarray's data structures makes it suitable for dealing
with multi-dimensional scientific data, and its use of dimension names
instead of axis labels (``dim='time'`` instead of ``axis=0``) makes such
arrays much more manageable than the raw numpy ndarray: with xarray, you don't
need to keep track of the order of arrays dimensions or insert dummy dimensions
(e.g., ``np.newaxis``) to align arrays."


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

    conda install -c conda-forge xarray
    conda install -c anaconda netcdf4
    conda install -c conda-forge pyresample
    conda install -c anaconda basemap 

If you don't install Anaconda then you can install these packages by following the instructions on each project's website.


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