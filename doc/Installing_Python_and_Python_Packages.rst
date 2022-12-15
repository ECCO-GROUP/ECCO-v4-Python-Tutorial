**************************
Python and Python Packages
**************************

The ECCO Python tutorial is compatible with Python 3.  It relies on several packages including **ecco_v4-py** which include codes to facilitate loading, plotting, and performing calculations on ECCOv4 state estimate fields.  

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

The *ecco_v4_py* package is a library of routines for analyzing the ECCO the Version 4 state estimate. The latest version can always be found on our `github repository`_ 


Below are three **options** or installing the *ecco_v4_py* Python package.

.. attention::

    Use only one of the options below!


Option 1: Clone into the repository using git (recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Cloning into the *ecco_v4_py* repository using `git` 
is recommended because 

a) you can easily see and modify source code
b) you share your improvements with the community.

To use `git` to clone into the project simply run the following commands
(in the example below the Python files will go into ~/ECCOv4-py/)

.. code-block:: bash

    > mkdir ~/ECCOv4-py
    > cd ~/ECCOv4-py
    > git clone https://github.com/ECCO-GROUP/ECCOv4-py.git


Option 2: Download the repository using git (less recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This method downloads the source code but if you make changes it is harder to share those changes with the community using git.

.. code-block:: bash
	
    > mkdir ~/ECCOv4-py
    > cd ~/ECCOv4-py
    > wget https://github.com/ECCO-GROUP/ECCOv4-py/archive/master.zip
    > unzip master.zip
    > rm master.zip

Option 3: Use the *conda* package manager (less recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*ecco_v4_py* is available via the *conda* package management system (see https://anaconda.org/conda-forge/ecco_v4_py ) 
Installing using conda should install all of the required dependencies.

.. code-block:: bash
	
    conda install ecco_v4_py
    
    
If for some reason the above command returns an error, include the ``-c`` option to point to the channel where the package is found (*conda-forge*).

.. code-block:: bash

    conda install -c conda-forge ecco_v4_py


Option 4: Use the *pip* package manager (not at all recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*ecco_v4_py* is available via the *pip* package manager (see https://pypi.org/project/ecco-v4-py/ ) Before using pip, you must first install the PROJ and GEOS libraries (see next section). 

.. code-block:: bash
	
    pip install ecco_v4_py
    
    

Installing Dependencies
-----------------------
   
.. DANGER::
    While conda is recommended because it automatically installs the required the GEOS (Geometry Engine) and PROJ (generic coordinate transformation software) binary libraries, you can install those libraries yourself. 

Instructions for installing the GEOS library can be found on the `geos website`_.  

Instructions for installing the PROJ library can be found on the `proj website`_.  

Some users have reported difficulties installing these libraries on their platforms.  For that reason, we recommend using Options 1-3.



Using the *ecco_v4_py* in your programs
---------------------------------------

Assuming you downloaded the *ecco_v4_py* routines to ``/home/username/ECCOv4-py`` then simply add these three lines to the top of your Python programs (or Jupyter Notebooks)

.. code-block:: python

    import sys
    sys.path.append('/home/username/ECCOv4-py')
    import ecco_v4_py as ecco


If you you installed the package using pip then the *ecco_v4_py* library will be automatically installed and will be ready to import into your Python program via the following commands:  

.. code-block:: python

    import ecco_v4_py as ecco


.. _proj website: https://proj.org/install.html
.. _geos website: https://libgeos.org/
.. _github repository: https://github.com/ECCO-GROUP/ECCOv4-py/tree/master/ecco_v4_py
