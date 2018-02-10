*****************
Tutorial Overview
*****************


What are the tutorials?
-----------------------


What are they *not*?
--------------------

Python lessons.



Where can I get the tutorials as Jupyter notebooks?
---------------------------------------------------

All tutorials are available as Jupyter notebooks here:
https://github.com/ECCO-GROUP/ECCO-v4-Python-Tutorial/tree/master/notebooks



Summary of Tutorials
--------------------




What if I don't like the way you do X?
--------------------------------------

That's ok.  


What if I find a mistake?
-------------------------

Definitely please tell us!  You can write Ian.Fenty at jpl.nasa.gov


Problems installing libraries?
------------------------------

Sometimes the pip Python package manager doesn't download the latest version of libraries from pypi.  On osx machines this could be because old versions of the libraries are stored in a local cache directory.  Deleting the pip cache directory and Using the "--no-cache-dir" option seems to solve this problem. 

::

    rm -fr ~/Library/Caches/pip/*
    pip install ecco_v4_py --no-cache-dir


Sometimes libraries get updated.  To update packages installed with conda use the the following command:

::

    conda update anaconda


